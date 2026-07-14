#!/usr/bin/env python
from statistics import mean, median
from io import TextIOWrapper
from zipfile import ZipFile
from glob import glob
import pandas as pd
import numpy as np
import requests
import argparse
import json
import csv
import sys
import re
import os


def sample_key(sample):
    """Normalize sample identifiers for dict lookups and ordering."""
    if sample is None or (isinstance(sample, float) and pd.isna(sample)):
        return ""
    return str(sample).strip()

def mqc_base_sample_name(sample_name_str):
    """Strip paired-end / raw suffixes from a MultiQC sample name."""
    for suffix in (' Read 1', ' Read 2', '_R1_raw', '_R2_raw', '_R1', '_R2'):
        if suffix in sample_name_str:
            return sample_name_str.split(suffix)[0].strip('_')
    return sample_name_str

def mqc_read_end(sample_name_str):
    """Return 'f' or 'r' for paired-end MultiQC sample names."""
    if ' Read 2' in sample_name_str or '_R2' in sample_name_str:
        return 'r'
    if ' Read 1' in sample_name_str or '_R1' in sample_name_str:
        return 'f'
    return 'f'

def fastqc_max_sequence_length(seq_len_value):
    """Parse FastQC Sequence length field (e.g. '35-76') into max read length."""
    if seq_len_value is None or (isinstance(seq_len_value, float) and pd.isna(seq_len_value)):
        return None
    seq_len_str = str(seq_len_value).strip()
    if not seq_len_str:
        return None
    if '-' in seq_len_str:
        try:
            return int(float(seq_len_str.split('-')[-1]))
        except (ValueError, TypeError):
            return None
    try:
        return int(float(seq_len_str))
    except (ValueError, TypeError):
        return None

def resolve_multiqc_read_length(all_fields, prefix='raw'):
    """Resolve read length from FastQC max sequence length."""
    max_lengths = []
    for suffix in ('_f', '_r'):
        key = f'{prefix}_max_sequence_length{suffix}'
        if all_fields.get(key) not in (None, ''):
            try:
                max_lengths.append(int(float(all_fields[key])))
            except (ValueError, TypeError):
                pass
    return max(max_lengths) if max_lengths else None

def lookup_raw_fastqc_sample(base_name, end, raw_fastqc, mqc_sample=None):
    """Find the saved FastQC sample key for a grouped read end."""
    candidates = []
    if mqc_sample:
        candidates.append(mqc_sample)
    if end == 'r':
        candidates.extend([f"{base_name}_R2", f"{base_name} Read 2"])
    else:
        candidates.extend([f"{base_name}_R1", f"{base_name} Read 1", base_name])
    for key in candidates:
        if key in raw_fastqc:
            return key
    return None

def iter_multiqc_plot_series(plot_section):
    """Yield (sample_name, xy_pairs) from MultiQC plot datasets (1.12 list/{data} or 1.26+ {lines}/{pairs})."""
    datasets = plot_section.get('datasets') if isinstance(plot_section, dict) else None
    if not datasets:
        return
    dataset0 = datasets[0]
    if isinstance(dataset0, list):
        series_list = dataset0
    elif isinstance(dataset0, dict):
        series_list = dataset0.get('lines') or dataset0.get('data') or []
    else:
        return
    for series in series_list:
        if not isinstance(series, dict):
            continue
        name = series.get('name')
        pairs = series.get('pairs', series.get('data'))
        if name is None or not pairs:
            continue
        yield name, pairs

def get_runsheet_order(runsheet_path):
    """Read runsheet and return ordered list of sample names"""
    if not runsheet_path or not os.path.exists(runsheet_path):
        return None
    try:
        df = pd.read_csv(runsheet_path, dtype=str)
        if 'Sample Name' in df.columns:
            return [sample_key(s) for s in df['Sample Name'].tolist()]
    except Exception as e:
        print(f"WARNING: Error reading runsheet: {str(e)}")
    return None

def get_runsheet_aliases(runsheet_path):
    """Map Original Sample Name (ISA) and other aliases to runsheet Sample Name."""
    aliases = {}
    if not runsheet_path or not os.path.exists(runsheet_path):
        return aliases
    try:
        df = pd.read_csv(runsheet_path, dtype=str)
    except Exception as e:
        print(f"WARNING: Error reading runsheet for sample aliases: {str(e)}")
        return aliases
    if 'Sample Name' not in df.columns:
        return aliases
    for _, row in df.iterrows():
        canonical = sample_key(row['Sample Name'])
        if not canonical:
            continue
        aliases[canonical] = canonical
        if 'Original Sample Name' in df.columns:
            orig = sample_key(row.get('Original Sample Name'))
            if orig and orig != canonical:
                aliases[orig] = canonical
    return aliases

def canonical_sample_name(sample, aliases):
    """Return the processing sample name for a sample identifier."""
    s = sample_key(sample)
    if not s:
        return s
    return aliases.get(s, s)

def remap_parser_dict(data, aliases):
    """Re-key parser output so ISA and processing sample names share one key."""
    if not aliases or not data:
        return data
    remapped = {}
    for key, value in data.items():
        canon = canonical_sample_name(key, aliases)
        if canon in remapped and isinstance(remapped[canon], dict) and isinstance(value, dict):
            merged = dict(remapped[canon])
            for field, field_value in value.items():
                if field_value not in (None, '') or field not in merged or merged[field] in (None, ''):
                    merged[field] = field_value
            remapped[canon] = merged
        else:
            remapped[canon] = value
    return remapped


def generate_validation_report(fieldnames, populated_fields, mode, assay_suffix, paired_end, runsheet=None, filled_fields=None, validation_mismatches=None):
    """Generate a validation report showing which columns are missing data"""
    
    # Check if this is an ERCC dataset by looking at the runsheet
    is_ercc_dataset = False
    try:
        if runsheet and os.path.exists(runsheet):
            df = pd.read_csv(runsheet)
            if 'has_ERCC' in df.columns:
                # Check if any row has has_ERCC = 1 (assuming it's a boolean column)
                is_ercc_dataset = df['has_ERCC'].any() if df['has_ERCC'].dtype in ['int64', 'bool'] else False
    except Exception:
        # If we can't determine, assume it's not ERCC
        is_ercc_dataset = False
    
    # Define field categories
    metadata_fields = [
        'osd_num', 'sample', 'organism', 'tissue', 'sequencing_instrument', 
        'library_selection', 'library_layout', 'strandedness', 'read_depth', 
        'read_length', 'rrna_contamination', 'rin', 'organism_part', 'cell_line', 
        'cell_type', 'secondary_organism', 'strain', 'animal_source', 'seed_source', 
        'source_accession'
    ]
    
    # Only include 'mix' field if this is an ERCC dataset
    if is_ercc_dataset:
        metadata_fields.append('mix')
    
    gene_count_fields = ['gene_detected_gt10', 'gene_total', 'gene_detected_gt10_pct']
    
    fastqc_raw_fields = [f for f in fieldnames if f.startswith('raw_')]
    fastqc_trimmed_fields = [f for f in fieldnames if f.startswith('trimmed_')]
    
    # For single-end data, exclude _r fields from FastQC sections
    if not paired_end:
        fastqc_raw_fields = [f for f in fastqc_raw_fields if not f.endswith('_r')]
        fastqc_trimmed_fields = [f for f in fastqc_trimmed_fields if not f.endswith('_r')]
    
    star_fields = [
        'uniquely_mapped_percent', 'multimapped_percent', 'multimapped_toomany_percent',
        'unmapped_tooshort_percent', 'unmapped_other_percent'
    ]
    
    bowtie2_fields = [
        'total_reads', 'overall_alignment_rate', 'aligned_none', 'aligned_one', 'aligned_multi'
    ]
    
    featurecounts_fields = [
        'total_count', 'num_assigned', 'pct_assigned', 'num_unassigned_nofeatures',
        'num_unassigned_ambiguity', 'pct_unassigned_nofeatures', 'pct_unassigned_ambiguity'
    ]
    
    rseqc_fields = [
        'mean_genebody_cov_5_20', 'mean_genebody_cov_40_60', 'mean_genebody_cov_80_95',
        'ratio_genebody_cov_3_to_5', 'pct_sense', 'pct_antisense', 'pct_undetermined',
        'cds_exons_pct', '5_utr_exons_pct', '3_utr_exons_pct', 'introns_pct', 
        'tss_up_1kb_pct', 'tss_up_1kb_5kb_pct', 'tss_up_5kb_10kb_pct', 'tes_down_1kb_pct', 
        'tes_down_1kb_5kb_pct', 'tes_down_5kb_10kb_pct', 'other_intergenic_pct'
    ]
    
    # Add inner distance fields only for paired-end data
    if paired_end:
        rseqc_fields.extend(['peak_inner_dist', 'peak_inner_dist_pct_reads'])
    
    rsem_fields = [
        'num_uniquely_aligned', 'pct_uniquely_aligned', 'pct_multi_aligned',
        'pct_filtered', 'pct_unalignable'
    ]
    
    def get_missing_fields(field_list):
        return [f for f in field_list if f in fieldnames and f not in populated_fields]
    
    report_filename = f'qc_validation{assay_suffix}.txt'
    with open(report_filename, 'w') as f:
        # Header
        f.write("QC metrics validation report\n\n")
        mode_display = "Microbes" if mode == 'microbes' else "Default"
        data_type = "Paired end" if paired_end else "Single end"
        f.write(f"Mode: {mode_display}\n")
        f.write(f"Data type: {data_type}\n\n")
        f.write("-" * 40 + "\n\n")
        
        # Calculate summary stats
        missing_metadata_count = len([f for f in metadata_fields if f in fieldnames and f not in populated_fields])
        total_metadata_count = len([f for f in metadata_fields if f in fieldnames])
        
        # Count expected MultiQC fields based on mode (excluding metadata)
        expected_multiqc_fields = gene_count_fields + fastqc_raw_fields + fastqc_trimmed_fields + rseqc_fields
        if mode == 'microbes':
            expected_multiqc_fields += bowtie2_fields + featurecounts_fields
        else:
            expected_multiqc_fields += star_fields + rsem_fields
        
        missing_multiqc_count = len([f for f in expected_multiqc_fields if f in fieldnames and f not in populated_fields])
        total_multiqc_count = len([f for f in expected_multiqc_fields if f in fieldnames])
        
        # Summary section
        f.write("Summary:\n")
        f.write(f"* {missing_metadata_count}/{total_metadata_count} Metadata fields are empty\n")
        f.write(f"* {missing_multiqc_count}/{total_multiqc_count} expected MultiQC fields are empty\n\n")
        
        f.write("Categories missing entries:\n\n")
        
        # Metadata
        missing_metadata = get_missing_fields(metadata_fields)
        if missing_metadata:
            f.write("* Metadata:\n")
            for field in missing_metadata:
                f.write(f"** {field}\n")
            f.write("\n")
        
        # Gene count
        missing_gene_count = get_missing_fields(gene_count_fields)
        if missing_gene_count:
            f.write("* Gene Count:\n")
            for field in missing_gene_count:
                f.write(f"** {field}\n")
            f.write("\n")
        
        # FastQC Raw
        missing_fastqc_raw = get_missing_fields(fastqc_raw_fields)
        if missing_fastqc_raw:
            f.write("* FastQC Raw:\n")
            for field in missing_fastqc_raw:
                f.write(f"** {field}\n")
            f.write("\n")
        
        # FastQC Trimmed
        missing_fastqc_trimmed = get_missing_fields(fastqc_trimmed_fields)
        if missing_fastqc_trimmed:
            f.write("* FastQC Trimmed:\n")
            for field in missing_fastqc_trimmed:
                f.write(f"** {field}\n")
            f.write("\n")
        
        # Mode-specific sections
        if mode == 'microbes':
            # For microbes mode, check Bowtie2 and FeatureCounts (STAR/RSEM expected to be empty)
            missing_bowtie2 = get_missing_fields(bowtie2_fields)
            if missing_bowtie2:
                f.write("* Bowtie2 (Prokaryotes):\n")
                for field in missing_bowtie2:
                    f.write(f"** {field}\n")
                f.write("\n")
            
            missing_featurecounts = get_missing_fields(featurecounts_fields)
            if missing_featurecounts:
                f.write("* FeatureCounts (Prokaryotes):\n")
                for field in missing_featurecounts:
                    f.write(f"** {field}\n")
                f.write("\n")
        else:
            # For eukaryotic mode, check STAR and RSEM (Bowtie2/FeatureCounts expected to be empty)
            missing_star = get_missing_fields(star_fields)
            if missing_star:
                f.write("* STAR (Eukaryotes):\n")
                for field in missing_star:
                    f.write(f"** {field}\n")
                f.write("\n")
            
            missing_rsem = get_missing_fields(rsem_fields)
            if missing_rsem:
                f.write("* RSEM (Eukaryotes):\n")
                for field in missing_rsem:
                    f.write(f"** {field}\n")
                f.write("\n")
        
        # RSeQC (relevant for both modes)
        missing_rseqc = get_missing_fields(rseqc_fields)
        if missing_rseqc:
            f.write("* RSeQC:\n")
            for field in missing_rseqc:
                f.write(f"** {field}\n")
            f.write("\n")
        
        # Auto-filled fields section
        if filled_fields:
            f.write("Auto-filled fields (from MultiQC data):\n")
            for field in sorted(filled_fields):
                if field == 'read_depth':
                    f.write(f"** {field} (filled from raw_total_sequences_f)\n")
                elif field == 'read_length':
                    f.write(f"** {field} (filled from raw_max_sequence_length_f/r)\n")
                else:
                    f.write(f"** {field}\n")
            f.write("\n")
        
        # Validation mismatches section
        if validation_mismatches:
            f.write("Validation mismatches (assay table vs MultiQC data):\n")
            for sample, field, assay_value, multiqc_value in sorted(validation_mismatches):
                f.write(f"** Sample: {sample}, Field: {field}\n")
                f.write(f"   Assay table value: {assay_value}\n")
                f.write(f"   MultiQC value: {multiqc_value}\n")
            f.write("\n")


def main(osd_num, paired_end, assay_suffix, mode, runsheet=None):

    # Handle OSD number: if empty/None, skip metadata fetch and use empty string for CSV
    if not osd_num or osd_num.strip() == '':
        osd_num = None
        had_osd_prefix = False
    else:
        # Handle OSD number format: if it starts with "OSD-", extract the number part
        # Otherwise, use it as-is (for custom identifiers like "LEAF_Brapa")
        had_osd_prefix = osd_num.startswith('OSD-')
        if had_osd_prefix:
            osd_num = osd_num.split('-')[1]
        # If no OSD prefix, use the string as-is

    # Create the multiqc_data list with conditionally selected parsers based on mode
    multiqc_data = [
        parse_isa(),
        parse_fastqc('raw', assay_suffix),
        parse_fastqc('trimmed', assay_suffix)
    ]
    
    # Add the appropriate parsers based on mode
    if mode == 'microbes':
        multiqc_data.append(parse_bowtie2(assay_suffix))
        #multiqc_data.append(parse_featurecounts(assay_suffix))
        multiqc_data.append(parse_featurecounts(assay_suffix))
    else:
        multiqc_data.append(parse_star(assay_suffix))
        multiqc_data.append(parse_rsem(assay_suffix))
    
    # Add the remaining common parsers
    multiqc_data.extend([
        parse_genebody_cov(assay_suffix),
        parse_infer_exp(assay_suffix),
        parse_read_dist(assay_suffix),
        get_genecount(assay_suffix, mode)
    ])

    if paired_end:
        multiqc_data.append(parse_inner_dist(assay_suffix))

    sample_aliases = get_runsheet_aliases(runsheet)
    if sample_aliases:
        multiqc_data = [remap_parser_dict(data_source, sample_aliases) for data_source in multiqc_data]

    samples = {canonical_sample_name(s, sample_aliases) for ss in multiqc_data for s in ss if sample_key(s)}

    # Order samples according to runsheet if provided
    runsheet_order = get_runsheet_order(runsheet)
    if runsheet_order:
        ordered_samples = [s for s in runsheet_order if s in samples]
        extra_samples = sorted(samples - set(runsheet_order))
        samples = ordered_samples + extra_samples
    else:
        samples = sorted(samples, key=lambda s: (not s.isdigit(), int(s) if s.isdigit() else s))

    # Only fetch metadata if osd_num is provided
    metadata = get_metadata(osd_num) if osd_num else {}

    fieldnames = [
        'osd_num', 'sample', 'organism', 'tissue', 'sequencing_instrument', 'library_selection', 'library_layout', 'strandedness', 'read_depth', 'read_length', 'rrna_contamination', 'rin', 'organism_part', 'cell_line', 'cell_type', 'secondary_organism', 'strain', 'animal_source', 'seed_source', 'source_accession', 'mix',

        # Gene count
        'gene_detected_gt10', 'gene_total', 'gene_detected_gt10_pct',

        # FastQC (raw)
        'raw_total_sequences_f', 'raw_avg_sequence_length_f', 'raw_median_sequence_length_f', 'raw_max_sequence_length_f', 'raw_quality_score_mean_f', 'raw_quality_score_median_f', 'raw_percent_duplicates_f',
        'raw_percent_gc_f', 'raw_gc_min_1pct_f', 'raw_gc_max_1pct_f', 'raw_gc_auc_25pct_f', 'raw_gc_auc_50pct_f', 'raw_gc_auc_75pct_f', 'raw_n_content_sum_f',
        'raw_total_sequences_r', 'raw_avg_sequence_length_r', 'raw_median_sequence_length_r', 'raw_max_sequence_length_r', 'raw_quality_score_mean_r', 'raw_quality_score_median_r', 'raw_percent_duplicates_r',
        'raw_percent_gc_r', 'raw_gc_min_1pct_r', 'raw_gc_max_1pct_r', 'raw_gc_auc_25pct_r', 'raw_gc_auc_50pct_r', 'raw_gc_auc_75pct_r', 'raw_n_content_sum_r',

        # FastQC (trimmed)
        'trimmed_total_sequences_f', 'trimmed_avg_sequence_length_f', 'trimmed_median_sequence_length_f', 'trimmed_quality_score_mean_f', 'trimmed_quality_score_median_f', 'trimmed_percent_duplicates_f',
        'trimmed_percent_gc_f', 'trimmed_gc_min_1pct_f', 'trimmed_gc_max_1pct_f', 'trimmed_gc_auc_25pct_f', 'trimmed_gc_auc_50pct_f', 'trimmed_gc_auc_75pct_f', 'trimmed_n_content_sum_f',
        'trimmed_total_sequences_r', 'trimmed_avg_sequence_length_r', 'trimmed_median_sequence_length_r', 'trimmed_quality_score_mean_r', 'trimmed_quality_score_median_r', 'trimmed_percent_duplicates_r',
        'trimmed_percent_gc_r', 'trimmed_gc_min_1pct_r', 'trimmed_gc_max_1pct_r', 'trimmed_gc_auc_25pct_r', 'trimmed_gc_auc_50pct_r', 'trimmed_gc_auc_75pct_r', 'trimmed_n_content_sum_r',

        # STAR
        'uniquely_mapped_percent', 'multimapped_percent', 'multimapped_toomany_percent', 'unmapped_tooshort_percent', 'unmapped_other_percent',

        # Bowtie2
        'total_reads', 'overall_alignment_rate', 'aligned_none', 'aligned_one', 'aligned_multi',

        # FeatureCounts
        'total_count', 'num_assigned', 'pct_assigned', 'num_unassigned_nofeatures', 'num_unassigned_ambiguity', 'pct_unassigned_nofeatures', 'pct_unassigned_ambiguity',
        
        # RSeQC
        'mean_genebody_cov_5_20', 'mean_genebody_cov_40_60', 'mean_genebody_cov_80_95', 'ratio_genebody_cov_3_to_5',
        'pct_sense', 'pct_antisense', 'pct_undetermined',
        'peak_inner_dist', 'peak_inner_dist_pct_reads',
        'cds_exons_pct', '5_utr_exons_pct', '3_utr_exons_pct', 'introns_pct', 'tss_up_1kb_pct', 'tss_up_1kb_5kb_pct', 'tss_up_5kb_10kb_pct', 'tes_down_1kb_pct', 'tes_down_1kb_5kb_pct', 'tes_down_5kb_10kb_pct', 'other_intergenic_pct',

        # RSEM
        'num_uniquely_aligned', 'pct_uniquely_aligned', 'pct_multi_aligned', 'pct_filtered', 'pct_unalignable'
    ]

    # Make a set of fieldnames for fast lookup
    fieldnames_set = set(fieldnames)

    output_filename = f'qc_metrics{assay_suffix}.csv'
    
    # Track which fields have data for validation report
    populated_fields = set()
    # Track which fields were auto-filled
    filled_fields = set()
    # Track validation mismatches (assay table vs MultiQC)
    validation_mismatches = []
    
    with open(output_filename, mode='w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for sample in samples:
            # Collect all fields for this sample
            all_fields = {}
            for data_source in multiqc_data:
                source_data = None
                if sample in data_source:
                    source_data = data_source[sample]
                else:
                    for source_key, source_value in data_source.items():
                        if canonical_sample_name(source_key, sample_aliases) == sample:
                            source_data = source_value
                            break
                if source_data:
                    for k, v in source_data.items():
                        # Only keep fields that are in the fieldnames list
                        if k in fieldnames_set:
                            all_fields[k] = v
                            if v is not None and v != '':  # Track populated fields
                                populated_fields.add(k)
                        else:
                            # Optionally add debug output to see which fields are being skipped
                            # print(f"Skipping field not in fieldnames: {k}")
                            pass
            
            # Track fields that are always populated
            populated_fields.add('osd_num')
            populated_fields.add('sample')
            
            # Track populated metadata fields and normalize if needed
            for k, v in metadata.items():
                if v is not None and v != '':
                    populated_fields.add(k)
                    # Normalize strandedness and library_selection in metadata
                    if k == 'strandedness':
                        normalized_v = str(v).upper()
                        if v != normalized_v:
                            validation_mismatches.append((sample, 'strandedness', v, normalized_v))
                            metadata[k] = normalized_v
                    elif k == 'library_selection':
                        lib_sel_lower = str(v).lower()
                        normalized_v = v
                        if 'ribo' in lib_sel_lower:
                            normalized_v = 'ribo-depletion'
                        elif 'poly' in lib_sel_lower:
                            normalized_v = 'polyA enrichment'
                        if v != normalized_v:
                            metadata[k] = normalized_v
            
            # Validate and fill in missing read_depth and read_length from raw FastQC data
            # If read_depth exists, validate it matches raw_total_sequences_f (total sequences = read depth)
            if all_fields.get('read_depth') and all_fields.get('read_depth') != '':
                if all_fields.get('raw_total_sequences_f'):
                    try:
                        assay_read_depth = int(float(all_fields['read_depth']))
                        multiqc_read_depth = int(float(all_fields['raw_total_sequences_f']))
                        if assay_read_depth != multiqc_read_depth:
                            validation_mismatches.append((sample, 'read_depth', assay_read_depth, multiqc_read_depth))
                    except (ValueError, TypeError):
                        pass
            # If read_depth is empty, fill from raw_total_sequences_f
            elif all_fields.get('raw_total_sequences_f'):
                try:
                    all_fields['read_depth'] = int(float(all_fields['raw_total_sequences_f']))
                    populated_fields.add('read_depth')
                    filled_fields.add('read_depth')
                except (ValueError, TypeError):
                    pass
            
            multiqc_read_length = resolve_multiqc_read_length(all_fields)
            # If read_length exists, validate against FastQC max length
            if all_fields.get('read_length') and all_fields.get('read_length') != '':
                if multiqc_read_length is not None:
                    try:
                        assay_read_length = int(float(all_fields['read_length']))
                        if assay_read_length != multiqc_read_length:
                            validation_mismatches.append((sample, 'read_length', assay_read_length, multiqc_read_length))
                    except (ValueError, TypeError):
                        pass
            # If read_length is empty, fill from FastQC max length
            elif multiqc_read_length is not None:
                all_fields['read_length'] = multiqc_read_length
                populated_fields.add('read_length')
                filled_fields.add('read_length')
            
            # Normalize strandedness (convert to uppercase)
            if all_fields.get('strandedness') and all_fields.get('strandedness') != '':
                original_strandedness = all_fields['strandedness']
                normalized_strandedness = str(original_strandedness).upper()
                if original_strandedness != normalized_strandedness:
                    validation_mismatches.append((sample, 'strandedness', original_strandedness, normalized_strandedness))
                    all_fields['strandedness'] = normalized_strandedness
            
            # Normalize library_selection (only if contains 'ribo' or 'poly', otherwise pass through unchanged)
            if all_fields.get('library_selection') and all_fields.get('library_selection') != '':
                original_lib_sel = all_fields['library_selection']
                lib_sel_lower = str(original_lib_sel).lower()
                normalized_lib_sel = original_lib_sel  # Default: pass through unchanged
                
                if 'ribo' in lib_sel_lower:
                    normalized_lib_sel = 'ribo-depletion'
                elif 'poly' in lib_sel_lower:
                    normalized_lib_sel = 'polyA enrichment'
                # If no match, normalized_lib_sel == original_lib_sel, so no change is made
                
                if original_lib_sel != normalized_lib_sel:
                    all_fields['library_selection'] = normalized_lib_sel
            
            # Write rows with osd_num and sample fields
            # Add OSD- prefix only if original input had it, otherwise use as-is
            # If osd_num is None/empty, use empty string
            osd_display = 'OSD-' + osd_num if had_osd_prefix and osd_num else (osd_num if osd_num else '')
            writer.writerow({'osd_num': osd_display, 'sample': sample, **metadata, **all_fields})
    
    # Generate validation report
    try:
        generate_validation_report(fieldnames, populated_fields, mode, assay_suffix, paired_end, runsheet, filled_fields, validation_mismatches)
    except Exception as e:
        print(f"WARNING: Failed to generate validation report: {str(e)}")


def get_metadata(osd_num):
    data = {}
    if not osd_num:
        return data  # Return empty dict if no osd_num provided
    try:
        r = requests.get('https://osdr.nasa.gov/osdr/data/osd/meta/' + osd_num)
        r.raise_for_status()  # Raise an exception for bad status codes
        
        # Source accession
        comments = r.json()['study']['OSD-' + osd_num]['studies'][0]['comments']
        source_accession = [c for c in comments if c['name'] == 'Data Source Accession'][0]

        if source_accession and source_accession['value']:
            data['source_accession'] = source_accession['value']
    except (requests.RequestException, KeyError, IndexError, ValueError, Exception) as e:
        print(f"WARNING: Error fetching metadata from OSD API: {str(e)}")
        # Return empty dict if API call fails
    
    return data


def parse_isa():
    try:
        with ZipFile(glob('*ISA.zip')[0]) as z:
            with z.open([i for i in z.namelist() if i.startswith('s_')][0]) as f:
                s = list(csv.DictReader(TextIOWrapper(f, 'utf-8'), delimiter='\t'))

            with z.open([i for i in z.namelist() if i.startswith('a_') and 'rna-seq' in i.lower().replace('_', '-') and 'mirna' not in i.lower() and 'single-cell' not in i.lower()][0]) as f:
                a = list(csv.DictReader(TextIOWrapper(f, 'utf-8'), delimiter='\t'))

        data = {}

        sample_fields = {
            'characteristics[organism]': 'organism',
            'characteristics[material type]': 'tissue',
            'material type': 'tissue',
            'characteristics[organism part]': 'organism_part',
            'factor value[organism part]': 'organism_part',
            'characteristics[cell line]': 'cell_line',
            'characteristics[cell line,http://purl.obolibrary.org/obo/clo_0000031,obi]': 'cell_line',
            'characteristics[cell type]': 'cell_type',
            'characteristics [cell type]': 'cell_type',
            'characteristics[secondary organism]': 'secondary_organism',
            'characteristics[strain]': 'strain',
            'characteristics[animal source]': 'animal_source',
            'characteristics[seed source]': 'seed_source'
        }
        assay_fields = {
            'parameter value[sequencing instrument]': 'sequencing_instrument',
            'parameter value[library selection]': 'library_selection',
            'parameter value[library layout]': 'library_layout',
            'parameter value[strandedness]': 'strandedness',
            'parameter value[stranded]': 'strandedness',
            'parameter value[read depth]': 'read_depth',
            'parameter value[read length]': 'read_length',
            'parameter value[rrna contamination]': 'rrna_contamination',
            'parameter value[rrna estimation]': 'rrna_contamination',
            'parameter value[qa score]': 'rin',
            'parameter value[spike-in mix number]': 'mix'
        }

        for row in a:
            key = sample_key(row['Sample Name'])
            data[key] = {assay_fields[k.lower()]: v.strip() for k, v in row.items() if k.lower() in assay_fields}

        for row in s:
            key = sample_key(row['Sample Name'])
            if key not in data:  # samples could be in other assays
                continue

            sample_data = {sample_fields[k.lower()]: v.strip() for k, v in row.items() if k.lower() in sample_fields}
            data[key] = {**data[key], **sample_data}

        return data
    except (FileNotFoundError, KeyError, IndexError, ValueError, Exception) as e:
        print(f"WARNING: Error processing ISA data: {str(e)}")
        return {}


def parse_fastqc(prefix, assay_suffix):
    try:
        with open(f'{prefix}_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        # Find FastQC section by looking for FastQC-specific fields
        fastqc_section = None
        for section in j['report_general_stats_data']:
            sample_data = next(iter(section.values()), {})
            if 'total_sequences' in sample_data and 'percent_gc' in sample_data:
                fastqc_section = section
                break

        if not fastqc_section:
            return {}

        raw_fastqc = j.get('report_saved_raw_data', {}).get('multiqc_fastqc', {})

        # Group the samples by base name for paired end data
        sample_groups = {}
        for sample in fastqc_section.keys():
            sample_name_str = sample_key(sample)
            base_name = sample_key(mqc_base_sample_name(sample_name_str))
            end = mqc_read_end(sample_name_str)
            sample_groups.setdefault(base_name, {'f': None, 'r': None})[end] = sample

        data = {}
        # Process each sample group
        for base_name, reads in sample_groups.items():
            data[base_name] = {}
            
            # Process forward read
            if reads['f']:
                for k, v in fastqc_section[reads['f']].items(): 
                    if k != 'percent_fails':
                        data[base_name][prefix + '_' + k + '_f'] = v
                raw_key = lookup_raw_fastqc_sample(base_name, 'f', raw_fastqc, reads['f'])
                max_len = fastqc_max_sequence_length(raw_fastqc.get(raw_key, {}).get('Sequence length') if raw_key else None)
                if max_len is not None:
                    data[base_name][prefix + '_max_sequence_length_f'] = max_len
                        
            # Process reverse read
            if reads['r']:
                for k, v in fastqc_section[reads['r']].items(): 
                    if k != 'percent_fails':
                        data[base_name][prefix + '_' + k + '_r'] = v
                raw_key = lookup_raw_fastqc_sample(base_name, 'r', raw_fastqc, reads['r'])
                max_len = fastqc_max_sequence_length(raw_fastqc.get(raw_key, {}).get('Sequence length') if raw_key else None)
                if max_len is not None:
                    data[base_name][prefix + '_max_sequence_length_r'] = max_len

        # Process other stats sections (quality, GC, etc)
        for section, suffix in [
            ('fastqc_per_base_sequence_quality_plot', 'quality_score'),
            ('fastqc_per_sequence_gc_content_plot', 'gc'),
            ('fastqc_per_base_n_content_plot', 'n_content')
        ]:
            if section in j['report_plot_data']:
                for series_name, pairs in iter_multiqc_plot_series(j['report_plot_data'][section]):
                    sample = sample_key(series_name)
                    base_name = sample_key(mqc_base_sample_name(sample))
                    read_suffix = '_r' if mqc_read_end(sample) == 'r' else '_f'

                    # Skip if we don't have this sample
                    if base_name not in data:
                        continue
                        
                    # Process based on the section
                    if suffix == 'quality_score':
                        data[base_name][prefix + '_quality_score_mean' + read_suffix] = mean([i[1] for i in pairs])
                        data[base_name][prefix + '_quality_score_median' + read_suffix] = median([i[1] for i in pairs])
                    elif suffix == 'gc':
                        gc_data_1pct = [i[0] for i in pairs if i[1] >= 1]
                        if gc_data_1pct:
                            data[base_name][prefix + '_gc_min_1pct' + read_suffix] = gc_data_1pct[0]
                            data[base_name][prefix + '_gc_max_1pct' + read_suffix] = gc_data_1pct[-1]
                            
                            gc_data_cum = list(np.cumsum([i[1] for i in pairs]))
                            data[base_name][prefix + '_gc_auc_25pct' + read_suffix] = list(i >= 25 for i in gc_data_cum).index(True)
                            data[base_name][prefix + '_gc_auc_50pct' + read_suffix] = list(i >= 50 for i in gc_data_cum).index(True)
                            data[base_name][prefix + '_gc_auc_75pct' + read_suffix] = list(i >= 75 for i in gc_data_cum).index(True)
                    elif suffix == 'n_content':
                        data[base_name][prefix + '_n_content_sum' + read_suffix] = sum([i[1] for i in pairs])

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError, ValueError, TypeError) as e:
        print(f"WARNING: Could not process {prefix} FastQC data: {e}")
        return {}


def parse_star(assay_suffix):
    try:
        with open(f'align_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        align_fields = ['uniquely_mapped_percent', 'multimapped_percent', 'multimapped_toomany_percent', 'unmapped_tooshort_percent', 'unmapped_other_percent']

        for sample in j['report_general_stats_data'][0].keys():
            data[sample] = {k:v for k, v in j['report_general_stats_data'][0][sample].items() if k in align_fields}

        return data
    except (FileNotFoundError, KeyError, IndexError):
        return {}


def parse_bowtie2(assay_suffix):
    """
    Parse Bowtie2 alignment statistics from MultiQC report
    """
    try:
        # Open the data file
        with open(f'align_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            data = json.load(f)
        
        # Create a dictionary to hold sample data
        samples_data = {}
        
        # Extract the stats for each sample
        for section in data['report_general_stats_data']:
            for sample, stats in section.items():
                # Clean up sample name to remove read identifiers
                base_sample = re.sub(r'_R[12]$', '', sample)
                
                if base_sample not in samples_data:
                    samples_data[base_sample] = {}
                
                # Report these 5 key metrics regardless of whether data is paired-end or single-end
                if 'total_reads' in stats:
                    samples_data[base_sample]['total_reads'] = stats['total_reads']
                
                if 'overall_alignment_rate' in stats:
                    samples_data[base_sample]['overall_alignment_rate'] = stats['overall_alignment_rate']
                
                # Handle aligned none, one, and multi for both paired and unpaired data
                # For paired data
                if 'paired_aligned_none' in stats:
                    samples_data[base_sample]['aligned_none'] = stats['paired_aligned_none']
                    samples_data[base_sample]['aligned_one'] = stats['paired_aligned_one']
                    samples_data[base_sample]['aligned_multi'] = stats['paired_aligned_multi']
                
                # For unpaired/single-end data
                elif 'unpaired_aligned_none' in stats:
                    samples_data[base_sample]['aligned_none'] = stats['unpaired_aligned_none']
                    samples_data[base_sample]['aligned_one'] = stats['unpaired_aligned_one']
                    samples_data[base_sample]['aligned_multi'] = stats['unpaired_aligned_multi']
        
        return samples_data
    
    except FileNotFoundError:
        print(f"WARNING: Could not find Bowtie2 MultiQC data file")
        return {}
    except Exception as e:
        print(f"ERROR parsing Bowtie2 data: {str(e)}")
        return {}


def parse_genebody_cov(assay_suffix):
    try:
        with open(f'geneBody_cov_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for cov_data in j['report_plot_data']['rseqc_gene_body_coverage_plot']['datasets'][0]['lines']:
            sample = cov_data['name']

            data[sample] = {
                'mean_genebody_cov_5_20': mean([i[1] for i in cov_data['pairs'] if 5 <= i[0] <= 20]),
                'mean_genebody_cov_40_60': mean([i[1] for i in cov_data['pairs'] if 40 <= i[0] <= 60]),
                'mean_genebody_cov_80_95': mean([i[1] for i in cov_data['pairs'] if 80 <= i[0] <= 95])
            }

            # Handle division by zero
            if data[sample]['mean_genebody_cov_5_20'] == 0:
                data[sample]['ratio_genebody_cov_3_to_5'] = None
            else:
                data[sample]['ratio_genebody_cov_3_to_5'] = data[sample]['mean_genebody_cov_80_95'] / data[sample]['mean_genebody_cov_5_20']

        return data
    except (FileNotFoundError, KeyError, IndexError):
        return {}


def parse_infer_exp(assay_suffix):
    try:
        with open(f'infer_exp_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        key_dict = {
            'se_sense': 'pct_sense',
            'se_antisense': 'pct_antisense',
            'pe_sense': 'pct_sense',
            'pe_antisense': 'pct_antisense',
            'failed': 'pct_undetermined'
        }

        data = {s.replace('_infer_expt', ''):{key_dict[k]:v * 100 for k, v in d.items() if k in key_dict} for s, d in j['report_saved_raw_data']['multiqc_rseqc_infer_experiment'].items()}

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process infer experiment data")
        return {}


def parse_inner_dist(assay_suffix):
    try:
        with open(f'inner_dist_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for dist_data in j['report_plot_data']['rseqc_inner_distance_plot']['datasets'][1]['lines']:
            sample = dist_data['name']

            max_dist = sorted(dist_data['pairs'], key=lambda i: i[1], reverse=True)[0]
            data[sample] = {
                'peak_inner_dist': max_dist[0],
                'peak_inner_dist_pct_reads': max_dist[1]
            }

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process inner distance data")
        return {}


def parse_read_dist(assay_suffix):
    try:
        with open(f'read_dist_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        read_fields = ['cds_exons_tag_pct', '5_utr_exons_tag_pct', '3_utr_exons_tag_pct', 'introns_tag_pct', 'tss_up_1kb_tag_pct', 'tes_down_1kb_tag_pct', 'other_intergenic_tag_pct']

        for sample, read_data in j['report_saved_raw_data']['multiqc_rseqc_read_distribution'].items():
            data[sample] = {k:v for k, v in read_data.items() if k in read_fields}

            data[sample]['tss_up_1kb_5kb_pct'] = read_data['tss_up_5kb_tag_pct'] - read_data['tss_up_1kb_tag_pct']
            data[sample]['tss_up_5kb_10kb_pct'] = read_data['tss_up_10kb_tag_pct'] - read_data['tss_up_5kb_tag_pct']
            data[sample]['tes_down_1kb_5kb_pct'] = read_data['tes_down_5kb_tag_pct'] - read_data['tes_down_1kb_tag_pct']
            data[sample]['tes_down_5kb_10kb_pct'] = read_data['tes_down_10kb_tag_pct'] - read_data['tes_down_5kb_tag_pct']

        data = {s.replace('_read_dist', ''):{k.replace('_tag', ''):v for k, v in d.items()} for s, d in data.items()}

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process read distribution data")
        return {}


def parse_rsem(assay_suffix):
    try:
        with open(f'RSEM_count_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for sample, count_data in j['report_saved_raw_data']['multiqc_rsem'].items():
            total_reads = count_data['Unique'] + count_data['Multi'] + count_data['Filtered'] + count_data['Unalignable']

            data[sample] = {
                'num_uniquely_aligned': count_data['Unique'],
                'pct_uniquely_aligned': count_data['Unique'] / total_reads * 100,
                'pct_multi_aligned': count_data['Multi'] / total_reads * 100,
                'pct_filtered': count_data['Filtered'] / total_reads * 100,
                'pct_unalignable': count_data['Unalignable'] / total_reads * 100
            }

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError, ZeroDivisionError):
        print(f"WARNING: Could not process RSEM data")
        return {}


def parse_featurecounts(assay_suffix):
    try:
        with open(f'FeatureCounts_multiqc{assay_suffix}_data/multiqc_data.json') as f:
            j = json.loads(f.read())

        data = {}

        for sample, count_data in j['report_saved_raw_data']['multiqc_featurecounts'].items():
            data[sample] = {
                'total_count': count_data['Total'],
                'num_assigned': count_data['Assigned'],
                'pct_assigned': count_data['percent_assigned'],
                'num_unassigned_nofeatures': count_data['Unassigned_NoFeatures'],
                'num_unassigned_ambiguity': count_data['Unassigned_Ambiguity'],
                'pct_unassigned_nofeatures': count_data['Unassigned_NoFeatures'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0,
                'pct_unassigned_ambiguity': count_data['Unassigned_Ambiguity'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0
            }

        return data
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError):
        print(f"WARNING: Could not process FeatureCounts data")
        return {}


def get_genecount(assay_suffix, mode):
    try:
        # For microbes mode, prioritize FeatureCounts, otherwise prioritize RSEM
        if mode == 'microbes':
            if os.path.exists(f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'
            elif os.path.exists(f'RSEM_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'RSEM_Unnormalized_Counts{assay_suffix}.csv'
            else:
                print(f"WARNING: Could not find any count files for gene count statistics")
                return {}  # No counts file found
        else:
            # Default behavior - RSEM first, then FeatureCounts
            if os.path.exists(f'RSEM_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'RSEM_Unnormalized_Counts{assay_suffix}.csv'
            elif os.path.exists(f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'):
                count_file = f'FeatureCounts_Unnormalized_Counts{assay_suffix}.csv'
            else:
                print(f"WARNING: Could not find any count files for gene count statistics")
                return {}  # No counts file found
            
        df = pd.read_csv(count_file, index_col=0)
        df = df[~df.index.str.contains('^ERCC-')]

        data = (df > 10).sum().to_frame(name='gene_detected_gt10')
        data['gene_total'] = df.shape[0]
        data['gene_detected_gt10_pct'] = data.gene_detected_gt10 / data.gene_total * 100

        return data.to_dict(orient='index')
    except (FileNotFoundError, pd.errors.EmptyDataError, ValueError, Exception) as e:
        print(f"WARNING: Error processing gene count data: {str(e)}")
        return {}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--osd-num')
    parser.add_argument('--paired', action='store_true')
    parser.add_argument('--assay_suffix', default='')
    parser.add_argument('--mode', default='default')
    parser.add_argument('--runsheet', default=None)
    args = parser.parse_args()

    main(args.osd_num, args.paired, args.assay_suffix, args.mode, args.runsheet)
