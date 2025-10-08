# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing bulk RNA sequence datasets through GeneLab's RNAseq consensus processing pipeline (RCP).


## Examples

1. [Runsheet for GLDS-48](single_end_runsheet/GLDS-48_bulkRNASeq_v1_runsheet.csv) (single end dataset)
2. [Runsheet for GLDS-194](paired_end_runsheet/GLDS-194_bulkRNASeq_v1_runsheet.csv) (paired end dataset)
3. [Runsheet for GLDS-605](paired_end_runsheet/GLDS-605_bulkRNASeq_v2_runsheet.csv) (paired end dataset with technical replicates)



## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |
| has_ERCC | bool | Set to True if ERCC spike-ins are included in the samples. This ensures ERCC normalized DGE is performed in addition to standard DGE. | True |
| paired_end | bool | Set to True if the samples were sequenced as paired-end. If set to False, samples are assumed to be single-end. | False |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| read1_path | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the forward reads fastq.gz file. | /my/data/sample_1.fastq.gz |
| read2_path | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the reverse reads fastq.gz file. For single-end data, this column should be omitted. | /my/data/sample_2.fastq.gz |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Space Flight |
| Original Sample Name | string | Used to map the sample name that will be used for processing to the original sample name. This is often identical except in cases where the original name includes spaces or weird characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |

## Entry Point Columns (for step-based workflows)

For entry points other than the default `raw_reads`, include the relevant columns below to specify the input files:

| Column Name | Entry Point | Type | Description | Example |
|:------------|:------------|:-----|:------------|:--------|
| trimmed_read1_path | trimmed_reads | string | Path to trimmed forward reads file | /my/data/sample_1_trimmed.fastq.gz |
| trimmed_read2_path | trimmed_reads | string | Path to trimmed reverse reads file (paired-end only) | /my/data/sample_2_trimmed.fastq.gz |
| bam_path | bam_files | string | Path to aligned BAM file | /my/data/sample_aligned.bam |
| genes_results_path | genes_results | string | Path to RSEM .genes.results file (`--mode default` only) | /my/data/sample.genes.results |
| counts_table_path | counts_table | string | Path to raw counts table file | /my/data/RSEM_Unnormalized_Counts.csv|
| dge_table_path | dge_table | string | Path to DGE table  | /my/data/differential_expression.csv |

## Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Source Name | string | Identifier linking samples. Used for handling technical replicates during differential gene expression analysis. Multiple samples with the same Source Name may be collapsed during analysis depending on the Has Tech Reps setting. | RR3_BSL_B7 |
| Has Tech Reps | bool | Indicates whether this sample is a technical replicate that should be collapsed with other samples sharing the same Source Name. Set to True for technical replicates that should be collapsed, False for distinct samples that should remain separate even if they share a Source Name. | False |
