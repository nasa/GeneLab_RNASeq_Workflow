#!/usr/bin/env python
""" This script takes as input a text file with software versions reported in format:
"TaskName":
    "SoftwareName": version

"TaskName2":
    "SoftwareName2": version2

and outputs a markdown table with the format:

|Program|Version|Relevant Links|
|---|---|---|
"""
from __future__ import annotations
import sys
import re
from pathlib import Path
import json
import yaml
import pandas as pd
from packaging import version

CONFIG = {
    "rnaseq": [
        ["NF_RCP", "https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq"],
        ["Nextflow", "https://github.com/nextflow-io/nextflow"],
        ["dp_tools", "https://github.com/J-81/dp_tools"],
        ["FastQC", "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"],
        ["MultiQC", "https://multiqc.info/"],
        ["Cutadapt", "https://cutadapt.readthedocs.io/en/stable/"],
        ["TrimGalore!", "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"],
        ["STAR", "https://github.com/alexdobin/STAR"],
        ["RSEM", "https://github.com/deweylab/RSEM"],
        ["Bowtie2", "https://github.com/BenLangmead/bowtie2"],
        ["Subread", "https://subread.sourceforge.net/"],
        ["SAMtools", "http://www.htslib.org/"],
        ["gtfToGenePred", "http://hgdownload.cse.ucsc.edu/admin/exe/"],
        ["genePredToBed", "http://hgdownload.cse.ucsc.edu/admin/exe/"],
        ["infer_experiment.py", "http://rseqc.sourceforge.net/#infer-experiment-py"],
        ["geneBody_coverage.py", "http://rseqc.sourceforge.net/#genebody-coverage-py"],
        ["inner_distance.py", "http://rseqc.sourceforge.net/#inner-distance-py"],
        ["read_distribution.py", "http://rseqc.sourceforge.net/#read-distribution-py"],
        ["R", "https://www.r-project.org/"],
        ["Bioconductor", "https://bioconductor.org"],
        ["BiocParallel", "https://bioconductor.org/packages/release/bioc/html/BiocParallel.html"],
        ["DESeq2", "https://bioconductor.org/packages/release/bioc/html/DESeq2.html"],
        ["tximport", "https://github.com/mikelove/tximport"],
        ["tidyverse", "https://www.tidyverse.org"],
        ["dplyr", "https://dplyr.tidyverse.org/"],
        ["knitr", "https://yihui.org/knitr/"],
        ["stringr", "https://github.com/tidyverse/stringr"],
        ["yaml", "https://github.com/yaml/yaml"],
        ["pandas", "https://github.com/pandas-dev/pandas"],
        ["seaborn", "https://seaborn.pydata.org/"],
        ["matplotlib", "https://matplotlib.org/stable"],
        ["numpy", "https://numpy.org/"],
        ["scipy", "https://scipy.org/"]
    ]
}

class NumericAsStringSafeLoader(yaml.SafeLoader):
    """YAML loader that does NOT coerce ints/floats; numeric-looking scalars stay as strings."""


# Remove implicit resolvers for ints and floats so they get loaded as strings
for ch, resolvers in list(NumericAsStringSafeLoader.yaml_implicit_resolvers.items()):
    NumericAsStringSafeLoader.yaml_implicit_resolvers[ch] = [
        (tag, regexp)
        for tag, regexp in resolvers
        if tag not in ("tag:yaml.org,2002:int", "tag:yaml.org,2002:float")
    ]


class QuotedStringDumper(yaml.SafeDumper):
    """YAML dumper that always quotes strings to avoid downstream numeric coercion."""


def _represent_str_always_quoted(dumper, data: str):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')


QuotedStringDumper.add_representer(str, _represent_str_always_quoted)

def normalize_name(name: str, known_names: list) -> str:
    """Match software name against known names, ignoring case and special chars"""
    name_clean = re.sub(r'[^a-zA-Z0-9]', '', name.lower())
    for known, _ in known_names:
        if re.sub(r'[^a-zA-Z0-9]', '', known.lower()) == name_clean:
            return known
    print(f"Warning: Unknown software detected: {name}")
    return name

def compare_versions(v1, v2):
    """Safely compare version strings"""
    try:
        return version.parse(str(v1)) > version.parse(str(v2))
    except version.InvalidVersion:
        # Fallback for non-standard version strings
        return str(v1) > str(v2)


def prefer_more_precise_representation(new_version: str, existing_version: str) -> bool:
    """When two versions are equal numerically, prefer the one with more precision (e.g., 1.30 vs 1.3).

    Returns True if new_version should replace existing_version based on precision of representation.
    """
    try:
        if version.parse(new_version) == version.parse(existing_version):
            # Compare count of digits to keep representations like 3.20 instead of 3.2
            new_digits = len(re.sub(r"\D", "", new_version))
            existing_digits = len(re.sub(r"\D", "", existing_version))
            if new_digits != existing_digits:
                return new_digits > existing_digits
            # If digit-count ties, prefer longer string (e.g., 1.0.0 vs 1.0)
            return len(new_version) > len(existing_version)
    except version.InvalidVersion:
        # If parsing fails for either, fall back to string length heuristic
        if new_version == existing_version:
            return False
        return len(new_version) > len(existing_version)
    return False

def main(versions_json_path: Path, output_path: Path, assay: str = 'rnaseq', workflow: str = None, workflow_version: str = None):
    software_urls = {name: url for name, url in CONFIG[assay]}
    known_names = CONFIG[assay]
    
    processed_versions = {}
    
    # Add workflow version if provided
    if workflow and workflow_version:
        processed_versions[workflow] = workflow_version
    
    with versions_json_path.open() as f:
        data = yaml.load(f, Loader=NumericAsStringSafeLoader)
        # Flatten nested structure
        for task_info in data.values():
            if isinstance(task_info, dict):
                for software, ver in task_info.items():
                    normalized_name = normalize_name(software, known_names)
                    # Strip quotes and convert to string
                    ver_str = str(ver).strip("'\"")
                    if normalized_name not in processed_versions:
                        processed_versions[normalized_name] = ver_str
                    else:
                        existing = processed_versions[normalized_name]
                        if compare_versions(ver_str, existing) or prefer_more_precise_representation(ver_str, existing):
                            processed_versions[normalized_name] = ver_str
    
    # Force gtfToGenePred and genePredToBed to version 469 IF they are present
    if "gtfToGenePred" in processed_versions:
        processed_versions["gtfToGenePred"] = "469"
    if "genePredToBed" in processed_versions:
        processed_versions["genePredToBed"] = "469"
    

    results = []
    for program, ver in processed_versions.items():
        results.append({
            "Program": program,
            "Version": ver,  # Already cleaned string
            "Relevant Links": software_urls.get(program, "")
        })

    df = pd.DataFrame(results)
    if not df.empty:
        # Create YAML dict before setting index (keep versions as strings)
        versions_dict = {row["Program"]: row["Version"] for _, row in df.iterrows()}
        
        # Now set index for markdown output
        df = df.set_index(keys="Program")
        
        # Split known and unknown software
        config_order = [name for name, _ in CONFIG[assay]]
        known_software = [x for x in config_order if x in df.index]
        unknown_software = [x for x in df.index if x not in config_order]
        unknown_software.sort()
        
        df = df.reindex(known_software + unknown_software)

        # Write outputs
        df.to_markdown(output_path, index=True)
        
        yaml_output = output_path.with_suffix('.yaml')
        with yaml_output.open('w') as f:
            yaml.dump(versions_dict, f, sort_keys=False, default_flow_style=False, Dumper=QuotedStringDumper)
        
        print(f"Wrote {output_path} and {yaml_output}")
    else:
        print("No software versions found to process")

if __name__ == "__main__":
    import click

    @click.command()
    @click.argument('input', type=click.Path(exists=True))
    @click.argument('output', type=click.Path())
    @click.option('--assay', type=click.Choice(['rnaseq']), default='rnaseq')
    @click.option('--workflow', type=str, help='Workflow name')
    @click.option('--workflow_version', type=str, help='Workflow version')
    def cli(input, output, assay, workflow, workflow_version):
        main(Path(input), Path(output), assay, workflow, workflow_version)

    cli()