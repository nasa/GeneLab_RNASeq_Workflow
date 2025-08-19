# GeneLab RNAseq Consensus Processing Workflow

> GeneLab, part of [NASA's Open Science Data Repository (OSDR)](https://www.nasa.gov/osdr), has wrapped each step of the RNASeq consensus processing pipeline ([RCP](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNASeq)), starting with version F, into a Nextflow workflow with validation and verification of output files built in after each step. This repository contains the Nextflow workflow code (NF_RCP) along with instructions for installation and usage starting with NF_RCP version 2.0.2. For previous versions, refer to the table in the GeneLab_Data_Processing repository [workflow documentation](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNASeq/Workflow_Documentation) which lists (and links to) each RCP version and the corresponding workflow code. Exact workflow run info and RCP version used to process specific datasets that have been released are available in the \*nextflow_processing_info.txt file on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/), which can be found under 'Files' -> 'GeneLab Processed RNA-Seq Files' -> 'Supplemental Materials'.

## General Workflow Information

### Implementation Tools

The current GeneLab RNAseq consensus processing pipelines (RCP) for eukaryotic organisms ([GL-DPPD-7101-G](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-G.md)) and prokaryotic organisms ([GL-DPPD-7115](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7115_Versions/GL-DPPD-7115.md)) are implemented as a single [Nextflow](https://nextflow.io/) DSL2 workflow that utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in containers. This workflow (NF_RCP) is run using the command line interface (CLI) of any unix-based system. While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow. See the [NF_RCP Workflow & Subworkflows](#nf_rcp-workflow--subworkflows) section below for more details on the NF_RCP workflow, including installation and execution information.

For datasets that include ERCC spike-ins, ERCC analysis is implemented as a separate, manual step using a [jupyter](https://jupyter.org) notebook running in a [conda](https://conda-forge.org/docs/) environment. See the [ERCC Analysis Workflow](#ercc-analysis-workflow) section below for more information on how to utilize the ERCC Analysis jupyter notebook.

<br>

# NF_RCP Workflow & Subworkflows

### NF_RCP Resource Requirements

The table below details the default maximum resource allocations for individual Nextflow processes.

| Mode      | Organism Type         | Default CPU Cores | Default Memory |
|-----------|-----------------------|-------------------|----------------|
| `default` | Eukaryotic organisms  | 16                | 72 GB          |
| `microbes`| Prokaryotic organisms | 8                 | 16 GB          |

> **Note:** These per-process resource allocations are defaults. They can be adjusted by modifying `cpus` and `memory`  directives in the configuration files: [`local.config`](workflow_code/conf/local.config) (local execution) and [`slurm.config`](workflow_code/conf/slurm.config) (SLURM clusters).

> **Click links below to show/hide workflow diagrams**

<details open>
<summary>NF_RCP workflow for GL-DPPD-7101-G (Eukaryotes)</summary>
<p align="center">
<a href="images/NF_RCP_euk_wf_diagram.png"><img src="images/NF_RCP_euk_wf_diagram.png"></a>
</p>
</details>

<details>
<summary>NF_RCP workflow for GL-DPPD-7115 (Prokaryotes)</summary>
<p align="center">
<a href="images/NF_RCP_prok_wf_diagram.png"><img src="images/NF_RCP_prok_wf_diagram.png"></a>
</p>
</details>

---
The NF_RCP workflow is composed of three subworkflows as shown in the workflow diagrams above.
Below is a description of each subworkflow and the additional output files generated that are not already indicated in the [GL-DPPD-7101-G](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-G.md) and [GL-DPPD-7115](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7115_Versions/GL-DPPD-7115.md) pipeline documents:

1. **Analysis Staging Subworkflow**

   - Description:
     - This subworkflow extracts the metadata parameters (e.g. organism, library layout) needed for processing from the OSD/GLDS ISA archive and retrieves the raw reads files hosted on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).
       > *OSD/GLDS ISA archive*: ISA directory containing Investigation, Study, and Assay (ISA) metadata files for a respective GLDS dataset - the *ISA.zip file is located under 'Files' -> 'Study Metadata Files' for any GeneLab Data Set (GLDS) in the [OSDR](https://osdr.nasa.gov/bio/repo/).

2. **RNAseq Consensus Pipeline Subworkflow**

   - Description:
     - This subworkflow uses the staged raw data and metadata parameters from the Analysis Staging Subworkflow to generate processed data using either:
       - [Version G of the GeneLab RCP](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-G.md) when the `--mode` parameter is omitted (default)
       - [The GeneLab Prokaryotic RCP](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7115_Versions/GL-DPPD-7115.md) when using `--mode microbes`
       
       The selection impacts the choice of aligner and read counter tools used in the workflow.

3. **V&V Pipeline Subworkflow**

   - Description:
     - This subworkflow performs validation and verification (V&V) on the raw and processed data files in real-time.  It performs a series of checks on the output files generated and flags the results, using the flag codes indicated in the table below, which are outputted as a set of log files. 
     
       **V&V Flags**:

       |Flag Codes|Flag Name|Interpretation|
       |:---------|:--------|:-------------|
       | 20    | GREEN   | Indicates the check passed all validation conditions |
       | 30    | YELLOW  | Indicates the check was flagged for minor issues (e.g. slight outliers) |
       | 50    | RED     | Indicates the check was flagged for moderate issues (e.g. major outliers) |
       | 80    | HALT    | Indicates the check was flagged for severe issues that trigger a processing halt (e.g. missing data) |

<br>

---
## Utilizing the Workflow

1. [Install Nextflow and Singularity](#1-install-nextflow-and-singularity)  
   1a. [Install Nextflow](#1a-install-nextflow)  
   1b. [Install Singularity](#1b-install-singularity)
2. [Download the Workflow Files](#2-download-the-workflow-files)  
3. [Fetch Singularity Images](#3-fetch-singularity-images)  
4. [Run the Workflow](#4-run-the-workflow)  
   4a. [Approach 1: Run the workflow on a GeneLab RNAseq dataset with automatic retrieval of reference fasta and gtf files](#4a-approach-1-run-the-workflow-on-a-genelab-rnaseq-dataset-with-automatic-retrieval-of-reference-fasta-and-gtf-files)  
   4b. [Approach 2: Run the workflow on a GeneLab RNAseq dataset with custom reference fasta and gtf files](#4b-approach-2-run-the-workflow-on-a-genelab-rnaseq-dataset-with-custom-reference-fasta-and-gtf-files)  
   4c. [Approach 3: Run the workflow on a non-GeneLab dataset using a user-created runsheet with automatic retrieval of reference fasta and gtf files](#4c-approach-3-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-automatic-retrieval-of-reference-fasta-and-gtf-files)  
   4d. [Approach 4: Run the workflow on a non-GeneLab dataset using a user-created runsheet with custom reference fasta and gtf files](#4d-approach-4-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-custom-reference-fasta-and-gtf-files)
5. [Additional Output Files](#5-additional-output-files)  

<br>

---

### 1. Install Nextflow and Singularity 

#### 1a. Install Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install Anaconda, we recommend installing a Miniforge version appropriate for your system, as documented on the [conda-forge website](https://conda-forge.org/download/), where you can find basic binaries for most systems. More detailed miniforge documentation is available in the [miniforge github repository](https://github.com/conda-forge/miniforge).
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab RCP workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. See the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_RCP GeneLab workflow for processing RNAseq data are in the [workflow_code](workflow_code) directory. To get a 
copy of latest NF_RCP version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_RNASeq_Workflow/releases/download/NF_RCP_2.0.2/NF_RCP_2.0.2.zip

unzip NF_RCP_2.0.2.zip
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_RCP workflow:
> Note: This command should be run in the location containing the `NF_RCP_2.0.2` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes. Approximately 8GB of RAM is needed to download and build the Singularity images.

```bash
bash NF_RCP_2.0.2/bin/prepull_singularity.sh NF_RCP_2.0.2/config/by_docker_image.config
```


Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow

While in the location containing the `NF_RCP_2.0.2` directory that was downloaded in [step 2](#2-download-the-workflow-files), you are now able to run the workflow.

Both workflows automatically load reference files and organism-specific gene annotation files from the [GeneLab annotations table](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv). For organisms not listed in the table or to use alternative reference files, additional workflow parameters can be specified.

 Below are four examples of how to run the NF_RCP workflow:
> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --reference_version) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

> Note: To use Docker instead of Singularity, use `-profile docker` in the Nextflow run command. Nextflow will automatically pull images as needed.

> Note: The `-resume` parameter can be used to resume a previously interrupted workflow from where it left off (see [Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html#modify-and-resume)) or to restart the workflow from a specific point by changing relevant parameters, which will re-execute that process and all downstream affected processes.

<br>

#### 4a. Approach 1: Run the workflow on a GeneLab RNAseq dataset with automatic retrieval of reference fasta and gtf files

```bash
nextflow run NF_RCP_2.0.2/main.nf \ 
   -profile singularity,local \
   --accession OSD-194 
```

> Note: For prokaryotic RNAseq datasets, add the parameter `--mode microbes` to run the workflow using the prokaryotic pipeline ([GL-DPPD-7115](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7115_Versions/GL-DPPD-7115.md)). The default value of this parameter is `default`, which will use the eukaryotic pipeline ([GL-DPPD-7101-G](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-G.md)).

<br>

#### 4b. Approach 2: Run the workflow on a GeneLab RNAseq dataset with custom reference fasta and gtf files

```bash
nextflow run NF_RCP_2.0.2/main.nf \ 
   -profile singularity,local \
   --accession OSD-194 \
   --reference_version 112 \
   --reference_source ensembl \ 
   --reference_fasta <url/or/path/to/fasta> \ 
   --reference_gtf <url/or/path/to/gtf>
```

> Note: The `--reference_source` and `--reference_version` parameters should match the reference source and version number of the reference fasta and gtf files used. 

> Note: For gene annotations in the differential expression output table, see the optional `--gene_annotations_file` parameter described in the [Optional Parameters](#optional-parameters) section.

<br>

#### 4c. Approach 3: Run the workflow on a non-GeneLab dataset using a user-created runsheet with automatic retrieval of reference fasta and gtf files

```bash
nextflow run NF_RCP_2.0.2/main.nf \ 
   -profile singularity,local \
   --runsheet_path </path/to/runsheet> 
```

> Note: Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

<br>

#### 4d. Approach 4: Run the workflow on a non-GeneLab dataset using a user-created runsheet with custom reference fasta and gtf files

```bash
nextflow run NF_RCP_2.0.2/main.nf \ 
   -profile singularity \
   --accession OSD-194 \
   --reference_version 112 \
   --reference_source ensembl \ 
   --reference_fasta <url/or/path/to/fasta> \ 
   --reference_gtf <url/or/path/to/gtf> 
```
> Note: This approach should be used for organisms not listed in the [GeneLab annotations table](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv). 

> Note: The `--reference_source` and `--reference_version` parameters should match the reference source and version number of the reference fasta and gtf files used. 

> Note: For gene annotations in the differential expression output table, see the optional `--gene_annotations_file` parameter described in the [Optional Parameters](#optional-parameters) section.

<br>

#### Required Parameters For All Approaches:

* `NF_RCP_2.0.2/main.nf` - Instructs Nextflow to run the NF_RCP workflow 

* `-profile` - Specifies the configuration profile(s) to load, `singularity` instructs Nextflow to setup and use singularity for all software called in the workflow; use `local` for local execution ([local.config](workflow_code/conf/local.config)) or `slurm` for SLURM cluster execution ([slurm.config](workflow_code/conf/slurm.config))
  > Note: The output directory will be named `GLDS-#` when using a OSD or GLDS accession as input, or `results` when running the workflow with only a runsheet as input.

<br>

**Additional Required Parameters For [Approach 1](#4a-approach-1-run-the-workflow-on-a-genelab-rnaseq-dataset-with-automatic-retrieval-of-reference-fasta-and-gtf-files):**

* `--accession` - The OSD or GLDS ID for the dataset to be processed, eg. `GLDS-194` or `OSD-194`

<br>

**Additional Required Parameters For [Approach 2](#4b-approach-2-run-the-workflow-on-a-genelab-rnaseq-dataset-with-custom-reference-fasta-and-gtf-files):**

* `--accession` - The OSD or GLDS ID for the dataset to be processed, eg. `GLDS-194` or `OSD-194`

* `--reference_version` - specifies the reference source version to use for the reference genome (Ensembl release `112` is used in this example); only needed when using Ensembl as the reference source

* `--reference_source` - specifies the source of the reference files used (the source indicated in the Approach 2 example is `ensembl`) 

* `--reference_fasta` - specifies the URL or path to a fasta file 

* `--reference_gtf` - specifies the URL or path to a gtf file

<br>

**Additional Required Parameters For [Approach 3](#4c-approach-3-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-automatic-retrieval-of-reference-fasta-and-gtf-files):**

* `--runsheet_path` - specifies the path to a local runsheet; if not provided, a runsheet is automatically generated using OSDR metadata (type: string, default: null)

<br>

**Additional Required Parameters For [Approach 4](#4d-approach-4-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-custom-reference-fasta-and-gtf-files):**

* `--runsheet_path` - specifies the path to a local runsheet; if not provided, a runsheet is automatically generated using OSDR metadata (type: string, default: null)

* `--reference_version` - specifies the reference source version to use for the reference genome (Ensembl release `112` is used in this example); only needed when using Ensembl as the reference source

* `--reference_source` - specifies the source of the reference files used (the source indicated in the Approach 2 example is `ensembl`) 

* `--reference_fasta` - specifies the URL or path to a fasta file 

* `--reference_gtf` - specifies the URL or path to a gtf file

<br>

#### Optional Parameters:

* `--gene_annotations_file` - Specifies the URL or path to a gene annotation file that adds additional gene annotation columns to the differential expression output table. This can be:

  - The file listed in the `genelab_annots_link` column of the [GeneLab annotations table](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)
  - A custom gene annotation file where:
    - For organisms listed in the GeneLab annotations table: gene IDs must be in a column with the same name as column 1 of the GeneLab organism-specific gene annotation file
    - For organisms not listed in the table: gene IDs must be in a column named `gene_id`
  
  Only genes included in the specified annotations file will receive additional annotations in the output.

* `--outdir` - specifies the base directory where the output directory will be created (type: string, default: ".")  

* `--force_single_end` - forces the analysis to use single end processing; for paired end datasets, this means only R1 is used; for single end datasets, this should have no effect (type: boolean, default: false)  

* `--reference_store_path` - specifies the directory to store the reference fasta and gtf files (type: string, default: "./References")  

* `--derived_store_path` - specifies the directory to store the tool-specific indices created during processing (type: string, default: "./DerivedReferences")

* `--mode` - specifies which pipeline to use: set to `default` to run GL-DPPD-7101-G pipeline or set to `microbes` for the GL-DPPD-7115 prokaryotic pipeline (type: string, default: "default")
  > Note: This allows the workflow to process either eukaryotic (default) or prokaryotic RNAseq data using the appropriate pipeline.

* **DGE Filtering Parameters** - Options for filtering genes prior to DGE analysis:

  * `--dge_filter_method` - Method for filtering genes based on raw counts prior to differential expression analysis (type: string, default: "sum_threshold"):
    - `sum_threshold`: Remove genes with total counts ≤ threshold 
    - `sample_percent`: Remove genes if >X% of total counts are concentrated in top Y highest-count samples
    - `min_samples`: Remove genes not expressed in at least X samples
    - `count_per_sample`: Remove genes with total counts less than X * sample_count
  
  * `--dge_filter_sum_threshold` - Threshold value for total raw counts (type: number, default: 10)
  
  * `--dge_filter_sample_percent_threshold` - Percentage threshold for concentration in top samples (type: integer, default: 90)
  
  * `--dge_filter_sample_percent_max_samples` - Number of top samples to check for concentration (type: integer, default: 1)
  
  * `--dge_filter_min_samples_threshold` - Minimum number of samples where gene must be expressed (type: integer, default: 1)
  
  * `--dge_filter_count_per_sample_threshold` - Multiplier for sample-scaled count threshold (type: number, default: 1)

<br>

**Additional Optional Parameters:**

All parameters listed above and additional optional arguments for the RCP workflow, including debug related options that may not be immediately useful for most users, can be viewed by running the following command:

```bash
nextflow run NF_RCP_2.0.2/main.nf --help
```

See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details common to all nextflow workflows.

<br>

---

### 5. Additional Output Files

The outputs from the Analysis Staging and V&V Pipeline Subworkflows are described below:
> Note: The outputs from the RNAseq Consensus Pipeline Subworkflow are documented in the [GL-DPPD-7101-G](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-G.md) processing protocol.

**Analysis Staging Subworkflow**

   - Output:
     - Metadata/\*_bulkRNASeq_v1_runsheet.csv (table containing metadata required for processing, including the raw reads files location)
     - Metadata/\*-ISA.zip (the ISA archive of the OSD datasets to be processed, downloaded from the OSDR)
   
   
**V&V Pipeline Subworkflow**

   - Output:
     - VV_Logs/VV_log_final_GLbulkRNAseq.csv (table containing V&V flags for all checks performed)
     - VV_Logs/VV_log_final_only_issues_GLbulkRNAseq.csv (table containing V&V flags ONLY for checks that produced a flag code >= 30)
     - VV_Logs/VV_log_VV_RAW_READS_GLbulkRNAseq.csv (table containing V&V flags ONLY for raw reads checks)
     - VV_Logs/VV_log_VV_TRIMMED_READS_GLbulkRNAseq.csv (table containing V&V flags for trimmed reads checks ONLY)
     - VV_Logs/VV_log_VV_ALIGNMENT_GLbulkRNAseq.csv (table containing V&V flags for alignment file checks ONLY)
     - VV_Logs/VV_log_VV_RSEQC_GLbulkRNAseq.csv (table containing V&V flags for RSeQC file checks ONLY)
     - VV_Logs/VV_log_VV_COUNTS_GLbulkRNAseq.csv (table containing V&V flags for gene quantification file checks ONLY) 
     - VV_Logs/VV_log_VV_DESEQ2_ANALYSIS_GLbulkRNAseq.csv (table containing V&V flags for DESeq2 Analysis output checks ONLY)

**Processing Information Archive**

   - Output:
     - GeneLab/processing_info_GLbulkRNAseq.zip (Archive containing workflow execution metadata)
       - processing_info/samples.txt (single column list of all sample names in the dataset)
       - processing_info/nextflow_log_GLbulkRNAseq.txt (Nextflow execution logs captured via `nextflow log`)
       - processing_info/nextflow_run_command_GLbulkRNAseq.txt (Exact command line used to initiate the workflow)

**QC metrics summary**

  - Output:
    - GeneLab/qc_metrics_GLbulkRNAseq.csv (comma-separated text file containing a summary of qc metrics and metadata for the dataset, see the [QC metrics README](./QC_metrics_README.md) for a complete list of field definitions)
<br>

Standard Nextflow resource usage logs are also produced as follows:
> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

**Nextflow Resource Usage Logs**

   - Output:
     - nextflow_info/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
     - nextflow_info/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
     - nextflow_info/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)
     - nextflow_info/pipeline_dag_{timestamp}.html (a visualization of the workflow process DAG)

<br>

---

<br>

# ERCC Analysis Workflow

## Running ERCC Analysis

1. [Install conda and the ERCC analysis environment](#1-install-conda-and-the-ercc-analysis-environment)
2. [Launch Jupyter environment](#2-launch-jupyter-environment)  
3. [Execute the ERCC Analysis](#3-execute-the-ercc-analysis)
<br>

---

### 1. Install conda and the ERCC analysis environment

We recommend installing a Miniforge version appropriate for your system, as documented on the [conda-forge website](https://conda-forge.org/download/), where you can find basic binaries for most systems. More detailed miniforge documentation is available in the [miniforge github repository](https://github.com/conda-forge/miniforge).

Once conda is installed on your system, create the ercc_analysis conda environment:
```bash
conda env create -f NF_RCP_2.0.2/envs/ercc_analysis.yml
```

### 2. Launch Jupyter environment

Create an  working folder for the ERCC Analysis in the main analysis directory created in [Step 4](#4-run-the-workflow) above and copy the jupyter notebook file into it (downloaded in [Step 2](#2-download-the-workflow-files)) and start a local jupyter notebook server.

```bash
# create working directory
mkdir ERCC_Analysis

# copy jupyter notebook from workflow_code folder to new working directory
cp NF_RCP_2.0.2/bin/combined_ercc_analysis.ipynb ERCC_Analysis

# activate the ercc_analysis environment
conda activate ercc_analysis

# start jupyter notebook:
jupyter notebook
```
The `jupyter notebook` command will open a page in the default web browser showing a file listing of the folder it was launched from. If the web page does not open automatically, navigate to [http://localhost:8888/tree](http://localhost:8888/tree).

### 3. Execute the ERCC Analysis

Open the combined_ercc_analysis.ipynb notebook found in the ERCC_Analysis folder created above and follow the instructions for running in the analysis in the notebook.

Once the last cell has been executed, export the notebook contents to HTML using the "Save and Export Notebook As" option in the "File" menu.

<br>

---

# Licenses

The software for the RNAseq pipeline and workflow is released under the [NASA Open Source Agreement (NOSA) Version 1.3](License/RNA_Sequencing_NOSA_License.pdf).


### 3rd Party Software Licenses

Licenses for the 3rd party open source software utilized in the RNAseq pipeline and workflow can be found in the [3rd_Party_Licenses sub-directory](License/3rd_Party_Licenses). 

<br>

---

## Notices

Copyright © 2025 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved. 

### Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT. 

The "GeneLab RNA Sequencing Processing Pipeline and Workflow" software also makes use of 3rd party Open Source software, released under the licenses indicated above.  A complete listing of 3rd Party software notices and licenses made use of in "GeneLab RNA Sequencing Processing Pipeline and Workflow" can be found in the [3rd Party Licenses README.md](License/3rd_Party_Licenses/README.md) file. 

<br>

---
**Developed by:**  
Amanda Saravia-Butler    
Jonathan Oribello  

**Maintained by:**  
Alexis Torres  

**Contributors:**
Crystal Han   
