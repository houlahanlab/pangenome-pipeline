# Pangenome:SomaticSNV pipeline

## Description
Align paired-end FASTQ files to a pangenome and perform somatic variant detection from tumour/normal paired samples.

## Input files
### Pipeline configuration
The pangenome-somatic pipeline requires a configuration file in YAML format. This YAML file provides paths to all input files (graph- and linear- reference files, target-regions, etc.), tool executables (to enable cross-system compatibility) and per-stage resource requirements.

See [example_pipeline_config.yaml](../examples/example_pipeline_config.yaml) . This file contains example filepaths and resource requirements, determined from testing a run with 40 WXS BAMs from TCGA.

### Reference Files
Download your desired pangenome reference [here](https://data.humanpangenome.org/alignments), including the .gbz, .dist and .min files.

Extract the linear-reference genome as included in the above pangenome:
<pre><code># note that this will result in non-traditional contig ordering
vg paths paths --extract-fasta -x graph_reference.gbz --paths-by GRCh38 > grch38_linear.fa
samtools dict -o grch38_linear.dict grch38_linear.fa
sed 's/GRCh38#0#chr/chr/g' grch38_linear.fa > grch38_linear_fixed.fa
samtools faidx grch38_linear_fixed.fa
</code></pre>

The 'ref_gbz', 'ref_dist' and 'ref_min' keys in your tool_config.yaml refer to the graph-based reference files (ie, graph_reference.gbz, graph_reference.dist and graph_reference.min). These are used for the VG alignment.

The 'ref_paths' key in your tool_config.yaml is used during VG SURJECT to select/output sequences mapped to the contigs of interest. This could be a .dict file (ie, grch38_linear.dict), or a separate text file of contig IDs for chromosomes of interest (ie, GRCh38#0#chr1-22, GRCh38#0#chrX, GRCh38#0#chrY) with 1 contig per line. The order of contigs in your final BAM will be the same as the order found in this file.

Note the contig stem (ie, GRCh38#0#) as this can be supplied via the tool_config.yaml (ref_pattern: 'GRCh38#0#'); doing so will remove it from the surjected BAM files allowing for compatibility with standard reference files used downstream (ie, GATK's bundle files / gnomad / 1000g / etc.). If not used, you must ensure all downstream reference-data files have contig names consistent with your reference fasta (see next).

The 'reference' key in your tool_config.yaml refers to a linear reference with standard contig names (ie, chr1-22, X, Y), such as 'grch38_linear_fixed.fa' in the above example. This will be used for BQSR and variant calling steps.


### Sample Data
The pangenome-somatic pipeline accepts either FASTQ files (paired-end only; 1 set per sample) in YAML format:

<pre><code>---
TCGA-FD-A3B5:
    TCGA-FD-A3B5-01A:
        type: tumour
        libraries:
            TCGA-FD-A3B5-01A-11D-A20D-08:
                runlanes:
                    L001:
                        fastq:
                            R1: /path/to/TCGA-FD-A3B5-01A-11D-A20D-08_R1.fastq.gz
                            R2: /path/to/TCGA-FD-A3B5-01A-11D-A20D-08_R2.fastq.gz
    TCGA-FD-A3B5-10A:
        type: normal
        libraries:
            TCGA-FD-A3B5-10A-01D-A20D-08:
                runlanes:
                    L001:
                        fastq:
                            R1: /path/to/TCGA-FD-A3B5-10A-01D-A20D-08_R1.fastq.gz
                            R2: /path/to/TCGA-FD-A3B5-10A-01D-A20D-08_R2.fastq.gz
</code></pre>

See [example_fastq_config.yaml](../examples/example_fastq_config.yaml) for an extended example.

## Workflow
The pangenome-somatic pipeline can be run using a single command, using the completed tool configuration YAML to determine all necessary steps/parameters:

<pre><code>cd /path/to/git/pangenome-pipeline

python scripts/run_pangenome_somatic_pipeline.py \
-t /path/to/project_tool_config.yaml \
-d /path/to/project_fastq_config.yaml \
-o /path/to/project/output_directory \
--alignment \
--strelka \
--mutect2 \
{optional: --remove, --dry-run }
</code></pre>

This will generate the directory structure and code needed. Specifically, the pipeline will run in 2 stages:
1. Data processing
  - alignment to the graph-reference using VG
  - extraction of reads aligned to the linear reference and conversion to BAM format
  - GATK's base-quality score recalibration
  - QC
2. Somatic SNV detection using T/N samples
  - using Strelka2 and Mutect2

### Pipeline arguments
<pre><code>python scripts/run_pangenome_somatic_pipeline.py -h </code></pre>
options:

| Option (short) | Option (full) | Description |
| :---           | :---:         | ---:          |
| -h             | --help        | show this help message and exit   |
| -t             | --tool_config | path to tool config file in YAML format (reference files and tool versions)   |
| -d             | --data_config | path to data config file in YAML format (patient/sample fastq files)   |
| -o             | --out_dir     | path to project output directory   |
|                | --alignment   | should the alignment step be run?   |
|                | --strelka     | should the strelka variant-calling step be run?   |
|                | --mutect2     | should the mutect2 variant-calling step be run?   |
|                | --remove      | if used, intermediate files will be removed   |
|                | --dry_run     | if used, directories/scripts will be created but no jobs will be submitted  |

### Step 1:
<pre><code># 1) Data processing
python scripts/run_alignments.py \
-t /path/to/project_tool_config.yaml \
-d /path/to/project_fastq_config.yaml \
-o /path/to/project/output_directory/alignments \
-b /path/to/project/output_directory/alignments/processed_bams.yaml \
{optional: --remove, --dry-run }
</code></pre>

From project_tool_config.yaml:
<pre><code># reference files
reference: /path/to/resources/references/linear_reference.fa
ref_gbz: /path/to/resources/references/graph_reference.gbz
ref_dist: /path/to/resources/references/graph_reference.dist
ref_min: /path/to/resources/references/graph_reference.min
ref_paths: /path/to/resources/references/primary_contigs.txt
ref_pattern: # the stem added to contig names in the graph (ie, "GRCh38#0#")

known_indels: /path/to/resources/gatk-resource-bundle/Homo_sapiens_assembly38.known_indels.vcf.gz
germline_resource: /path/to/resources/germline_resource/af-only-gnomad.hg38.vcf.gz
target_regions: /path/to/resources/capture_panels/Agilent_SureSelect_All_Exon_V5.bed.gz

# enter the executable for the following tools, as found on your specific system
vg_path: 'apptainer exec /workspace/lab/houlahanlab/software/apptainers/vg_v1.61.0.sif vg'
gatk_path: 'gatk'

# stage parameters for slurm
giraffe:
  time: '8:0:0'
  mem: '48G'
  cpus: 2
gatk:
  time: '8:0:0'
  mem: '8G'
  cpus: 2
</code></pre>

For each sample in project_fastq_config.yaml, this will create a call-ready BAM file: /path/to/project/output_directory/alignments/PATIENT/SAMPLE/SAMPLE_recalibrated.bam

This will output the paths for these BAM files here: /path/to/project/output_directory/alignments/processed_bams.yaml

QC metrics from all samples will also collected and output here: 
- /path/to/project/output_directory/alignments/date_PROJECTID_vg_alignment_statistics.tsv
- /path/to/project/output_directory/alignments/date_PROJECTID_gatk_alignment_statistics.tsv
- /path/to/project/output_directory/alignments/date_PROJECTID_contamination_estimates.tsv

### Step 2A Strelka2:
<pre><code># 2A) Variant Calling with Strelka2
python scripts/run_strelka.py \
-t /path/to/project_tool_config.yaml \
-d /path/to/project/alignments/processed_bams.yaml \
-o /path/to/project/output_directory/variant-calling/Strelka \
{optional: --remove, --dry-run }
</code></pre>

From project_tool_config.yaml:
<pre><code># reference files
reference: /path/to/resources/references/linear_reference.fa
target_regions: /path/to/resources/capture_panels/Agilent_SureSelect_All_Exon_V5.bed.gz

# enter the executable for the following tools, as found on your specific system
conda_path: /workspace/lab/houlahanlab/software/micromamba/bin/conda # used for bcftools/bgzip/tabix
strelka_path: 'apptainer exec /workspace/lab/houlahanlab/software/apptainers/strelka_2.9.10--h9ee0642_1.sif'

# stage parameters for slurm
# strelka requires the reference fasta and test BAMs to contain the same contig information
# if using a subset of contigs (ie, primary only), provide the path to a subsetted fasta here
# otherwise the above reference will be used
strelka:
  ref: /path/to/resources/references/linear_reference_subset.fa
  time: '4:0:0'
  mem: '12G'
  cpus: 1
</code></pre>

For each tumour sample in processed_bams.yaml, this will create a filtered VCF (somatic SNVs with PASS filter): /path/to/project/output_directory/variant-calling/Strelka/PATIENT/SAMPLE/SAMPLE_Strelka_pass_variants.vcf.gz

### Step 2B Mutect2:
<pre><code># 2B) Variant Calling with Mutect2
python scripts/run_mutect2.py \
-t /path/to/project_tool_config.yaml \
-d /path/to/project/alignments/processed_bams.yaml \
-o /path/to/project/output_directory/variant-calling/Mutect2 \
{optional: --remove, --dry-run }
</code></pre>

From project_tool_config.yaml:
<pre><code># reference files
reference: /path/to/resources/references/linear_reference.fa
target_regions: /path/to/resources/capture_panels/Agilent_SureSelect_All_Exon_V5.bed.gz
pon: /path/to/resources/germline_resource/1000g_pon.hg38.vcf.gz
germline_resource: /path/to/resources/germline_resource/af-only-gnomad.hg38.vcf.gz

# enter the executable for the following tools, as found on your specific system
conda_path: /workspace/lab/houlahanlab/software/micromamba/bin/conda # used for bcftools/bgzip/tabix
gatk_path: 'gatk'

# stage parameters for slurm
mutect2:
  time: '4:0:0'
  mem: '4G'
  cpus: 1
</code></pre>

For each tumour sample in processed_bams.yaml, this will create a filtered VCF (somatic SNVs and INDELs with PASS filter): /path/to/project/output_directory/variant-calling/Mutect2/PATIENT/SAMPLE/SAMPLE_Mutect2_pass_variants.vcf.gz

