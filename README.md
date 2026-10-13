# pangenome-pipeline (v1.0.0)
Pipelines for pangenome-based analyses.

## Introduction
This is a collection of pipelines to be used for NGS analyses using the pangenome reference.
- run_alignment.py: this will align reads from fastq files (paired end) to a graph-based reference file
- run_strelka.py: this will run the Strelka2 variant-caller on paired T/N samples
- run_somaticsniper.py: this will run the SomaticSniper variant-caller on paired T/N samples
- run_mutect2.py: this will run the Mutect2 variant-caller on paired T/N samples or T-only samples
- run_deepsomatic.py: this will run the DeepSomatic variant-caller on paired T/N samples or T-only samples

Start by creating a clone of the repository:

<pre><code>mkdir /workspace/lab/houlahanlab/your_macid/git
cd /workspace/lab/houlahanlab/your_macid/git
git clone git@github.com:houlahanlab/pangenome-pipeline.git
</code></pre>

## Dependencies
The pangenome-pipeline was developed using python. Jobs are managed either using the Slurm workload manager (parallization across samples) or using a local-only method (each step runs in sequence).

Required tools/packages:
* python3
 * python packages: os, subprocess, yaml, argparse, datetime, simple_slurm
* R (>= 4.5.3)
 * R packages: argparse

* VG, samtools, bcftools, bgzip, tabix, GATK (>= v4), Strelka, SomaticSniper (with bam-readcount), DeepSomatic

Recommended/optional:
* apptainer

## Example
To set up your config files, see includes examples:
* example_data.yaml: this lists the paths to your fastq files while providing the donor/sample/library/lane information
* example_tool.yaml: this lists the paths to your desired reference genome and provides the compute resource requirements for each slurm job

<pre><code>cd /workspace/lab/houlahanlab/your_macid/git/pangenome-pipeline

# test run the pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python scripts/run_pangenome_somatic_pipeline.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/output --alignment --mutect2 --strelka --somaticsniper --deepsomatic -m slurm --remove --dry_run

# now, submit the jobs to slurm!
python scripts/run_pangenome_somatic_pipeline.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/output --alignment --mutect2 --strelka --somaticsniper --deepsomatic -m slurm --remove
</code></pre>

Each step can be run separately:
<pre><code>cd /workspace/lab/houlahanlab/your_macid/git/pangenome-pipeline

# test run the alignment pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python scripts/run_alignments.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/output/alignments -b /path/to/output/alignments/recalibrated_bam_files.yaml -m slurm --remove --dry_run

# now, submit the jobs to slurm!
python scripts/run_alignments.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/output/alignments -b /path/to/output/alignments/recalibrated_bam_files.yaml -m slurm --remove
</code></pre>

<pre><code> # test run the Strelka variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python scripts/run_strelka.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Strelka -m slurm --remove --dry_run

# now, submit the jobs to slurm!
python scripts/run_strelka.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Strelka -m slurm --remove

# test run the SomaticSniper variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python scripts/run_somaticsniper.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/SomaticSniper -m slurm --remove --dry_run

# now, submit the jobs to slurm!
python scripts/run_somaticsniper.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/SomaticSniper -m slurm --remove

# test run the Mutect2 variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python scripts/run_mutect2.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Mutect2 -m slurm --remove --dry_run

# now, submit the jobs to slurm!
python scripts/run_mutect2.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Mutect2 -m slurm --remove

# test run the DeepSomatic variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python scripts/run_deepsomatic.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/DeepSomatic -m slurm --remove --dry_run

# now, submit the jobs to slurm!
python scripts/run_deepsomatic.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/DeepSomatic -m slurm --remove
</code></pre>

