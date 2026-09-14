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
The pangenome-pipeline was developed using python and is currently compatible using the Slurm workload manager.
Required tools/packages (these should all be available in the shared lab environment):
- python3
- python packages: os, subprocess, yaml, argparse, datetime, simple_slurm
- VG, GATK (>= v4), Strelka, SomaticSniper (with bam-readcount), DeepSomatic

## Example
To set up your config files, see includes examples:
- example_data.yaml: this lists the paths to your fastq files while providing the donor/sample/library/lane information
- example_tool.yaml: this lists the paths to your desired reference genome and provides the compute resource requirements for each slurm job

<pre><code>cd /workspace/lab/houlahanlab/your_macid/git/pangenome-pipeline

# test run the alignment pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python run_alignments.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/output/alignments -b /path/to/output/alignments/recalibrated_bam_files.yaml --remove --dry_run

# now, submit the jobs to slurm!
python run_alignments.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/output/alignments -b /path/to/output/alignments/recalibrated_bam_files.yaml --remove
</code></pre>

<pre><code> # test run the Strelka variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python run_strelka.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Strelka --remove --dry_run

# now, submit the jobs to slurm!
python run_strelka.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Strelka --remove

# test run the SomaticSniper variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python run_somaticsniper.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/SomaticSniper --remove --dry_run

# now, submit the jobs to slurm!
python run_somaticsniper.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/SomaticSniper --remove

# test run the Mutect2 variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python run_mutect2.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Mutect2 --remove --dry_run

# now, submit the jobs to slurm!
python run_mutect2.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/Mutect2 --remove

# test run the DeepSomatic variant-calling pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python run_deepsomatic.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/DeepSomatic --remove --dry_run

# now, submit the jobs to slurm!
python run_deepsomatic.py -t /path/to/your/tool_config.yaml -d /path/to/output/alignments/recalibrated_bam_files.yaml -o /path/to/variant_callers/DeepSomatic --remove
</code></pre>

