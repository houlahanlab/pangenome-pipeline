# pangenome-pipeline
Pipelines for pangenome-based analyses.

## Introduction
This is a collection of pipelines to be used for NGS analyses using the pangenome reference.
- run_alignment.py: this will align reads from fastq files (paired end) to a graph-based reference file

Start by creating a clone of the repository:

<pre><code>mkdir /workspace/lab/houlahanlab/<user>/git
cd /workspace/lab/houlahanlab/<user>/git
git clone git@github.com:houlahanlab/pangenome-pipeline.git
</code></pre>

## Dependencies
The pangenome-pipeline was developed using python and is currently compatible using the Slurm workload manager.
- python packages os, subprocess, yaml, argparse, datetime, simple_slurm


## Example
To set up your config files, see includes examples:
- example_data.yaml: this lists the paths to your fastq files while providing the donor/sample/library/lane information
- example_tool.yaml: this lists the paths to your desired reference genome and provides the compute resource requirements for each slurm job

<pre><code>cd /workspace/lab/houlahanlab/<user>/git/pangenome-pipeline

# test run the pipeline; this will generate the output directory structure and scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python run_alignments.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/some/output/directory -b /path/to/some/output/directory/output_bam_files.yaml --remove --dry_run

# now, submit the jobs to slurm!
python run_alignments.py -t /path/to/your/tool_config.yaml -d /path/to/your/data_config.yaml -o /path/to/some/output/direct
ory -b /path/to/some/output/directory/output_bam_files.yaml --remove
</code></pre>

