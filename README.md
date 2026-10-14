# pangenome-pipeline (v1.0.0)
Pipelines for pangenome-based analyses.

## Introduction
This is a collection of pipelines to be used for NGS analyses using the pangenome reference.
* run_alignment.py: this will align reads from fastq files (paired end) to a graph-based reference file
* run_strelka.py: this will run the Strelka2 variant-caller on paired T/N samples
* run_somaticsniper.py: this will run the SomaticSniper variant-caller on paired T/N samples
* run_mutect2.py: this will run the Mutect2 variant-caller on paired T/N samples or T-only samples
* run_deepsomatic.py: this will run the DeepSomatic variant-caller on paired T/N samples or T-only samples

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
    * note that all of these are conventiently packaged using docker [here](/path/to/some/repo)

Recommended/optional:
* apptainer


## Examples
To run the benchmarking example using GIAB's HG008 T/N sample, see [here](./examples/benchmarking/README).

To run your own project, set up your config files (see those in examples for templates):
* [example_data.yaml](./examples/example_data.yaml): this lists the paths to your fastq files while providing the donor/sample/library/lane information for use in the SAM header
* [example_tool.yaml](./examples/example_tool.yaml): this lists the paths to your desired reference genome and provides the compute resource requirements for each slurm job

See [documentation](./docs/pipeline-how-to.md) for further examples.

