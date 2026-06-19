"""This pipeline will run the Mutect2 variant-caller on VG-aligned, 
GATK-recalibrated BAM files."""

import os
import subprocess
import yaml
import argparse
#import pysam

from datetime import datetime
from simple_slurm import Slurm

from utilities import write_script,collect_job_stats,is_file_missing_or_empty,get_java_mem

### subroutines
def get_sm_tag(bam_file):
    """Extract SM tag from BAM header (this will be the VCF sample ID)."""

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Retrieve the 'RG' (Read Group) lines from the header dictionary
        read_groups = bam.header.get('RG', [])
        
        # Extract the 'SM' tag from each read group if it exists
        sample_name = {rg['SM'] for rg in read_groups if 'SM' in rg}
        
    return list(sample_name)[0]


def run_mutect2(reference, tumour, vcf, normal, normal_id, tmp_dir, java_mem, n_cpus=None, gnomad=None, pon=None, targets=None):
    """Run Mutect2 on either T/N pair or in tumour-only mode"""

    java = f"--java-options '{java_mem} -XX:ParallelGCThreads=2 -XX:ConcGCThreads=1 -Djava.io.tmpdir={tmp_dir}'"

    tool_comment = "# run mutect2"
    tool_command = [
        "gatk", java, "Mutect2",
        "-R", reference,
        "-I", tumour,
        "-O", vcf
    ]

    if normal:
        tool_command.append(f"-I {normal} -normal {normal_id}")

    if n_cpus:
        tool_command.append(f"--native-pair-hmm-threads {n_cpus}")

    if targets:
        tool_command.append(f"--interval-padding 10 --interval-set-rule INTERSECTION --intervals {targets}")
    
    if gnomad:
        tool_command.append(f"--germline-resource {gnomad}")

    if pon:
        tool_command.append(f"--panel-of-normals {pon}")


    final_command = tool_comment + "\n" + " ".join(map(str, tool_command))
    return final_command


def run_filter_vcf(conda_path, reference, vcf, filtered_vcf, java_mem, tmp_dir):
    """Filter raw Mutect2 variant calls."""

    tmp_vcf = vcf.replace('raw','variants')

    java = f"--java-options '{java_mem} -XX:ParallelGCThreads=2 -XX:ConcGCThreads=1 -Djava.io.tmpdir={tmp_dir}'"

    tool_comment1 = "# annotate mutect2 variants with filters"
    tool_command1 = [
        "gatk", java, "FilterMutectCalls",
        "-R", reference,
        "-V", vcf,
        "-O", tmp_vcf
    ]

    tool_comment2 = "# filter variants to only PASS"
    tool_command2 = [
        f'eval "$({conda_path} shell.bash hook)"',
        f"conda activate samtools_env",
        f"bcftools view -f PASS {tmp_vcf} > {filtered_vcf}",
        f"bgzip {filtered_vcf}",
        f"tabix -p vcf {filtered_vcf}.gz",
        "conda deactivate"
    ]

    final_command = tool_comment1 + "\n" + " ".join(map(str, tool_command1))
    final_command += "\n" + tool_comment2 + "\n" + "\n".join(map(str, tool_command2))
    return final_command


def main(tool_config:str, data_config:str, out_dir:str, remove:bool = False, dry_run:bool = False):

    # set up logging
    log_directory = os.path.join(out_dir, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    current_run = datetime.now().strftime("%Y-%m-%d-%H-%M")
    if dry_run:
        log_file = f"{log_directory}/mutect2_pipeline__dryrun.log"
    else:
        log_file = f"{log_directory}/mutect2_pipeline__{current_run}.log"
    
    log = open(log_file, 'w')

    print("---", file = log)
    print("Running the Mutect2 variant-calling pipeline.", file = log)
    print(f"  Tool config used: {tool_config}", file = log)
    print(f"  Data config used: {data_config}", file = log)
    print(f"  Outputting to: {out_dir}", file = log)
    if remove:
        print("  >> clean mode selected (--remove); intermediate files will be removed upon completion", file = log)
    if dry_run:
        print("  >> trial mode selected (--dry_run); jobs will not be submitted", file = log)

    # load tool/paramter information
    with open(tool_config, 'r') as file:
        tool_data = yaml.safe_load(file)

    reference = tool_data['reference']

    # load sample information
    with open(data_config, 'r') as file:
        sample_data = yaml.safe_load(file)

    n_patients = len(sample_data.keys())

    # add details to log file
    print(f"\nLinear reference file used is {reference}", file = log)
    print("Using gatk version: 4.3.0.0 [with Picard version 2.27.5]", file = log)
    print(f"\n----------\nProcessing {n_patients} patients...", file = log)

    # initiate object to hold job_ids
    all_jobs = []

    # work through each sample
    for donor in sample_data.keys():

        # check for 'normal' sample
        normal_samples = list(sample_data[donor]['normal'].keys())

        normal_id, normal_bam = None, None
        if len(normal_samples) > 0:
            normal_id = normal_samples[0]
            normal_bam = sample_data[donor]['normal'][normal_id]
            
            if not os.path.exists(normal_bam):
                print(f"<<Error: Normal BAM file {normal_bam} not found. Running tumour-only mode for {donor}...>>", file = log)

        patient_directory = os.path.join(out_dir, donor)
        tmp_directory = os.path.join(patient_directory, 'TEMP')
        os.makedirs(tmp_directory, exist_ok=True)

        patient_jobs = []
        final_files = []
        cleanup_files = []

        n_samples = len(sample_data[donor].keys())
        print(f"\nProcessing {n_samples} tumour samples for donor: {donor}", file = log)

        for sample in (sample_data[donor]['tumour'].keys()):

            print(f"  Processing {sample}", file = log)

            tumour_bam = sample_data[donor]['tumour'][sample]

            # process all fastqs for each sample
            if not os.path.exists(tumour_bam):
                print(f"<<Error: BAM file {tumour_bam} not found. Skipping {sample}...>>", file = log)
                continue
 
            smp_directory = os.path.join(patient_directory, sample)
            os.makedirs(smp_directory, exist_ok=True)

            # indicate output files
            mutect2_vcf = os.path.join(smp_directory, f"{sample}_Mutect2_raw.vcf.gz")
            filtered_vcf = os.path.join(smp_directory, f"{sample}_Mutect2_pass_variants.vcf")

            job_id = []

            ## Step 1: Mutect2
            mutect2_command, filter_command = None, None
            java_mem = get_java_mem(tool_data['mutect2']['mem'])

            # run alignment
            if is_file_missing_or_empty(mutect2_vcf):
                mutect2_command = run_mutect2(
                    reference = reference,
                    targets = tool_data['target_regions'],
                    tumour = tumour_bam, 
                    normal = normal_bam,
                    normal_id = normal_id,
                    vcf = mutect2_vcf,
                    tmp_dir = tmp_directory,
                    java_mem = java_mem,
                    gnomad = tool_data['germline_resource'],
                    pon = tool_data['pon']
                )

            cleanup_files.append(mutect2_vcf)
            cleanup_files.append(f"{mutect2_vcf}.tbi")

            # collect metrics for final bam
            if is_file_missing_or_empty(filtered_vcf):
                filter_command = run_filter_vcf(
                    conda_path = tool_data['conda_path'],
                    reference = reference,
                    vcf = mutect2_vcf, 
                    filtered_vcf = filtered_vcf,
                    tmp_dir = tmp_directory,
                    java_mem = java_mem,
                )
            
            call_commands = list(filter(None, [mutect2_command, filter_command]))

            # write slurm script and sbatch (unless dry-run)
            print(f"    Submitting job for Mutect2...", file = log)

            if len(call_commands) > 0 and is_file_missing_or_empty(f"{filtered_vcf}.gz"):

                job_id, script = write_script(
                    args,
                    job_name = f"run_mutect2_{sample}",
                    command = call_commands,
                    job_parameters = tool_data['mutect2']
                    )

                patient_jobs.append(job_id)

                print(f"      >> sbatch {script}", file = log)
                if not dry_run:
                    print(f"      >> job id is {job_id}", file = log)

            else:
                print(f"      >> Mutect2 step already completed; no jobs to submit.", file = log)

            final_files.append(f"{filtered_vcf}.gz")

        # do cleanup once per donor
        if remove:
            cleanup_command = f"if [[ -s " + " && -s ".join(final_files) + " ]]; then"
            cleanup_command += f"\n  rm -rf {tmp_directory}"
            cleanup_command += "\n  rm " + "\n  rm ".join(cleanup_files)
            cleanup_command += "\nelse"
            cleanup_command += "\n  echo 'Final output VCF(s) is missing; not removing intermediate files!'"
            cleanup_command += "\n  exit 1"
            cleanup_command += "\nfi"
        
            print(f"\n  Running clean-up step for {donor}...", file = log)

            if len(patient_jobs) > 0:
                job_id, script = write_script(
                    args,
                    job_name = f"run_cleanup_{donor}",
                    command = [cleanup_command],
                    job_parameters = {"mem": "256M", "time": "1:0:0"},
                    depends = patient_jobs
                    )

                patient_jobs.append(job_id)

                print(f"    >> sbatch {script}", file = log)
                if not dry_run:
                    print(f"    >> job id is {job_id}", file = log)

        all_jobs += patient_jobs

    # write slurm script and sbatch (unless dry-run)
    if not dry_run:
        sacct_out = f"{out_dir}/logs/slurm_job_metrics_{current_run}.out"
        sacct_command = collect_job_stats(all_jobs, sacct_out)
        job_id, script = write_script(
            args,
            job_name = "collect_job_metrics",
            command = [sacct_command],
            job_parameters = {"mem": "256M", "time": "1:0:0"},
            depends = all_jobs 
            )

    # add details to log file
    print("\n----------\nAll patients processed.", file = log)

    n_jobs = len(all_jobs)
    print(f"\nTotal number of jobs to run: {n_jobs}", file = log)

    if not dry_run:
        print(f"  >> Final job id is: {job_id}", file = log)

    log.close()


### main
if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description = "Run Mutect2 on T/N or T-only BAMs")
    parser.add_argument("-t", "--tool_config", help = "tool config (reference files and tool versions)")
    parser.add_argument("-d", "--data_config", help = "data config (bams listed in yaml format)")
    parser.add_argument("-o", "--out_dir", help = "path to output directory")
    parser.add_argument("--remove", action = "store_true", help = "should intermediates be removed?")
    parser.add_argument("--dry_run", action = "store_true", help = "is this a test?")
    args = parser.parse_args()

    main(
        tool_config = args.tool_config, 
        data_config = args.data_config, 
        out_dir = args.out_dir,
        remove = args.remove, 
        dry_run = args.dry_run
    )
