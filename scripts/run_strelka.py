"""This pipeline will run the Strelka2 variant-caller on VG-aligned, 
GATK-recalibrated BAM files."""

import os
import yaml
import argparse

from datetime import datetime
from utilities import write_script,collect_job_stats,is_file_missing_or_empty

### subroutines
def get_sm_tag(bam_file:str) -> str:
    """Extract SM tag from BAM header (this will be the VCF sample ID)."""

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Retrieve the 'RG' (Read Group) lines from the header dictionary
        read_groups = bam.header.get('RG', [])
        
        # Extract the 'SM' tag from each read group if it exists
        sample_name = {rg['SM'] for rg in read_groups if 'SM' in rg}
        
    return list(sample_name)[0]


def run_strelka(strelka_executable:str, reference:str, targets:str, tumour:str, normal:str, rundir:str) -> str:
    """Run Strelka on the tumour/normal pair."""

    tool_comment = "# run strelka"
    tool_command1 = [
        strelka_executable, "configureStrelkaSomaticWorkflow.py",
        "--tumourBam", tumour,
        "--normalBam", normal,
        "--referenceFasta", reference,
        "--runDir", rundir
    ]

    if targets:
        tool_command1.append(f"--exome --callRegions {targets}")

    tool_command2 = f"{strelka_executable} {rundir}/runWorkflow.py -m local -j 8"

    final_command = tool_comment + "\n" + " ".join(map(str, tool_command1)) + "\n\n" + tool_command2
    return final_command


def run_filter_vcf(bcftools:str, tabix:str, vcf:str, filtered_vcf:str) -> str:
    """Filter raw Strelka variant calls."""

    tool_comment = "# filter variants to only PASS"
    tool_command = [
        f"{bcftools} view -f PASS -Oz -o {filtered_vcf} {vcf}",
        f"{tabix} {filtered_vcf}",
    ]

    final_command = tool_comment + "\n" + "\n".join(map(str, tool_command))
    return final_command


def strelka_pipeline(tool_config:str, data_config:str, out_dir:str, mode:str = 'slurm', remove:bool = False, dry_run:bool = False, wait_for:str = None) -> str:

    # load tool/paramter information
    with open(tool_config, 'r') as file:
        tool_data = yaml.safe_load(file)

    # set up logging
    log_directory = os.path.join(out_dir, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    current_run = datetime.now().strftime("%Y-%m-%d-%H-%M")
    if dry_run:
        log_file = f"{log_directory}/strelka_pipeline__dryrun.log"
    else:
        log_file = f"{log_directory}/strelka_pipeline__{current_run}.log"
    
    with open(log_file, 'w') as log:
        log.write("---\n")
        log.write("Running the Strelka variant-calling pipeline.\n")
        log.write(f"  Tool config used: {tool_config}\n")
        log.write(f"  Data config used: {data_config}\n")
        log.write(f"  Outputting to: {out_dir}\n")
        if remove:
            log.write("  >> clean mode selected (--remove); intermediate files will be removed upon completion\n")
        if mode == 'local':
            log.write("  >> local mode selected (--mode local); jobs will be run locally if able\n")
        elif mode == 'slurm':
            log.write("  >> slurm mode selected (--mode slurm); jobs will be submitted using Slurm\n")
        if dry_run:
            log.write("  >> trial mode selected (--dry_run); jobs will not be submitted\n")

    # extract tool paths from config
    bcftools_executable, strelka_executable, tabix_executable = None, None, None

    if 'container_path' in tool_data:
        bcftools_executable = f"{tool_data['container_path']} bcftools"
        strelka_executable = f"{tool_data['container_path']}"
        tabix_executable = f"{tool_data['container_path']} tabix"
    else:
        bcftools_executable = tool_data['bcftools_path']
        strelka_executable = tool_data['strelka_path']
        tabix_executable = tool_data['tabix_path']

    # extract reference files from config
    reference = tool_data['reference']
    if tool_data['strelka']['ref']:
        reference = tool_data['strelka']['ref']

    # add details to log file
    with open(log_file, 'a') as log:
        log.write(f"\nLinear reference file used is {reference}\n")
        log.write(f"\nUsing Strelka version: {strelka_executable}\n")
        log.write("Using samtools/bcftools version: 1.23.1\n")

    # load sample information
    with open(data_config, 'r') as file:
        sample_data = yaml.safe_load(file)

    n_patients = len(sample_data.keys())
    
    with open(log_file, 'a') as log:
        log.write(f"\n----------\nProcessing {n_patients} patients...\n")

    # is this pipeline waiting for any prior jobs?
    prior_jobs = []
    if wait_for:
        prior_jobs = str(wait_for).split(",")
        
    # initiate object to hold job_ids
    all_jobs = []

    # work through each sample
    for donor in sample_data.keys():

        with open(log_file, 'a') as log:
            log.write(f"\nProcessing donor: {donor}\n")

        patient_directory = os.path.join(out_dir, donor)
        os.makedirs(patient_directory, exist_ok=True)

        patient_jobs = []
        final_files = []
        cleanup_files = []

        # check for 'normal' sample
        normal_samples = list(sample_data[donor]['normal'].keys())

        normal_id, normal_bam = None, None
        if len(normal_samples) > 0:
            normal_id = normal_samples[0]
            normal_bam = sample_data[donor]['normal'][normal_id]
            
            with open(log_file, 'a') as log:
                log.write(f"\n  Using normal BAM: {normal_bam}\n")

        else:
            with open(log_file, 'a') as log:
                log.write(f"  <<Error: No normal BAM provided. Skipping {donor}...>>\n")
            continue

        # process tumour sample(s)
        n_samples = len(sample_data[donor]['tumour'].keys())

        with open(log_file, 'a') as log:
            log.write(f"\n  Processing {n_samples} tumour samples for donor: {donor}\n")

        for sample in (sample_data[donor]['tumour'].keys()):

            with open(log_file, 'a') as log:
                log.write(f"    Processing {sample}\n")

            tumour_bam = sample_data[donor]['tumour'][sample]

            smp_directory = os.path.join(patient_directory, sample)
            os.makedirs(smp_directory, exist_ok=True)

            # indicate output files
            strelka_vcf = os.path.join(smp_directory, "results/variants/somatic.snvs.vcf.gz")
            filtered_vcf = os.path.join(smp_directory, f"{sample}_Strelka_pass_variants.vcf.gz")

            ## Step 1: Strelka
            strelka_command, filter_command = None, None

            # run variant caller
            if is_file_missing_or_empty(strelka_vcf):
                strelka_command = run_strelka(
                    strelka_executable,
                    reference = reference,
                    targets = tool_data['target_regions'],
                    tumour = tumour_bam, 
                    normal = normal_bam, 
                    rundir = smp_directory
                )

            cleanup_files.append(os.path.join(smp_directory, 'workspace'))

            # filter variant calls
            if is_file_missing_or_empty(filtered_vcf):
                filter_command = run_filter_vcf(
                    bcftools = bcftools_executable,
                    tabix = tabix_executable,
                    vcf = strelka_vcf, 
                    filtered_vcf = filtered_vcf
                )
            
            call_commands = list(filter(None, [strelka_command, filter_command]))

            # write slurm script and sbatch (unless dry-run)
            with open(log_file, 'a') as log:
                log.write(f"      Submitting job for Strelka...\n")

                if len(call_commands) > 0 and is_file_missing_or_empty(filtered_vcf):   

                    strelka_job_id, strelka_script = write_script(
                        log_dir = log_directory,
                        dry_run = dry_run,
                        mode = mode,
                        job_name = f"run_strelka_{sample}",
                        command = call_commands,
                        depends = prior_jobs,
                        job_parameters = tool_data['strelka']
                    )

                    patient_jobs.append(strelka_job_id)

                    log.write(f"       >> RUNSCRIPT: {strelka_script}\n")
                    if not dry_run:
                        log.write(f"       >> job id is {strelka_job_id}\n")

                else:
                    log.write(f"       >> Strelka step already completed; no jobs to submit.\n")

            final_files.append(filtered_vcf)

        # do cleanup once per donor
        if remove:
            cleanup_command = "# check for final VCF(s), then cleanup intermediates"
            cleanup_command += f"\nif [[ -s " + " && -s ".join(final_files) + " ]]; then"
            cleanup_command += "\n  rm -rf " + "\n  rm -rf ".join(cleanup_files)
            cleanup_command += "\nelse"
            cleanup_command += "\n  echo 'Final output VCF(s) is missing; not removing intermediate files!'"
            cleanup_command += "\n  exit 1"
            cleanup_command += "\nfi"
        
            with open(log_file, 'a') as log:
                log.write(f"\n  Running clean-up step for {donor}...\n")

                if len(patient_jobs) > 0:
                    clean_job_id, clean_script = write_script(
                        log_dir = log_directory,
                        dry_run = dry_run,
                        mode = mode,
                        job_name = f"run_cleanup_{donor}",
                        command = [cleanup_command],
                        job_parameters = {"mem": "256M", "time": "1:0:0"},
                        depends = patient_jobs
                    )

                    patient_jobs.append(clean_job_id)

                    log.write(f"    >> RUNSCRIPT: {clean_script}\n")
                    if not dry_run:
                        log.write(f"    >> job id is {clean_job_id}\n")

        all_jobs += patient_jobs

    # add details to log file
    with open(log_file, 'a') as log:
        log.write("\n----------\nAll patients processed.\n")

        n_jobs = len(all_jobs)
        log.write(f"\nTotal number of jobs processed: {n_jobs}\n")

    # write slurm script and sbatch (unless dry-run)
    sacct_job_id, sacct_script = None, None
    
    if mode == 'slurm' and not dry_run:
        sacct_out = f"{out_dir}/logs/slurm_job_metrics_{current_run}.out"
        sacct_command = collect_job_stats(all_jobs, sacct_out)
        sacct_job_id, sacct_script = write_script(
            log_dir = log_directory,
            dry_run = dry_run,
            mode = mode,
            job_name = "collect_job_metrics",
            command = [sacct_command],
            job_parameters = {"mem": "256M", "time": "1:0:0"},
            depends = all_jobs 
        )

        with open(log_file, 'a') as log:
            log.write(f"  >> Final job id is: {sacct_job_id}\n")

    return (sacct_job_id)


### main
if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description = "Run Strelka on T/N BAMs")
    parser.add_argument("-t", "--tool_config", help = "tool config (reference files and tool versions)")
    parser.add_argument("-d", "--data_config", help = "data config (bams listed in yaml format)")
    parser.add_argument("-o", "--out_dir", help = "path to output directory")
    parser.add_argument("-m", "--mode", help = "run mode (local or slurm); defaults to slurm", default = 'slurm')
    parser.add_argument("--remove", action = "store_true", help = "should intermediates be removed?")
    parser.add_argument("--dry_run", action = "store_true", help = "is this a test?")
    parser.add_argument("--wait_for", help = "prior job ids (as a comma-separated list) that this function is dependent on")
    args = parser.parse_args()

    strelka_pipeline(
        tool_config = args.tool_config,
        data_config = args.data_config,
        out_dir = args.out_dir,
        mode = args.mode,
        remove = args.remove,
        dry_run = args.dry_run,
        wait_for = args.wait_for
    )
