"""This pipeline will run the SomaticSniper variant-caller on VG-aligned, 
GATK-recalibrated BAM files."""

import os
import yaml
import argparse

from datetime import datetime
from utilities import write_script,collect_job_stats,is_file_missing_or_empty

### subroutines
def run_somaticsniper(somsnip_executable:str, reference:str, tumour:str, normal:str, vcf:str) -> str:
    """Run SomaticSniper on the tumour/normal pair."""

    tool_comment = "# run somatic-sniper"
    tool_command = [
        f"{somsnip_executable}/build/bin/bam-somaticsniper",
        "-Q 40 -G -L -F vcf -f", reference,
        tumour,
        normal,
        vcf
    ]

    final_command = tool_comment + "\n" + " ".join(map(str, tool_command))
    return final_command


def run_sniper_filters(somsnip_executable:str, rc_executable:str, samtools:str, bcftools:str, bgzip:str, tabix:str, reference:str, tumour:str, vcf:str, stem:str, targets:str) -> str:
    """Apply variant filters."""

    script_path = f"{somsnip_executable}/src/scripts"
    hc_vcf = f"{stem}_SomaticSniper_hc.vcf"
    out_vcf = f"{stem}_SomaticSniper_hc_filtered.vcf.gz"

    these_regions = ""
    if targets:
        these_regions = f"-R {targets}"

    tool_comment1 = "# apply somaticsniper filters"
    tool_command1 = [
        f"{samtools} calmd -b {tumour} {reference} > {stem}_NM.bam",
        f"{samtools} index {stem}_NM.bam",
        f"{samtools} mpileup -f {reference} {stem}_NM.bam > {stem}.pileup",
        f"perl {script_path}/snpfilter.pl --snp-file {vcf} --indel-file {stem}.pileup --out-file {stem}.SNPfilter",
        f"perl {script_path}/prepare_for_readcount.pl --snp-file {stem}.SNPfilter --out-file {stem}.SNPfilter.pos",
        f"{rc_executable} bam-readcount -w 1 -b 15 -f {reference} -l {stem}.SNPfilter.pos {stem}_NM.bam > {stem}.readcounts",
        f"perl {script_path}/fpfilter.pl --snp-file {stem}.SNPfilter --readcount-file {stem}.readcounts",
        f"perl {script_path}/highconfidence.pl --snp-file {stem}.SNPfilter.fp_pass --out-file {hc_vcf}",
        f"{bgzip} {hc_vcf}",
        f"{tabix} -p vcf {hc_vcf}.gz"
    ]

    tool_comment2 = "# filter variants by depth and capture regions (if provided)"
    tool_command2 = [
        f"{bcftools} view -e 'FMT/DP < 5' {these_regions} -Oz -o {out_vcf} {hc_vcf}.gz",
        f"{tabix} -p vcf {out_vcf}"
    ]

    final_command = tool_comment1 + "\n" + "\n".join(map(str, tool_command1))
    final_command += "\n" + tool_comment2 + "\n" + "\n".join(map(str, tool_command2))
    return final_command


def somaticsniper_pipeline(tool_config:str, data_config:str, out_dir:str, mode:str = 'slurm', remove:bool = False, dry_run:bool = False, wait_for:str = None) -> str:

    # load tool/paramter information
    with open(tool_config, 'r') as file:
        tool_data = yaml.safe_load(file)
    
    # set up logging
    log_directory = os.path.join(out_dir, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    current_run = datetime.now().strftime("%Y-%m-%d-%H-%M")
    if dry_run:
        log_file = f"{log_directory}/somaticsniper_pipeline__dryrun.log"
    else:
        log_file = f"{log_directory}/somaticsniper_pipeline__{current_run}.log"
    
    with open(log_file, 'w') as log:
        log.write("---\n")
        log.write("Running the SomaticSniper variant-calling pipeline.\n")
        log.write(f"  Tool config used: {tool_config}\n")
        log.write(f"  Data config used: {data_config}\n")
        log.write(f"  Outputting to: {out_dir}\n")
        if remove:
            log.write("  >> clean mode selected (--remove); intermediate files will be removed upon completion\n")
        if dry_run:
            log.write("  >> trial mode selected (--dry_run); jobs will not be submitted\n")

    # extract tool paths from config
    somaticsniper_executable = tool_data['somaticsniper_path']
    readcount_executable = tool_data['readcount_path']

    # extract tool paths from config
    somaticsniper_executable = tool_data['somaticsniper_path']
    readcount_executable = tool_data['readcount_path']

    samtools_executable, bcftools_executable, bgzip_executable, tabix_executable = None, None, None, None

    if 'container_path' in tool_data:
        samtools_executable = f"{tool_data['container_path']} samtools"
        bcftools_executable = f"{tool_data['container_path']} bcftools"
        bgzip_executable = f"{tool_data['container_path']} bgzip"
        tabix_executable = f"{tool_data['container_path']} tabix"
    else:
        samtools_executable = tool_data['samtools_path']
        bcftools_executable = tool_data['bcftools_path']
        bgzip_executable = tool_data['bgzip_path']
        tabix_executable = tool_data['tabix_path']

    # extract reference files from config
    reference = tool_data['reference']

     # add details to log file
    with open(log_file, 'a') as log:
        log.write(f"\nLinear reference file used is {reference}\n")
        log.write(f"\nUsing SomaticSniper version: {somaticsniper_executable}\n")
        log.write(f"\nUsing bam-readcount version: {readcount_executable}\n")
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

        patient_directory = os.path.join(out_dir, donor)
        os.makedirs(patient_directory, exist_ok=True)

        patient_jobs = []
        final_files = []
        cleanup_files = []

        n_samples = len(sample_data[donor].keys())
        with open(log_file, 'a') as log:
            log.write(f"\nProcessing {n_samples} tumour samples for donor: {donor}\n")

        for sample in (sample_data[donor]['tumour'].keys()):

            with open(log_file, 'a') as log:
                log.write(f"  Processing {sample}\n")

            tumour_bam = sample_data[donor]['tumour'][sample]

            smp_directory = os.path.join(patient_directory, sample)
            os.makedirs(smp_directory, exist_ok=True)

            # indicate output files
            somaticsniper_vcf = os.path.join(smp_directory, f"{sample}_somaticsniper_variants.vcf")
            nm_bam = f"{sample}_NM.bam"
            tumour_pileup = f"{sample}.pileup"
            tumour_readcounts = f"{sample}.readcounts"
            filtered_vcf = os.path.join(smp_directory, f"{sample}_SomaticSniper_hc_filtered.vcf.gz")

            ## Step 1: SomaticSniper
            cd_command = f"# move to working directory\ncd {smp_directory}"
            somsnip_command, filter_command = None, None

            # run variant caller
            if is_file_missing_or_empty(somaticsniper_vcf):
                somsnip_command = run_somaticsniper(
                    somsnip_executable = somaticsniper_executable,
                    reference = reference,
                    tumour = tumour_bam, 
                    normal = normal_bam, 
                    vcf = os.path.basename(somaticsniper_vcf)
                )

            cleanup_files.append(f"{smp_directory}/{somaticsniper_vcf}")

            # filter variant calls
            if is_file_missing_or_empty(filtered_vcf):
                filter_command = run_sniper_filters(
                    somsnip_executable = somaticsniper_executable,
                    rc_executable = readcount_executable,
                    samtools = samtools_executable,
                    bcftools = bcftools_executable,
                    bgzip = bgzip_executable,
                    tabix = tabix_executable,
                    reference = reference,
                    targets = tool_data['target_regions'],
                    tumour = tumour_bam,
                    vcf = os.path.basename(somaticsniper_vcf),
                    stem = sample
                )
            
            cleanup_files.append(f"{smp_directory}/{nm_bam}")
            cleanup_files.append(f"{smp_directory}/{nm_bam}.bai")
            cleanup_files.append(f"{smp_directory}/{tumour_pileup}")
            cleanup_files.append(f"{smp_directory}/{tumour_readcounts}")
            cleanup_files.append(f"{smp_directory}/{sample}.SNPfilter*")
        
            call_commands = list(filter(None, [somsnip_command, filter_command]))

            # write slurm script and sbatch (unless dry-run)
            with open(log_file, 'a') as log:
                log.write(f"    Submitting job for SomaticSniper...\n")

                if len(call_commands) > 0 and is_file_missing_or_empty(filtered_vcf):   

                    sniper_job_id, sniper_script = write_script(
                        log_dir = log_directory,
                        dry_run = dry_run,
                        mode = mode,
                        job_name = f"run_somaticsniper_{sample}",
                        command = [cd_command] + call_commands,
                        depends = prior_jobs,
                        job_parameters = tool_data['somaticsniper']
                    )

                    patient_jobs.append(sniper_job_id)

                    log.write(f"      >> RUNSCRIPT: {sniper_script}\n")
                    if not dry_run:
                        log.write(f"      >> job id is {sniper_job_id}\n")

                else:
                    log.write(f"      >> SomaticSniper step already completed; no jobs to submit.\n")

            final_files.append(filtered_vcf)

        # do cleanup once per donor
        if remove:
            cleanup_command = "# check for final VCF(s), then cleanup intermediates"
            cleanup_command += f"\nif [[ -s " + " && -s ".join(final_files) + " ]]; then"
            cleanup_command += "\n  rm " + "\n  rm ".join(cleanup_files)
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

            log.write(f"  >> Final job id is: {sacct_job_id}\n")

    return (sacct_job_id)


### main
if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description = "Run SomaticSniper on T/N BAMs")
    parser.add_argument("-t", "--tool_config", help = "tool config (reference files and tool versions)")
    parser.add_argument("-d", "--data_config", help = "data config (bams listed in yaml format)")
    parser.add_argument("-o", "--out_dir", help = "path to output directory")
    parser.add_argument("-m", "--mode", help = "run mode (local or slurm); defaults to slurm", default = 'slurm')
    parser.add_argument("--remove", action = "store_true", help = "should intermediates be removed?")
    parser.add_argument("--dry_run", action = "store_true", help = "is this a test?")
    parser.add_argument("--wait_for", help = "prior job ids (as a comma-separated list) that this function is dependent on")
    args = parser.parse_args()

    somaticsniper_pipeline(
        tool_config = args.tool_config, 
        data_config = args.data_config, 
        out_dir = args.out_dir,
        mode = args.mode,
        remove = args.remove, 
        dry_run = args.dry_run,
        wait_for = args.wait_for
    )
