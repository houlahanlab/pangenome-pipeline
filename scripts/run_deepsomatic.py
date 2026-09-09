"""This pipeline will run the DeepSomatic variant-caller on VG-aligned, 
GATK-recalibrated BAM files."""

import os
import yaml
import argparse

from datetime import datetime
from utilities import write_script,collect_job_stats,is_file_missing_or_empty

### subroutines
def run_deepsomatic(tool:str, model:str, reference:str, tumour:str, tumour_id:str, vcf:str, tmp_dir:str, normal:str = None, normal_id:str = None, pon:str = None, targets:str = None, n_cpus=1) -> str:
    """Run DeepSomatic on the tumour/normal pair."""

    if not normal:
        model += "_TUMOR_ONLY"

    tool_comment = "# run deep somatic"
    tool_command = [
        tool, "run_deepsomatic",
        f"--model_type={model}",
        "--ref", reference,
        "--reads_tumor", tumour,
        "--sample_name_tumor", tumour_id,
        "--output_vcf", vcf,
        "--intermediate_results_dir", tmp_dir
    ]

    if n_cpus:
        tool_command.append(f"--num_shards {n_cpus}")
    else:
        tool_command.append(f"--num_shards 1")

    if targets:
        tool_command.append(f"--regions {targets}")

    if normal:
        tool_command.append(f"--reads_normal {normal} --sample_name_normal {normal_id}")

    if pon:
        tool_command.append(f"--pon_filtering {pon} --use_default_pon_filtering true")

    final_command = tool_comment + "\n" + " ".join(map(str, tool_command))
    return final_command


def run_filter_vcf(conda_path:str, vcf:str, filtered_vcf:str) -> str:
    """Filter raw DeepSomatic variant calls."""

    tool_comment = "# filter variants to only PASS"
    tool_command = [
        f'eval "$({conda_path} shell.bash hook)"',
        f"conda activate samtools_env",
        f"bcftools view -f PASS {vcf} > {filtered_vcf}",
        f"bgzip {filtered_vcf}",
        f"tabix -p vcf {filtered_vcf}.gz",
        "conda deactivate"
    ]

    final_command = tool_comment + "\n" + "\n".join(map(str, tool_command))
    return final_command


def deepsomatic_pipeline(tool_config:str, data_config:str, out_dir:str, remove:bool = False, dry_run:bool = False, wait_for:bool = False) -> str:

   # load tool/paramter information
    with open(tool_config, 'r') as file:
        tool_data = yaml.safe_load(file)

    # set up logging
    log_directory = os.path.join(out_dir, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    current_run = datetime.now().strftime("%Y-%m-%d-%H-%M")
    if dry_run:
        log_file = f"{log_directory}/deepsomatic_pipeline__dryrun.log"
    else:
        log_file = f"{log_directory}/deepsomatic_pipeline__{current_run}.log"
    
    log = open(log_file, 'w')

    print("---", file = log)
    print("Running the DeepSomatic variant-calling pipeline.", file = log)
    print(f"  Tool config used: {tool_config}", file = log)
    print(f"  Data config used: {data_config}", file = log)
    print(f"  Outputting to: {out_dir}", file = log)
    if remove:
        print("  >> clean mode selected (--remove); intermediate files will be removed upon completion", file = log)
    if dry_run:
        print("  >> trial mode selected (--dry_run); jobs will not be submitted", file = log)

    # extract tool paths from config
    deepsomatic_executable = tool_data['deepsomatic_path']

    # extract reference files from config
    reference = tool_data['reference']

    # add details to log file
    print(f"\nLinear reference file used is {reference}", file = log)
    print(f"\nUsing DeepSomatic version: {deepsomatic_executable}", file = log)
    print("Using samtools/bcftools version: 1.23.1", file = log)

    # load sample information
    with open(data_config, 'r') as file:
        sample_data = yaml.safe_load(file)

    n_patients = len(sample_data.keys())

    print(f"\n----------\nProcessing {n_patients} patients...", file = log)

    # is this pipeline waiting for any prior jobs?
    prior_jobs = []
    if wait_for:
        prior_jobs = wait_for.split(",")
        
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
                print(f"<<Error: Normal BAM file {normal_bam} not found. Skipping {donor}...>>", file = log)
                continue

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
            deepsomatic_vcf = os.path.join(smp_directory, f"{sample}_DeepSomatic_raw.vcf.gz")
            filtered_vcf = os.path.join(smp_directory, f"{sample}_DeepSomatic_pass_variants.vcf")

            ## Step 1: DeepSomatic
            deepsomatic_command, filter_command = None, None

            # run variant caller
            if is_file_missing_or_empty(deepsomatic_vcf):
                deepsomatic_command = run_deepsomatic(
                    tool = deepsomatic_executable,
                    model = tool_data['deepsomatic']['model'],
                    reference = reference,
                    targets = tool_data['target_regions'],
                    pon = tool_data['pon'],
                    tumour = tumour_bam,
                    tumour_id = sample, 
                    normal = normal_bam, 
                    normal_id = normal_id, 
                    vcf = deepsomatic_vcf,
                    tmp_dir = tmp_directory,
                    n_cpus = tool_data['deepsomatic']['cpus']
                )

            cleanup_files.append(deepsomatic_vcf)
            cleanup_files.append(f"{deepsomatic_vcf}.tbi")

            # filter variant calls
            if is_file_missing_or_empty(filtered_vcf):
                filter_command = run_filter_vcf(
                    conda_path = tool_data['conda_path'],
                    vcf = deepsomatic_vcf, 
                    filtered_vcf = filtered_vcf
                )
            
            call_commands = list(filter(None, [deepsomatic_command, filter_command]))

            # write slurm script and sbatch (unless dry-run)
            print(f"    Submitting job for DeepSomatic...", file = log)

            if len(call_commands) > 0 and is_file_missing_or_empty(f"{filtered_vcf}.gz"):

                ds_job_id, ds_script = write_script(
                    args,
                    job_name = f"run_deepsomatic_{sample}",
                    command = call_commands,
                    depends = prior_jobs,
                    job_parameters = tool_data['deepsomatic']
                )

                patient_jobs.append(ds_job_id)

                print(f"      >> sbatch {ds_script}", file = log)
                if not dry_run:
                    print(f"      >> job id is {ds_job_id}", file = log)

            else:
                print(f"      >> DeepSomatic step already completed; no jobs to submit.", file = log)

            final_files.append(f"{filtered_vcf}.gz")

        # do cleanup once per donor
        if remove:
            cleanup_command = "# check for final VCF(s), then cleanup intermediates"
            cleanup_command += f"\nif [[ -s " + " && -s ".join(final_files) + " ]]; then"
            cleanup_command += f"\n  rm -rf {tmp_directory}"
            cleanup_command += "\n  rm " + "\n  rm ".join(cleanup_files)
            cleanup_command += "\nelse"
            cleanup_command += "\n  echo 'Final output VCF(s) is missing; not removing intermediate files!'"
            cleanup_command += "\n  exit 1"
            cleanup_command += "\nfi"
        
            print(f"\n  Running clean-up step for {donor}...", file = log)

            if len(patient_jobs) > 0:
                clean_job_id, clean_script = write_script(
                    args,
                    job_name = f"run_cleanup_{donor}",
                    command = [cleanup_command],
                    job_parameters = {"mem": "256M", "time": "1:0:0"},
                    depends = patient_jobs
                )

                patient_jobs.append(clean_job_id)

                print(f"    >> sbatch {clean_script}", file = log)
                if not dry_run:
                    print(f"    >> job id is {clean_job_id}", file = log)

        all_jobs += patient_jobs

    # write slurm script and sbatch (unless dry-run)
    if not dry_run:
        sacct_out = f"{out_dir}/logs/slurm_job_metrics_{current_run}.out"
        sacct_command = collect_job_stats(all_jobs, sacct_out)
        sacct_job_id, sacct_script = write_script(
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
        print(f"  >> Final job id is: {sacct_job_id}", file = log)

    log.close()

    return (sacct_job_id)


### main
if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description = "Run DeepSomatic on T/N BAMs")
    parser.add_argument("-t", "--tool_config", help = "tool config (reference files and tool versions)")
    parser.add_argument("-d", "--data_config", help = "data config (bams listed in yaml format)")
    parser.add_argument("-o", "--out_dir", help = "path to output directory")
    parser.add_argument("--remove", action = "store_true", help = "should intermediates be removed?")
    parser.add_argument("--dry_run", action = "store_true", help = "is this a test?")
    parser.add_argument("--wait_for", help = "prior job ids (as a comma-separated list) that this function is dependent on")
    args = parser.parse_args()

    deepsomatic_pipeline(
        tool_config = args.tool_config, 
        data_config = args.data_config, 
        out_dir = args.out_dir,
        remove = args.remove, 
        dry_run = args.dry_run,
        wait_for = args.wait_for
    )
