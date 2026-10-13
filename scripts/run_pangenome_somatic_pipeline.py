"""This pipeline will submit jobs for the Pangenome:Somatic pipeline."""

import os
import yaml
import argparse

from datetime import datetime
from run_alignments import graph_alignment_pipeline
from run_mutect2 import mutect2_pipeline
from run_strelka import strelka_pipeline
from run_somaticsniper import somaticsniper_pipeline
from run_deepsomatic import deepsomatic_pipeline

SCRIPT_PATH = os.path.join(os.path.abspath(__file__), 'scripts')

### subroutines
def main(tool_config:str, data_config:str, out_dir:str, mode:str = 'slurm', run_alignment:bool = False, run_mutect2:bool = False, run_strelka:bool = False, run_somaticsniper:bool = False, run_deepsomatic:bool = False, remove:bool = False, dry_run:bool = False) -> None:

    # load tool/paramter information
    #with open(tool_config, 'r') as file:
        #tool_data = yaml.safe_load(file)

    # get config inputs
    #project_id = tool_data['project_name']
    #out_dir = tool_data['project_directory']

    # set up logging
    log_directory = os.path.join(out_dir, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    current_run = datetime.now().strftime("%Y-%m-%d-%H-%M")
    if dry_run:
        log_file = f"{log_directory}/pangenome-somatic-pipeline__dryrun.log"
    else:
        log_file = f"{log_directory}/pangenome-somatic-pipeline__{current_run}.log"
    
    with open(log_file, 'w') as log:
        log.write("---\n")
        log.write("Running the pangenome-somatic-pipeline.\n")
        log.write(f"  Tool config used: {tool_config}\n")
        log.write(f"  Outputting to: {out_dir}\n")
        if remove:
            log.write("  >> clean mode selected (--remove); intermediate files will be removed upon completion\n")
        if mode == 'local':
            log.write("  >> local mode selected (--mode local); jobs will be run locally if able\n")
        elif mode == 'slurm':
            log.write("  >> slurm mode selected (--mode slurm); jobs will be submitted using Slurm\n")
        if dry_run:
            log.write("  >> trial mode selected (--dry_run); jobs will not be submitted\n")


    # indicate potential directories
    alignment_directory = os.path.join(out_dir, 'alignments')
    mutect2_directory = os.path.join(out_dir, 'Mutect2')
    strelka_directory = os.path.join(out_dir, 'Strelka')
    somsnip_directory = os.path.join(out_dir, 'SomaticSniper')
    deepsom_directory = os.path.join(out_dir, 'DeepSomatic')


    # initiate object to hold job_ids
    alignment_job_id = None
    mutect2_job_id, strelka_job_id, somaticsniper_job_id, deepsomatic_job_id = None, None, None, None

    all_jobs = []

    with open(log_file, 'a') as log:
        log.write("\n----------\n\n")

   
        ### ALIGNMENT STEP ###
        log.write("Checking if alignment pipeline should be run...\n")

        if run_alignment:

            log.write("\n  Submitting alignment pipeline...\n")

            # run alignment pipeline
            bam_output = os.path.join(alignment_directory, 'processed_bams.yaml')
            log.write(f"    Using FASTQs in {data_config}\n")
            log.write(f"    Output directory: {alignment_directory}\n")
            log.write(f"    Output BAM config: {bam_output}\n")

            alignment_job_id = graph_alignment_pipeline(
                tool_config = tool_config,
                data_config = data_config,
                out_dir = alignment_directory,
                out_bams_filename = bam_output,
                mode = mode,
                remove = remove,
                dry_run = dry_run
            )

            all_jobs.append(alignment_job_id)

        else:
            log.write(f"\n  Alignment step not selected; will assume {data_config} lists BAMs instead of FASTQs...\n")
            bam_output = data_config

        log.write("\n----------\n\n")

        ### VARIANT CALLING STEPS ###
        log.write("Checking which variant-calling tool(s) should be run...\n")

        if run_mutect2:

            log.write("\n  Submitting Mutect2 pipeline...\n")

            # run mutect2 pipeline
            log.write(f"    Using BAMs in {bam_output}\n")
            log.write(f"    Output directory: {mutect2_directory}\n")

            mutect2_job_id = mutect2_pipeline(
                tool_config = tool_config,
                data_config = bam_output,
                out_dir = mutect2_directory,
                mode = mode,
                wait_for = alignment_job_id,
                remove = remove,
                dry_run = dry_run
            )

            all_jobs.append(mutect2_job_id)

        if run_strelka:

            log.write("\n  Submitting Strelka pipeline...\n")

            # run strelka pipeline
            log.write(f"    Using BAMs in {bam_output}\n")
            log.write(f"    Output directory: {strelka_directory}\n")

            strelka_job_id = strelka_pipeline(
                tool_config = tool_config,
                data_config = bam_output,
                out_dir = strelka_directory,
                mode = mode,
                wait_for = alignment_job_id,
                remove = remove,
                dry_run = dry_run
            )

            all_jobs.append(strelka_job_id)

        if run_somaticsniper:

            log.write("\n  Submitting SomaticSniper pipeline...\n")

            # run somaticsniper pipeline
            log.write(f"    Using BAMs in {bam_output}\n")
            log.write(f"    Output directory: {somsnip_directory}\n")

            somaticsniper_job_id = somaticsniper_pipeline(
                tool_config = tool_config,
                data_config = bam_output,
                out_dir = somsnip_directory,
                mode = mode,
                wait_for = alignment_job_id,
                remove = remove,
                dry_run = dry_run
            )

            all_jobs.append(somaticsniper_job_id)

        if run_deepsomatic:

            log.write("\n  Submitting DeepSomatic pipeline...\n")

            # run deepsomatic pipeline
            log.write(f"    Using BAMs in {bam_output}\n")
            log.write(f"    Output directory: {deepsom_directory}\n")

            deepsomatic_job_id = deepsomatic_pipeline(
                tool_config = tool_config,
                data_config = bam_output,
                out_dir = deepsom_directory,
                mode = mode,
                wait_for = alignment_job_id,
                remove = remove,
                dry_run = dry_run
            )

            all_jobs.append(deepsomatic_job_id)
    
        if not dry_run:
            final_job_ids = ",".join(map(str, all_jobs))
            log.write(f"  >> Final job id(s): {final_job_ids}\n")


### main
if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description = "Run Pangenome:Somatic pipeline")
    parser.add_argument("-t", "--tool_config", help = "tool config (reference files, input files, etc.)")
    parser.add_argument("-d", "--data_config", help = "data config (FASTQs or pre-processed BAMs, listed in yaml format)")
    parser.add_argument("-o", "--out_dir", help = "path to project directory")
    parser.add_argument("-m", "--mode", help = "run mode (local or slurm); defaults to slurm", default = 'slurm')
    parser.add_argument("--alignment", action = "store_true", help = "run the alignment step (input is fastqs)")
    parser.add_argument("--mutect2", action = "store_true", help = "run mutect2 variant calling")
    parser.add_argument("--strelka", action = "store_true", help = "run strelka variant calling")
    parser.add_argument("--somaticsniper", action = "store_true", help = "run somaticsniper variant calling")
    parser.add_argument("--deepsomatic", action = "store_true", help = "run deepsomatic variant calling")
    parser.add_argument("--remove", action = "store_true", help = "should intermediates be removed?")
    parser.add_argument("--dry_run", action = "store_true", help = "is this a test?")
    args = parser.parse_args()

    main(
        tool_config = args.tool_config,
        data_config = args.data_config,
        out_dir = args.out_dir,
        mode = args.mode,
        run_alignment = args.alignment,
        run_mutect2 = args.mutect2,
        run_strelka = args.strelka,
        run_somaticsniper = args.somaticsniper,
        run_deepsomatic = args.deepsomatic,
        remove = args.remove,
        dry_run = args.dry_run
    )
