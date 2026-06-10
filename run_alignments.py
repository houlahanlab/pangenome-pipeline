"""This pipeline will perform graph-based alignment from fastq files and output a
GATK-recalibrated BAM that can be used for variant-calling."""

import os
import subprocess
import yaml
import argparse

from datetime import datetime
from simple_slurm import Slurm

from utilities import write_script,collect_job_stats,is_file_missing_or_empty,get_java_mem

### subroutines
def verify_ref(gbz, dist_file=None, min_file=None):
    """Verify that all necessary reference files are available."""

    if dist_file is None:
        dist_file = gbz.with_suffix(".dist")

    if min_file is None:
        min_file = gbz.with_suffix(".min")

    if not os.path.exists(gbz):
        print(f"Error: gbz file not found.")

    if not os.path.exists(dist_file):
        print(f"Error: dist file not found.")

    if not os.path.exists(min_file):
        print(f"Error: min file not found.")
    
    return gbz, dist_file, min_file


def run_giraffe(vg, reference, dist_file, min_file, output_gam, fastq_files):
    """Perform graph-based alignment using VG Giraffe."""

    vg_comment = "# run alignment using vg giraffe"
    vg_command = [
        vg, "giraffe",
        "-p -t 16 -D 150",
        "-Z", reference,
        "-d", dist_file,
        "-m", min_file,
        "-f", fastq_files[0],
        "-f", fastq_files[1],
        ">", output_gam
    ]

    final_command = vg_comment + "\n" + " ".join(map(str, vg_command))
    return final_command


def run_vg_stats(vg, gam, stats):
    """Collect alignment stats from graph-aligned reads."""

    vg_comment = "# collect alignment statistics from graph"
    vg_command = [
        vg, "stats",
        "-a", gam,
        ">", stats
    ]

    final_command = vg_comment + "\n" + " ".join(map(str, vg_command))
    return final_command


def run_surject(vg, reference, gam, paths, bam):
    """Surject graph-based alignments to linear alignments (gam to bam)."""

    vg_comment = "# convert gam to bam"
    vg_command = [
        vg, "surject",
        "-x", reference,
        "-i -b", gam,
        "-F", paths,
        ">", bam
    ]

    final_command = vg_comment + "\n" + " ".join(map(str, vg_command))
    return final_command


def run_gatk_sort(reference, in_bam, out_bam, tmp_dir, lib, donor, sample, java_mem):
    """Add readgroups and output coordinate-sorted BAM."""

    java = f"--java-options '{java_mem} -XX:ParallelGCThreads=2 -XX:ConcGCThreads=1'"

    gatk_comment = "# sort and add readgroups to bam header"
    gatk_command = [
        "gatk", java, "AddOrReplaceReadGroups",
        "-R", reference,
        "-I", in_bam,
        "-O", out_bam,
        "--RGLB", lib,
        "--RGPU", donor,
        "--RGSM", sample,
        "--RGPL Illumina --SORT_ORDER coordinate --CREATE_INDEX",
        "--TMP_DIR", tmp_dir
    ]

    final_command = gatk_comment + "\n" + " ".join(map(str, gatk_command))
    return final_command


def run_gatk_markdup(in_bam, out_bam, metrics, tmp_dir, java_mem):
    """Mark duplicate reads and output an indexed BAM."""

    java = f"--java-options '{java_mem} -XX:ParallelGCThreads=2 -XX:ConcGCThreads=1'"

    gatk_comment = "# mark duplicate reads"
    gatk_command = [
        "gatk", java, "MarkDuplicates",
        "-I", in_bam,
        "-O", out_bam,
        "-M", metrics,
        "--ASSUME_SORT_ORDER coordinate --CREATE_INDEX",
        "--TMP_DIR", tmp_dir
    ]

    final_command = gatk_comment + "\n" + " ".join(map(str, gatk_command))
    return final_command


def run_gatk_bqsr(reference, in_bam, out_table, indels, out_bam, tmp_dir, java_mem):
    """Perform GATK's base-quality-score-recalibration."""

    java = f"--java-options '{java_mem} -XX:ParallelGCThreads=2 -XX:ConcGCThreads=1'"

    gatk_comment1 = "# create recalibration table"
    gatk_command1 = [
        "gatk", java, "BaseRecalibrator",
        "-R", reference,
        "-I", in_bam,
        "-O", out_table,
        "--known-sites", indels,
        "--tmp-dir", tmp_dir
    ]

    gatk_comment2 = "# perform recalibration"
    gatk_command2 = [
        "gatk", java, "ApplyBQSR",
        "-R", reference,
        "-I", in_bam,
        "-O", out_bam,
        "--bqsr-recal-file", out_table,
        "--tmp-dir", tmp_dir
    ]

    final_command = gatk_comment1 + "\n" + " ".join(map(str, gatk_command1)) 
    final_command += "\n" + gatk_comment2 + "\n" + " ".join(map(str, gatk_command2))
    return final_command


def run_gatk_qc(bam, metrics, tmp_dir):
    """Collect summary alignment metrics from recalibrated BAM."""

    gatk_comment = "# collect alignment summary metrics"
    gatk_command = [
        "gatk CollectAlignmentSummaryMetrics",
        "-I", bam,
        "-O", metrics,
        "--TMP_DIR", tmp_dir
    ]

    final_command = gatk_comment + "\n" + " ".join(map(str, gatk_command))
    return final_command


def main(tool_config:str, data_config:str, out_dir:str, out_bams_filename:str, remove:bool = False, dry_run:bool = False):

    # set up logging
    log_directory = os.path.join(out_dir, 'logs')
    os.makedirs(log_directory, exist_ok=True)

    current_run = datetime.now().strftime("%Y-%m-%d-%H-%M")
    if dry_run:
        log_file = f"{log_directory}/graph_alignment_pipeline__dryrun.log"
    else:
        log_file = f"{log_directory}/graph_alignment_pipeline__{current_run}.log"
    
    log = open(log_file, 'w')

    print("---", file = log)
    print("Running the VG/GATK graph-alignment and recalibration pipeline.", file = log)
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

    vg_executable = tool_data['vg_path']
    reference = tool_data['reference']

    gbz_file, dist_file, min_file = verify_ref(
        tool_data['ref_gbz'],
        tool_data['ref_dist'],
        tool_data['ref_min']
        )

    # load sample information
    with open(data_config, 'r') as file:
        sample_data = yaml.safe_load(file)

    n_patients = len(sample_data.keys())

    # add details to log file
    print(f"\nGraph-based reference file used is {gbz_file}", file = log)
    print(f"Linear reference file used is {reference}", file = log)
    print(f"\nUsing vg version: {vg_executable}", file = log)
    print("Using gatk version: 4.3.0.0 [with Picard version 2.27.5]", file = log)
    print(f"\n----------\nProcessing {n_patients} patients...", file = log)

    # initiate object to hold job_ids
    all_jobs = []

    # initiate object to hold final output bams
    out_bams = {}

    # work through each sample
    for donor in sample_data.keys():

        patient_directory = os.path.join(out_dir, donor)
        tmp_directory = os.path.join(patient_directory, 'TEMP')
        os.makedirs(tmp_directory, exist_ok=True)

        patient_jobs = []
        final_bams = []
        cleanup_files = []

        n_samples = len(sample_data[donor].keys())
        print(f"\nProcessing {n_samples} for donor: {donor}", file = log)

        for sample in (sample_data[donor].keys()):
 
            smp_directory = os.path.join(patient_directory, sample)
            os.makedirs(smp_directory, exist_ok=True)

            tissue = sample_data[donor][sample]['type']

            print(f"  Processing {sample}", file = log)

            for library in (sample_data[donor][sample]['libraries'].keys()):
            
                for lane in (sample_data[donor][sample]['libraries'][library]['runlanes'].keys()):

                    # build in command to concatenate fastqs if multiple lanes/libraries are present
                    #
                    #

                    fastq_files = [
                        sample_data[donor][sample]['libraries'][library]['runlanes'][lane]['fastq']['R1'],
                        sample_data[donor][sample]['libraries'][library]['runlanes'][lane]['fastq']['R2']
                    ]

            # process all fastqs for each sample
            if not os.path.exists(fastq_files[0]):
                print(f"<<Error: FASTQ file {fastq_files[0]} not found. Skipping {sample}...>>", file = log)
                continue
            if not os.path.exists(fastq_files[1]):
                print(f"<<Error: FASTQ file {fastq_files[1]} not found. Skipping {sample}...>>", file = log)
                continue

            # indicate output files
            output_gam = os.path.join(smp_directory, f"{sample}_vg_aligned.gam")
            stats_gam = os.path.join(smp_directory, f"{sample}_vg_aligned_stats.txt")
            surject_bam = os.path.join(smp_directory, f"{sample}_vg_aligned_surjected.bam")
            sorted_bam = os.path.join(smp_directory, f"{sample}_vg_aligned_surjected_sorted.bam")
            markdup_bam = os.path.join(smp_directory, f"{sample}_vg_aligned_surjected_sorted_markdup.bam")
            markdup_stats = os.path.join(smp_directory, f"{sample}_vg_aligned_surjected_sorted_markdup_metrics.txt")
            recal_table = os.path.join(smp_directory, f"{sample}_recal.table")
            recal_bam = os.path.join(smp_directory, f"{sample}_recalibrated.bam")
            bam_qc = os.path.join(smp_directory, f"{sample}_final_alignment_metrics.txt")

            cleanup_files.append(surject_bam)
            cleanup_files.append(sorted_bam)
            cleanup_files.append(sorted_bam.replace('bam','bai'))
            cleanup_files.append(markdup_bam)
            cleanup_files.append(markdup_bam.replace('bam','bai'))

            job_id = []

            ## Step 1: VG giraffe + surject
            giraffe_command, stats_command, surject_command = None, None, None 

            # run alignment
            if is_file_missing_or_empty(output_gam):
                giraffe_command = run_giraffe(
                    vg_executable,
                    gbz_file,
                    dist_file,
                    min_file,
                    output_gam,
                    fastq_files
                    )

            # collect alignment stats
            if is_file_missing_or_empty(stats_gam):
                stats_command = run_vg_stats(
                    vg_executable,
                    output_gam,
                    stats_gam
                    )

            # convert gam to bam
            if is_file_missing_or_empty(surject_bam):
                surject_command = run_surject(
                    vg_executable,
                    gbz_file,
                    output_gam,
                    tool_data['ref_paths'],
                    surject_bam
                    )

            align_commands = list(filter(None, [giraffe_command, stats_command, surject_command]))

            # write slurm script and sbatch (unless dry-run)
            print(f"    Submitting job for alignment step...", file = log)

            if len(align_commands) > 0 and is_file_missing_or_empty(recal_bam):        

                job_id, script = write_script(
                    args,
                    job_name = f"run_alignment_{sample}",
                    command = align_commands,
                    job_parameters = tool_data['giraffe']
                    )

                patient_jobs.append(job_id)

                print(f"      >> sbatch {script}", file = log)
                if not dry_run:
                    print(f"      >> job id is {job_id}", file = log)

            else:
                print(f"      >> Alignment step already completed; no jobs to submit.", file = log)


            ## Step 2: GATK markdups and recalibration
            sort_command, markdup_command, bqsr_command = None, None, None
            java_mem = get_java_mem(tool_data['gatk']['mem'])

            # update readgroup info and sort bam
            if is_file_missing_or_empty(sorted_bam):
                sort_command = run_gatk_sort(
                    reference,
                    in_bam = surject_bam,
                    out_bam = sorted_bam,
                    tmp_dir = tmp_directory,
                    lib = ";".join(sample_data[donor][sample]['libraries'].keys()),
                    donor = donor,
                    sample = sample,
                    java_mem = java_mem
                    )

            # mark duplicates
            if is_file_missing_or_empty(markdup_bam):
                markdup_command = run_gatk_markdup(
                    in_bam = sorted_bam,
                    out_bam = markdup_bam,
                    metrics = markdup_stats,
                    tmp_dir = tmp_directory,
                    java_mem = java_mem
                    )

            # perform base-quality score recalibration
            if is_file_missing_or_empty(recal_bam):
                bqsr_command = run_gatk_bqsr(
                    reference,
                    in_bam = markdup_bam,
                    out_table = recal_table,
                    indels = tool_data['known_indels'],
                    out_bam = recal_bam,
                    tmp_dir = tmp_directory,
                    java_mem = java_mem
                    )

            gatk_commands = list(filter(None, [sort_command, markdup_command, bqsr_command]))

            # write slurm script and sbatch (unless dry-run)
            print(f"    Submitting job for postprocessing step...", file = log)

            if len(gatk_commands) > 0 and is_file_missing_or_empty(recal_bam):

                job_id, script = write_script(
                    args,
                    job_name = f"run_postprocessing_{sample}",
                    command = gatk_commands,
                    job_parameters = tool_data['gatk'],
                    depends = [job_id]
                    )

                patient_jobs.append(job_id)

                print(f"      >> sbatch {script}", file = log)
                if not dry_run:
                    print(f"      >> job id is {job_id}", file = log)

            else:
                print(f"      >> Postprocessing step already completed; no jobs to submit.", file = log)

            # add final BAM to yaml file
            if donor not in out_bams:
                out_bams[donor] = {}
            if tissue not in out_bams[donor]:
                out_bams[donor][tissue] = {}

            out_bams[donor][tissue][sample] = recal_bam
            final_bams.append(recal_bam)

            ## Step 3: BAM QC
            qc_command = None

            # collect metrics for final bam
            if is_file_missing_or_empty(bam_qc):
                qc_command = run_gatk_qc(
                    recal_bam,
                    bam_qc,
                    tmp_directory
                    )

            # write slurm script and sbatch (unless dry-run)
            print(f"    Submitting job for QC step...", file = log)

            if qc_command is not None:

                job_id, script = write_script(
                    args,
                    job_name = f"run_bamqc_{sample}",
                    command = [qc_command],
                    job_parameters = {"mem": "4G", "time": "8:0:0"},
                    depends = [job_id]
                    )

                patient_jobs.append(job_id)

                print(f"      >> sbatch {script}", file = log)
                if not dry_run:
                    print(f"      >> job id is {job_id}", file = log)

            else:
                print(f"      >> QC step already completed; no jobs to submit.", file = log)


        # do cleanup once per donor
        if remove:
            cleanup_command = f"if [[ -s " + " && -s ".join(final_bams) + " ]]; then\n"
            cleanup_command += f"  rm -rf {tmp_directory}"
            cleanup_command += "\n  rm " + "\n  rm ".join(cleanup_files)
            cleanup_command += "else"
            cleanup_command += "\n  echo 'Final output bam is missing; not removing intermediate files!'"
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

            else:
                subprocess.run(f"rm -rf {tmp_directory}")

        all_jobs += patient_jobs

    # write final output yaml
    with open(out_bams_filename, 'w') as file:
        yaml.dump(out_bams, file, default_flow_style = False)
 
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
    print(f"YAML containing paths to recalibrated BAM files is: {out_bams_filename}", file = log)

    n_jobs = len(all_jobs)
    print(f"\nTotal number of jobs to run: {n_jobs}", file = log)

    if not dry_run:
        print(f"  >> Final job id is: {job_id}", file = log)

    log.close()


### main
if __name__ == "__main__":

    # get command line arguments
    parser = argparse.ArgumentParser(description = "Align fastq files")
    parser.add_argument("-t", "--tool_config", help = "tool config (reference files and tool versions)")
    parser.add_argument("-d", "--data_config", help = "data config (fastqs listed in yaml format)")
    parser.add_argument("-o", "--out_dir", help = "path to output directory")
    parser.add_argument("-b", "--out_bams_filename", help = "path to output yaml (listing aligned BAMs)")
    parser.add_argument("--remove", action = "store_true", help = "should intermediates be removed?")
    parser.add_argument("--dry_run", action = "store_true", help = "is this a test?")
    args = parser.parse_args()

    main(
        tool_config=args.tool_config, 
        data_config=args.data_config, 
        out_dir=args.out_dir, 
        out_bams_filename=args.out_bams_filename, 
        remove=args.remove, 
        dry_run=args.dry_run
    )

