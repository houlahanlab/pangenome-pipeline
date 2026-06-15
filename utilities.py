import os

from simple_slurm import Slurm

def export_slurm_script(slurm_obj, filename: str):
    """Safely extracts the script string, unescapes all heredoc modifications,
    and writes a clean, runnable copy to disk.
    """
    script_str = str(slurm_obj)
    
    # Need to unescape backslashes first.
    # Otherwise,  escaping on the other characters will be corrupted.
    script_str = script_str.replace("\\\\", "\\")
    script_str = script_str.replace("\\`", "`")
    script_str = script_str.replace("\\$", "$")
    script_str = script_str.replace("# ", "\n# ")
    
    with open(filename, "w") as f:
        f.write(script_str)


def get_java_mem(mem):
    """Format total memory requested into a java-appropriate argument."""

    if "G" in mem:
        less_mem = int(mem.replace("G","")) - 1
        java_mem = f"-Xmx{less_mem}g"
    else:
        java_mem = ""

    return java_mem


def write_script(args, job_name, command, job_parameters=[], depends=None, array=None):
    """Write slurm commands to file and submit jobs."""

    job_dir = f"{args.out_dir}/logs/{job_name}"
    script = f"{job_dir}/script.sh"
    os.makedirs(job_dir, exist_ok=True)

    job_mem = '1G'
    if 'mem' in job_parameters:
        job_mem = job_parameters['mem']

    job_time = '1:0:0'
    if 'time' in job_parameters:
        job_time = job_parameters['time']
 
    job_cpus = '1'
    if 'cpus' in job_parameters:
        job_cpus = job_parameters['cpus']
   
    job_queue = 'general'
    if 'G' in job_mem:
        memory = int(job_mem.replace('G',''))
        if memory > 100:
            job_queue = 'massive'

    slurm = Slurm(
        job_name = job_name,
        chdir = job_dir,
        time = job_time,
        mem = job_mem,
        cpus_per_task = job_cpus,
        partition = job_queue
    )

    if depends is not None:
        clean_depends = [x for x in depends if x]
        if clean_depends:
            str_depends = ':'.join(map(str, clean_depends))
            if job_name == "collect_job_metrics":
                slurm.set_dependency(f"afterany:{str_depends}")
            else:
                slurm.set_dependency(f"afterok:{str_depends}")
            slurm.add_arguments(kill_on_invalid_dep="yes")

    if array is not None:
        slurm.set_array(array)

    slurm.add_cmd("echo $HOSTNAME")
    slurm.add_cmd("echo $SLURM_NODELIST") # $SLURMD_NODENAME?
    slurm.add_cmd("set -e")

    for subcommand in command:
        slurm.add_cmd(subcommand)

    export_slurm_script(slurm, script)

    job_id = job_name
    if not args.dry_run:
        job_id = slurm.sbatch('' + Slurm.SLURM_ARRAY_TASK_ID);

    return job_id,script


def collect_job_stats(job_ids, outfile):
    """Collect slurm job metrics and convert to user-friendly table."""

    sacct_command = [
        'sacct -P --delimiter=","',
        '--format="User,JobID,JobName,Start,End,AllocCPUS,CPUTime,MaxRSS,State,ExitCode"',
        "-j", ",".join(map(str, job_ids)),
        '| sed "s/,/\\t/g"',
        ">", outfile + ";"
    ]

    final_command = " ".join(sacct_command)
    return final_command


def is_file_missing_or_empty(file: str) -> bool:
    """Determine whether a file already exists and/or is empty."""

    if not os.path.isfile(file):
        return True
    return os.path.getsize(file) == 0
