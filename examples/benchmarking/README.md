# Pangenome-somaticSNV-pipeline

## Benchmarking Data
HG008 is a matched tumour-normal dataset developed by the Genome in a Bottle Consortium[1][1]: HG008-T is a pancreatic ductal adenocarcinoma cell line, while HG008-N-D is matched normal duodenal tissue. Both samples have been sequenced using multiple technologies and platforms, making it ideal for use in benchmarking studies.

Sequencing data for HG008-T and HG008-N-D were downloaded from GIAB (WGS, in BAM format), along with a set of high-confidence somatic SNVs:

* https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam
* https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam
* https://42basepairs.com/download/web/giab/data_somatic/HG008/Liss_lab/analysis/Hartwig_Oncoanalyzer_BCM-ILMN_20260204/pave/HG008-T.sage.somatic.pave.vcf.gz

For consistency with out current study utilising whole-exome sequencing from TCGA, the above BAMs were filtered down to a set of targeted regions, focusing on chromosome 17, and subsampled to 50% coverage to reduce filesize. Finally, reads from these chr17-subsetted-and-subsampled BAM files were extracted back to FASTQ format.

<pre><code># example commands
samtools view -b -P -L [chr17_targeted_regions.bed.gz][bed] --subsample 0.5 HG008-T_Illumina_195x_GRCh38-GIABv3.bam \
    | samtools collate -O - \
    | samtools fastq -1 HG008T_R1.fastq.gz -2 HG008T_R2.fastq.gz -0 /dev/null -s /dev/null -n

samtools view -b -P -L [chr17_targeted_regions.bed.gz][bed] --subsample 0.5 HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam \
    | samtools collate -O - \
    | samtools fastq -1 HG008N_R1.fastq.gz -2 HG008N_R2.fastq.gz -0 /dev/null -s /dev/null -n
</code></pre>


## Reference Files
Download your desired pangenome reference [here](https://data.humanpangenome.org/alignments), including the .gbz, .dist and .min files.

* https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.d9.dist
* https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.d9.gbz
* https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.d9.min

Extract the linear-reference genome as included in the above pangenome:
<pre><code># note that this will result in non-traditional contig ordering
vg paths paths --extract-fasta -x hprc-v1.1-mc-grch38.d9.gbz --paths-by GRCh38 > grch38_linear.fa
samtools dict -o grch38_linear.dict grch38_linear.fa
sed 's/GRCh38#0#chr/chr/g' grch38_linear.fa > grch38_linear_fixed.fa
samtools faidx grch38_linear_fixed.fa
</code></pre>

Download the required germline-resource files:
    * https://storage.cloud.google.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
    * https://storage.cloud.google.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
    * https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
    * https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
    * https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
    * https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi


## Running the pipeline
To run the pangenome-somaticSNV pipeline using the above input files:
1. Fill in the paths to your FASTQ files in the provided [data_config.yaml][data_config]
2. Fill in the paths to your reference files in the provided [tool_config.yaml][tool_config]
3. Download the prepacked docker/apptainer image [here][docker]
    * then fill in the container path/call command in [tool_config.yaml][tool_config]
4. Run the pipeline!
    * if running locally, replace 'slurm' with 'local' in the below command

<pre><code>cd /path/to/git/pangenome-pipeline

# test run the pipeline; this will generate the output directory structure and \
#     scripts that will be submitted to slurm
# check that all of these seem reasonable before submitting the next command
python scripts/run_pangenome_somatic_pipeline.py \
    -t /path/to/your/tool_config.yaml \
    -d /path/to/your/data_config.yaml \
    -o /path/to/output \
    --alignment --mutect2 --strelka -m slurm --remove --dry_run

# now, submit the jobs to slurm!
python scripts/run_pangenome_somatic_pipeline.py \
    -t /path/to/your/tool_config.yaml \
    -d /path/to/your/data_config.yaml \
    -o /path/to/output \
    --alignment --mutect2 --strelka -m slurm --remove
</code></pre>

The expected runtimes and memory usage can be found in [expected_runtimes.txt][expected_resources].


## Benchmarking Results
A set of 'true positive' variants were identified from GIAB's high-confidence somatic variant callset by extracting SNV's with the PASS filter on chr17.

These variants, along with the expected variant allele frequency (VAF) and the observed VAFs from our pipeline for Strelka2 and Mutect2 can be found in [expected_results.tsv][expected_results]


[1]: doi.org/10.1038/s41597-025-05438-2
[bed]: chr17_targeted_regions.bed.gz
[data_config]: data_config.yaml
[tool_config]: tool_config.yaml
[docker]: /path/to/dockerhub
[expected_resources]: expected_runtimes.txt
[expected_results]: expected_results.tsv
