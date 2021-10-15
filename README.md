# seclip
A pipeline for processing SE eCLIP data to identify genomic locations of RNA binding proteins (RBPs).

# Usage
From the terminal of taco cluster, directly call `/storage/vannostrand/software/eclip/seclip` 
with 3 required parameters along with their values to submit your job for processing:

```shell
$ /storage/vannostrand/software/eclip/seclip \
  --ip_fastqs ip_1_fastq ip_2_fastq ip_3_fastq \
  --input_fastqs input_1_fastq input_2_fastq input_3_fastq \
  --names S1 S2 S3 
```

If you want to fine tune your analysis, check out and add other optional parameters 
when calling the pipeline:
```shell
$ /storage/vannostrand/software/eclip/seclip -h

usage: eclip [-h] 
             --input_fastqs INPUT_FASTQS [INPUT_FASTQS ...] 
             --ip_fastqs IP_FASTQS [IP_FASTQS ...] 
             --names NAMES [NAMES ...] 
             --adapters_fasta ADAPTERS_FASTA 
             [--barcodes_pattern BARCODES_PATTERN]
             [--dataset DATASET] [--outdir OUTDIR] 
             [--genome GENOME] [--repeat REPEAT] 
             [--species SPECIES] [--blacklist_bed BLACKLIST_BED] 
             [--track TRACK] [--track_genome TRACK_GENOME] 
             [--l2fc L2FC] [--l10p L10P]
             [--enrichment_filter ENRICHMENT_FILTER] 
             [--job JOB] [--email EMAIL] [--time TIME] 
             [--memory MEMORY] [--cpus CPUS] 
             [--scheduler {pbs,qsub,slurm,sbatch}] 
             [--hold_submit] [--debug] [--dry_run]

A pipeline for processing eCLIP data to identify genomic locations of RNA binding proteins (RBPs).

optional arguments:
  -h, --help            show this help message and exit
  --input_fastqs INPUT_FASTQS [INPUT_FASTQS ...]
                        Path(s) to INPUT FASTQ file(s).
  --ip_fastqs IP_FASTQS [IP_FASTQS ...]
                        Path(s) to IP FASTQ file(s).
  --names NAMES [NAMES ...]
                        Shortnames for each sample, e.g., rep1, rep2.
  --adapters_fasta ADAPTERS_FASTA
                        Path to the fasta file contains adapters and their sequences.
  --barcodes_pattern BARCODES_PATTERN
                        Pattern of barcodes for umi-tools using to extract UMIs (for single-end dataset only, 
                        default: NNNNNNNNNN
  --dataset DATASET     Name of the eCLIP dataset, default: eCLIP.
  --outdir OUTDIR       Path to the output directory, an pre-existing directory or current directory.
  --genome GENOME       Path to STAR reference genome index directory.
  --repeat REPEAT       Path to STAR repeat elements index directory.
  --species SPECIES     Species name (short name code) the dataset associated with, e.g., hg19, mm10.
  --blacklist_bed BLACKLIST_BED
                        Path to the bed file contains blacklist.
  --track TRACK         Name for the UCSC Genome Browser track, default: eCLIP
  --track_genome TRACK_GENOME
                        Genome name (a short name code) for the track.
  --l2fc L2FC           Only consider peaks at or above this log2 fold change cutoff.
  --l10p L10P           Only consider peaks at or above this log10 p value cutoff.
  --enrichment_filter ENRICHMENT_FILTER
                        Pre-filter peaks that are enriched over input.
  --job JOB             Name of your job, default: eCLIP.
  --email EMAIL         Email address for notifying you the start, end, and abort of you job.
  --time TIME           Time (in integer hours) for running your job.
  --memory MEMORY       Amount of memory (in GB) for all CPU cores.
  --cpus CPUS           Maximum number of CPU cores can be used for your job.
  --scheduler {pbs,qsub,slurm,sbatch}
                        Name (case insensitive) of the scheduler on your cluster.
  --hold_submit         Generate the submit script but hold it without submitting to the job scheduler.
  --debug               Invoke debug mode (for development use only).
  --dry_run             Print out steps and files involved in each step without actually running the pipeline.
```

# Example
## Understand your dataset
Assume we have a SE eCLIP dataset contains the following raw fastq files for read 1:

| Sample |  IP FASTQ             |  INPUT FASTQ            |
|--------|-----------------------|-------------------------|
|S1     |/path/to/ip1/fastq.gz   |/path/to/input1/fastq.gz |
|S2     |/path/to/ip2/fastq.gz   |/path/to/input2/fastq.gz |
|S3     |/path/to/ip3/fastq.gz   |/path/to/input3/fastq.gz |

And we are interested in identifying RBPs based on human genome assembly hg19.

## Kick off analysis
Based on the structure of our SE eCLIP dataset and our purpose, 
we can issue the following command from a terminal on taco cluster:
```shell
/storage/vannostrand/software/eclip/seclip \
    --ip_fastqs /path/to/ip1/fastq.gz /path/to/ip2/fastq.gz /path/to/ip3/fastq.gz \
    --input_fastqs /path/to/input1/fastq.gz /path/to/input2/fastq.gz /path/to/input3/fastq.gz \
    --names S1 S2 S3 \
    --outdir /path/to/output/directory
```

Notice that we added one additional parameter `outdir`, which will 
save all outputs into a designed directory instead of current work
directory (this is highly recommended). 

If there is no error regarding paths to fastq  files, the above 
command should be succeeded with a message like this:
```
Job submit script was saved to:
    /path/to/output/directory/submit.sh
Job eCLIP was successfully submitted with the following resources:
           Job name: eCLIP
   Output directory: /path/to/output/directory
    Number of cores: 16
         Job memory: 32
        Job runtime: 1-12:00 (D-HH:MM)
```

This tells you your job was successfully submitted to job scheduler on 
taco cluster and waiting for resources to start processing.

When `cd` into `/path/to/output/directory`, you will notice that a file named 
`submit.sh` was created with details for kicking off the pipeline. Except the 
lines start with `#`, which are directives for taco cluster job scheduler, all 
other lines should be self-explanatory. Refer to the usage section or contact 
me if you have any questions.

```shell
#!/usr/bin/env bash

#SBATCH -n 16                       # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t 1-12:00                  # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=32G                   # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name=eCLIP            # Short name for the job

#SBATCH --mail-user=fei.yuan@bcm.edu
#SBATCH --mail-type=ALL

export TMPDIR=/storage/vannostrand/tmp
export TEMP=/storage/vannostrand/tmp
export TMP=/storage/vannostrand/tmp

/storage/vannostrand/software/eclip/seclip \
  --ip_fastqs /path/to/ip1/fastq.gz /path/to/ip2/fastq.gz /path/to/ip3/fastq.gz \
  --input_fastqs /path/to/input1/fastq.gz /path/to/input2/fastq.gz /path/to/input3/fastq.gz \
  --names IP1 IP2 IP3 \
  --adapters_fasta /storage/vannostrand/software/eclip/tests/seRBFOX2/se.adapters.fasta \
  --barcodes_pattern NNNNNNNNNN \
  --species hg19 \
  --repeat /storage/vannostrand/reference_data/hg19/repbase_v2_star_index \
  --genome /storage/vannostrand/reference_data/hg19/genome_star_index \
  --outdir /path/to/output/directory \
  --dataset eCLIP \
  --track eCLIP \
  --track_genome hg19 \
  --l2fc 3 \
  --l10p 3 \
  --cpus 16
```

## Monitor your analysis job
By default, `seclip` will submit your job with name `eCLIP` or a name passed to 
parameter `job`. By issue the following command (replace the username with your 
own username):
```shell
$ sequeue -u [username]
```
This will give you something like this:
```
 JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
    77     mhgcp   PhyloP      username  R 1-00:59:41      1 mhgcp-m00
   208     mhgcp    eCLIP      username  R    1:01:37      1 mhgcp-c01
```
The first column lists job IDs for all of your submitted jobs. You can ignore the 
second column here. The third column lists job names for all of your submitted 
jobs. You should see the job name of your analysis here before it complete or 
gets killed due to some error. The fifth column ST tells you the status of you job:
R for running, PD or pending or waiting for resources to be allocated. You may 
also see some other status other than R or PD, you can ignore them for a while 
and check it again in few minutes, if it stays the same and not turns to R or PD, 
please let me know.

## Modify your analysis job
The `submit.sh` script was automatically generated using info provided by you and some 
default settings. If any of the parameter does not look right, you can modify it 
accordingly before your job start (see Monitor your analysis job section for how 
to get the status of your job). 

In case your job started immediately after you issue the command, you can cancel the 
running job, deleted anything created in the output directory and re-submit the job by 
using the modified `submit.sh` script or re-run the command with correct parameters.

To cancel a job, find its job id and issue the following command:
```shell
$ scancel [job_id]
```

After the job has been cancelled (check all your jobs to make sure it has been 
cancelled), you can modify the `submit.sh` and re-submit it:
```shell
$ cd /path/to/output/directory
$ sbatch submit.sh
```

Or you can start over by modifying the command you used to kick off yur job, then 
the `submit.sh` will be overwritten and automatically submitted.

## Understand output
### Run log
By default, the pipeline will have taco cluster write run log into the output 
directory with name `slurm-*`, where the `*` will be your job id. If you kicked 
off the job multiple times, you will see multiple `slurm-*` log files. Each log 
file contains details of the processing and any error. Check the latest run log 
to make sure there is no error for the analysis before you head to the results.

### Results
You can expect outputs of the following file structure:
```
Extract UMI: *.umi.fastq.gz
Cut adapt: *.trim.fastq.gz
           *.trim.metrics
           *.trim.trim.metrics
           *.trim.trim.trim.metrics
Repeat element map: *.re.map.bam
                    *.re.map.log
                    *.re.unmap.fastq.gz
Reference genome map: *.genome.map.bam
                      *.genome.map.log
                      *.genome.unmap.fastq.gz
Deduped and position sorted BAM: *.bam
RPM-normalized BigWig files: *.plus.bw, .minus.bw
Track hub: hub.txt
Clipper peaks: *.peak.clusters.bed
pureClip peaks: *.crosslink.sites.bed, *.binding.regions.bed
Input normalize annotated peaks: *.peak.clusters.normalized.compressed.annotated.bed
Reproducible peaks: *.reproducible.peaks.bed
```