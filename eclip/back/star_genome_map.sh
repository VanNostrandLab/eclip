#!/usr/bin/env bash

set -e

function usage() {
      cat << HELP

Map reads to reference genome using STAR.

Usage: star_genome_map FASTQ <-x index> [-o bam] [-c cpus] [-h] [-q]

Required:
  FASTQ                     Path to FASTQ file.
  -x, --index               STAR reference genome index directory.

Options:
  -o, --output_bam          Path to output BAM file.
  -c, --cpus                Maximum number of CPUS can be used for mapping.
  -h, --help                Print out help message and exit.
	-q, --quiet              Sort fastq file quietly without printing out message.

HELP
  exit 1
}

(( "$#" )) || usage

while (( "$#" )); do
  case "$1" in
    -x|--index)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        index=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -o|--output_bam)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        output_bam=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -c|--cpus)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        cpus=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -q|--quiet)
			quiet=1
			shift 1
			;;
    -h|--help)
      usage
      ;;
    *)
      fastq="$1"
      shift 1
      ;;
  esac
done

fastq="${fastq:-""}"
index="${index:-""}"
output_bam="${output_bam:-""}"
cpus="${cpu:-1}"
quiet="${quiet:-0}"

IFS=, read -r fastq1 fastq2 <<< "${fastq}"
[[ -f "${fastq1}" ]] || { echo "Input FASTQ file ${fastq1} does not exist or is not a file!"; exit 1; }
[[ -d "${index}" ]] || { echo "Reference genome index ${index} does not exist or is not a directory!"; exit 1; }
if [[ "${output_bam}" == "" ]]
then
	case "${fastq1}" in
			*.fq)
					name="${fastq1%.fq}";;
			*.fastq)
					name="${fastq1%.fastq}";;
			*.fq.gz)
					name="${fastq1%.fq.gz}";;
			*.fastq.gz)
					name="${fastq1%.fastq.gz}";;
			*)
					name="${fastq1}";;
	esac
else
	[[ "${output_bam}" == *.bam ]] && name="${output_bam%.bam}" || name="${output_bam}"
fi
[[ "${fastq1}" == *.gz ]] && read_file_command="zcat" || read_file_command="-"

folder=$(mktemp -d "${name}".repbase.map.XXXXXXXXXX)
trap 'rm -rf ${folder}' EXIT

IFS='' read -r cmd << EOF
STAR \
--runMode alignReads \
--runThreadN ${cpus} \
--alignEndsType EndToEnd \
--genomeDir ${index} \
--outBAMcompression 10 \
--outFileNamePrefix ${folder}/ \
--outFilterMultimapNmax 1 \
--outFilterScoreMin 10 \
--outReadsUnmapped Fastx \
--outSAMattrRGline ID:foo \
--outSAMattributes All \
--outSAMtype BAM Unsorted \
--readFilesCommand ${read_file_command} \
--readFilesIn ${fastq1}
EOF

if [[ "${fastq2}" != "" ]]
then
	[[ -f "${fastq2}" ]] || { echo "Input FASTQ file ${fastq2} does not exist or is not a file!"; exit 1; }
	cmd="${cmd}  ${fastq2}"
	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Mapping paired-end reads $fastq1 and"
	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")]                          $fastq2 ..."
else
	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Mapping single-end reads $fastq1 ..."
fi
cmd="${cmd} > /dev/null"

echo "${cmd}" | sed 's/  / \\\n                /g' | sed 's/ -/ \\\n  -/g' | sed 's/ >/ \\\n  >/g'
eval "${cmd}"

mv "${folder}"/Aligned.out.bam "${name}.bam"
mv "${folder}"/Log.final.out "${name}.log"

mate1="${folder}"/Unmapped.out.mate1
pigz -c -p "${cpus}" "${mate1}" > "${name}.mate1.gz"
mate2="${folder}"/Unmapped.out.mate2
[[ -f "${mate2}" ]] && pigz -c -p "${cpus}" "${mate2}" > "${name}.mate2.gz"

if [[ "${fastq2}" != "" ]]
then
	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Mapping paired-end reads $fastq1 and"
	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")]                          $fastq2 complete."
else
	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Mapping single-end reads $fastq1 complete."
fi