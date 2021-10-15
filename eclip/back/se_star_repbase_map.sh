#!/usr/bin/env bash

set -e

function usage() {
      cat << HELP

Map reads in a single FASTQ file to repbase repeat elements using STAR.

Usage: star_repbase_map -i fastq -o bam -x index [-c cpus] [-h]

Required:
  -i, --input_fastq               Input FASTQ file.
  -o, --output_bam                Output BAM file.
  -x, --index                     STAR repbase repeat elements index directory.

Options:
  -c, --cpus                      Maximum number of CPUS can be used for mapping.
  -h, --help                      Print out help message and exit.


HELP
  exit 1
}

(( "$#" )) || usage

while (( "$#" )); do
  case "$1" in
    -i|--input_fastq)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        fastq=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -o|--output_bam)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        bam=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -x|--index)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        index=$2
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
    -h|--help)
      usage
      ;;
    *)
      echo "Error: Unsupported argument $1" >&2
      usage
      ;;
  esac
done

[[ -f "${fastq}" ]] || { echo "Input FASTQ file $fastq does not exist or is not a file!"; exit 1; }
[[ -d "${index}" ]] || { echo "Repbase index $index does not exist or is not a directory!"; exit 1; }
[[ "${bam}" == *.bam ]] || { echo "Output BAM $bam does not end with .bam extension!"; exit 1; }
[[ "${fastq}" == *.gz ]] && read_file_command="zcat" || read_file_command="-"
cpus="${cpu:-1}"

folder=$(mktemp -d "${bam/.bam/}".XXXXXXXXXX.repbase.map)
trap 'rm -rf ${folder}' EXIT

cat <<- HELP
STAR \
  --runMode alignReads \
  --runThreadN ${cpus} \
  --alignEndsType EndToEnd \
  --genomeDir ${index} \
  --outBAMcompression 10 \
  --outFileNamePrefix ${folder}/ \
  --outFilterMultimapNmax 100 \
  --outFilterScoreMin 10 \
  --outReadsUnmapped Fastx \
  --outSAMattrRGline ID:foo \
  --outSAMattributes All \
  --outSAMtype BAM Unsorted \
  --readFilesCommand ${read_file_command} \
  --readFilesIn ${fastq}
HELP

STAR \
   --runMode alignReads \
   --runThreadN "${cpus}" \
   --alignEndsType EndToEnd \
   --genomeDir "${index}" \
   --outBAMcompression 10 \
   --outFileNamePrefix "${folder}"/ \
   --outFilterMultimapNmax 100 \
   --outFilterScoreMin 10 \
   --outReadsUnmapped Fastx \
   --outSAMmode Full \
   --outSAMattrRGline ID:foo \
   --outSAMattributes All \
   --outSAMtype BAM Unsorted \
   --readFilesCommand "${read_file_command}" \
   --readFilesIn "${fastq}"

mv "${folder}"/Aligned.out.bam "${bam}"
mv "${folder}"/Log.final.out "${bam/.bam/.log}"
pigz -c -p "${cpus}" "${folder}"/Unmapped.out.mate1 > "${bam/.bam/.mate1.gz}"