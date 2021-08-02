#!/usr/bin/env bash

set -e

input_fastq=","
output_fastq=","
output_metric=","
adapters=""
barcodes=","
barcodes_fasta=""
randomer_length=10
cpus=1

function usage() {
      cat << HELP

Cut adapters for eCLIP fastq file(s).

Usage: eclip_cut_adapt -i input_fastq -o output_fastq -m output_metric
                       [-a adapters] [-b barcodes] [-f barcodes_fasta]
                       [-r randomer_length] [-c cpus]

Required:
  -i, --input_fastq        Input fastq file(s).
  -o, --output_fastq       Output fastq file(s).
  -m, --output_metric      Metric output files.

Required for single-end dataset:
  -a, --adapters           Adapters fastq file.

Required for paired-end dataset:
  -b, --barcodes           Name of barcodes, e.g., A01, B06, NIL, ... .
  -f, --barcodes_fasta     Path to FASTA file of barcodes.
  -r, --randomer_length    Length of the randomer, default: ${randomer_length}.

Options:
  -c, --cups               Number of CPUs can be used, default: ${cpus}.

Notes:
  For paired-end dataset, comma (,) separated two elements need to be provided for input_fastq,
  output_fastq, and barcodes.
  For output_metric, comma (,) separated two elements need to be provided no matter single-end
  paired-end dataset, for save the first and the second round trimming metrics.

HELP
  exit 1
}

function error() {
    echo "$1"
    exit 1
}

function parse_barcodes() {
    local barcodes
    barcodes=$(grep ">" -v "${1}" | tr '\n', " ${2}")
    echo "${2} ${barcodes}"
}

[[ $# -eq 0 ]] && usage

while (( "$#" )); do
  case "$1" in
    -i|--input-fastq)
        if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
          input_fastq=$2
          shift 2
        else
          echo "Error: Argument for $1 is missing" >&2
          usage
        fi
        ;;
    -o|--output-fastq)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        output_fastq=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
        -m|--output_metric)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        output_metric=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -a|--adapters)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        adapters=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -b|--barcodes)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        barcodes=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -f|--barcodes_fasta)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        barcodes_fasta=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -r|--randomer_length)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        randomer_length=$2
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
    -*)
      echo "Error: Unsupported argument $1" >&2
      exit 1
      ;;
  esac
done

IFS=, read -r input_fastq1 input_fastq2 <<< "${input_fastq}"
IFS=, read -r output_fastq1 output_fastq2 <<< "${output_fastq}"
IFS=, read -r output_metric1 output_metric2 <<< "${output_metric}"

temp_fastq1=$(mktemp "${input_fastq1}".XXXXXX.fastq.gz)
temp_fastq2=$(mktemp "${input_fastq2}".XXXXXX.fastq.gz)
trap 'rm -f ${temp_fastq1} ${temp_fastq2}' EXIT

if [[ "${input_fastq2}" == "" ]]
then
  [ -f "${adapters}" ] || error "Single-end dataset without adapters provided!"
#  adapters=$(parse_barcodes "${adapters}" "-a")
  overlap=1
  echo cutadapt -O "${overlap}" -j "${cpus}" --match-read-wildcards --times 1 \
  -e 0.1 --quality-cutoff 6 -m 18 -a file:"${adapters}" -o "${temp_fastq1}" \
  "${input_fastq1}" > "${output_metric1}"
  cutadapt -O "${overlap}" -j "${cpus}" --match-read-wildcards --times 1 \
  -e 0.1 --quality-cutoff 6 -m 18 -a file:"${adapters}" -o "${temp_fastq1}" \
  "${input_fastq1}" > "${output_metric1}"

  overlap=5
  cutadapt -O "${overlap}" -j "${cpus}" --match-read-wildcards --times 1 \
  -e 0.1 --quality-cutoff 6 -m 18 -a file:"${adapters}" -o "${output_fastq1}" \
  "${temp_fastq1}" > "${output_metric2}"
else
  IFS=, read -r barcode1 barcode2 <<< "${barcodes}"
  [ -f "${barcodes_fasta}" ] || error "Paired-end dataset without barcodes FASTA file provided!"
  eclip_parse_barcodes.sh "${randomer_length}" "${barcodes_fasta}" "${barcode1}" "${barcode2}"
#  g=$(parse_barcodes "g_adapters.fasta" "-g")
#  a=$(parse_barcodes "a_adapters.fasta" "-a")
#  A=$(parse_barcodes "A_adapters.fasta" "-A")
  overlap=$(tr -d '\n' < "trim_first_overlap_length.txt")
  cutadapt -O "${overlap}" -j "${cpus}" --match-read-wildcards --times 1 \
  -e 0.1 --quality-cutoff 6 -m 18 \
  -g file:g_adapters.fasta -a file:a_adapters.fasta -A file:A_adapters.fasta \
  -o "${temp_fastq1}" -p "${temp_fastq2}" \
  "${input_fastq1}" "${input_fastq2}" > "${output_metric1}"

  overlap=$(tr -d '\n' < "trim_again_overlap_length.txt")
  cutadapt -O "${overlap}" -j "${cpus}" --match-read-wildcards --times 1 \
  -e 0.1 --quality-cutoff 6 -m 18 -A file:A_adapters.fasta \
  -o "${output_fastq1}" -p "${output_fastq2}" \
  "${temp_fastq1}" "${temp_fastq2}" > "${output_metric2}"
  rm g_adapters.fasta a_adapters.fasta A_adapters.fasta
  rm trim_first_overlap_length.txt trim_again_overlap_length.txt
fi