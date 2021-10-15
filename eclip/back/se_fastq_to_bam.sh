#!/usr/bin/env bash

set -e

input_fastq=""
output_fastq=""
memory="2G"

function usage() {
      cat << HELP

Sort FASTQ file by IDs.

Usage: eclip_fastq_sort -i input_fastq -o output_fastq -m memory

Required:
  -i, --input_fastq        Input fastq file.
  -o, --output_fastq       Output fastq file.

Options:
  -m, --memory             Maximum memory for sorting; suffix K/M/G recognized,
                           default: $memory.

Notes:
  For paired-end dataset, comma (,) separated two FASTQ files can be provided for input_fastq
  and output_fastq, respectively.
  If you want the output_fastq being gzip compressed, append .gz extension to
  output_fastq file(s).

HELP
  exit 1
}

[[ $# -eq 0 ]] && usage

function fast_sort() {
    if [[ "$1" == *.gz ]]
    then
      if [[ "$2" == *.gz ]]
      then
        zcat "$1" | fastq-sort --id --temporary-directory "$folder" -S "$memory" | gzip > "$2"
      else
        zcat "$1" | fastq-sort --id --temporary-directory "$folder" -S "$memory" > "$2"
      fi
    else
      if [[ "$2" == *.gz ]]
      then
        fastq-sort --id --temporary-directory "$folder" -S "$memory" "$1" | gzip > "$2"
      else
        fastq-sort --id --temporary-directory "$folder" -S "$memory" "$1" > "$2"
      fi
    fi
}

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
    -m|--memory)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        memory=$2
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

folder=$(mktemp "${output_fastq1}".XXXXXX.sort)
trap 'rm -rf ${folder}' EXIT

fast_sort "${input_fastq1}" "${output_fastq1}"
if [[ "${input_fastq2}" != "" ]]
then
  fast_sort "${input_fastq2}" "${output_fastq2}"
fi