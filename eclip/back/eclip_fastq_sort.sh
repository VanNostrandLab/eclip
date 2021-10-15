#!/usr/bin/env bash

set -e

fastq=""
memory="2G"
quiet=0

function usage() {
      cat << HELP

Sort FASTQ file by IDs.

Usage: eclip_fastq_sort FASTQ [-o output_fastq] [-m memory] [-q]

Required:
  FASTQ                    Input fastq file.

Options:
  -o, --output_fastq       Output fastq file.
  -m, --memory             Maximum memory for sorting; suffix K/M/G recognized, default: $memory.
  -q, --quiet              Sort fastq file quietly without printing out message.

Notes:
  - For paired-end dataset, comma (,) separated two FASTQ files can be provided for input_fastq
    and output_fastq, respectively.
  - Without providing output_fastq, 4 known fastq file extensions (.fq, .fastq, .fq.gz., and fastq.gz)
    will have .sort pre-append to one of them as default output extension. If known extensions were not
    found, .sort.fastq.gz will be append to input fastq as default output extension,
    e.g., a.fastq.gz --> a.sort.fastq.gz; abc.def --> abc.def.sort.fastq.gz
  - If you want the output_fastq being gzip compressed, append .gz extension to
    output_fastq file(s).

HELP
  exit 1
}

(( $# )) || usage

function fast_sort() {
		[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Sorting fastq $1 ..."
    if [[ "$1" == *.gz ]]
    then
      if [[ "$2" == *.gz ]]
      then
        cmd="zcat $1 | fastq-sort --id --temporary-directory $folder -S $memory | gzip > $2"
      else
        cmd="zcat $1 | fastq-sort --id --temporary-directory $folder -S $memory > $2"
      fi
    else
      if [[ "$2" == *.gz ]]
      then
        cmd="fastq-sort --id --temporary-directory $folder -S $memory $1 | gzip > $2"
      else
        cmd="fastq-sort --id --temporary-directory $folder -S $memory $1 > $2"
      fi
    fi
    [[ "${quiet}" == 1 ]] || echo "${cmd}"
    eval "${cmd}"
    [[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Sorting fastq $1 complete."
}

while (( "$#" )); do
  case "$1" in
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
    -q|--quiet)
    	quiet=1
    	shift 1
    	;;
    -h|--help)
      usage
      ;;
    *)
    	if [[ "${fastq}" == "" ]]
    	then
      	fastq="$1"
      	shift 1
      else
      	echo "Error: unsupported positional argument $1, only one positional argument allowed." >&2
				usage
			fi
      ;;
  esac
done

function output_name() {
    [[ "$1" != "" ]] || { echo ""; exit 0; }
    if [[ "$2" == "" ]]
    then
    	case $1 in
    			*.fq)
    					echo "${1%.fq}".sort.fq;;
    			*.fastq)
    					echo "${1%.fastq}".sort.fastq;;
    			*.fq.gz)
    					echo "${1%.fq.gz}".sort.fq.gz;;
    			*.fastq.gz)
    					echo "${1%.fastq.gz}".sort.fastq.gz;;
    			*)
    					echo "${1}".sort.fastq.gz;;
    	esac
    else
    	echo "$2"
    fi
}

IFS=, read -r fastq1 fastq2 <<< "${fastq}"
IFS=, read -r output1 output2 <<< "${output_fastq}"
output1=$(output_name "${fastq1}" "${output1}")
output2=$(output_name "${fastq2}" "${output2}")

folder=$(mktemp -d "${output1}".XXXXXX.sort)
trap 'rm -rf ${folder}' EXIT

fast_sort "${fastq1}" "${output1}"
[[ "${fastq2}" == "" ]] || fast_sort "${fastq2}" "${output2}"
