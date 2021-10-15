#!/usr/bin/env bash

set -e

fastq=""
output=""
memory="2G"
quiet=0
bc_pattern="NNNNNNNNNN"

function usage() {
      cat << HELP

Sort FASTQ file by IDs.

Usage: eclip_umi_extract FASTQ [-o output_fastq] [-q]

Required:
  FASTQ                    Input fastq file.

Options:
  -o, --output_fastq       Output fastq file.
  -b, --bc_pattern         Barcode pattern, default: ${bc_pattern}.
  -q, --quiet              Extract UMIs quietly without printing out message.

Notes:
  - Without providing output_fastq, 4 known fastq file extensions (.fq, .fastq, .fq.gz., and fastq.gz)
    will have .umi pre-append to one of them as default output extension. If known extensions were not
    found, .umi.fastq.gz will be append to input fastq as default output extension,
    e.g., a.fastq.gz --> a.umi.fastq.gz; abc.def --> abc.def.umi.fastq.gz
  - If you want the output_fastq being gzip compressed, append .gz extension to
    output_fastq file.

HELP
  exit 1
}

(( $# )) || usage

while (( "$#" )); do
  case "$1" in
    -o|--output-fastq)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        output=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage
      fi
      ;;
    -b|--bc_pattern)
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
    					echo "${1%.fq}".umi.fq;;
    			*.fastq)
    					echo "${1%.fastq}".umi.fastq;;
    			*.fq.gz)
    					echo "${1%.fq.gz}".umi.fq.gz;;
    			*.fastq.gz)
    					echo "${1%.fastq.gz}".umi.fastq.gz;;
    			*)
    					echo "${1}".umi.fastq.gz;;
    	esac
    else
    	echo "$2"
    fi
}

function fmt_cmd() {
    echo "$1" | sed 's/ -/ \\\n  -/g' | sed 's/ >/ \\\n  >/g' | sed 's/ 2>/ \\\n  2>/g'
}

[[ -f "${fastq}" ]] || { echo "Input FASTQ file ${fastq} does not exist or is not a file!"; exit 1; }
output=$(output_name "${fastq}" "${output}")

[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Extract UMIs in $fastq ..."

IFS='' read -r cmd << EOF
umi_tools extract \
--random-seed 1 \
--stdin $fastq \
--bc-pattern ${bc_pattern} \
--stdout $output \
> /dev/null
EOF

[[ "${quiet}" == 1 ]] || fmt_cmd "${cmd}"
eval "${cmd}"
[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Extract UMIs in $fastq complete."