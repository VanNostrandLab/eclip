#!/usr/bin/env bash

set -e

function usage() {
      cat << HELP

Cut adapters for eCLIP reads in fastq file(s).

Usage: eclip_cut_adapt FASTQ [-o output] [-c cpus] [-q]
                       [-a adapters] [-b barcodes] [-f barcodes_fasta] [-r randomer_length]

Required:
  FASTQ                     Path to input fastq file.

Required for single-end dataset:
  -a, --adapters_fasta      Adapters fastq file.

Required for paired-end dataset:
  -b, --barcodes            Comma separated names of barcodes, e.g., A01,B06 or NIL,NIL, ... .
  -f, --barcodes_fasta      Path to FASTA file of barcodes.
  -r, --randomer_length     Length of the randomer, default: ${randomer_length}.

Options:
	-o, --output              Path to output fastq file.
  -c, --cups                Number of CPUs can be used, default: ${cpus}.
	-q, --quiet              Sort fastq file quietly without printing out message.

Notes:
	For paired-end dataset, two fastq files and two output files can be provided as a comma
	separated element, e.g., fastq1,fastq2, or -o output1,output2.
HELP
  exit 1
}


function parse_barcodes() {
    local barcodes
    barcodes=$(grep ">" -v "${1}" | tr '\n', " ${2}")
    echo "${2} ${barcodes}"
}

(( "$#" )) || usage

while (( "$#" )); do
  case "$1" in
    -a|--adapters_fasta)
      if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
        adapters_fasta=$2
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
    -o|--output)
			if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
				output=$2
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

fastq="${fastq:-""}"
IFS=, read -r fastq1 fastq2 <<< "${fastq}"
[[ -f "${fastq1}" ]] || { echo "Input FASTQ file ${fastq1} does not exist or is not a file!"; exit 1; }
cpus="${cpus:-1}"
quiet="${quiet:-0}"
randomer_length="${randomer_length:-10}"
adapters_fasta="${adapters_fasta:-""}"
barcodes="${barcodes:-""}"
barcodes_fasta="${barcodes_fasta:-""}"
output="${output:-""}"

if [ -n "${output}" ]
then
	IFS=, read -r output1 output2 <<< "${output}"
	metric1="${output1%.fastq.gz}"
	metric2="${output1%.fastq.gz}"
	metric1="${metric1%.trim}.trim.metrics"
	metric2="${metric2%.trim}.trim.trim.metrics"
else
	output1="${fastq1%.fastq.gz}.trim.fastq.gz"
  output2="${fastq2%.fastq.gz}.trim.fastq.gz"
  metric1="${output1%.fastq.gz}.metrics"
  metric2="${output1%.fastq.gz}.trim.metrics"
fi

temp_fastq1=$(mktemp "${output1}".XXXXXXXXXX)
temp_fastq2=$(mktemp "${output2}".XXXXXXXXXX)
trap 'rm -f ${temp_fastq1} ${temp_fastq2}' EXIT

function fmt_cmd() {
    echo "$1" | sed 's/  / \\\n  /g' | sed 's/ -/ \\\n  -/g' | sed 's/ >/ \\\n  >/g' | sed 's/ 2>/ \\\n  2>/g'
}

if [[ "${fastq2}" == "" ]]
then
  [ -f "${adapters_fasta}" ] || { echo "Single-end dataset without adapters_fasta provided!"; exit 1; }
  [[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Cut adapters in single-end reads $fastq1 ..."
  overlap=1

  IFS='' read -r cmd << EOF
cutadapt \
--quiet \
--overlap ${overlap} \
--cores ${cpus} \
--match-read-wildcards \
--times 1 \
--error-rate 0.1 \
--quality-cutoff 6 \
--minimum-length 18 \
--adapter file:${adapters_fasta} \
--output ${temp_fastq1}  ${fastq1} \
> ${metric1}
EOF
  [[ "${quiet}" == 1 ]] || fmt_cmd "${cmd}"
  eval "${cmd}"

  overlap=5
  IFS='' read -r cmd << EOF
cutadapt \
--quiet \
--overlap ${overlap} \
--cores ${cpus} \
--match-read-wildcards \
--times 1 \
--error-rate 0.1 \
--quality-cutoff 6 \
--minimum-length 18 \
--adapter file:${adapters_fasta} \
--output ${output1}  ${temp_fastq1} \
> ${metric2}
EOF
	[[ "${quiet}" == 1 ]] || fmt_cmd "${cmd}"
  eval "${cmd}"
  [[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Cut adapters in $fastq1 complete."
else
	[[ -f "${fastq2}" ]] || { echo "Input FASTQ file ${fastq2} does not exist or is not a file!"; exit 1; }
	[[ ${barcodes} != "" ]] || { echo "Paired-end dataset without comma separated barcodes provided!"; exit 1; }
  IFS=, read -r barcode1 barcode2 <<< "${barcodes}"
  [ -f "${barcodes_fasta}" ] || { echo "Paired-end dataset without barcodes FASTA file provided!"; exit 1; }

  [[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Cut adapters in paired-end reads: $fastq1 and"
  [[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")]                                   $fastq2 ..."

  cmd="eclip_parse_barcodes.sh ${randomer_length} ${barcodes_fasta} ${barcode1} ${barcode2} > /dev/null"
	[[ "${quiet}" == 1 ]] || echo "${cmd}"
	eval "${cmd}"

  overlap=$(tr -d '\n' < "trim_first_overlap_length.txt")
  IFS='' read -r cmd << EOF
cutadapt \
--quiet \
--overlap ${overlap} \
--cores ${cpus} \
--match-read-wildcards \
--times 1 \
--error-rate 0.1 \
--quality-cutoff 6 \
--minimum-length 18 \
--front file:g_adapters.fasta \
--adapter file:a_adapters.fasta \
-A file:A_adapters.fasta \
--output ${temp_fastq1} \
--paired-output ${temp_fastq2}  ${fastq1}  ${fastq2} \
> ${metric1}
EOF
	[[ "${quiet}" == 1 ]] || fmt_cmd "${cmd}"
	eval "${cmd}"

  overlap=$(tr -d '\n' < "trim_again_overlap_length.txt")
  IFS='' read -r cmd << EOF
cutadapt \
--quiet \
--overlap ${overlap} \
--cores ${cpus} \
--match-read-wildcards \
--times 1 \
--error-rate 0.1 \
--quality-cutoff 6 \
--minimum-length 18 \
-A file:A_adapters.fasta \
--output ${output1} \
--paired-output ${output2}  ${temp_fastq1}  ${temp_fastq2} \
> ${metric2}
EOF
	[[ "${quiet}" == 1 ]] || fmt_cmd "${cmd}"
	eval "${cmd}"

	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")] Cut adapters in paired-end reads: $fastq1 and"
	[[ "${quiet}" == 1 ]] || echo "[$(date +"%H:%M:%S")]                                   $fastq2 complete."

  rm "${barcode1}" "${barcode2}" barcodea.fasta barcodeb.fasta
  rm g_adapters_default.fasta a_adapters_default.fasta
  rm g_adapters.fasta a_adapters.fasta A_adapters.fasta
  rm trim_first_overlap_length.txt trim_again_overlap_length.txt
fi