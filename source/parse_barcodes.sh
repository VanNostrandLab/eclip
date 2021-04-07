#!/usr/bin/env bash

adapter5prime=ACACGACGCTCTTCCGATCT
adapter3prime=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
trim_first_overlap_length=1

randommer_length=$1
barcodes_fasta=$2
barcodeida=$3
barcodeidb=$4

echo $barcodeida > $barcodeida
echo $barcodeidb > $barcodeidb

if [ "${barcodeida}" == NIL ] || [ "${barcodeida}" == RIL ]; then
  echo RIL > barcodea.fasta
  echo RIL > barcodeb.fasta
  titlea=RIL
  barcodea=
  titleb=RIL
  barcodeb=
  sizebarcodea=${#barcodea}
  sizebarcodeb=${#barcodeb}
else
  grep -A 1 $barcodeida $barcodes_fasta > barcodea.fasta
  grep -A 1 $barcodeidb $barcodes_fasta > barcodeb.fasta
  {
    read -r titlea
    read -r barcodea
  }  < barcodea.fasta
  {
    read -r titleb
    read -r barcodeb
  }  < barcodeb.fasta
  titlea=${titlea:1}
  titleb=${titleb:1}
fi

sizebarcodea=${#barcodea}
sizebarcodeb=${#barcodeb}

echo
printf "%b\n" "titlea : ${titlea}"
printf "%b\n" "titleb : ${titleb}"
printf "%b\n" "barcodea : ${barcodea}"
printf "%b\n" "barcodeb : ${barcodeb}"
printf "%b\n" "sizebarcodea : ${sizebarcodea}"
printf "%b\n" "sizebarcodeb : ${sizebarcodeb}"

longestbarcodesize=$(( sizebarcodea > sizebarcodeb ? sizebarcodea : sizebarcodeb ))
trim_again_overlap_length=$((longestbarcodesize - 2))
trim_again_overlap_length=$(( trim_again_overlap_length > 0 ? trim_again_overlap_length : 5 ))
printf "%b\n" "trim_again_overlap_length : ${trim_again_overlap_length}"
printf "%b\n" "randommer_length : ${randommer_length}"
printf "%b\n" "adapter5prime : ${adapter5prime}"
printf "%b\n" "adapter3prime : ${adapter3prime}"

revtr_barcodea=`echo ${barcodea} | rev | tr ATGC TACG`
revtr_barcodeb=`echo ${barcodeb} | rev | tr ATGC TACG`
echo
printf "%b\n" "barcodea : ${barcodea}"
printf "%b\n" "barcodeb : ${barcodeb} "
printf "%b\n" "revtr_barcodea : ${revtr_barcodea}"
printf "%b\n" "revtr_barcodeb : ${revtr_barcodeb}"

revtr_adapter5prime=`echo ${adapter5prime} | rev | tr ATGC TACG`
revtr_adapter3prime=`echo ${adapter3prime} | rev | tr ATGC TACG`
echo
printf "%b\n" "adapter5prime : ${adapter5prime}"
printf "%b\n" "adapter3prime : ${adapter3prime}"
printf "%b\n" "revtr_adapter5prime : ${revtr_adapter5prime}"
printf "%b\n" "revtr_adapter3prime : ${revtr_adapter3prime}"

if [ "${barcodeida}" == NIL ] || [ "${barcodeida}" == RIL ]; then
  g_adapters_fasta=">adapter5prime\n${adapter5prime:10:20}\n"
else
  g_adapters_fasta=">adapter5prime_barcodea\n${adapter5prime:10:20}${barcodea}\n>adapter5prime_barcodeb\n${adapter5prime:10:20}${barcodeb}\n"
fi
echo
printf "%b\n" "g_adapters_fasta : \n${g_adapters_fasta}"

a_adapters_fasta=">randommer_adapter3prime\n"$(printf %${randommer_length}s | tr " " "N")${adapter3prime}
if [ "${randommers_length}" = "0" ]; then
  a_adapters_fasta=""
fi
echo
printf "%b\n" "a_adapters_fasta : \n${a_adapters_fasta}"

revtr_barcodea__revtr_adapter5prime=${revtr_barcodea}${revtr_adapter5prime}
echo "revtr_barcodea__revtr_adapter5prime : ${revtr_barcodea__revtr_adapter5prime}"
revtr_barcodeb__revtr_adapter5prime=${revtr_barcodeb}${revtr_adapter5prime}
echo "revtr_barcodeb__revtr_adapter5prime : ${revtr_barcodeb__revtr_adapter5prime}"
echo

if [ $sizebarcodea -eq 0 ]; then
  A_adapters_seed="AGATCGGAAGAGCGT;GATCGGAAGAGCGTC;ATCGGAAGAGCGTCG;TCGGAAGAGCGTCGT"
  A_adapters_seed="${A_adapters_seed};CGGAAGAGCGTCGTG;GGAAGAGCGTCGTGT;GAAGAGCGTCGTGTA;AAGAGCGTCGTGTAG"
else
  for i in $( seq 0 $(( $sizebarcodea -1 )) )
  do
    A_adapters_seed="${A_adapters_seed}${revtr_barcodea__revtr_adapter5prime:${i}:15};"
    echo "${i} - a: ${A_adapters_seed}"
    A_adapters_seed="${A_adapters_seed}${revtr_barcodeb__revtr_adapter5prime:${i}:15};"
    echo "${i} - b: ${A_adapters_seed}"
  done
  for i in {0..5}
  do
    echo "${i}: ${revtr_adapter5prime:${i}:15}"
    A_adapters_seed="${A_adapters_seed}${revtr_adapter5prime:${i}:15};"
  done
fi

A_adapters_seed=`echo ${A_adapters_seed} | tr ";" "\n"`
A_adapters=""
count=0
OLD_IFS=$IFS
while IFS=$'\t' read A_adapter; do
  A_adapters_fasta="${A_adapters_fasta}>A_adapter_${count}\n${A_adapter}\n"
  count=$(($count+1))
done  <<<  "$A_adapters_seed"
echo
IFS=$OLD_IFS
printf "A_adapters_seed :\n${A_adapters_seed}"
echo
printf "A_adapters_fasta : \n${A_adapters_fasta}"

echo -n "$trim_first_overlap_length" > trim_first_overlap_length.txt
echo -n "$trim_again_overlap_length" > trim_again_overlap_length.txt

echo -e "\n\n" > g_adapters_default.fasta
echo -e "\n\n" > a_adapters_default.fasta

echo -e $g_adapters_fasta > g_adapters.fasta
echo -e $a_adapters_fasta > a_adapters.fasta
echo -e $A_adapters_fasta > A_adapters.fasta

set -ex
