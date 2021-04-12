#!/usr/bin/env bash
set -e
unset PYTHONPATH

conda_install () {
  echo "    Installing $3 ...";
  conda install --channel="$1" --prefix="$2" --yes --quiet "$3" >/dev/null
  echo "    Successfully Installed $3.";
}

echo "Checking environment ..."
for package in conda wget git
do
  if command -v "${package}" >/dev/null; then
    echo "    ${package} ... installed";
  else
    echo "    Installing CLIP pipeline requires ${package} but it's not installed. Aborted!"; exit 1;
  fi
done
echo "Checking environment complete."

eclip="$( cd .; pwd -P )"
venv="${eclip}/venv"

echo "Set up virtual environment for clip ...";

echo "    Installing Python (3.8) ..."
conda create --prefix="${venv}" --yes --quiet python=3.8 >/dev/null
echo "    Successful installed Python (3.8)."

conda_install "conda-forge" "${venv}" cython

echo "    Installing Perl and Perl packages ... "
conda_install "bioconda" "${venv}" perl
conda_install "bioconda" "${venv}" perl-app-cpanminus
"${venv}/bin/cpanm" Statistics::Basic --quiet >/dev/null;
"${venv}/bin/cpanm" Statistics::Distributions --quiet >/dev/null;
"${venv}/bin/cpanm" Statistics::R --quiet >/dev/null;
echo "    Successful installed 3 Perl packages."

for package in bedtools cutadapt umi_tools fastq-tools samtools=1.9 star=2.4.0j pysam ruffus pandas loguru pureclip
do
  conda_install "bioconda" "${venv}" "${package}"
done

for package in ucsc-bedgraphtobigwig ucsc-bedsort
do
  conda_install "bioconda" "${venv}" "${package}"
done

"${venv}/bin/pip" install future --prefix="${venv}" --quiet

echo "    Installing eclipdemux ..."
git clone --quiet https://github.com/VanNostrandLab/eclipdemux.git "${eclip}/eclipdemux";
"${venv}/bin/pip" install "${eclip}/eclipdemux" --prefix="${venv}" --quiet;
rm -rf "${eclip}/eclipdemux";
echo "    Successfully Installed eclipdemux."

echo "    Installing clipper ..."
git clone --quiet https://github.com/VanNostrandLab/clipper.git "${eclip}/clipper";
"${venv}/bin/pip" install "${eclip}/clipper" --prefix="${venv}" --quiet;
rm -rf "${eclip}/clipper"
echo "    Successfully Installed clipper."