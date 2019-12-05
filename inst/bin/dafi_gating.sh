#!/usr/bin/env bash

usage="$(basename "$0") [--help] [-s n] -- DAFi gating on single txt format FCM sample. Please use FCSTrans to first preprocess and convert the sample FCS format to txt format.

where:
    --help  show this help text
    --nclusters specify number of initial clusters to use for the k-means clustering on the whole sample (default: 100)
    --nreclusters specify number of clusters to use for the re-clustering after applying filtering on the parent population (default: same as nclusters)"

if grep -q avx2 /proc/cpuinfo; then optimize=true; else optimize=false; fi

numOfClusters=100
numOfReClusters=$numOfClusters
numOfCores=24
seed=2

while [ $# -gt 0 ]; do case $1 in
    --help) echo $usage; exit;;
    --inputfile) inputDir="$2"; echo "Input directory = "$inputDir;;
    --inclusioncfg) inclusionCfg="$2"; echo "inclusion config = "$inclusionCfg;;
    --exclusioncfg) exclusionCfg="$2"; echo "exclusion config = "$exclusionCfg;;
    --nclusters) numOfClusters="$2"; echo "number of clusters = "$numOfClusters;;
    --nreclusters) numOfReClusters="$2"; echo "number of reclusters = "$numOfReClusters;;
    *) break;;
  esac; shift; shift
done

cwd=$(pwd)
Label=${PWD##*/}

mkdir -p config
echo "Copying inclusion config "$inclusionCfg
cp $inclusionCfg $cwd/config/inclusion.config
echo "Copying exclusion config "$exclusionCfg
cp $exclusionCfg $cwd/config/exclusion.config

echo "copying FCS files..."
mkdir -p FCS
if [[ -d $inputDir ]]; then
	cp -a $inputDir/. FCS/
elif [[ -f $inputDir ]]; then
	cp $inputDir FCS/
else 
	while read file; do 
		cp "${file}" FCS 
	done < $inputDir
fi

cp $cwd/config/inclusion.config pipeline.config

for file in $cwd/FCS/*
	do
	filepath=$(basename $file )
	extension="${filepath##*.}"
	filename="${filepath%.*}"
	echo Processing File:  $filename
	mkdir $filename
	cd $filename
	dafi_gating $file $cwd/config/inclusion.config $cwd/config/exclusion.config $numOfClusters $numOfReClusters $numOfCores $seed
	cd ..
done

echo "done!"
