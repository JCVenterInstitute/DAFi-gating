#!/usr/bin/env bash

usage="$(basename "$0") [--help] [-s n] -- DAFi pipeline script to transform, gate, and generate reports on the flowcytometry FCS files

where:
    --help  show this help text
    --no_opt  disables intel AVX2 level optimization in case of compatibility issues (default: AVX2 enabled if detected)
    --nclusters specify number of initial clusters to use for the k-means clustering on the whole sample (default: 100)
    --nreclusters specify number of clusters to use for the re-clustering after applying filtering on the parent population (default: same as nclusters)
    --ncores specify number of cpu cores to use for the clustering algorithm (default: autodetected)
    --seed specify seed number for the random number generator (default: 2)"

if grep -q avx2 /proc/cpuinfo; then optimize=true; else optimize=false; fi

numOfClusters=100
numOfReClusters=$numOfClusters
numOfCores=$(grep -c ^processor /proc/cpuinfo)
seed=2

while [[ "$#" > 1 ]]; do case $1 in
    --help) echo $usage; exit;;
    --no_opt) optimize=false;;
    --nclusters) numOfClusters="$2";;
    --nreclusters) numOfReClusters="$2";;
    --ncores) numOfCores="$2";;
    --seed) seed="$2";;
    *) break;;
  esac; shift; shift
done

if [ -d "work" ]; then
        cd ~/work
        cwd=$(pwd)
else
        cwd=$(pwd)
fi
cwd=$(pwd)
Label=${PWD##*/}

#Transformation and compensation using FCSTrans
#Convert FCS binary to txt format 
if [ -d "TXT" ]; then
        rm -Rf TXT
        mkdir TXT
else
        mkdir TXT
fi
cd TXT
run_FCSTrans2TXT.R $cwd/FCS
cd ..

#Arrange columns and replace labels
if [ -d "Preprocessed" ]; then
        rm -Rf Preprocessed
        mkdir Preprocessed
else
        mkdir Preprocessed
fi
cd Preprocessed
count=1
for f in ../TXT/*.txt
do
filename=$(basename "$f");
(ArrangeHeader.awk $cwd/config/header.lst $cwd/config/header.txt $f > $filename)&
if [ $(($count % $numOfCores)) -eq 0 ]; then wait; fi
count=$((count+1));
done
wait
cd ..
cp $cwd/config/inclusion.config pipeline.config

#Gating and filtering via DAFi framework
if [ -d "Gated" ]; then
        rm -Rf Gated
        mkdir Gated
else
        mkdir Gated
fi
cd Gated
for file in $cwd/Preprocessed/*
	do
	filepath=$(basename $file )
	extension="${filepath##*.}"
	filename="${filepath%.*}"
	echo Processing File:  $filename
	mkdir $filename
	cd $filename
	if [ "$optimize" = true ] ; then
		dafi_intel $file $cwd/config/inclusion.config $cwd/config/exclusion.config $numOfClusters $numOfReClusters $numOfCores $seed
	else
		dafi_gating $file $cwd/config/inclusion.config $cwd/config/exclusion.config $numOfClusters $numOfReClusters $numOfCores $seed
	fi
	cd ..
done

#Create composite population percentage/events table across data set
parsePercentage.awk */DU*percentage.txt > Batch_percentages.txt
parsePercentage.awk */DU*events.txt > Batch_events.txt
cd ..
time jupyter nbconvert --ExecutePreprocessor.timeout=10000 --to html_embed --execute AutoReport.ipynb
