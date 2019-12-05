#!/usr/bin/env bash

usage="$(basename "$0") [--help] [-s n] -- DAFi pipeline script to transform, gate, and generate reports on the flowcytometry FCS files

where:
    --help  show this help text
    --nclusters specify number of initial clusters to use for the k-means clustering on the whole sample (default: 200)
    --nreclusters specify number of clusters to use for the re-clustering after applying filtering on the parent population (default: 500)
    --sizeparam Specify size of individual dots (default = 1)
    --sizeofaxislabels Specify size of axis labels (default = 14)
    --titleshow Set whether to show filename as title
    --sizeofpoplabels Specify the size of population label (min = 5, max = 20, default 10)
    --showmultigates Set whether to show multiple gates on same 2D plot
    --showparent Set whether to show parent population or not
    --sort Sort the population drawing sequence
    --showgrid Set grid visible or not
    --hidegatelines Set whether to hide gate lines or not"

if grep -q avx2 /proc/cpuinfo; then optimize=true; else optimize=false; fi

numOfClusters=200
numOfReClusters=500
numOfCores=$(grep -c ^processor /proc/cpuinfo)
seed=2

while [ $# -gt 0 ]; do case $1 in
    --help) echo $usage; exit;;
    --inputdir) inputDir="$2"; echo "Input directory = "$inputDir;;
    --inclusioncfg) inclusionCfg="$2"; echo "inclusion config = "$inclusionCfg;;
    --exclusioncfg) exclusionCfg="$2"; echo "exclusion config = "$exclusionCfg;;
    --nclusters) nclusters="$2"; echo "number of clusters = "$nclusters;;
    --nreclusters) nreclusters="$2"; echo "number of clusters for reclustering = "$nreclusters;;
    --sizeparam) sizeparam="$2"; echo "size of individual dots = "$sizeparam;;
    --sizeofaxislabels) sizeofaxislabels="$2"; echo "size of axis labels = "$sizeofaxislabels;;
    --sizeofpoplabels) sizeofpoplabels="$2"; echo "size of population label  = "$sizeofpoplabels;;
    --titleshow) 
        if [ "$2" == 1 ]; then 
            titleshow="--titleshow" 
        fi; echo "show filename as title = "$titleshow;;
    --showmultigates) 
        if [ "$2" == 1 ]; then 
            showmultigates="--showmultigates"
        fi; echo "show multiple gates  = "$showmultigates;;
    --showparent) 
        if [ "$2" == 1 ]; then 
            showparent="--showparent" 
        fi; echo "show parent population = "$showparent;;
    --sort) 
        if [ "$2" == 1 ]; then 
            sortp="--sort" 
        fi; echo "Sort the population drawing sequence = "$sortp;;
    --showgrid) 
        if [ "$2" == 1 ]; then 
            showgrid="--gridshow" 
        fi; echo "grid visible = "$showgrid;;
    --hidegatelines) 
        if [ "$2" == 1 ]; then 
            hidegatelines="--hidegatelines"
        fi; echo "hide gate lines = "$hidegatelines;;
    *) break;;
  esac; shift; shift
done



cwd=$(pwd)
Label=${PWD##*/}

mkdir -p config
echo "Copying inclusion config "$inclusionCfg
cp $inclusionCfg $cwd/config/inclusion.config
echo "Copyong exclusion config "$exclusionCfg
cp $exclusionCfg $cwd/config/exclusion.config

#[ -z "$inputComp" ] && echo "No input compensation file specified" || cp $inputComp $cwd/config
echo "copying FCS files..."
mkdir -p FCS
while read file; do
        cp "${file}" FCS
done < $inputDir

cp $cwd/config/inclusion.config pipeline.config

#Gating and filtering via DAFi framework
echo "Gating...."
if [ -d "Gated" ]; then
        rm -Rf Gated
        mkdir Gated
else
        mkdir Gated
fi
cd Gated
for file in $cwd/FCS/*
        do
        filepath=$(basename $file )
        extension="${filepath##*.}"
        filename="${filepath%.*}"
        echo Processing File:  $filename
        mkdir $filename
        cd $filename
        if [ "$optimize" = true ] ; then
                dafi_intel $file $cwd/config/inclusion.config $cwd/config/exclusion.config $nclusters $nreclusters $numOfCores $seed
        else
                dafi_gating $file $cwd/config/inclusion.config $cwd/config/exclusion.config $nclusters $nreclusters $numOfCores $seed
        fi
        cd ..
done

#Create composite population percentage/events table across data set
echo "generating tables..."
parsePercentage.awk */DU*percentage.txt > Batch_percentages.txt
parsePercentage.awk */DU*events.txt > Batch_events.txt
cd ..
echo "generating html report..."

cp /var/DAFi-gating/Notebooks/StaticCompositePlots.ipynb .
time jupyter nbconvert --ExecutePreprocessor.timeout=10000 --to html_embed --execute StaticCompositePlots.ipynb > jpy.out 2>jpyerror.txt
mkdir Reports
mv *.html Reports
source activate py27
mkdir Plots
cd Plots
/var/DAFi-gating/Python/DAFi-plotting.py -f $cwd"/Gated/*/flock*.txt" --config ../pipeline.config --sizeofdot $sizeparam --pthread 24 --sizeofaxislabels $sizeofaxislabels --sizeofpoplabels $sizeofpoplabels $titleshow $showmultigates $showparent $sortp $showgrid $hidegatelines
source deactivate
echo "done!"
