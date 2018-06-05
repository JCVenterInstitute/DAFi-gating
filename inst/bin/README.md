# User-directed unsupervised identification of cell populations - Scripts
	
## Description: 
Utility scripts to be used within docker/singularity containers to run the DAFi pipeline

## Requirement:
Build and run jpynotebook docker/singularity container
(https://github.com/JCVenterInstitute/DAFi-gating/tree/master/docker/JupyterNotebook)

## dafi_pipeline.sh 
Script for running custom panel configuration

```
usage="$(basename "$0") [--help] [-s n] -- DAFi pipeline script to transform, gate, and generate reports on the flowcytometry FCS 
files

where:
    --help  show this help text
    --no_opt  disables intel AVX2 level optimization in case of compatibility issues (default: AVX2 enabled if detected)
    --nclusters specify number of initial clusters to use for the k-means clustering on the whole sample (default: 100)
    --nreclusters specify number of clusters to use for the re-clustering after applying filtering on the parent population (default: 
same as nclusters)
    --ncores specify number of cpu cores to use for the clustering algorithm (default: autodetected)
    --seed specify seed number for the random number generator (default: 2)"
```

### Sample docker run command

```bash
docker run --rm -p 8888:8888 -v $host_dir:/home/jovyan/work dafi/jpynote start.sh dafi_pipeline.sh --nclusters 200
```

where dafi/jpynote is the tag to the dafi jpynotebook container in your system and $host_dir is the directory with the dataset to be 
analyzed. 

### Requirements for the host directory ($host_dir): 

Within the host directory to be mapped to the working directory of the docker container

1) Place all fcs files within a subfolder named FCS ($host_dir/FCS)
	
2) Place inclusion.config: a 11-column tab delimited file, for recursive data filtering within $host_dir/config subfolder ($host_dir/config/inclusion.config)

|Pop_ID|DimensionX|DimensionY|Min_X|Max_X|Min_Y|Max_Y|Parent_ID|Cluster_Type(0: Clustering; 1: Bisecting; 2: Slope-based)|Visualize_or_Not|Recluster_or_Not|Cell_Phenotype(optional)|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|1|1|4|20|70|5|55|0|0|0|0|Lymphocyte|
|2|1|2|30|90|0|110|1|1|0|1|Singlets|
|3|4|5|100|150|80|140|2|2|1|1|LiveSinglets|
|4|19|17|76|140|106|200|3|1|0|1|CD4T|
|5|19|17|76|140|55|105|3|1|0|0|CD8T|
|6|8|7|81|140|50|120|3|1|0|0|CD4Treg|
|7|8|7|20|80|25|90|3|1|0|0|CD4Tnonreg|
	
3) Place exclusion.config: a 11-column tab delimited file with the same format, but for reversed filtering within $host_dir/config ($host_dir/config/exclusion.config)
subfolder

|Pop_ID|DimensionX|DimensionY|Min_X|Max_X|Min_Y|Max_Y|Parent_ID|Cluster_Type(0: Clustering; 1: Bisecting; 2: Slope-based)|Visualize_or_Not|Recluster_or_Not|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|1|1|4|0|85|100|200|0|0|1|0|

4) Place a header list text file named header.lst in the $host_dir/config folder to specify the marker/column order in comma delimited format ($host_dir/config/header.lst)

Example content of header.lst

```
FSC-A,FSC-H,FSC-W,SSC-A,SSC-H,SSC-W,CD3,CD4,CD8,CCR7,CD45RA
```

5) Place a header label replacement text file named header.txt in the $host_dir/config folder to specify change in column name from flurorencence dye name to marker name 
$host_dir/config/header.txt)

Example content of header.txt

First column = fluorescence name
Second column = marker name
(can input multiple variations)

```
Alexa Fluor 488-A	LIVE
Am Cyan-A	IgD
AmCyan-A	IgD
AmCyan-A	IgD
APC-A	CD38
APC-Cy7-A	CD20
APC-H7-A	CD20
B515-A	LIVE
B710-A	CD19
BD Horizon V450-A	CD3
BD Horizon V500-A	IgD
BD Horizon CD3	CD3
FITC-A	LIVE
G560-A	CD24
G780-A	CD27
Pacific Blue-A	CD3
PE-A	CD24
PE-Cy7-A	CD27
PE Cy7 YG-A	CD27
PE-Cy7 YG-A	CD27
PerCP-Cy5-5-A	CD19
PE YG-A	CD24
R660-A	CD38
R780-A	CD20
V450-A	CD3
```

## LJI132_pipeline.sh
Simplified script for running La Jolla Institute for the 132 batch dataset. Configuration for the three panels are already embedded 

```
usage="$(basename "$0") [--help] [-s n] -- DAFi pipeline script to transform, gate, and generate reports on the flowcytometry FCS
files

where:
    --help  show this help text
    --no_opt  disables intel AVX2 level optimization in case of compatibility issues (default: AVX2 enabled if detected)
    --nclusters specify number of initial clusters to use for the k-means clustering on the whole sample (default: 100)
    --nreclusters specify number of clusters to use for the re-clustering after applying filtering on the parent population (default:
same as nclusters)
    --ncores specify number of cpu cores to use for the clustering algorithm (default: autodetected)
    --seed specify seed number for the random number generator (default: 2)"
    --panel specify the panel number (default: 1)
```

### Sample docker run command

```bash
docker run --rm -p 8888:8888 -v $host_dir:/home/jovyan/work dafi/jpynote start.sh LJI132_pipeline.sh --panel 1 --nclusters 200
```

where dafi/jpynote is the tag to the dafi jpynotebook container in your system and $host_dir is the directory with the dataset to be
analyzed.


### Requirements for the host directory ($host_dir):

Within the host directory to be mapped to the working directory of the docker container

1) Place all fcs files within a subfolder named FCS ($host_dir/FCS)

