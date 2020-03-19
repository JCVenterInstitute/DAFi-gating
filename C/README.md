# User-directed unsupervised identification of cell populations
## Using the DAFi C
### Author: Yu "Max" Qian, Ph.D., mqian@jcvi.org or qianyu_cs@yahoo.com, Ivan Chang, Ph.D., ichang@jcvi.org, and Bob Sinkovits, Ph.D., sinkovit@sdsc.edu
### Copyright: Author and J. Craig Venter Institute
### Paper link: https://www.biorxiv.org/content/early/2017/09/26/193912

### Usage: 
```
dafi_gating data_file gating_spec_file filter_spec_file initial_max_num_clusters re_cluster_max_num_clusters num_threads seed
```

data_file is a tab delimited file with compensated and transformed values (please use FCSTrans to convert FCS files to tab delimited file https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3932304/) 
	
	gating_spec_file: a 12-column tab delimited file, for recursive data filtering
	Pop_ID	DimensionX	DimensionY	Min_X	Max_X	Min_Y	Max_Y	Parent_ID	Cluster_Type(0: Clustering; 1: Bisecting; 2: Slope-based)	Visualize_or_Not	Recluster_or_Not	Cell_type(optional)
	1       1       4       20      70      5       55      0       0       0	0	Lymph
	2       1       2       30      90      0       110     1       1       0	1	Singlets
	3       4	5	100     150     80      140     2       2       1	1	LiveSinglets
	4       19      17      76      140     106     200     3       1       0	1	CD4T
	5       19      17      76      140     55      105     3       1       0	0	CD8T
	6       8       7       81      140     50      120     3       1       0	0	CD4Treg
	7       8       7       20      80      25      90      3       1       0	0	CD4Tnonreg
	
	filter_spec file: a 11-column tab delimited file with the same format, but for reversed filtering:
	1       1       4       0       85      100     200     0       0       1	0

1) initial_max_num_clusters is the maximum number of clusters for the initial clustering before filtering (default: 500)
2) re_cluster_max_num_clusters is the maximum number of clusters for all of the subsequent re-clustering steps (default: 1000)
3) num_threads is the number of threads available to partition the events for parallelization (default: 8)
4) seed is the initial seed number for the random number generator used in the cluster seeding (default: 2)
	
## Example:
```
dafi_gating preprocessed_fcs.txt forward.config reverse.config 500 1000 24 2
```
Please compile with intel optimization flags:
```
icc -O3 -xHost -o dafi_gating DAFi-gating_omp.c -lm
```
