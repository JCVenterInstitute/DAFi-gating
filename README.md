# DAFi-gating
# User-directed unsupervised identification of cell populations
	
	Author: Yu "Max" Qian, Ph.D., mqian@jcvi.org or qianyu_cs@yahoo.com
	
	Copyright: Author and J. Craig Venter Institute
	
	Usage: DAFi-gating data_file gating_spec_file filter_spec_file max_num_clusters (by default 100) 
	
	data_file is a tab delimited file with compensated and transformed values.
	
	gating_spec_file: a 10-column tab delimited file, events inside will be kept
	Pop_ID	DimensionX	DimensionY	Min_X	Max_X	Min_Y	Max_Y	Parent_ID	Cluster_Type	Visualize_or_Not
	1       1       4       20      70      5       55      0       0       0
	2       6       9       30      90      0       110     1       1       0
	3       10      19      100     150     80      140     2       0       1
	4       19      17      76      140     106     200     3       1       0
	5       19      17      76      140     55      105     3       1       0
	6       8       7       81      140     50      120     3       1       0
	7       8       7       20      80      25      90      3       1       0
	
	filter_spec file: a 10-column tab delimited file with the same format, but the events inside will be removed instead of being kept:
	1       1       4       0       85      100     200     0       0       1
