# Pertussis CyTOF data analysis

### Authors: Aishwarya Mandava, M.S. (amandava@jcvi.org), and Yu "Max" Qian, Ph.D (mqian@jcvi.org or qianyu_cs@yahoo.com)

## Data Transformation and Preprocessing:
The preprocessing of FCS files was done in R. The pre-gating channels and lineage channels of expression matrix from FCS file were arcsinh transformed using *asinh* function in R with a cofactor of *2*. The instrument channels however were not transformed using arcsinh but linearly transformed and re-scaled. The co-factor was empirically determined based on the data distributions as in ***Arcsinh_Co-factor_LJI_CyTOF.png***. The matrix was then linearly rescaled to a range of *0 to 4096* for DAFi visualization, followed by renaming the marker names in R. The preprocessed files, in TXT format, were then used to apply the command line version of DAFi Gating. The input configuration files have co-ordinates of the gates used for each subject. 

### Channel	Channel names
Pre-gating channels:	DNA-1, DNA-2, Viability;
Lineage channels:	CD45, CD3, RGCC, CD19, IL31RA, CD38, CD20-CD4, CD64, CD11c, IGHG1, CX3CR1, CD86, CD123, GBP2, CCL4, CCL2, CD45RA, CD274, OX40, CD161, CD1c, CD80, BCMA, CD33, TNF, CD40, CXCL10, CCR7, IFIT2, CD25, CD54, CCL7, IgM-CD8a, CD14, HLA-DR, CD56, CD16;
Instrument channels:	Time, Event length, Center, Offset, Width, Residual;

### Arcsinh in R:
```
expr <- exprs(FCS)
expr <- cbind(asinh(expr[,c(pregating_channels,lineage_channels)]/cofactor), expr[,instrument_channels])
```
### Linear rescaling: 
```
a*(data-min(data))/(max(data)-min(data)
a=4096
```
### DAFi-Gating:
Inputs are inclusion configuration and exclusion configuration gating files, which can be found in folder *Config*.
Command line:
```
./DAFi_Gating.sh <Input FCS file> <Inclusion config file> <exclusion config file>  <initial cluster size> <re-cluster size>
```
The DAFi command line including number of clusters and seed number can be found in file *DAFi_Gating.sh*

### DAFi configuration files:
Configuration files for each subject can be found in folder *Config*.

### Gating for singlets (***SingletGateforDAFi.png***)
The DAFi gating boundaries for the singlet gate on DNA-1 vs DNA-2 were defined based on kernel density estimation for each subject using ***kde2d*** from MASS package in R with 50 grid points in each direction. The events with low levels on DNA are debris and high level on DNA are doublets. For each subject, we used the sample of the first visit to identify the initial singlet gate location. It is important to note that this gating location is used for the DAFi analysis to identify the cluster of singlet cells in high-dimensional space, instead of outputing the cells defined by the density contour lines. 

### Composites of 2D dot plots
2D dot plots showing the cell populations identified by sequential DAFi-gating analysis can be found in folder *Composites". 

### Reference for data preprocessing:
http://biosurf.org/cytof_data_scientist.html



