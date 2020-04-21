Pertussis CyTOF data analysis

Data Transformation and Preprocessing:
The preprocessing of FCS files was done in R. The pre-gating channels and lineage channels of expression matrix from FCS file were arcsinh transformed using asinh function in R with a cofactor of 2. The instrument channels however were not transformed using arcsinh. The co-factor was empirically determined based on the distributions. The matrix was then linearly rescaled to a range of 0 to 4096 for DAFi visualization, followed by this the marker names were renamed in R. The preprocessed “.txt” files were then used to apply the command line version of DAFi Gating. The input configuration files have co-ordinates of the gates used for each subject. 

Channel	Channel names
Pre-gating channels	DNA-1, DNA-2, Viability
Lineage channels	CD45, CD3, RGCC, CD19, IL31RA, CD38, CD20-CD4, CD64, CD11c, IGHG1, CX3CR1, CD86, CD123, GBP2, CCL4, CCL2, CD45RA, CD274, OX40, CD161, CD1c, CD80, BCMA, CD33, TNF, CD40, CXCL10, CCR7, IFIT2, CD25, CD54, CCL7, IgM-CD8a, CD14, HLA-DR, CD56, CD16
Instrument channels	Time, Event length, Center, Offset, Width, Residual

Arcsinh in R:
expr <- exprs(FCS)
expr <- cbind(asinh(expr[,c(pregating_channels,lineage_channels)]/cofactor), expr[,instrument_channels])

Linear rescaling: 
a*(data-min(data))/(max(data)-min(data)
a=4096

DAFi-Gating:
Inputs are inclusion configuration and exclusion configuration gating files. 
Implementation: ./DAFi_Gating.sh <Input FCS file> <Inclusion config file> <exclusion config file>  <initial cluster size> <re-cluster size>

Gating for singlets:
The gating boundaries for DNA-1 vs DNA-2 were defined based on kernel density estimation for each subject using kde2d from MASS package in R with 50 grid points in each direction. The events with low levels on DNA are debris and high level on DNA are doublets, so, the singlets gate in DAFi inclusion file was chosen based on the kernel density estimation of the first visit for each subject. 
