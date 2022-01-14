### Authors: 
Yu "Max" Qian, Ph.D., mqian@jcvi.org or qianyu.cs@gmail.com


These are C code for extracting keywords from the TXT results after running FCSTrans to convert and transform the binary FCS files. The output keywords can be used for checking whether there are any typo in the FCS file header or any unstained channels. It is important to note that any typo or inconsistent column names in the header need to be corrected before applying automated analysis and comparison of all the samples.
The header can be updated with user-specified marker names by using ReplaceHeaderRemoveTime, which also removes the time column because it is not useful for identification of cell populations.
