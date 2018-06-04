# User-directed unsupervised identification of cell populations
	
## DAFi: User-directed unsupervised filtering and identification of cell populations from flow cytometry data

### Author: 
Yu "Max" Qian, Ph.D., mqian@jcvi.org or qianyu_cs@yahoo.com, Ivan Chang, Ph.D., ichang@jcvi.org, and Bob Sinkovits, Ph.D., sinkovit@sdsc.edu
	
### Copyright: Author and J. Craig Venter Institute
	
Paper link: https://www.biorxiv.org/content/early/2017/09/26/193912

## Description: 
The DAFi package provides a new framework for cell population identification for flow cytometry (FCM) data. The framework is compatible with many existing clustering algorithms such as Kmeans, Kmeans++, mini-batch Kmeans, gaussian mixture models, k-medoids, self-organizing map, etc, and allows user to input user defined gating hierarchy to convert the aforementioned unsupervised algorithms into powerful semi-supervised and automatic cell population identification approach. First, the data is read and preprocessed, then configuration files are parsed and implemented as user defined gating directions. Clustering algorithm choosen in the initiation is then applied to the data events of the whole FCM data at the start of the iterative gating/filtering loop, and optionally again at each specified population subset (reclustering). At each gating step, the user can specify to apply bisecting (events filtering based on user defined boundaries just like in manual gating), slope-based (events filtering based on user defined slopes), or cluster centroids based filtering (filtering all events members of a cluster based on their centroid's inclusion or exclusion by the user defined gates). Outputs include population events and percentages table, as well as an events printout table consisting of all event's transformed channel values and population membership info for external analysis and plotting. Several built-in plotting options are also available, e.g. 2D dot plots of the user specified gating channels, and centroids overlay to the 2D dot plots.


### DAFi R docker image
If you have docker containerization system enabled, you can download the pre-configured Dockerbuild and build the R-DAFi dockerized container that will allow you to run a local R-studio server with all necessary packages.

Move the Dockerbuild file to a directory, then run with:

```
docker build -t "r_dafi:r_dafi"
```

after the image is built, instantiate it with (change password to your own):

```
docker run -d -p 8787:8787 -e PASSWORD=dafi -e ROOT=TRUE r_dafi:r_dafi
``` 

then open the Rstudio on localhost via a webbrowser:

http://localhost:8787

and log in using user name rstudio and password dafi
