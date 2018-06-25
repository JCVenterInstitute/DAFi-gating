
# DAFi Automated Jupyter Container (based on Jupyter Notebook Scientific Python Stack and GenePattern Server)

https://hub.docker.com/r/rarsenal/dafi-gp_jupyter/

docker pull rarsenal/dafi-gp_jupyter


## What it Gives You

* C based DAFi-gating pipeline with automated reports
* Jupyter Notebook 5.2.x
* GenePattern Notebook Plugin for Jupyter
* Conda Python 3.x environment
* bokeh, datashader, dask, pandas, matplotlib, scipy, seaborn, scikit-learn, scikit-image, sympy, cython, patsy, statsmodel, cloudpickle, dill, numba, bokeh, vincent, beautifulsoup, xlrd pre-installed
* Unprivileged user `jovyan` (uid=1000, configurable, see options) in group `users` (gid=100) with ownership over `/home/jovyan` and `/opt/conda`
* [tini](https://github.com/krallin/tini) as the container entrypoint and [start-notebook.sh](../base-notebook/start-notebook.sh) as the default command
* `/usr/local/bin/start-notebook.d` directory for custom init scripts that you can add in derived images
* A [start-singleuser.sh](../base-notebook/start-singleuser.sh) script useful for running a single-user instance of the Notebook server, as required by JupyterHub
* A [start.sh](../base-notebook/start.sh) script useful for running alternative commands in the container (e.g. `ipython`, `jupyter kernelgateway`, `jupyter lab`)
* Options for HTTPS, password auth, and passwordless `sudo`
* GenePattern workflow Engine with DAFi modules

## Basic Use

The following command starts a container with the Notebook server listening for HTTP connections on port 8888 with a randomly generated authentication token configured.

```
docker run -it --rm -p 8888:8888 -p 8883:8080 rarsenal/dafi-gp_jupyter
```

Take note of the authentication token included in the notebook startup log messages. Include it in the URL you visit to access the Notebook server or enter it in the Notebook login form.

## Notebook Options

The Docker container executes a [`start-notebook.sh` script](../base-notebook/start-notebook.sh) script by default. The `start-notebook.sh` script handles the `NB_UID`, `NB_GID` and `GRANT_SUDO` features documented in the next section, and then executes the `jupyter notebook`.

You can pass [Jupyter command line options](https://jupyter.readthedocs.io/en/latest/projects/jupyter-command.html) through the `start-notebook.sh` script when launching the container. For example, to secure the Notebook server with a custom password hashed using `IPython.lib.passwd()` instead of the default token, run the following:

```
docker run -d -p 8888:8888 -p 8883:8080 rarsenal/dafi-gp_jupyter start-notebook.sh --NotebookApp.password='sha1:74ba40f8a388:c913541b7ee99d15d5ed31d4226bf7838f83a50e'
```

For example, to set the base URL of the notebook server, run the following:

```
docker run -d -p 8888:8888 -p 8883:8080 rarsenal/dafi-gp_jupyter start-notebook.sh --NotebookApp.base_url=/some/path
```

For example, to disable all authentication mechanisms (not a recommended practice):

```
docker run -d -p 8888:8888 -p 8883:8080 rarsenal/dafi-gp_jupyter start-notebook.sh --NotebookApp.token=''
```

You can sidestep the `start-notebook.sh` script and run your own commands in the container. See the *Alternative Commands* section later in this document for more information.


## Command-line Options

## dafi_pipeline.sh 
Script for running custom panel configuration in the container as a command line service

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
docker run --rm -p 8888:8888 -v $host_dir:/home/jovyan/work -p 8883:8080 rarsenal/dafi-gp_jupyter start.sh dafi_pipeline.sh --nclusters 200
```

where rarsenal/dafi-gp_jupyter is the tag to the dafi-gp_jupyter container in your system and $host_dir is the directory with the dataset to be 
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
docker run --rm -p 8888:8888 -v $host_dir:/home/jovyan/work -p 8883:8888 rarsenal/dafi-gp_jupyter start.sh LJI132_pipeline.sh --panel 1 --nclusters 200
```

where rarsenal/dafi-gp_jupyter is the tag to the dafi jpynotebook container in your system and $host_dir is the directory with the dataset to be
analyzed.


### Requirements for the host directory ($host_dir):

Within the host directory to be mapped to the working directory of the docker container

* Place all fcs files within a subfolder named FCS ($host_dir/FCS)




## Docker Options

You may customize the execution of the Docker container and the Notebook server it contains with the following optional arguments.

* `-e GEN_CERT=yes` - Generates a self-signed SSL certificate and configures Jupyter Notebook to use it to accept encrypted HTTPS connections.
* `-e NB_UID=1000` - Specify the uid of the `jovyan` user. Useful to mount host volumes with specific file ownership. For this option to take effect, you must run the container with `--user root`. (The `start-notebook.sh` script will `su jovyan` after adjusting the user id.)
* `-e NB_GID=100` - Specify the gid of the `jovyan` user. Useful to mount host volumes with specific file ownership. For this option to take effect, you must run the container with `--user root`. (The `start-notebook.sh` script will `su jovyan` after adjusting the group id.)
* `-e GRANT_SUDO=yes` - Gives the `jovyan` user passwordless `sudo` capability. Useful for installing OS packages. For this option to take effect, you must run the container with `--user root`. (The `start-notebook.sh` script will `su jovyan` after adding `jovyan` to sudoers.) **You should only enable `sudo` if you trust the user or if the container is running on an isolated host.**
* `-v /some/host/folder/for/work:/home/jovyan/work` - Mounts a host machine directory as folder in the container. Useful when you want to preserve notebooks and other work even after the container is destroyed. **You must grant the within-container notebook user or group (`NB_UID` or `NB_GID`) write access to the host directory (e.g., `sudo chown 1000 /some/host/folder/for/work`).**

## SSL Certificates

You may mount SSL key and certificate files into a container and configure Jupyter Notebook to use them to accept HTTPS connections. For example, to mount a host folder containing a `notebook.key` and `notebook.crt`:

```
docker run -d -p 8888:8888 \
    -v /some/host/folder:/etc/ssl/notebook \
    DAFi/jpyenv start-notebook.sh \
    --NotebookApp.keyfile=/etc/ssl/notebook/notebook.key
    --NotebookApp.certfile=/etc/ssl/notebook/notebook.crt
```

Alternatively, you may mount a single PEM file containing both the key and certificate. For example:

```
docker run -d -p 8888:8888 \
    -v /some/host/folder/notebook.pem:/etc/ssl/notebook.pem \
    DAFi/jpyenv start-notebook.sh \
    --NotebookApp.certfile=/etc/ssl/notebook.pem
```

In either case, Jupyter Notebook expects the key and certificate to be a base64 encoded text file. The certificate file or PEM may contain one or more certificates (e.g., server, intermediate, and root).

For additional information about using SSL, see the following:

* The [docker-stacks/examples](https://github.com/jupyter/docker-stacks/tree/master/examples) for information about how to use [Let's Encrypt](https://letsencrypt.org/) certificates when you run these stacks on a publicly visible domain.
* The [jupyter_notebook_config.py](jupyter_notebook_config.py) file for how this Docker image generates a self-signed certificate.
* The [Jupyter Notebook documentation](https://jupyter-notebook.readthedocs.io/en/latest/public_server.html#using-ssl-for-encrypted-communication) for best practices about running a public notebook server in general, most of which are encoded in this image.


## Conda Environments

The default Python 3.x [Conda environment](http://conda.pydata.org/docs/using/envs.html) resides in `/opt/conda`. 

The commands `jupyter`, `ipython`, `python`, `pip`, and `conda` (among others) are available in both environments. For convenience, you can install packages into either environment regardless of what environment is currently active using commands like the following:

```
# install a package into the default (python 3.x) environment
pip install some-package
conda install some-package
```


## Alternative Commands


### start.sh

The `start.sh` script supports the same features as the default `start-notebook.sh` script (e.g., `GRANT_SUDO`), but allows you to specify an arbitrary command to execute. For example, to run the text-based `ipython` console in a container, do the following:

```
docker run -it --rm DAFi/jpyenv start.sh ipython
```

Or, to run JupyterLab instead of the classic notebook, run the following:

```
docker run -it --rm -p 8888:8888 DAFi/jpyenv start.sh jupyter lab
```

This script is particularly useful when you derive a new Dockerfile from this image and install additional Jupyter applications with subcommands like `jupyter console`, `jupyter kernelgateway`, etc.

### Others

You can bypass the provided scripts and specify your an arbitrary start command. If you do, keep in mind that certain features documented above will not function (e.g., `GRANT_SUDO`).