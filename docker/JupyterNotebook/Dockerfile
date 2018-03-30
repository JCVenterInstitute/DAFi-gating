FROM jupyter/minimal-notebook

USER root

MAINTAINER Ivan Chang (ichang@jcvi.org)

# libav-tools for R and matplotlib anim
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libav-tools \
    fonts-dejavu \
    tzdata \
    gfortran \
    gcc &&  apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN conda install --quiet --yes \
    'nomkl' \
    'cmake' \
    'ipywidgets=7.1*' \
    'pandas=0.19*' \
    'numexpr=2.6*' \
    'matplotlib=2.0*' \
    'scipy=0.19*' \
    'seaborn=0.7*' \
    'scikit-learn=0.18*' \
    'scikit-image=0.12*' \
    'sympy=1.0*' \
    'dask' \
    'statsmodels=0.8*' \
    'cloudpickle=0.2*' \
    'numba=0.31*' \
    'bokeh=0.12*' \
    'sqlalchemy=1.1*' \
    'hdf5=1.8.17'  \
    'r-base=3.4.1' \
    'r-irkernel=0.8*' \
    'r-plyr=1.8*' \
    'r-devtools=1.13*' \
    'r-tidyverse=1.1*' \
    'r-shiny=1.0*' \
    'r-rmarkdown=1.8*' \
    'r-forecast=8.2*' \
    'r-rsqlite=2.0*' \
    'r-reshape2=1.4*' \
    'r-caret=6.0*' \
    'r-rcurl=1.95*' \
    'r-crayon=1.3*' \
    'r-htmltools=0.3*' \
    'r-htmlwidgets=1.0*' \
    'rstudio' \
    'r-hexbin=1.27*' && \
    conda install --quiet --yes -c conda-forge jupyter_contrib_nbextensions && \
    conda install --quiet --yes -c conda-forge nbrsessionproxy && \
    conda install --quiet --yes -c bokeh datashader && \
    conda install --quiet --yes -c ioam holoviews bokeh && \
    conda remove --quiet --yes --force jupyterlab jupyterlab_launcher qt pyqt && \
    conda clean -tipsy && \
    # Activate ipywidgets extension in the environment that runs the notebook server
    jupyter nbextension enable --py widgetsnbextension --sys-prefix && \
    npm cache clean && \
    rm -rf /home/$NB_USER/.cache/yarn && \
    rm -rf /home/$NB_USER/.node-gyp && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

# Install Multicore-TSNE (no conda/pip yet), compile a general non-optimized version of C dafi, and copy binary files
RUN mkdir /home/$NB_USER/bin && \
    mkdir /home/$NB_USER/Notebooks && \
    cd /tmp && \
    git clone https://github.com/DmitryUlyanov/Multicore-TSNE.git && \
    cd Multicore-TSNE/ && \
    yes | pip install . && \
    cd /tmp && \
    rm -rf Multicore-TSNE && \
    git clone https://github.com/JCVenterInstitute/DAFi-gating.git && \
    cd DAFi-gating/C
    gcc -o -o dafi_gating DAFi-gating_omp.c -lm && \
    mv dafi_gating /home/$NB_USER/bin && \
    cd .. && \
    cp FCSTrans/* /home/$NB_USER/bin && \
    cp Notebooks/* /home/$NB_USER/Notebooks && \
    cp inst/bin/* /hoome/$NB_USER/bin && \
    cd .. && \
    rm -rf DAFi-gating && \
    fix-permissions /home/$NB_USER/bin && \
    fix-permissions /home/$NB_USER/Notebooks && \
    fix-permissions /home/$NB_USER
    

# Import matplotlib the first time to build the font cache.
ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot" && \
    fix-permissions /home/$NB_USER

USER $NB_UID

