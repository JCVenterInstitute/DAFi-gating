FROM bioconductor/release_base2

MAINTAINER ichang@jcvi.org

#setup R configs
#install libs needed for ClusterR and tex needed for building vignette 
RUN apt-get update && \
    apt-get -y install --fix-missing libtiff5-dev libgmp-dev libpng-dev texinfo texlive texlive-fonts-extra


#install R packages
RUN Rscript -e "source('https://bioconductor.org/biocLite.R')" && \
    Rscript -e "biocLite('ClusterR')" && \
    Rscript -e "biocLite('flowCore')" && \
    Rscript -e "biocLite('flowViz')" && \
    Rscript -e "biocLite('FlowSOM')" && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('devtools')" && \
    Rscript -e "devtools::install_github('JCVenterInstitute/DAFi-gating', build_vignettes = FALSE)"

#RUN export ADD=shiny && bash /etc/cont-init.d/add
