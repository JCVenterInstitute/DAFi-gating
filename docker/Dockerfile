FROM bioconductor/release_base2
#setup R configs
RUN apt-get update 
RUN apt-get install -y libtiff5-dev
RUN apt-get install -y libgmp-dev
RUN apt-get install -y libpng-dev

#for building vignettes
RUN apt-get install -y texinfo
RUN apt-get install -y texlive
RUN apt-get install -y texlive-fonts-extra


RUN Rscript -e "source('https://bioconductor.org/biocLite.R')"
RUN Rscript -e "biocLite('ClusterR')"

RUN Rscript -e "biocLite('flowCore')"
RUN Rscript -e "biocLite('flowViz')"
RUN Rscript -e "biocLite('FlowSOM')"
#RUN export ADD=shiny && bash /etc/cont-init.d/add
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "devtools::install_github('JCVenterInstitute/DAFi-gating', build_vignettes = FALSE)"
