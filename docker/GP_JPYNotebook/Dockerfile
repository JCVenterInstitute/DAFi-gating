FROM rarsenal/dafi-jupyter

USER root

MAINTAINER Ivan Chang (ichang@jcvi.org)

# libav-tools for R and matplotlib anim
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openjdk-8-jre &&  apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ADD start.sh /usr/local/bin/start.sh
ADD start-notebook.sh /usr/local/bin/start-notebook.sh

RUN mkdir /usr/src/myapp

COPY install/* /usr/src/myapp/

RUN cd /usr/src/myapp && \
	wget http://software.broadinstitute.org/cancer/software/genepattern/download_installer/installers/gp-3.9.10-153/Linux/NoVM/GPserver.bin && \
	chmod +x GPserver.bin && \
	PATH=$PATH:/usr/src/myapp && \
	./GPserver.bin -f install-gp.properties && \
	rm GPserver.bin && \
	fix-permissions /opt/genepattern

EXPOSE 8080
#EXPOSE 9001
VOLUME /home/$NB_USER/work
USER $NB_UID
