FROM rarsenal/dafi-jupyter

USER root

MAINTAINER Ivan Chang (ichang@jcvi.org)

# libav-tools for R and matplotlib anim
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openjdk-8-jre &&  apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /usr/src/myapp

COPY install/* /usr/src/myapp/

# wget http://software.broadinstitute.org/cancer/software/genepattern/download_installer/installers/gp-3.9.10-153/Linux/NoVM/GPserver.bin && \
        
RUN cd /usr/src/myapp && \
	wget https://github.com/genepattern/genepattern-server/releases/download/v3.9.11-rc.4%2Bb228/GPserver-v3.9.11-rc.4+b228.bin && \
        chmod +x GPserver-v3.9.11-rc.4+b228.bin && \
	PATH=$PATH:/usr/src/myapp && \
	./GPserver-v3.9.11-rc.4+b228.bin -f install-gp.properties && \
	rm GPserver-v3.9.11-rc.4+b228.bin && \
	fix-permissions /opt/genepattern

EXPOSE 8080
#EXPOSE 9001
VOLUME /home/$NB_USER/work

ADD start-notebook.sh /usr/local/bin/start-notebook.sh
ADD run_with_path.sh /opt/genepattern/resources/wrapper_scripts/run_with_path.sh
ADD start.sh /usr/local/bin/start.sh

RUN cd /var/DAFi-gating && \
	git pull

# add details to  the .profile so GP can see them
ADD profile /home/jovyan/profile
RUN cat /home/jovyan/profile >> /home/jovyan/.profile


USER $NB_UID
