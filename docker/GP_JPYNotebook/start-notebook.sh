#!/bin/bash
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

set -e

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.
GPONLY=0	 # Whether to run GP only or start notebook server
# Sync with Github repo for latest DAFi scripts and notebooks

while getopts "u?f?g?:" opt; do
    case "$opt" in
    u)  cd /var/DAFi-gating
        git pull

        cp /var/DAFi-gating/GP/resources/* /opt/genepattern/resources
	echo "copying resources from github"
        ;;
    f)  cd /var/DAFi-gating
        git pull

        cp /var/DAFi-gating/GP/fiona_resources/* /opt/genepattern/resources
	echo "copying FIONA resources from github"
        ;;
    g)	GPONLY=1
    esac
done

declare -a ARGS
for var in "$@"; do
    # Ignore known bad arguments
    if [ "$var" == '-f' ] || [ "$var" == '-u' ] || [ "$var" == '-g' ]; then
        continue
    fi
    ARGS[${#ARGS[@]}]="$var"
done

# Syncing DAFi scripts and notebooks from Github repo
cd $HOME
if [ -d "work" ]; then
        cp /var/DAFi-gating/Notebooks/*.ipynb ~/work
else
        cp /var/DAFi-gating/Notebooks/*.ipynb $HOME
fi

if [[ ! -z "${JUPYTERHUB_API_TOKEN}" ]]; then
  # launched by JupyterHub, use single-user entrypoint
	if [[ "$GPONLY" -eq "0" ]]; then
		nohup /opt/genepattern/StartGenePatternServer &

		exec /usr/local/bin/start-singleuser.sh ${ARGS[@]}
	else
		echo "Running in GenePattern only mode"
		/opt/genepattern/StartGenePatternServer
	fi


else
  if [[ ! -z "${JUPYTER_ENABLE_LAB}" ]]; then
    . /usr/local/bin/start.sh jupyter lab ${ARGS[@]}
  else
	if [[ "$GPONLY" -eq "0" ]]; then
		nohup /opt/genepattern/StartGenePatternServer &
		. /usr/local/bin/start.sh jupyter notebook ${ARGS[@]}
	else
		echo "Running in GenePattern only mode"
		/opt/genepattern/StartGenePatternServer
	fi
  fi
fi

cd $HOME
