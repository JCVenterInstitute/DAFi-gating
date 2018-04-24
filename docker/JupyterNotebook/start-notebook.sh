#!/bin/bash
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

set -e

# Syncing DAFi scripts and notebooks from Github repo
cd /var/DAFi-gating
git pull
cd $HOME
if [ -d "work" ]; then
        cp /var/DAFi-gating/Notebooks/*.ipynb ~/work
else
        cp /var/DAFi-gating/Notebooks/*.ipynb $HOME
fi

if [[ ! -z "${JUPYTERHUB_API_TOKEN}" ]]; then
  # launched by JupyterHub, use single-user entrypoint
  exec /usr/local/bin/start-singleuser.sh $*
else
  if [[ ! -z "${JUPYTER_ENABLE_LAB}" ]]; then
    . /usr/local/bin/start.sh jupyter lab $*
  else
    . /usr/local/bin/start.sh jupyter notebook $*
  fi
fi

cd $HOME
