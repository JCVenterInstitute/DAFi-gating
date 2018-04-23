#!/bin/bash
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

set -e

cd /tmp/DAFi-gating
git pull
mkdir ~/Notebooks
cp Notebooks/*.ipynb ~/Notebooks

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

