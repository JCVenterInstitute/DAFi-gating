#!/bin/bash
export XDG_CACHE_HOME=/tmp/.cache/
export LD_LIBRARY_PATH=/opt/conda/lib
export PATH=/var/DAFi-gating/inst/bin:/var/DAFi-gating/FCSTrans:/var/DAFi-gating/Python:/var/DAFi-gating/inst/bin:/var/DAFi-gating/FCSTrans:/var/DA
echo PATH IS $PATH

$@



