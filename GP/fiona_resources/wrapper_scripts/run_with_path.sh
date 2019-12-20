#!/bin/bash
export XDG_CACHE_HOME=/tmp/.cache/
export LD_LIBRARY_PATH=/opt/conda/lib
export PATH=/var/DAFi-gating/inst/bin:/var/DAFi-gating/FCSTrans:/var/DAFi-gating/Python:/var/DAFi-gating/inst/bin:/var/DAFi-gating/FCSTrans:/var/DAFi-gating/Python:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

$@
