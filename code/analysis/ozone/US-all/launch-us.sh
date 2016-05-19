#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
OZONE=$HOME/repos-git/spatial-skew-t/code/analysis/ozone/US-all
L=$1

$RBIN CMD BATCH --vanilla --no-save $OZONE/us-all-$L.R $OZONE/us-all-$L.out 2>&1
exit 0
