#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
SIM=$HOME/repos-git/spatial-skew-t/code/analysis/simstudy

$RBIN CMD BATCH --vanilla --no-save $SIM/pot2-c.R $SIM/pot2-c.out 2>&1
exit 0
