#!/bin/bash

export VERBOSITY=${VERBOSITY:-2}
export LOGFILE=/tmp/iip.out

PORT=7000
COUNTER=0
echo $NB_IIP_PROCESS
while [[ $COUNTER -lt $NB_IIP_PROCESS ]]; do
    echo "spawn process"
    spawn-fcgi -f /opt/iipsrv/src/iipsrv.fcgi -p $(($PORT+$COUNTER))
    let COUNTER=COUNTER+1
done
