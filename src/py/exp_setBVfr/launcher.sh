#!/bin/bash

mode=$1
catchment=$2
scenario=$3

echo $mode $catchment $scenario

source ~/python_venv/smash/bin/activate

python /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/exp_setbvfr.py -m  $mode -c $catchment -s $scenario
