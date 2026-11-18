#!/bin/bash

catchment=$1
echo $catchment
source ~/python_venv/smash/bin/activate
#python /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/exp_setbvfr.py -m "init" -c $catchment -s "s1"

#exec python /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/exp_setbvfr.py -m "init" -c "$catchment" -s "s1"
exec python /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/exp_setbvfr.py -m "init" -c "$catchment" -s "s1"
