#!/bin/bash

#catchment=$1

until [ $# = 0 ]
    do
        case $1 in
            -catchment)
            shift
            catchment=$1
            shift
            ;;
            -bin)
            shift
            bin=$1
            shift
            ;;
            -script)
            shift
            script=$1
            shift
            ;;
            *)
            echo "Bad input parameters $1."
            shift
            ;;
        esac
done

echo $catchment

source $HOME/python_venv/smash/bin/activate

exec $bin $script -m "init" -c "$catchment" -s "s1"


#~ echo $catchment >> succeed_process.txt
