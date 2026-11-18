#!/bin/bash
memory_limit=0.75

until [ $# = 0 ]
    do
        case $1 in
            -stop)
            shift
            if [ -e ./status_task_killer.txt ] ; then
                task_killer_pid=$(grep "PID:" ./status_task_killer.txt | cut -d : -f 2)
                kill $task_killer_pid
                rm ./status_task_killer.txt
                echo "task_killer killed"
                exit
            else
                echo "task_killer is not running."
            fi
            ;;
            -memory_limit)
            shift
            memory_limit=$1
            shift
            ;;
            -launcher_path)
            shift
            launcher_path=$1
            shift
            ;;
            *)
            echo "Bad input parameters $1."
            shift
            ;;
        esac
done

if [ ! -e ./status_task_killer.txt ] ; then
    echo "PID:"$$ > ./status_task_killer.txt
else
    echo "task_killer already running."
    exit
fi

     
while [ ! -e ./stop_task_killer ] ;
do

    free_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $7}')
    total_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $2}')
    #syst_load=$(echo "scale=2;$free_mem/$total_mem" | bc)
    
    #process_mem_total=$(ps aux | awk -v proc_pat="$process_pattern" '$11 ~ proc_pat {sum += $6} END {printf "%.2f", sum/1024}')

    mapfile -t PIDS_LIST < <(awk 'NF {print $1}' ./running_process.txt)
    if [ "${#PIDS_LIST[@]}" -gt 0 ] ; then
        process_mem_total=$(ps -o rss= -p "$(IFS=,; echo "${PIDS_LIST[*]}")" | awk '{sum += $1} END {printf "%.2f", sum/1024}')
    fi
    if [ -z $process_mem_total ] ; then
        process_mem_total=0
    fi
    echo "Process total memory $process_mem_total."
    lim=$(echo "$memory_limit*$free_mem" | bc)
    echo $lim
    
    if [ "$( echo "$process_mem_total > $memory_limit*$free_mem" | bc)" -eq 1 ] ; then
    #if [ 1 -eq '$(echo "$syst_load < $mem_threshold" | bc)' ] ;
        
        echo "Total process memory load upper the memory limit. We kill the last process to free memory... "
        
        last_pid=$(tail -n 1 ./running_process.txt | awk '{print $1}')
        last_catchment=$(tail -n 1 ./running_process.txt | awk '{print $2}')
        
        echo $last_pid $last_catchment > ./killed_process.txt
        
        # tue un processus et actualise running process
        kill $last_pid
        sed -i "/${last_catchment}/d" ./running_process.txt
        #cat running_process.txt killed_process.txt | uniq -u > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/running_process.txt
        
        #echo $last_catchment >> /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt
        
        #echo "" > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/killed_process.txt
        
        rm "$launcher_path/smash_gamma_$last_catchment"
        
        #update proc lim for the task manager.
        #nb_cur_proc=$(cat ./running_process.txt | wc -l)
        nb_cur_proc=$(ps -p "$(awk '{printf "%s%s", sep, $1; sep=","}' ./running_process.txt)" -o pid=,comm= | wc -l)
        sed -i "s/^nb_proc_max:.*/nb_proc_max:$nb_cur_proc" ./status_task_manager.txt
        echo "Limit the number of process to $nb_cur_proc"
        
        echo -e "Process smash_gamma_$last_catchment with $last_pid killed.\n"
    fi
    
    sleep 2
    
done

