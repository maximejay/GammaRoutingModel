#!/bin/bash
memory_threshold=0.7
memory_limit=0.75
process_pattern="smash_manager"
launcher_prog="/home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/launcher.sh"
user_nb_proc_max=10
queued_process=./queued_process.txt

until [ $# = 0 ]
    do
        case $1 in
            -stop)
            shift
            if [ -e ./status_task_manager.txt ] ; then
                task_manager_pid=$(grep "PID:" ./status_task_manager.txt | cut -d : -f 2)
                kill $task_manager_pid
                rm ./status_task_manager.txt
                echo "task_manager killed"
                
                if [ -e ./status_task_killer.txt ] ; then
                    task_killer_pid=$(grep "PID:" ./status_task_killer.txt | cut -d : -f 2)
                    kill $task_killer_pid
                    rm ./status_task_killer.txt
                    echo "task_killer killed"
                    exit
                else
                    echo "task_killer is not running."
                fi
                
                exit
            else
                echo "task_manager is not running."
            fi
            ;;
            -script)
            shift
            launcher_prog=$1
            shift
            ;;
            -memory_threshold)
            shift
            memory_threshold=$1
            shift
            ;;
            -memory_limit)
            shift
            memory_limit=$1
            shift
            ;;
            -process_pattern)
            shift
            process_pattern=$1
            shift
            ;;
            -nb_proc_max)
            shift
            user_nb_proc_max=$1
            shift
            ;;
            -args_list)
            rm ./queued_process.txt
            shift
            until [ $# = 0 ]
            do
                if [ ${1::1} -eq "-" ] ; then
                    break
                fi
                echo $1 >> ./queued_process.txt
                shift
            done
            ;;
            -args_filelist)
            shift
            queued_process=$1
            shift
            ;;
            *)
            echo "Bad input parameters $1."
            shift
            ;;
        esac
done

if [ ! -e ./status_task_manager.txt ] ; then
    echo "PID:"$$ > ./status_task_manager.txt
    echo "nb_proc_max:$user_nb_proc_max" >> ./status_task_manager.txt
else
    echo "task_manager already running."
    exit
fi

launcher_path=$(dirname launcher_prog)

#
if [ ! -e $queued_process ] ; then 
    echo "$queued_process does not exist"
    #ls /home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh | grep -E -o "[A-Z][0-9]{1,7}" > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt
fi

#launch the task killer
nohup bash ./task_killer.sh -launcher_path $launcher_path -memory_limit $memory_limit > log_task_killer.log &

next_catchment=$(cat $queued_process <(awk '{print $2}' ./running_process.txt) | sort | uniq -u | head -n 1)

while [ ! -z $next_catchment ] ;
do
    echo "process args $next_catchment"
    free_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $7}')
    total_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $2}')
    syst_load=$(echo "scale=2;1-$free_mem/$total_mem" | bc)
    
    nb_cur_proc=$(ps -p "$(awk '{printf "%s%s", sep, $1; sep=","}' ./running_process.txt)" -o pid=,comm= | wc -l)
    #nb_cur_proc=$(cat ./running_process.txt | wc -l)
    nb_proc_max=$(grep "nb_proc_max:" ./status_task_manager.txt | cut -d : -f 2)
    process_mem_total=$(ps aux | awk -v proc_pat="$process_pattern" '$11 ~ proc_pat {sum += $6} END {printf "%.2f", sum/1024}')
    if [ -z $process_mem_total ] ; then
        process_mem_total=0
    fi
    
    if [ "$nb_cur_proc" -gt 0 ] ; then
        mem_per_proc=$(echo "$process_mem_total/$nb_cur_proc" | bc)
    else
        mem_per_proc=0
    fi
    if [ "$(echo "$mem_per_proc > 0" | bc)" -eq 1 ] ; then
        new_nb_proc_max=$(echo "scale=0; $memory_limit*$total_mem/$mem_per_proc" | bc)
    else
        new_nb_proc_max=0
    fi
    
    if [ "$new_nb_proc_max" -gt "$nb_proc_max" ] && [ "$new_nb_proc_max" -lt "$user_nb_proc_max" ] ; then
        echo "Raise the number of process limitation to $new_nb_proc_max"
        nb_proc_max=$new_nb_proc_max
        sed -i "s/nb_proc_max=.*/nb_proc_max=$nb_proc_max" ./status_task_manager.txt
    fi

    if [ "$(echo "$syst_load < $memory_threshold" | bc)" -eq 1 ] && [ $nb_cur_proc -lt $nb_proc_max ]; then
        
        echo "System memory load lower than the memory limit threshold. We launch a new process ... "
        
        alias_launcher_prog="${launcher_path}/${process_pattern}_${next_catchment}"

        if [[ -e "${alias_launcher_prog}" ]];
        then
            rm "${alias_launcher_prog}"
        fi
        ln -s $launcher_prog $alias_launcher_prog

        nohup ${alias_launcher_prog} $next_catchment > ${launcher_path}/nohup_${process_pattern}_${next_catchment}.log & PID=$!

        echo $PID ${next_catchment} >> ./running_process.txt
        
        #sed -i "/${next_catchment}/d" ./queued_process.txt
        #cat queued_process.txt running_process.txt | uniq -u > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt
        
        echo -e "Process $alias_launcher_prog with PID '$PID' launched.\n"
        
    fi
    
    sleep 2
    
    next_catchment=$(cat $queued_process <(awk '{print $2}' ./running_process.txt) | sort | uniq -u | head -n 1)
    
done

