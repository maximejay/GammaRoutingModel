#!/bin/bash
memory_threshold=0.5 # fraction x of the current free mem; free mem = available mem + current managed process mem. At the begening, there is no managed process. Above this threshold no more task will be launched until the memory lower again.
memory_limit=0.75 # fraction x of the current free mem ; free mem = available mem + current managed process mem. At the begening, there is no managed process. Above this threshold, the last task will be killed to free memory.
process_pattern="smash_manager"
launcher_prog="/home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/launcher.sh"
user_nb_proc_max=5 #Max number of parrallel task to run
freq=10 #Frequency duration in second before launching a new task.
queued_process=./queued_process.txt #File with input args in coluumn

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
            -freq)
            shift
            freq=$1
            shift
            ;;
            -process_alive_info)
            shift
            if [ -e ./running_process.txt ] & [ $(cat ./running_process.txt | wc -l) -gt 0 ] ; then
                ps -p "$(awk '{printf "%s%s", sep, $1; sep=","}' ./running_process.txt)" -o pid=,comm=
            else
                echo "No running process."
            fi
            exit
            ;;
            -kill_all)
            shift
            while read PID TAG
            do
                MEM=$(ps -p "$PID" -o rss= 2>/dev/null)
                if [ -n "$MEM" ]; then
                    kill $PID
                    echo "Kill $PID"
                    sed -i "/${PID}/d" ./running_process.txt
                fi
            done < ./running_process.txt
            exit
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
    echo "nb_cur_proc:0" >> ./status_task_manager.txt
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

next_catchment=$(cat $queued_process <(awk '{print $2}' ./running_process.txt) | sort | uniq -u | sed '/^$/d' | head -n 1)

while [ ! -z $next_catchment ] ;
do
    echo "Process args $next_catchment"
    TOTAL_MEM=$(free -m | head -n +2 | tail -n 1 | awk '{print $2}')
    AVAILABLE_MEM=$(free -m | head -n +2 | tail -n 1 | awk '{print $4}')
    ALL_MEM=$(ps -eo rss= | awk '{sum+=$1} END {print sum/1024}')
    
    SELECTED_MEM=0
    while read PID TAG
    do
        MEM=$(ps -p "$PID" -o rss= 2>/dev/null)
        
        if [ -n "$MEM" ]; then
            SELECTED_MEM=$(echo "scale=2;$SELECTED_MEM + $MEM/1024" | bc)
        fi
    done < ./running_process.txt
    
    #RESERVED_MEM=$(echo "$TOTAL_MEM - $ALL_MEM + $SELECTED_MEM" | bc)
    RESERVED_MEM=$(echo "$AVAILABLE_MEM + $SELECTED_MEM" | bc)
    
    syst_load=$(echo "scale=2;$SELECTED_MEM/$RESERVED_MEM" | bc)
    
    echo "Total memory                        : $TOTAL_MEM Mo"
    echo "Consumed memory                     : $ALL_MEM Mo"
    echo "Free memory                         : $AVAILABLE_MEM Mo"
    echo "Current memory of monitored process : $SELECTED_MEM Mo"
    echo "Current claimed memory              : $RESERVED_MEM Mo"
    echo "Claimed memory load                 : $( echo "$syst_load*100" | bc)% / $( echo "$memory_threshold*100" | bc )%"
    
    nb_cur_proc=0
    if [ -e ./running_process.txt ] && [ $(cat ./running_process.txt | wc -l) -gt 0 ] ; then
        nb_cur_proc=$(ps -p "$(awk '{printf "%s%s", sep, $1; sep=","}' ./running_process.txt)" -o pid=,comm= | wc -l)
        sed -i "s/^nb_cur_proc:.*/nb_cur_proc:$nb_cur_proc/g" ./status_task_manager.txt
    fi
    
    nb_proc_max=$(grep "nb_proc_max:" ./status_task_manager.txt | cut -d : -f 2)
    
    mem_per_proc=0
    if [ "$nb_cur_proc" -gt 0 ] ; then
        mem_per_proc=$(echo "$SELECTED_MEM/$nb_cur_proc" | bc)
    fi
    
    new_nb_proc_max=0
    if [ "$(echo "$mem_per_proc > 0" | bc)" -eq 1 ] ; then
        new_nb_proc_max=$(echo "scale=0; $memory_limit*$RESERVED_MEM/$mem_per_proc" | bc)
    fi
    
    if [ "$new_nb_proc_max" -gt "$nb_proc_max" ] && [ "$new_nb_proc_max" -lt "$user_nb_proc_max" ] ; then
        echo "Raise the number of process limitation to $new_nb_proc_max"
        nb_proc_max=$new_nb_proc_max
        sed -i "s/^nb_proc_max=.*/nb_proc_max=$nb_proc_max/g" ./status_task_manager.txt
    fi

    if [ "$(echo "$syst_load < $memory_threshold" | bc)" -eq 1 ] && [ $nb_cur_proc -lt $nb_proc_max ]; then
        
        echo "System memory load lower than the memory limit threshold and number of process lower than the limitation. We launch a new process ... "
        
        alias_launcher_prog="${launcher_path}/${process_pattern}_${next_catchment}"

        if [[ -e "${alias_launcher_prog}" ]];
        then
            rm "${alias_launcher_prog}"
        fi
        ln -s $launcher_prog $alias_launcher_prog

        nohup ${alias_launcher_prog} $next_catchment > ${launcher_path}/nohup_${process_pattern}_${next_catchment}.log & PID=$!

        echo $PID ${next_catchment} >> ./running_process.txt
        
        echo -e "Process $alias_launcher_prog with PID '$PID' launched."
    else
        echo "Memory limit / limit number of process reached. Wait for some process to end."
    fi
    
    sleep $freq
    
    next_catchment=$(cat $queued_process <(awk '{print $2}' ./running_process.txt) | sort | uniq -u | head -n 1)
    
    echo -e "\n"
    
done

