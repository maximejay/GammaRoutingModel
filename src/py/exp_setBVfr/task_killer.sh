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

    TOTAL_MEM=$(free -m | head -n +2 | tail -n 1 | awk '{print $2}')
    ALL_MEM=$(ps -eo rss= | awk '{sum+=$1} END {print sum/1024}')
    AVAILABLE_MEM=$(free -m | head -n +2 | tail -n 1 | awk '{print $7}')
    
    #process_mem_total=$(ps aux | awk -v proc_pat="$process_pattern" '$11 ~ proc_pat {sum += $6} END {printf "%.2f", sum/1024}')
    SELECTED_MEM=0
    while read PID TAG
    do
        MEM=$(ps -p "$PID" -o rss= 2>/dev/null)
        
        if [ -n "$MEM" ]; then
            SELECTED_MEM=$(echo "scale=2;$SELECTED_MEM + $MEM/1024" | bc)
        fi
    done < ./running_process.txt
    
    RESERVED_MEM=$(echo "$AVAILABLE_MEM + $SELECTED_MEM" | bc)
    MEM_LIM=$(echo "$memory_limit*$RESERVED_MEM" | bc)
    
    echo "Total memory                        : $TOTAL_MEM Mo"
    echo "Consumed memory                     : $ALL_MEM Mo"
    echo "Free memory                         : $AVAILABLE_MEM Mo"
    echo "Current memory of monitored process : $SELECTED_MEM Mo"
    echo "Current claimed memory              : $RESERVED_MEM Mo"
    echo "Max claimed memory threshold        : $MEM_LIM Mo"
    
    if [ "$( echo "$SELECTED_MEM > $MEM_LIM" | bc)" -eq 1 ] ; then
    #if [ 1 -eq '$(echo "$syst_load < $mem_threshold" | bc)' ] ;
        
        echo "Monitored process memory upper the max memory threshold. We kill the last process to free memory ... "
        
        last_pid=$(tail -n 1 ./running_process.txt | awk '{print $1}')
        last_catchment=$(tail -n 1 ./running_process.txt | awk '{print $2}')
        
        echo $last_pid $last_catchment >> ./killed_process.txt
        
        # tue un processus et actualise running process
        kill $last_pid
        sed -i "/${last_catchment}/d" ./running_process.txt
        
        rm "$launcher_path/alias_bin/smash_gamma_$last_catchment"
        
        #update proc lim for the task manager.
        if [ -e ./running_process.txt ] && [ $(cat ./running_process.txt | wc -l) -gt 0 ] ; then
            nb_cur_proc=$(ps -p "$(awk '{printf "%s%s", sep, $1; sep=","}' ./running_process.txt)" -o pid=,comm= | wc -l)
            sed -i "s/^nb_proc_max:.*/nb_proc_max:$nb_cur_proc/g" ./status_task_manager.txt
            echo "Limit the number of process to $nb_cur_proc"
        fi
        
        echo -e "Process smash_gamma_$last_catchment with $last_pid killed."
    fi
    
    echo -e "\n"
    
    sleep 2
    
done

