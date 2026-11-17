mem_threshold=0.2

free_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $7}')
total_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $2}')
syst_load=$(echo "scale=2;$free_mem/$total_mem" | bc)


     
while [ ! -e /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/stop_task_killer ] :
do

    if [ 1 -eq '$(echo "$syst_load < $mem_threshold" | bc)' ] ;
        
        echo "System memory load upper than the memory limit threshold. We kill the last process to free memory... "
        
        last_pid=$(tail -n 1 /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/running_process.txt | awk '{print $1}')
        last_catchment=$(tail -n 1 /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/running_process.txt | awk '{print $2}')
        
        echo $last_pid $last_catchment > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/killed_process.txt
        
        # tue un processus et actualise running process
        kill $last_pid & echo $last_pid >> /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/killed_process.txt 
        cat running_process.txt killed_process.txt | uniq -u > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/running_process.txt
        
        echo $last_catchment >> /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt
        
        echo "" > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/killed_process.txt
        
        rm "/home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/smash_gamma_$last_catchment"
        
        echo "Process $last_alias_launcher with $last_pid killed."
    fi
    
    sleep 2
    
done

ps -u maxime | grep "$next_catchment"
