mem_threshold=0.2
launcher_prog="/home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/launcher.sh"
ls /home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh | grep -E -o "[A-Z][0-9]{1,7}" > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt

free_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $7}')
total_mem=$(free -m | head -n +2 | tail -n 1 | awk '{print $2}')
syst_load=$(echo "scale=2;$free_mem/$total_mem" | bc)

#launch the task killer
nohup bash /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/task_killer.sh > log_task_killer.log &

next_catchment=$(cat /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt <(awk '{print $1}' /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/running_process.txt) | uniq -u | head -n 1)

while [ -z $next_catchment ] :
do

    if [ 1 -eq '$(echo "$syst_load > $mem_threshold" | bc)' ] ;
        
        echo "System memory load lower than the memory limit threshold. We launch a new process ... "
        
        alias_launcher_prog="/home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/smash_gamma_${next_catchment}"

        if [[ -e "${alias_launcher_prog}" ]];
        then
            rm "${alias_launcher_prog}"
        fi
        ln -s $launcher_prog $alias_launcher_prog

        nohup ${alias_launcher_prog} "init" $catchment "s1" > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/nohup_log.log & PID=$!

        echo $PID ${next_catchment} >> /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/running_process.txt
        
        sed -i "/${catchment}/d" /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt
        #cat queued_process.txt running_process.txt | uniq -u > /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt
        
        echo "Process $alias_launcher_prog with $PID launched."
        
    fi
    
    sleep 60
    
    next_catchment=$(cat /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/queued_process.txt <(awk '{print $1}' /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/running_process.txt) | uniq -u | head -n 1)
    
done

touch /home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/stop_task_killer
