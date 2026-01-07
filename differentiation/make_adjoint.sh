#!/bin/bash
#This simple shell script will help you to differentiate the current GRD code ang get the adjoin

SOURCES=../src/f90/

#Create directories
mkdir backward/
mkdir forward/

#Copy the required code for diferentiation
cp $SOURCES/mod_routing_mesh.f90  $SOURCES/mod_routing_setup.f90 $SOURCES/mod_routing_parameters.f90 $SOURCES/mod_routing_states.f90 $SOURCES/mod_routing_results.f90 $SOURCES/mod_gamma_routing.f90 $SOURCES/run_forward.f90 $SOURCES/cost_function.f90 .

#$SOURCES/mod_gamma_function.f90

#Create the adjoin code
#make all
make dJ_dQin_adj
make dJ_dQin_tlm

#Modification du code
#sed -i "s/COMMON\_SIMU\_DIFF/COMMON\_SIMU\_DIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/MOD\_GAMMA\_ROUTING\_SETUP\_DIFF/MOD\_GAMMA\_ROUTING\_SETUP\_DIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/MOD\_GAMMA\_ROUTING\_MESH\_DIFF/MOD\_GAMMA\_ROUTING\_MESH\_DIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/MOD\_GAMMA\_ROUTING\_PARAMETERS\_DIFF/MOD\_GAMMA\_ROUTING\_PARAMETERS\_DIFF\_D/g" ./forward/TLM_d.f90
#sed -i "s/MODULE\_GAMMA\_FUNCTION\_DIFF/MODULE\_GAMMA\_FUNCTIONS\_DIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/MOD\_GAMMA\_ROUTING\_STATES\_DIFF/MOD\_GAMMA\_ROUTING\_STATES\_DIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/MOD\_GAMMA\_ROUTING\_DIFF/MOD\_GAMMA\_ROUTING\_DIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/ROUTING\_HYDROGRAM\_FORWARD\_NODIFF/ROUTING\_HYDROGRAM\_FORWARD\_NODIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/COST\_FUNCTION\_NODIFF/COST\_FUNCTION\_NODIFF\_D/g" ./forward/TLM_d.f90
sed -i "s/REGULARIZATION\_NODIFF/REGULARIZATION\_NODIFF\_D/g" ./forward/TLM_d.f90

sed -i "s/.*USE.*MOD\_GAMMA\_ROUTING\_PARAMETERS\_DIFF\_D/USE MOD\_GAMMA\_ROUTING\_PARAMETERS/g" ./forward/TLM_d.f90
sed -i "s/.*USE.*MOD\_GAMMA\_ROUTING\_PARAMETERS\_DIFF/USE  MOD\_GAMMA\_ROUTING\_PARAMETERS/g" ./backward/ADJ_b.f90


#Ajout le calul de la fonciton cout dans l adjoint : insertion d'une ligne avant l'autre (peu aussi se faire avec awk '/pattern/{print "somthing"}1' monfich)
sed -i 's/.*CALL.*COST\_FUNCTION\_B.*/ CALL COST\_FUNCTION\(routing_setup,routing_mesh,routing_parameter,observations,qnetwork,tab_cost,cost\)\n&/' ./backward/ADJ_b.f90

#Copy the adjoin code in the main directory 
cp ./backward/ADJ_b.f90 $SOURCES/AADJ_b.f90
cp ./forward/TLM_d.f90 $SOURCES/ATLM_d.f90


#remove copies of sources
rm mod_routing_mesh.f90 mod_routing_setup.f90 mod_routing_parameters.f90 mod_routing_states.f90 mod_routing_results.f90 mod_gamma_routing.f90 run_forward.f90 cost_function.f90

#Final message
echo "Differentiation done !"
echo "Check errors or warnings in the log files before to use the generated code !"
echo "Re-build your code  before use !"

