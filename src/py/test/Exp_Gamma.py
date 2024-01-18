import numpy as np
import matplotlib.pyplot as plt
import gamma

import smash
import functions_smash_plot
import functions_smash_stats

if 'model' in locals():
    del(model)

#Simple jonction test case for the Gamma routing model
#Here is the meshing o1-o4 are the nodes

    # ~ !Confluence test case
    # ~ !       o2
    # ~ !  o1  /
    # ~ !  \  /
    # ~ !   o3
    # ~ !   |
    # ~ !   o4
    # ~ !

model=gamma.Model()

#parametre du model
model.routing_setup_init(npdt=100,dt=900.,vmin=0.1,vmax=10.,mode_discretization_step=0.1,spreading_discretization_step=0.2,ponderation_regul=10000.0,velocity_computation="qm3",varying_spread=1)

#creation du maillage (Manual): Gamma have no meshing capabilities, you can use the interface to smash to build a mesh from a regular grid

#4 nodes and 2 maximum upstream nodes
model.routing_mesh_init(nb_nodes=4,nb_upstream_nodes=2)
#Nodes computation order (node number)
model.routing_mesh.upstream_to_downstream_nodes=np.array([1,2,3,4])
#Linker between one node with the upstream ones
model.routing_mesh.nodes_linker=0
model.routing_mesh.nodes_linker[:,2]=np.array([1,2]) #the third node (index=2) has upstream node 1 and 2 (index are different than node number ! This is due to the communication between fortran and Python)
model.routing_mesh.nodes_linker[1,3]=3 #the downstream node o4 targets o3

model.routing_mesh.surface=np.array([1,1,1,1]) #the drained surface in km2 à set to 1km (no influence). It is only important for the velocity computation vs Qspe
model.routing_mesh.dx=np.array([1000.0,7000.0,2000.0,1000.]) #Distance vers le noeud aval. A long bief slows the model (Tabulation of the Gamma routing coefficients is long because the delay is big...). For better performances, we should split the bief.

#The controlloed nodes where obeervations are available : only one downstream control : node 4
model.routing_mesh.controlled_nodes[0]=4

#Create the mesh
model.routing_mesh_update()

#Testing combinaison of differents routing parameters
Xi=[0.5,1.5] #hydraulics coef
S=[0.5,2] #spreading coef

for i in range(2):
    
    for j in range(2):
        
        param_xi=Xi[i]
        param_S=S[j]
        #Initilaise the parameters
        model.routing_parameters_init(hydraulics_coefficient=param_xi,spreading=param_S)
        
        #creating an array of inflow
        inflows=np.zeros(shape=(model.routing_setup.npdt,model.routing_mesh.nb_nodes))
        inflows[0:2,0]=2.0
        inflows[0:2,1]=4.0
        
        #run the model (direct run)
        model.run(inflows)
        
        fig,ax=plt.subplots()
        plt.rc('font', size=18) 
        ax.plot(model.routing_results.discharges[:,0],label="Gamma at node 01",lw=2) #nodes are in Fortran index: in python array we need -1: how to fix that ?
        ax.plot(model.routing_results.discharges[:,1],label="Gamma at node 02",lw=2) #nodes are in Fortran index: in python array we need -1: how to fix that ?
        ax.plot(model.routing_results.discharges[:,2],label="Gamma at node 03",lw=2) #nodes are in Fortran index: in python array we need -1: how to fix that ?
        ax.plot(model.routing_results.discharges[:,model.routing_mesh.controlled_nodes[0]-1],label="Gamma at node 04",lw=2) #nodes are in Fortran index: in python array we need -1: how to fix that ?
        
        ax.legend(loc='upper left')
        ax.axes.grid(True,alpha=.7, ls="--")
        ax.axes.set_xlabel("Time-Step (hours)")
        ax.axes.set_ylabel("Discharges (m^3/s)")
        ax.set_title(f"Hydraulics Coef={Xi[i]}, Spreading={S[j]}")
        fig.show()
        plot=(fig,ax)
        functions_smash_plot.save_figure(plot,figname=f"Exp_Gamma_{i}{j}.pdf",xlim=[0,40],xsize=12,ysize=10)



#Simple twin experiment:
#Store the simulated discharge in an array: consider that is the true discharges, i.e our observation vector
observations=np.zeros(shape=model.routing_results.discharges.shape)
observations[:,:]=model.routing_results.discharges[:,:]


fig,ax=plt.subplots()
ax.plot(observations)
fig.show()

fig,ax=plt.subplots()
ax.plot(model.routing_results.discharges)
fig.show()

fig,ax=plt.subplots()
ax.plot(model.routing_results.velocities)
fig.show()


#model.routing_parameters_change(hydraulics_coefficient=1.0,spreading=900.)
#changing parameters
model.routing_parameters.hydraulics_coefficient=0.7
model.routing_parameters.spreading=1.4

#run the model
model.run(inflows,states_init=0)

#compute the cost function between the "observation" and the new simulated discharges
model.cost_function(observations,model.routing_results.discharges)

cost_initial=model.routing_results.costs.copy()

#calibrate the parameters to fit the "observed" discharges
model.calibration(inflows,observations,states_init=0)

cost_final=model.routing_results.costs.copy()

fig,ax=plt.subplots()
ax.plot(model.routing_results.discharges)
fig.show()

fig,ax=plt.subplots()
ax.plot(model.routing_results.velocities)
fig.show()
