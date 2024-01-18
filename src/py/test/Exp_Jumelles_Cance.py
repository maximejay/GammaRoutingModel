
import numpy as np
import matplotlib.pyplot as plt

import gamma
import smash
import scipy
import functions_smash_plot
import functions_smash_stats



#EXPERIENCE JUMELLE
optimize_parameter=True
optimize_states=False

#Smash Gamma model structure gr-g 
setup_cance,mesh_cance=smash.load_dataset('Cance')
setup_cance.update({"save_qsim_domain":True})
setup_cance["structure"]="gr-g"
smash_model = smash.Model(setup_cance, mesh_cance)

#Set states and parameters
smash_model.states.hp=0.1
smash_model.parameters.cp=200.
smash_model.parameters.cft=100.

#direct run of SMash model
smash_model.run(inplace=True)

#Build the coupled Gamma model depending the smash configuration. Add any optionnal arguments
model_gamma=gamma.smashplug.ConfigureGammaWithSmash(smash_model,dt=900,vmin=0.1,vmax=10.,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=0.0,velocity_computation="qm3",varying_spread=1,criteria='rmse',ponderation_cost=1.,spreading_uniform=1)

#Set default parameter. This is not the best method. We must use the method model_gamma.routing_parameters_change() for that ! After setting parameters we must update the states of the model, but model_gamma.run() do it for you by default.
model_gamma.routing_parameters.hydraulics_coefficient=1.0
model_gamma.routing_parameters.spreading=1.

#get the inflows
Inflows=gamma.smashplug.GetGammaInflowFromSmash(smash_model,dt=900)

#run the direct model (update states before)
model_gamma.run(Inflows)


#Get observation for twin experiment
observations=np.zeros(shape=model_gamma.routing_results.discharges.shape)
observations[:,:]=model_gamma.routing_results.discharges[:,:]

#Set new parameter
smash_model.states.hp=0.1
smash_model.parameters.cp=100.
smash_model.parameters.cft=100.
model_gamma.routing_parameters.hydraulics_coefficient=0.5
model_gamma.routing_parameters.spreading=1.5

#runt the direct model
smash_model.run(inplace=True)
#Get the output of Smash => input of gamma
New_inflows=gamma.smashplug.GetGammaInflowFromSmash(smash_model,dt=900)
#run the direct gamma model with new inflows
model_gamma.run(New_inflows,states_init=0)
#Compute the cost
model_gamma.cost_function(observations,model_gamma.routing_results.discharges)

#Or use a direct call method
gamma.smashplug.RunCoupledModel(smash_model,model_gamma)
model_gamma.cost_function(observations,model_gamma.routing_results.discharges)

#define a control vector
if optimize_parameter:
    ControlVector=gamma.smashplug.VectorizeModelParameters(smash_model,model_gamma,control_parameters_list=["cp","hydraulics_coefficient","spreading"])

#define a control vector
if optimize_states:
    ControlVector=gamma.smashplug.VectorizeModelParameters(smash_model,model_gamma,control_parameters_list=["hp"])

#Compute the gradient
gradient=gamma.smashplug.ComputeModelGradients(ControlVector, smash_model, model_gamma, New_inflows, observations)

#Run the the direct model, compute the cost and the gradients
cost,gradient=gamma.smashplug.ComputeCostAndGradients(ControlVector['X'],ControlVector,smash_model,model_gamma,observations)

# ~ #Optimize the distributed parameters
if optimize_parameter:
    BestControlVector,optimized_smash,optimized_gamma=gamma.smashplug.OptimizeCoupledModel(
        smash_model,
        model_gamma,
        observations,
        control_parameters_list=["cp","hydraulics_coefficient","spreading"],
        bounds={"cp":[0.1,1000.],"hydraulics_coefficient":[0.3,5.],"spreading":[0.5,3.]},
        maxiter=20,
        tol=0.00001,
        ScaleGradientsByBounds=True,
        ScaleGammaGradients=True
        )

#Optimize the distributed states of Smash
if optimize_states:
    BestControlVector,optimized_smash,optimized_gamma=gamma.smashplug.OptimizeCoupledModel(
        smash_model,
        model_gamma,
        observations,
        control_parameters_list=["hp"],
        bounds={"hp":[0.01,1.]},
        maxiter=20,
        tol=0.00001
        )

abcisses=gamma.smashplug.GetCorrespondingTimeStepRange(smash_model,model_gamma)
fig,ax=plt.subplots()
plt.rc('font', size=18) 
ax.plot(abcisses["ts_gamma"],observations[:,model_gamma.routing_mesh.controlled_nodes[0]-1],label="Reference SMASH-gr-g-Gamma (Reference)",color='0',lw=2)
ax.plot(abcisses["ts_gamma"],optimized_gamma.routing_results.discharges[:,model_gamma.routing_mesh.controlled_nodes[0]-1],label="SMASH-gr-g-Gamma (Calibrated)",lw=2) #nodes are in Fortran index: in python array we need -1: how to fix that ?
ax.plot(abcisses["ts_gamma"],model_gamma.routing_results.discharges[:,model_gamma.routing_mesh.controlled_nodes[0]-1],label="SMASH-gr-g-Gamma (Background)",lw=2)
ax.legend(loc='upper left')
ax.axes.grid(True,alpha=.7, ls="--")
ax.axes.set_xlabel("Time-Step (hours)")
ax.axes.set_xlabel("Discharges (m^3/s)")
fig.show()
plot=(fig,ax)

functions_smash_plot.save_figure(plot,figname="ExpJumelle.pdf",xlim=[0,1500],ylim=[0,600],xsize=24,ysize=10)
functions_smash_plot.save_figure(plot,figname="ExpJumelle_zoom1.pdf",xlim=[550,725],ylim=[0,450])
functions_smash_plot.save_figure(plot,figname="ExpJumelle_zoom2.pdf",xlim=[1150,1400],ylim=[0,600])


#analyse des paramètres obtenus
np.mean(optimized_smash.parameters.cp,where=(optimized_smash.mesh.active_cell==1))
if optimize_parameter:
    np.mean(optimized_gamma.routing_parameters.hydraulics_coefficient)
    np.mean(optimized_gamma.routing_parameters.spreading)

#Convert Gamma vector parameter to gridded data
Grid_Hydraulic_Coef=gamma.smashplug.GammaVectorsToSmashGrid(optimized_gamma.routing_parameters.hydraulics_coefficient,smash_model,model_gamma)
Spreading=gamma.smashplug.GammaVectorsToSmashGrid(optimized_gamma.routing_parameters.spreading,smash_model,model_gamma)

functions_smash_plot.plot_image(Grid_Hydraulic_Coef,figname="ExpJumelle_CoefHydraulics.pdf",mask=smash_model.mesh.active_cell,title="Hydraulic coefficients",title_font_size=14)
functions_smash_plot.plot_image(Spreading,figname="ExpJumelle_Spreading.pdf",mask=smash_model.mesh.active_cell,title="Spreading coefficients",title_font_size=14)
functions_smash_plot.plot_image(optimized_smash.parameters.cp,figname="ExpJumelle_Cp.pdf",mask=smash_model.mesh.active_cell,title="Capacities of the production reservoir",title_font_size=14)



#Other useful Functions for managing obs/sim discharges between smash and gamma:

#convert gridded discharges from Gamma to smash discharges vector 
SmashDischargesVectors=gamma.smashplug.GriddedDataToSmashVectors(model_gamma.routing_results.discharges,model_gamma,smash_dt=3600)
#convert Vector discharges from Smash to Gamma gridded discharges 
GammaGriddedObservation=gamma.smashplug.SmashDataVectorsToGrid(smash_model.input_data.qobs,model_gamma,smash_dt=3600)
#Interpolate observations from smash
InterpolatedObservationsAsVector=gamma.smashplug.InterpolatedObservationsAtGauge(smash_model,900)
#interpolate smash qobs to gamma dt
InterpolateSmashQobsFor1GaugeToGammaDT=gamma.smashplug.FitDataToNewDt(smash_model.input_data.qobs[0,:],3600,900)
#interpolate Gamma qsim to smash dt
InterpolateQsimAt1CellToSmashDT=gamma.smashplug.FitDataToNewDt(model_gamma.routing_results.discharges[:,model_gamma.routing_mesh.controlled_nodes[0]-1],900,3600)





