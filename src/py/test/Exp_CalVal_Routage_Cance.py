
import numpy as np
import matplotlib.pyplot as plt

import gamma
import smash
import scipy
import functions_smash_plot
import functions_smash_stats

#we can optimize parameters or initial states
optimize_parameter=True
optimize_states=False

#EXPERIENCE DE BASE CANCE :
#Calibration SMASH gra
#Setting SMASH grg
#Calibrated only parameter gamma

#Smash reference model structure gr-a 
setup_cance,mesh_cance=smash.load_dataset('Cance')
setup_cance.update({"save_qsim_domain":True})
setup_cance["structure"]="gr-a"
model_smash_gr = smash.Model(setup_cance, mesh_cance)
model_smash_gr.states.hp=0.1
model_smash_gr.run(inplace=True)

#Here smash parameter are normalized
model_smash_gr_calibrated = model_smash_gr.optimize(
        mapping="distributed",
        algorithm="l-bfgs-b",
        options={"maxiter": 15},
        control_vector=['cp','cft','lr'],
        bounds={"cp":[0.1,1000.],"cft":[0.1,1000],"lr":[1.,100.]}
    )



#Smash reference model structure gr-ag
setup_cance,mesh_cance=smash.load_dataset('Cance')
setup_cance.update({"save_qsim_domain":True})
setup_cance["structure"]="gr-g"
model_smash_gr_g = smash.Model(setup_cance, mesh_cance)
model_smash_gr.states.hp=0.1
model_smash_gr_g.parameters=model_smash_gr_calibrated.parameters.copy()
model_smash_gr_g.run(inplace=True)


#Smash Gamma model structure gr-g (manual settings)
model_gamma_manual_settings=gamma.smashplug.ConfigureGammaWithSmash(model_smash_gr_g,dt=900,vmin=0.1,vmax=10.,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=0.0,velocity_computation="qm3",varying_spread=1,criteria="nse",ponderation_cost=1000.,spreading_uniform=1,pdt_start_optim=1600)

#Set default parameter
model_gamma_manual_settings.routing_parameters.hydraulics_coefficient=0.5
model_gamma_manual_settings.routing_parameters.spreading=1.0
model_gamma_manual_settings.routing_mesh.controlled_nodes[1:3]=0

#Get the inflow from smash
Inflows=gamma.smashplug.GetGammaInflowFromSmash(model_smash_gr_g,dt=900)
model_gamma_manual_settings.run(Inflows)

#get observed discharges:
GammaGriddedObservation=gamma.smashplug.SmashDataVectorsToGrid(model_smash_gr_g.input_data.qobs,model_gamma_manual_settings,smash_dt=model_smash_gr_g.setup.dt)

#optimize parameters
#Here smash parameter are not normalized, we use ScaleGradientsByBounds and ScaleGammaGradients instead !
BestControlVector,optimized_smash_model_cance,optimized_gamma_model_cance=gamma.smashplug.OptimizeCoupledModel(
        model_smash_gr_g,
        model_gamma_manual_settings,
        GammaGriddedObservation,
        control_parameters_list=["hydraulics_coefficient","spreading"],
        bounds={"hydraulics_coefficient":[0.3,5.0],"spreading":[0.5,3.]},
        maxiter=15,
        tol=0.00001,
        ScaleGradientsByBounds=True,
        ScaleGammaGradients=True
        )

# Get the abscisse for easy plotting
abcisses=gamma.smashplug.GetCorrespondingTimeStepRange(optimized_smash_model_cance,optimized_gamma_model_cance)
fig,ax=plt.subplots()
plt.rc('font', size=18) 
ax.plot(abcisses["ts_gamma"],model_gamma_manual_settings.routing_results.discharges[:,model_gamma_manual_settings.routing_mesh.controlled_nodes[0]-1],label="SMASH-grg-Gamma (Ebauche)",lw=2)

ax.plot(abcisses["ts_gamma"],optimized_gamma_model_cance.routing_results.discharges[:,optimized_gamma_model_cance.routing_mesh.controlled_nodes[0]-1],label="SMASH-grg-Gamma (Calibrated)",lw=2) #nodes are in Fortran index: in python array we need -1: how to fix that ?

ax.plot(abcisses["ts_smash"],model_smash_gr_calibrated.output.qsim[0,:],label="SMASH-gra (Calibrated)",lw=2)

ax.plot(abcisses["ts_smash"],model_smash_gr_calibrated.input_data.qobs[0,:],label="Observations",color='0',lw=2)

ax.legend(loc='upper left')
ax.axes.grid(True,alpha=.7, ls="--")
ax.axes.set_xlabel("Time-Step (hours)")
ax.axes.set_ylabel("Discharges (m^3/s)")
fig.show()
plot=(fig,ax)

functions_smash_plot.save_figure(plot,figname="ExpCalValRoutage.pdf",xlim=[0,1500],ylim=[0,350],xsize=24,ysize=10)
functions_smash_plot.save_figure(plot,figname="ExpCalValRoutage_zoom1.pdf",xlim=[550,725],ylim=[0,250],xsize=12,ysize=10)
functions_smash_plot.save_figure(plot,figname="ExpCalValRoutage_zoom2.pdf",xlim=[1150,1400],ylim=[0,350],xsize=12,ysize=10)


#Convert Gazmma vector parameter to gridded data
Grid_Hydraulic_Coef=gamma.smashplug.GammaVectorsToSmashGrid(optimized_gamma_model_cance.routing_parameters.hydraulics_coefficient,optimized_smash_model_cance,optimized_gamma_model_cance)
Spreading=gamma.smashplug.GammaVectorsToSmashGrid(optimized_gamma_model_cance.routing_parameters.spreading,optimized_smash_model_cance,optimized_gamma_model_cance)


functions_smash_plot.plot_image(Grid_Hydraulic_Coef,figname="ExpCalValRoutage_CoefHydraulics.pdf",mask=optimized_smash_model_cance.mesh.active_cell,title="Hydraulic coefficients",title_font_size=14)

functions_smash_plot.plot_image(Spreading,figname="ExpCalValRoutage_Spreading.pdf",mask=optimized_smash_model_cance.mesh.active_cell,title="Spreading coefficients",title_font_size=14)

functions_smash_plot.plot_image(model_smash_gr_calibrated.parameters.cp,figname="ExpCalValRoutage_Cp.pdf",mask=optimized_smash_model_cance.mesh.active_cell,title="Capacities of the production reservoir",title_font_size=14)

functions_smash_plot.plot_image(model_smash_gr_calibrated.parameters.cft,figname="ExpCalValRoutage_Cft.pdf",mask=optimized_smash_model_cance.mesh.active_cell,title="Capacities of the transfert reservoir",title_font_size=14)


#get observed discharges:
GammaGriddedObservation=gamma.smashplug.SmashDataVectorsToGrid(optimized_smash_model_cance.input_data.qobs,optimized_gamma_model_cance,smash_dt=optimized_smash_model_cance.setup.dt)

#Compute the nse
nse_gra=functions_smash_stats.nse(model_smash_gr_calibrated.input_data.qobs[0,:],model_smash_gr_calibrated.output.qsim[0,:],lacuna=-99.)
nse_grg=functions_smash_stats.nse(GammaGriddedObservation[:,optimized_gamma_model_cance.routing_mesh.controlled_nodes[0]-1],optimized_gamma_model_cance.routing_results.discharges[:,optimized_gamma_model_cance.routing_mesh.controlled_nodes[0]-1],lacuna=-99.)

print("Calage NSE_gra=",nse_gra)
print("Calage NSE_grg=",nse_grg)




#validation:
#New Smash Setup
setup_cance,mesh_cance=smash.load_dataset('Cance')
setup_cance.update({"save_qsim_domain":True})
setup_cance["structure"]="gr-g"
setup_cance["start_time"]='2013-09-01 00:00'
setup_cance["end_time"]='2014-04-01 00:00'
setup_cance["pet_directory"]="/home/maxime/DATA/ETP-SFR-FRA-INTERA_L93/"
setup_cance["prcp_directory"]="/home/maxime/DATA/PLUIE/"
model_smash_grg_validation = smash.Model(setup_cance, mesh_cance)

#New Coupling setup
setup_cance["structure"]="gr-a"
setup_cance.update({"save_qsim_domain":False})
model_smash_gra_validation = smash.Model(setup_cance, mesh_cance)

#set parameters and states
model_smash_grg_validation.parameters=model_smash_gr_calibrated.parameters.copy()
model_smash_gra_validation.parameters=model_smash_gr_calibrated.parameters.copy()
#Set previous states for restart at last time step
#model_smash_grg_validation.states=optimized_smash_model_cance.outputs.fstates.copy()
#model_smash_gra_validation.states=optimized_smash_model_cance.outputs.fstates.copy()
model_smash_grg_validation.states.hp=0.1
model_smash_gra_validation.states.hp=0.1

model_gamma_validation=gamma.smashplug.ConfigureGammaWithSmash(model_smash_grg_validation,dt=900,vmin=0.1,vmax=10.,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=0.0,velocity_computation="qm3")

#Set optimized parameter
model_gamma_validation.routing_parameters_change(hydraulics_coefficient=optimized_gamma_model_cance.routing_parameters.hydraulics_coefficient.copy(),spreading=optimized_gamma_model_cance.routing_parameters.spreading.copy())

#set previous states
# ~ #Attention au shape des states qui sont différents, prendre les states correspondant aux shape du nouveau modèle. Les shapes des 2 modèle sont différents car l'un provient d'un calage et les la longueur de la fenêtre mémoire est calculé à partir des bornes ou de la valeur des paramètres si celle-ci n'est pas spécifiée
# ~ shapeOfstates=model_gamma_validation.routing_states.states.shape
# ~ model_gamma_validation.routing_states.states=optimized_gamma_model_cance.routing_states.states[0:shapeOfstates[0],:].copy()
# ~ model_gamma_validation.routing_states.remainder=optimized_gamma_model_cance.routing_states.remainder[0:shapeOfstates[0],:].copy()

#run the coupled model
gamma.smashplug.RunCoupledModel(model_smash_grg_validation,model_gamma_validation)

model_smash_gra_validation.run(inplace=True)

abcisses=gamma.smashplug.GetCorrespondingTimeStepRange(model_smash_grg_validation,model_gamma_validation)
fig,ax=plt.subplots()
plt.rc('font', size=18) 
ax.plot(abcisses["ts_gamma"],model_gamma_validation.routing_results.discharges[:,model_gamma_validation.routing_mesh.controlled_nodes[0]-1],label="SMASH-grg-Gamma (Validation)",lw=2) #nodes are in Fortran index: in python array we need -1: how to fix that ?

ax.plot(abcisses["ts_smash"],model_smash_gra_validation.output.qsim[0,:],label="SMASH-gra (Validation)",lw=2)

ax.plot(abcisses["ts_smash"],model_smash_grg_validation.input_data.qobs[0,:],label="Observations",color='0',lw=2)

ax.legend(loc='upper left')
ax.axes.grid(True,alpha=.7, ls="--")
ax.axes.set_xlabel("Time-Step (hours)")
ax.axes.set_ylabel("Discharges (m^3/s)")
fig.show()
plot=(fig,ax)

functions_smash_plot.save_figure(plot,figname="ExpCalValRoutage_validation_period.pdf",xlim=[0,5000],ylim=[0,100],xsize=24,ysize=10)
functions_smash_plot.save_figure(plot,figname="ExpCalValRoutage_validation_period_zoom1.pdf",xlim=[150,250],ylim=[0,20])
functions_smash_plot.save_figure(plot,figname="ExpCalValRoutage_validation_period_zoom2.pdf",xlim=[2760,2860],ylim=[0,100])



#get observed discharges:
GammaGriddedObservation=gamma.smashplug.SmashDataVectorsToGrid(model_smash_grg_validation.input_data.qobs,model_gamma_validation,smash_dt=model_smash_gr_g.setup.dt)

#compute the nse
nse_gra=functions_smash_stats.nse(model_smash_grg_validation.input_data.qobs[0,:],model_smash_gra_validation.output.qsim[0,:],lacuna=-99.)
nse_grg=functions_smash_stats.nse(GammaGriddedObservation[:,model_gamma_validation.routing_mesh.controlled_nodes[0]-1],model_gamma_validation.routing_results.discharges[:,model_gamma_validation.routing_mesh.controlled_nodes[0]-1],lacuna=-99.)

print("NSE_gra=",nse_gra)
print("NSE_grg=",nse_grg)
