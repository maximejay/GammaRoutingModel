import numpy as np
import matplotlib.pyplot as plt

import gamma
import smash
import scipy
import functions_smash_plot
import functions_smash_stats


# Hydrological Experiment : Calibration + Validation, Comparison between Gr-G and Gr-a


# Smash reference model structure gr-a
setup_cance, mesh_cance = smash.load_dataset("Cance")
setup_cance.update({"save_qsim_domain": True})
setup_cance["structure"] = "gr-a"
model_smash_gr = smash.Model(setup_cance, mesh_cance)
model_smash_gr.states.hp = 0.1
model_smash_gr.run(inplace=True)

# Here smash parameter are normalized
model_smash_gr_calibrated = model_smash_gr.optimize(
    mapping="distributed",
    algorithm="l-bfgs-b",
    options={"maxiter": 15},
    control_vector=["cp", "cft", "lr"],
    bounds={"cp": [0.1, 1000.0], "cft": [0.1, 1000], "lr": [1.0, 100.0]},
)


# Smash Gamma model structure gr-g
setup_cance, mesh_cance = smash.load_dataset("Cance")
setup_cance.update({"save_qsim_domain": True})
setup_cance["structure"] = "gr-g"
smash_model = smash.Model(setup_cance, mesh_cance)
smash_model.states.hp = 0.1

# configure the Gamma model
model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
    smash_model,
    dt=900,
    vmin=0.1,
    vmax=10.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.1,
    ponderation_regul=0.0,
    velocity_computation="qm3",
    varying_spread=1,
    spreading_uniform=1,
    criteria="nse",
    ponderation_cost=10000.0,
    pdt_start_optim=1600,
)

# Set parameter
model_gamma.routing_parameters.hydraulics_coefficient = 0.5
model_gamma.routing_parameters.spreading = 1.0
model_gamma.routing_mesh.controlled_nodes[1:3] = 0

# direct run of the coupled model
gamma.smashplug.RunCoupledModel(smash_model, model_gamma)

# get observed discharges:
GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
    smash_model.input_data.qobs, model_gamma, smash_dt=smash_model.setup.dt
)

# compute the initial cost
initial_cost = model_gamma.cost_function(
    GammaGriddedObservation, model_gamma.routing_results.discharges
)

# Get the control vector
ControlVector = gamma.smashplug.VectorizeModelParameters(
    smash_model,
    model_gamma,
    control_parameters_list=["cp", "cft", "hydraulics_coefficient", "spreading"],
)
# Adding bounds to the control vector triggers the Normalisations of the gradients with respect to the bounds (It is added automatically in the optimize routine and it can be controlled with the flag ScaleGradientsByBounds)
ControlVector.update(
    {
        "bounds": {
            "cp": [0.1, 1000.0],
            "cft": [0.1, 1000.0],
            "hydraulics_coefficient": [0.3, 5.0],
            "spreading": [0.5, 3.0],
        }
    }
)

# get the inflows
GammaInflows = gamma.smashplug.GetGammaInflowFromSmash(smash_model, dt=900)

# Compute the gradient
gradient = gamma.smashplug.ComputeModelGradients(
    ControlVector, smash_model, model_gamma, GammaInflows, GammaGriddedObservation
)


# optimize parameters
# Here smash parameter are not normalized. That cause slow convergence. We try some gradient normalisation technics : ScaleGradientsByBounds and ScaleGammaGradients
# the spreading coefficient will uniformly calibrated => the gradient of this variable is scaled by the nb of nodes !
BestControlVector, optimized_smash_model, optimized_gamma_model = (
    gamma.smashplug.OptimizeCoupledModel(
        smash_model,
        model_gamma,
        GammaGriddedObservation,
        control_parameters_list=["cp", "cft", "hydraulics_coefficient", "spreading"],
        bounds={
            "cp": [0.1, 1000.0],
            "cft": [0.1, 1000.0],
            "hydraulics_coefficient": [0.3, 5.0],
            "spreading": [0.5, 3.0],
        },
        maxiter=15,
        tol=0.00001,
        ScaleGradientsByBounds=True,
        ScaleGammaGradients=True,
    )
)


# check the final cost
# convert Vector discharges from Smash to Gamma gridded discharges
GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
    optimized_smash_model.input_data.qobs,
    optimized_gamma_model,
    smash_dt=optimized_smash_model.setup.dt,
)
# compute the final cost
model_gamma.cost_function(
    GammaGriddedObservation, optimized_gamma_model.routing_results.discharges
)

# get the abscisse in time step for smash and gamma (easy plotting capabilities)
abcisses = gamma.smashplug.GetCorrespondingTimeStepRange(smash_model, model_gamma)

# make a plot
fig, ax = plt.subplots()
plt.rc("font", size=18)

ax.plot(
    abcisses["ts_gamma"],
    model_gamma.routing_results.discharges[
        :, model_gamma.routing_mesh.controlled_nodes[0] - 1
    ],
    label="SMASH-grg-Gamma, Qsim Initial",
    lw=2,
)  # nodes are in Fortran index: in python array we need -1: how to fix that ?

ax.plot(
    abcisses["ts_gamma"],
    optimized_gamma_model.routing_results.discharges[
        :, model_gamma.routing_mesh.controlled_nodes[0] - 1
    ],
    label="SMASH-grg-Gamma, Qsim Optimal",
    lw=2,
)  # nodes are in Fortran index: in python array we need -1: how to fix that ?

ax.plot(
    abcisses["ts_smash"],
    model_smash_gr_calibrated.output.qsim[0, :],
    label="Reference SMASH-gra, Qsim Optimal",
    lw=2,
)

ax.plot(
    abcisses["ts_gamma"],
    GammaGriddedObservation[:, model_gamma.routing_mesh.controlled_nodes[0] - 1],
    label="Observed discharges (OBS)",
    color="0",
    lw=2,
)

ax.axes.set_xlabel("Time-Step (hours)")
ax.axes.set_ylabel("Discharges (m^3/s)")

ax.legend(loc="upper left", fontsize=18)
ax.axes.grid(True, alpha=0.7, ls="--")

fig.show()
plot = (fig, ax)

functions_smash_plot.save_figure(
    plot, figname="ExpCalVal.pdf", xlim=[0, 1500], xsize=24, ysize=10
)
functions_smash_plot.save_figure(
    plot,
    figname="ExpCalVal_zoom1.pdf",
    xlim=[550, 725],
    ylim=[0, 250],
    xsize=12,
    ysize=10,
)
functions_smash_plot.save_figure(
    plot,
    figname="ExpCalVal_zoom2.pdf",
    xlim=[1150, 1400],
    ylim=[0, 350],
    xsize=12,
    ysize=10,
)

# convert gamma vector parameters to gridded data
Grid_Hydraulic_Coef = gamma.smashplug.GammaVectorsToSmashGrid(
    optimized_gamma_model.routing_parameters.hydraulics_coefficient,
    smash_model,
    model_gamma,
)

Spreading = gamma.smashplug.GammaVectorsToSmashGrid(
    optimized_gamma_model.routing_parameters.spreading, smash_model, model_gamma
)

# plot that parameters
functions_smash_plot.plot_image(
    Grid_Hydraulic_Coef,
    figname="ExpCalVal_CoefHydraulics.pdf",
    mask=smash_model.mesh.active_cell,
    title="Hydraulic coefficients",
    title_font_size=14,
)

functions_smash_plot.plot_image(
    Spreading,
    figname="ExpCalVal_Spreading.pdf",
    mask=smash_model.mesh.active_cell,
    title="Spreading coefficients",
    title_font_size=14,
)

functions_smash_plot.plot_image(
    optimized_smash_model.parameters.cp,
    figname="ExpCalVal_Cp.pdf",
    mask=smash_model.mesh.active_cell,
    title="Capacities of the production reservoir",
    title_font_size=14,
)

functions_smash_plot.plot_image(
    optimized_smash_model.parameters.cft,
    figname="ExpCalVal_cft.pdf",
    mask=smash_model.mesh.active_cell,
    title="Capacities of the transfert reservoir",
    title_font_size=14,
)


# get observed discharges:
GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
    model_smash_gr_calibrated.input_data.qobs,
    optimized_gamma_model,
    smash_dt=model_smash_gr_calibrated.setup.dt,
)

# Compute the nse
nse_gra = functions_smash_stats.nse(
    model_smash_gr_calibrated.input_data.qobs[0, :],
    model_smash_gr_calibrated.output.qsim[0, :],
    lacuna=-99.0,
)
nse_grg = functions_smash_stats.nse(
    GammaGriddedObservation[
        :, optimized_gamma_model.routing_mesh.controlled_nodes[0] - 1
    ],
    optimized_gamma_model.routing_results.discharges[
        :, optimized_gamma_model.routing_mesh.controlled_nodes[0] - 1
    ],
    lacuna=-99.0,
)

print("Calage NSE_gra=", nse_gra)
print("Calage NSE_grg=", nse_grg)


# Validation here

# New setup
setup_cance, mesh_cance = smash.load_dataset("Cance")
setup_cance.update({"save_qsim_domain": True})
setup_cance["structure"] = "gr-g"
setup_cance["start_time"] = "2013-09-01 00:00"
setup_cance["end_time"] = "2014-04-01 00:00"
setup_cance["pet_directory"] = "/home/maxime/DATA/ETP-SFR-FRA-INTERA_L93/"
setup_cance["prcp_directory"] = "/home/maxime/DATA/PLUIE/"
model_smash_grg_validation = smash.Model(setup_cance, mesh_cance)

# New SMASH Gamma
setup_cance["structure"] = "gr-a"
setup_cance.update({"save_qsim_domain": False})
model_smash_gra_validation = smash.Model(setup_cance, mesh_cance)

# set parameters and states
model_smash_grg_validation.parameters = optimized_smash_model.parameters.copy()
model_smash_gra_validation.parameters = model_smash_gr_calibrated.parameters.copy()
# in case of restart at last time step
# model_smash_grg_validation.states=optimized_smash_model_cance.outputs.fstates.copy()
# model_smash_gra_validation.states=optimized_smash_model_cance.outputs.fstates.copy()
model_smash_grg_validation.states.hp = 0.1
model_smash_gra_validation.states.hp = 0.1

model_gamma_validation = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash_grg_validation,
    dt=900,
    vmin=0.1,
    vmax=10.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.1,
    ponderation_regul=0.0,
    velocity_computation="qm3",
)

# Set optimized parameter
model_gamma_validation.routing_parameters_change(
    hydraulics_coefficient=optimized_gamma_model.routing_parameters.hydraulics_coefficient.copy(),
    spreading=optimized_gamma_model.routing_parameters.spreading.copy(),
)

# set previous states #in case of restart at last time step
# ~ #Attention au shape des states qui sont différents, prendre les states correspondant aux shape du nouveau modèle. Les shapes des 2 modèle sont différents car l'un provient d'un calage et les la longueur de la fenêtre mémoire est calculé à partir des bornes ou de la valeur des paramètres si celle-ci n'est pas spécifiée
# TODO: We should provide a function for that
# ~ shapeOfstates=model_gamma_validation.routing_states.states.shape
# ~ model_gamma_validation.routing_states.states=optimized_gamma_model_cance.routing_states.states[0:shapeOfstates[0],:].copy()
# ~ model_gamma_validation.routing_states.remainder=optimized_gamma_model_cance.routing_states.remainder[0:shapeOfstates[0],:].copy()

# run the coupled model SMASH-Gamma
gamma.smashplug.RunCoupledModel(model_smash_grg_validation, model_gamma_validation)

# run the Smash reference model
model_smash_gra_validation.run(inplace=True)

# Get the abscisse for plotting
abcisses = gamma.smashplug.GetCorrespondingTimeStepRange(
    model_smash_grg_validation, model_gamma_validation
)
fig, ax = plt.subplots()
plt.rc("font", size=18)
ax.plot(
    abcisses["ts_gamma"],
    model_gamma_validation.routing_results.discharges[
        :, model_gamma_validation.routing_mesh.controlled_nodes[0] - 1
    ],
    label="SMASH-grg-Gamma (Validation)",
    lw=2,
)  # nodes are in Fortran index: in python array we need -1: how to fix that ?

ax.plot(
    abcisses["ts_smash"],
    model_smash_gra_validation.output.qsim[0, :],
    label="SMASH-gra (Validation)",
    lw=2,
)

ax.plot(
    abcisses["ts_smash"],
    model_smash_grg_validation.input_data.qobs[0, :],
    label="Observations",
    color="0",
    lw=2,
)

ax.legend(loc="upper left")
ax.axes.grid(True, alpha=0.7, ls="--")
ax.axes.set_xlabel("Time-Step (hours)")
ax.axes.set_ylabel("Discharges (m^3/s)")
fig.show()
plot = (fig, ax)

functions_smash_plot.save_figure(
    plot, figname="ExpCalVal_validation_period.pdf", xlim=[0, 5000], xsize=24, ysize=10
)
functions_smash_plot.save_figure(
    plot,
    figname="ExpCalVal_validation_period_zoom1.pdf",
    xlim=[150, 250],
    ylim=[0, 50],
    xsize=12,
    ysize=10,
)
functions_smash_plot.save_figure(
    plot,
    figname="ExpCalVal_validation_period_zoom2.pdf",
    xlim=[2760, 2860],
    ylim=[0, 100],
    xsize=12,
    ysize=10,
)

# get observed discharges:
GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
    model_smash_grg_validation.input_data.qobs,
    model_gamma_validation,
    smash_dt=model_smash_grg_validation.setup.dt,
)

# Compute the NSE
nse_gra = functions_smash_stats.nse(
    model_smash_gra_validation.input_data.qobs[0, :],
    model_smash_gra_validation.output.qsim[0, :],
    lacuna=-99.0,
)
nse_grg = functions_smash_stats.nse(
    GammaGriddedObservation[
        :, model_gamma_validation.routing_mesh.controlled_nodes[0] - 1
    ],
    model_gamma_validation.routing_results.discharges[
        :, model_gamma_validation.routing_mesh.controlled_nodes[0] - 1
    ],
    lacuna=-99.0,
)

print("NSE_gra=", nse_gra)
print("NSE_grg=", nse_grg)
