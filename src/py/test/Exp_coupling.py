import numpy as np
import matplotlib.pyplot as plt
import smash

import gamma_routing_model as gamma
import scipy


# import functions_smash_plot
# import functions_smash_stats

# Smash reference model structure gr-a
setup_cance, mesh_cance = smash.factory.load_dataset("Cance")
setup_cance["routing_module"] = "zero"
smash_model = smash.Model(setup_cance, mesh_cance)
smash_model.forward_run()

# smash_model.response.qt.shape
setup_cance, mesh_cance = smash.factory.load_dataset("Cance")
smash_optimize = smash.Model(setup_cance, mesh_cance)
optimize_options = {
    "parameters": ["cp", "ct", "llr"],
    "bounds": {"cp": (1, 1000), "ct": (1, 1000), "llr": (1, 1000)},
    "termination_crit": {
        "maxiter": 15,
        "factr": 1e6,
    },
}
cost_options = {
    "gauge": "dws",
    "end_warmup": "2014-10-02 00:00",
}
smash_optimize.optimize(
    mapping="distributed",
    optimize_options=optimize_options,
    cost_options=cost_options,
)

# smash compilé avec les options de debuggage


# configure the Gamma model
model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
    smash_model,
    dt=900,
    vmin=0.1,
    vmax=5.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.1,
    ponderation_regul=0.0,
    velocity_computation="qm3",
    varying_spread=0,
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
    smash_model.response_data.q, model_gamma, smash_dt=smash_model.setup.dt
)

# compute the initial cost
model_gamma.cost_function(GammaGriddedObservation, model_gamma.routing_results.discharges)

boundaries = gamma.smashplug.get_boundaries(
    model_gamma,
    smash_model,
    control_parameters_list=["cp", "ct", "hydraulics_coefficient", "spreading"],
    bounds={
        "cp": [0.1, 1000.0],
        "ct": [0.1, 1000.0],
        "hydraulics_coefficient": [0.3, 5.0],
        "spreading": [0.5, 3.0],
    },
)

# Get the control vector
ControlVector = gamma.smashplug.VectorizeModelParameters(
    smash_model,
    model_gamma,
    control_parameters_list=["cp", "ct", "hydraulics_coefficient", "spreading"],
    bounds={
        "cp": [0.1, 1000.0],
        "ct": [0.1, 1000.0],
        "hydraulics_coefficient": [0.3, 5.0],
        "spreading": [0.5, 3.0],
    },
)

# ControlVector.update(
#     {
#         "bounds": {
#             "cp": [0.1, 1000.0],
#             "ct": [0.1, 1000.0],
#             "hydraulics_coefficient": [0.3, 5.0],
#             "spreading": [0.5, 3.0],
#         }
#     }
# )

control_vector = gamma.smashplug.normalize_control_vector(ControlVector)

# get the inflows
GammaInflows = gamma.smashplug.GetGammaInflowFromSmash(smash_model, dt=900)


cost, grad = gamma.smashplug.ComputeCostAndGradients(
    ControlVector["X"],
    ControlVector,
    smash_model,
    model_gamma,
    GammaGriddedObservation,
    True,
    True,
)

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
        control_parameters_list=["cp", "ct", "hydraulics_coefficient", "spreading"],
        bounds={
            "cp": [0.1, 1000.0],
            "ct": [0.1, 1000.0],
            "hydraulics_coefficient": [0.3, 5.0],
            "spreading": [0.5, 3],
        },
        maxiter=30,
        tol=0.00001,
        ScaleGradientsByBounds=False,
        ScaleGammaGradients=False,
    )
)

import matplotlib.pyplot as plt
from gamma_routing_model.smashplug.coupling import functions_smash_plot


def save_figure(
    plot,
    figname="myfigure",
    xsize=8,
    ysize=6,
    transparent=False,
    dpi=80,
    xlim=[None, None],
    ylim=[None, None],
):

    if isinstance(plot, tuple):
        fig = plot[0]
        ax = plot[1]
    else:
        fig = plot
        ax = None

    if ax is not None:
        if ylim[0] != None:
            ax.set_ylim(bottom=ylim[0])
        if ylim[1] != None:
            ax.set_ylim(top=ylim[1])
        if xlim[0] != None:
            ax.set_xlim(left=xlim[0])
        if xlim[1] != None:
            ax.set_xlim(right=xlim[1])

    fig.set_size_inches(xsize, ysize, forward=False)
    fig.savefig(figname, transparent=transparent, dpi=dpi, bbox_inches="tight")


fig, ax = plt.subplots()
ax.plot(smash_optimize.response_data.q.transpose()[:, 0], label="Qobs", color="0", lw=2)
ax.plot(smash_optimize.response.q.transpose()[:, 0], label="Qsmash", lw=2)
x = np.arange(0, 1440, 0.25)
ax.plot(
    x,
    model_gamma.routing_results.discharges[
        :, np.argmax(model_gamma.routing_mesh.cumulated_surface)
    ],
    lw=2,
    label="QGamma_initial",
)
ax.plot(
    x,
    optimized_gamma_model.routing_results.discharges[
        :, np.argmax(model_gamma.routing_mesh.cumulated_surface)
    ],
    lw=2,
    label="QGamma_final",
)

ax.legend(loc="upper left")
ax.axes.grid(True, alpha=0.7, ls="--")
ax.axes.set_xlabel("Time-Step (hours)")
ax.axes.set_ylabel("Discharges (m^3/s)")
fig.show()

save_figure(
    (fig, ax),
    figname="ExpCoupling_calib_uniform_spread.pdf",
    xlim=[0, 1500],
    ylim=[0, 350],
    xsize=48,
    ysize=10,
)

# Convert Gazmma vector parameter to gridded data
Grid_Hydraulic_Coef = gamma.smashplug.GammaVectorsToSmashGrid(
    optimized_gamma_model.routing_parameters.hydraulics_coefficient,
    optimized_smash_model,
    optimized_gamma_model,
)
Spreading = gamma.smashplug.GammaVectorsToSmashGrid(
    optimized_gamma_model.routing_parameters.spreading,
    optimized_smash_model,
    optimized_gamma_model,
)


functions_smash_plot.plot_image(
    Grid_Hydraulic_Coef,
    figname="ExpCalValRoutage_CoefHydraulics.pdf",
    mask=optimized_smash_model.mesh.active_cell,
    title="Hydraulic coefficients",
    title_font_size=14,
)

functions_smash_plot.plot_image(
    Spreading,
    figname="ExpCalValRoutage_Spreading.pdf",
    mask=optimized_smash_model.mesh.active_cell,
    title="Spreading coefficients",
    title_font_size=14,
)

functions_smash_plot.plot_image(
    optimized_smash_model.parameters.cp,
    figname="ExpCalValRoutage_Cp.pdf",
    mask=optimized_smash_model.mesh.active_cell,
    title="Capacities of the production reservoir",
    title_font_size=14,
)

functions_smash_plot.plot_image(
    optimized_smash_model.parameters.ct,
    figname="ExpCalValRoutage_Cft.pdf",
    mask=optimized_smash_model.mesh.active_cell,
    title="Capacities of the transfert reservoir",
    title_font_size=14,
)
