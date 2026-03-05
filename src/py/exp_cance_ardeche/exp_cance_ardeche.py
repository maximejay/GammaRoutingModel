# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import smash

import gamma_routing_model as gamma
import scipy
import smashbox as sb


setup_cance, mesh_cance = smash.factory.load_dataset("Cance")


sbc = sb.SmashBox()

sbc.myparam.list_param()
sbc.myparam.set_param("outletsID", list(mesh_cance["code"]))
sbc.myparam.set_param(
    "outlets_shapefile", "/home/maxime/DATA/BNBVlight/BNBV_light_bassins.shp"
)
sbc.myparam.set_param(
    "outlets_database_fields",
    {
        "coord_x": "X_L93",
        "coord_y": "Y_L93",
        "area": "SURF",
        "id": "CODE_SITE",
        "id_shapefile": "CODE_BNBV",
    },
)
warmup_time = 150

sbc.newmodel("Cance")
sbc.Cance.generate_mesh()
sbc.Cance.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
sbc.Cance.forward_run(warmup=warmup_time)

sbc.Cance.myplot.plot_mesh()
sbc.Cance.myplot.plot_hydrograph(columns=[0, 1, 2])
sbc.Cance.myplot.multiplot_parameters()


sbc.newmodel("Cance_Gamma")
sbc.Cance_Gamma.generate_mesh()
sbc.Cance_Gamma.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
sbc.Cance_Gamma.mysetup.update_setup({"routing_module": "zero", "return_opt_grad": "qe"})
sbc.Cance_Gamma.forward_run(warmup=warmup_time)


# Todo build a copy function of model in smashbox
sbc.newmodel("Cance_Gamma_calibrated")
sbc.Cance_Gamma_calibrated.generate_mesh()
sbc.Cance_Gamma_calibrated.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
sbc.Cance_Gamma_calibrated.mysetup.update_setup(
    {"routing_module": "zero", "return_opt_grad": "qe"}
)
sbc.Cance_Gamma_calibrated.model()

# configure the Gamma model
model_gamma_warm = gamma.smashplug.ConfigureGammaWithSmash(
    sbc.Cance_Gamma.warmup_model,
    dt=900,
    vmin=0.1,
    vmax=10.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.1,
    ponderation_regul=0.0,
    velocity_computation="qmm",
    varying_spread=1,
    spreading_uniform=1,
    criteria="nse",
    ponderation_cost=10000.0,
    pdt_start_optim=1600,
)
model_gamma_warm.routing_parameters.hydraulics_coefficient = 1.0
model_gamma_warm.routing_parameters.spreading = 2.0
gamma.smashplug.RunCoupledModel(sbc.Cance_Gamma.warmup_model, model_gamma_warm)


# configure the Gamma model
model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
    sbc.Cance_Gamma.mysmashmodel,
    dt=900,
    vmin=0.1,
    vmax=10.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.1,
    ponderation_regul=10.0,
    velocity_computation="qmm",
    varying_spread=0,
    spreading_uniform=1,
    criteria="nse",
    ponderation_cost=10000.0,
    pdt_start_optim=1600,
)
# Set parameter
model_gamma.routing_parameters.hydraulics_coefficient = 1.0
model_gamma.routing_parameters.spreading = 2.0

# Set initial_states : not useful, since the calibration involve a jump in discharge at the first time steps. may be better to start at 0.
# model_gamma.routing_memory.remainder_init = model_gamma_warm.routing_memory.remainder
# model_gamma.routing_memory.states_init = model_gamma_warm.routing_memory.states
# model_gamma.routing_memory = (
#     model_gamma_warm.routing_memory.copy()
# )  # does not work ! pb _handle not copied !

gamma.smashplug.RunCoupledModel(sbc.Cance_Gamma.mysmashmodel, model_gamma)

sbc.Cance_Gamma.myplot.plot_hydrograph(columns=[0, 1, 2])

model_gamma.routing_mesh.controlled_nodes = 0
model_gamma.routing_mesh.controlled_nodes[0] = model_gamma.routing_mesh.gauge_nodes[
    np.argmax(
        sbc.Cance_Gamma.mysmashmodel.mesh.area[
            model_gamma.routing_mesh.gauge_name_index - 1
        ]
    )
]


# Calage
# control_parameters_list = ["cp", "ct", "hydraulics_coefficient", "spreading"]
# control_parameters_list = ["hydraulics_coefficient", "spreading"]
control_parameters_list = ["hydraulics_coefficient"]
# get observed discharges:
GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
    sbc.Cance_Gamma.mysmashmodel.response_data.q,
    model_gamma,
    smash_dt=sbc.Cance_Gamma.mysmashmodel.setup.dt,
)


BestControlVector, optimized_smash_model, optimized_gamma_model = (
    gamma.smashplug.OptimizeCoupledModel(
        sbc.Cance_Gamma.mysmashmodel,
        model_gamma,
        GammaGriddedObservation,
        control_parameters_list=control_parameters_list,
        bounds={
            # "cp": [1.0, 2000.0],
            # "ct": [1.0, 1000.0],
            "hydraulics_coefficient": [0.3, 5.0],
            # "spreading": [0.5, 5.0],
        },
        maxiter=10,
        tol=None,
        ScaleGradients=False,
        ScaleGammaGradientsBySurface=False,
        optim_type="local",
        local_optimizer="L-BFGS-B",  # L-BFGS-B | trust-constr | SLSQP | TNC
    )
)

sbc.Cance_Gamma_calibrated.mysmashmodel = optimized_smash_model
sbc.Cance_Gamma_calibrated.myplot.plot_hydrograph(columns=[0, 1, 2])
sbc.Cance_Gamma.myplot.plot_hydrograph(columns=[0, 1, 2])
sbc.Cance_Gamma_calibrated.myplot.multiplot_parameters()
sbc.Cance.myplot.multiplot_parameters()


ctrl = ""
for c in control_parameters_list:
    ctrl = ctrl + "-" + c

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


gamma.smashplug.functions_smash_plot.plot_image(
    Grid_Hydraulic_Coef,
    figname=f"CoefHydraulics_VS{model_gamma.routing_setup.varying_spread}_SU{model_gamma.routing_setup.spreading_uniform}_CAL{ctrl}.pdf",
    mask=optimized_smash_model.mesh.active_cell,
    title="Hydraulic coefficients",
    title_font_size=14,
)

gamma.smashplug.functions_smash_plot.plot_image(
    Spreading,
    figname=f"ExpCalValRoutage_Spreading_VS{model_gamma.routing_setup.varying_spread}_SU{model_gamma.routing_setup.spreading_uniform}_CAL{ctrl}.pdf",
    mask=optimized_smash_model.mesh.active_cell,
    title="Spreading coefficients",
    title_font_size=14,
)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

g_node = 380
qobs = np.where(
    GammaGriddedObservation[:, g_node - 1] < 0,
    np.nan,
    GammaGriddedObservation[:, g_node - 1],
)
qgamma_init = model_gamma.routing_results.discharges[:, g_node - 1]
ax.plot(
    qobs,
    label="Qobs",
    color="0",
    lw=2,
)
ax.plot(
    qgamma_init,
    label="qsim",
    color="red",
    lw=2,
)
fig.show()


def plot_q(model_smash, model_gamma, g_node):

    col_gamma = g_node - 1
    col_smash = (
        model_gamma.routing_mesh.gauge_name_index[
            model_gamma.routing_mesh.gauge_nodes == g_node
        ]
        - 1
    )

    fig, ax = plt.subplots()

    qobs = np.where(
        model_smash.response_data.q.transpose() < 0.0,
        np.nan,
        model_smash.response_data.q.transpose(),
    )

    ax.plot(
        qobs[:, col_smash],
        label="Qobs",
        color="0",
        lw=2,
    )

    ax.plot(model_smash.response.q.transpose()[:, col_smash], label="Qsmash", lw=2)

    x = np.arange(0, sbc.Cance_Gamma.mysmashmodel.response.q.shape[1], 0.25)
    ax.plot(
        x,
        model_gamma.routing_results.discharges[:, col_gamma],
        lw=2,
        label="QGamma_final",
    )

    ax.legend(loc="upper left")
    ax.axes.grid(True, alpha=0.7, ls="--")
    ax.axes.set_xlabel("Time-Step (hours)")
    ax.axes.set_ylabel("Discharges (m^3/s)")
    fig.show()


plot_q(sbc.Cance_Gamma.mysmashmodel, model_gamma, g_node=380)
