# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import smash
import gamma_routing_model as gamma
import scipy
import smashbox as sb
import pyhdf5_handler
import os
import datetime

bv = "Cance"

scenarios_gamma = {
    "s1": {
        "control_vector": ["cp", "ct", "hydraulics_coefficient"],
        "varying_spread": 0,
        "spreading_uniform": 1,
    },
    "s2": {
        "control_vector": ["cp", "ct", "hydraulics_coefficient", "spreading"],
        "varying_spread": 1,
        "spreading_uniform": 1,
    },
    "s3": {
        "control_vector": ["cp", "ct", "hydraulics_coefficient", "spreading"],
        "varying_spread": 1,
        "spreading_uniform": 0,
    },
    "s4": {
        "control_vector": ["hydraulics_coefficient"],
        "varying_spread": 0,
        "spreading_uniform": 1,
    },
    "s5": {
        "control_vector": ["hydraulics_coefficient", "spreading"],
        "varying_spread": 1,
        "spreading_uniform": 1,
    },
    "s6": {
        "control_vector": ["hydraulics_coefficient", "spreading"],
        "varying_spread": 1,
        "spreading_uniform": 0,
    },
}

scenarios_smash = {
    "s1": {
        "control_vector": ["cp", "ct", "llr"],
    },
    "s2": {
        "control_vector": [
            "llr",
        ],
    },
}

bounds = {
    "cp": [1, 2000],
    "ct": [1, 1000],
    "llr": [1, 200],
    "hydraulics_coefficient": [0.3, 5],
    "spreading": [0.5, 5],
}


start_time = "2008-01-01 00:00"
end_time = "2018-01-01 00:00"

start_time = "2014-01-01 00:00"
end_time = "2014-04-01 00:00"

val_start_time = "2018-01-01 00:00"
val_end_time = "2025-01-01 00:00"

gamma_dt = 900.0

warmup_time = 150
pdt_start_optim = int(warmup_time * 24 * 3600 / gamma_dt)
end_warmup = datetime.datetime.fromisoformat(start_time) - datetime.timedelta(
    days=warmup_time
)


gamma_settings = {
    "vmin": 0.1,
    "vmax": 10.0,
    "mode_discretization_step": 0.1,
    "spreading_discretization_step": 0.1,
    "ponderation_regul": 10.0,
    "velocity_computation": "qmm",
    "varying_spread": 1,
    "spreading_uniform": 1,
    "criteria": "nse",
    "ponderation_cost": 10000.0,
    "pdt_start_optim": pdt_start_optim,
}

exu = {
    "Cance": ["V3524010", "V3515010", "V3517010"],
    "Ardeche": ["V5014010", "V5015210", "V5004030"],
}


sbc = sb.SmashBox()

sbc.myparam.list_param()
sbc.myparam.set_param("outletsID", exu[bv])
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

os.makedirs("./figures", exist_ok=True)

# create required model
sbc.newmodel("Cance")
sbc.Cance.generate_mesh()
sbc.Cance.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
sbc.Cance.mysetup.update_setup({"start_time": val_start_time, "end_time": val_end_time})
sbc.Cance.forward_run(warmup=warmup_time)

sbc.Cance.save_model_container(f"{bv}.hdf5")
sbc.Cance.save_as_smash_model(f"./{bv}")


sbc.Cance.myplot.plot_mesh(fig_settings={"figname": "./figures/mesh_{bv}.pdf"})
# sbc.Cance.myplot.plot_hydrograph(columns=[0, 1, 2])
# sbc.Cance.myplot.multiplot_parameters()

sbc.copymodel("Cance", "Cance_Gamma")
sbc.Cance_Gamma.mysetup.update_setup({"routing_module": "zero", "return_opt_grad": "qe"})
sbc.Cance_Gamma.model()  # rebuild model because setup has changed but without reading the data
sbc.copymodeldata("Cance", "Cance_Gamma")  # copy the data
sbc.Cance_Gamma.forward_run()
sbc.Cance_Gamma.save_model_container(f"{bv}.hdf5")
sbc.Cance_Gamma.save_as_smash_model(f"./{bv}")


# optimize smash model
for s in scenarios_smash:
    sbc.copymodel(bv, f"{bv}_{s}")
    new_bounds = {}
    for p in scenarios_smash[s]["control_vector"]:
        new_bounds.update({p: bounds[p]})

    optimize_options = {
        "parameters": scenarios_smash[s]["control_vector"],
        "bounds": new_bounds,
        "termination_crit": {
            "maxiter": 30,
            "factr": 1e6,
        },
    }
    cost_options = {
        "gauge": "dws",
        "end_warmup": end_warmup,
    }
    model = getattr(sb, f"{bv}_{s}")
    model.optimize(
        mapping="distributed",
        optimize_options=optimize_options,
        cost_options=cost_options,
    )
    model.save_model_container(f"{bv}.hdf5")
    model.save_as_smash_model(f"./{bv}")


for s in scenarios_gamma:
    sbc.copymodel("Cance_Gamma", f"{bv}_Gamma_calibrated_{s}")
    new_bounds = {}
    for p in scenarios_smash[s]["control_vector"]:
        new_bounds.update({p: bounds[p]})

    # configure the Gamma model
    model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
        sbc.Cance_Gamma.mysmashmodel, dt=gamma_dt, **gamma_settings
    )

    # Set parameter
    model_gamma.routing_parameters.hydraulics_coefficient = 1.0
    model_gamma.routing_parameters.spreading = 2.0

    model = getattr(sb, f"{bv}_Gamma_calibrated_{s}")

    gamma.smashplug.RunCoupledModel(model.mysmashmodel, model_gamma)

    model_gamma.routing_mesh.controlled_nodes = 0
    model_gamma.routing_mesh.controlled_nodes[0] = model_gamma.routing_mesh.gauge_nodes[
        np.argmax(
            sbc.Cance_Gamma.mysmashmodel.mesh.area[
                model_gamma.routing_mesh.gauge_name_index - 1
            ]
        )
    ]

    BestControlVector, optimized_smash_model, optimized_gamma_model = (
        gamma.smashplug.OptimizeCoupledModel(
            sbc.Cance_Gamma.mysmashmodel,
            model_gamma,
            # GammaGriddedObservation,
            control_parameters_list=scenarios_gamma[s]["control_vector"],
            bounds=new_bounds,
            maxiter=20,
            maxfun=20,
            tol=1e-6,
            ScaleGradients=False,
            ScaleGammaGradientsBySurface=True,
            optim_type="local",
            local_optimizer="L-BFGS-B",  # L-BFGS-B | trust-constr | SLSQP | TNC
        )
    )

    model.mysmashmodel = optimized_smash_model.copy()
    model.save_model_container(f"{bv}.hdf5")
    model.save_as_smash_model(f"./{bv}")

    fig, ax = gamma.smashplug.plot.multiplot_parameters(
        sbc.Cance_Gamma_calibrated, optimized_gamma_model
    )

    sbc.plot.plot.save_figure(f"{bv}_parameters_{s}.pdf", xsize=12, ysize=10)

    fig, ax = gamma.smashplug.plot.plot_hydrograph(
        (
            {
                "model": sbc.Cance.mysmashmodel,
                "fig_settings": {"color": "blue", "label": "Smash Regional"},
            },
            {
                "model": sbc.Cance_Gamma.mysmashmodel,
                "fig_settings": {"color": "green", "label": "Smash Regional + Gamma"},
            },
            {
                "model": model.mysmashmodel,
                "fig_settings": {
                    "color": "red",
                    "label": "Smash Regional + Gamma-Calibrated",
                },
            },
            # {
            #     "model": sbc.Cance.optimize_model,
            #     "fig_settings": {
            #         "color": "orange",
            #         "label": "Smash Rgional + llr-Calibrated",
            #     },
            # },
        )
    )

    sbc.plot.plot.save_figure(f"{bv}_hydrograph_{s}.pdf", xsize=12, ysize=10)


# save figure


# validation temporelle
sbc.newmodel("Cance_validation")
sbc.Cance.generate_mesh()
sbc.Cance.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
sbc.Cance.mysetup.update_setup({"date_start": val_start_time, "date_end": val_end_time})
sbc.Cance.forward_run(warmup=warmup_time)

sbc.copymodel("Cance_validation", "Cance_Gamma_validation")
sbc.Cance_Gamma.mysetup.update_setup({"routing_module": "zero", "return_opt_grad": "qe"})
sbc.Cance_Gamma.model()  # rebuild model because setup has changed but without reading the data
sbc.copymodeldata("Cance", "Cance_Gamma")  # copy the data
sbc.Cance_Gamma.forward_run()


# configure the Gamma model for warmup (not used here)
model_gamma_validation_warm = gamma.smashplug.ConfigureGammaWithSmash(
    sbc.Cance_Gamma.warmup_model, dt=gamma_dt, **gamma_settings
)
model_gamma_validation_warm.routing_parameters.hydraulics_coefficient = (
    optimized_gamma_model.routing_parameters.hydraulics_coefficient
)
model_gamma_validation_warm.routing_parameters.spreading = (
    optimized_gamma_model.routing_parameters.spreading
)

gamma.smashplug.RunCoupledModel(sbc.Cance_Gamma.warmup_model, model_gamma_validation_warm)


# configure the Gamma model
model_gamma_validation = gamma.smashplug.ConfigureGammaWithSmash(
    sbc.Cance_Gamma.mysmashmodel, dt=gamma_dt, **gamma_settings
)
model_gamma_validation.routing_parameters.hydraulics_coefficient = (
    optimized_gamma_model.routing_parameters.hydraulics_coefficient
)
model_gamma_validation.routing_parameters.spreading = (
    optimized_gamma_model.routing_parameters.spreading
)
# Set initial_states : not useful, since the calibration involve a jump in discharge at the first time steps. may be better to start at 0.
model_gamma_validation.routing_memory.remainder_init = (
    model_gamma_validation_warm.routing_memory.remainder
)
model_gamma_validation.routing_memory.states_init = (
    model_gamma_validation_warm.routing_memory.states
)

gamma.smashplug.RunCoupledModel(sbc.Cance_Gamma.warmup_model, model_gamma_validation)

# sbc.Cance_Gamma_calibrated.myplot.multiplot_parameters()

# Grid_Hydraulic_Coef = gamma.smashplug.GammaVectorsToSmashGrid(
#     optimized_gamma_model.routing_parameters.hydraulics_coefficient,
#     optimized_smash_model,
#     optimized_gamma_model,
# )
# Spreading = gamma.smashplug.GammaVectorsToSmashGrid(
#     optimized_gamma_model.routing_parameters.spreading,
#     optimized_smash_model,
#     optimized_gamma_model,
# )


# fig, ax = sb.plot.plot.plot_image(
#     Spreading,
#     mask=optimized_smash_model.mesh.active_cell,
#     ax_settings={"title": "Spreading coefficients", "title_fontsize": 14},
#     fig_settings={
#         "figname": f"ExpCalValRoutage_Spreading_VS{model_gamma.routing_setup.varying_spread}_SU{model_gamma.routing_setup.spreading_uniform}_CAL{ctrl}.pdf"
#     },
# )
# fig.show()

# fig, ax = sb.plot.plot.plot_image(
#     Grid_Hydraulic_Coef,
#     mask=optimized_smash_model.mesh.active_cell,
#     ax_settings={"title": "Hydraulic coefficients", "title_fontsize": 14},
#     fig_settings={
#         "figname": f"CoefHydraulics_VS{model_gamma.routing_setup.varying_spread}_SU{model_gamma.routing_setup.spreading_uniform}_CAL{ctrl}.pdf"
#     },
# )
# fig.show()


# fig, ax = sb.plot.plot.plot_hydrograph(
#     sbc.Cance.mysmashmodel,
#     plot_rainfall=True,
#     plot_settings_sim={"color": "blue", "label": "Smash Regional"},
# )
# fig, ax = sb.plot.plot.plot_hydrograph(
#     sbc.Cance_Gamma.mysmashmodel,
#     figure=(fig, *ax),
#     plot_rainfall=False,
#     plot_qobs=False,
#     plot_settings_sim={"color": "green", "label": "Smash Regional + Gamma"},
# )
# fig, ax = sb.plot.plot.plot_hydrograph(
#     sbc.Cance_Gamma_calibrated.mysmashmodel,
#     figure=(fig, *ax),
#     plot_rainfall=False,
#     plot_qobs=False,
#     plot_settings_sim={"color": "red", "label": "Smash Regional + Gamma-Calibrated"},
# )
# fig, ax = sb.plot.plot.plot_hydrograph(
#     sbc.Cance.optimize_model,
#     figure=(fig, *ax),
#     plot_rainfall=False,
#     plot_qobs=False,
#     plot_settings_sim={"color": "orange", "label": "Smash Rgional + llr-Calibrated"},
# )
# fig.show()


# import matplotlib.pyplot as plt

# fig, ax = plt.subplots()

# g_node = 380
# qobs = np.where(
#     GammaGriddedObservation[:, g_node - 1] < 0,
#     np.nan,
#     GammaGriddedObservation[:, g_node - 1],
# )
# qgamma_init = model_gamma.routing_results.discharges[:, g_node - 1]
# ax.plot(
#     qobs,
#     label="Qobs",
#     color="0",
#     lw=2,
# )
# ax.plot(
#     qgamma_init,
#     label="qsim",
#     color="red",
#     lw=2,
# )
# fig.show()


# def plot_q(model_smash, model_gamma, g_node):

#     col_gamma = g_node - 1
#     col_smash = (
#         model_gamma.routing_mesh.gauge_name_index[
#             model_gamma.routing_mesh.gauge_nodes == g_node
#         ]
#         - 1
#     )

#     fig, ax = plt.subplots()

#     qobs = np.where(
#         model_smash.response_data.q.transpose() < 0.0,
#         np.nan,
#         model_smash.response_data.q.transpose(),
#     )

#     ax.plot(
#         qobs[:, col_smash],
#         label="Qobs",
#         color="0",
#         lw=2,
#     )

#     ax.plot(model_smash.response.q.transpose()[:, col_smash], label="Qsmash", lw=2)

#     x = np.arange(0, model_smash.response.q.shape[1], 0.25)
#     ax.plot(
#         x,
#         model_gamma.routing_results.discharges[:, col_gamma],
#         lw=2,
#         label="QGamma_final",
#     )

#     ax.legend(loc="upper left")
#     ax.axes.grid(True, alpha=0.7, ls="--")
#     ax.axes.set_xlabel("Time-Step (hours)")
#     ax.axes.set_ylabel("Discharges (m^3/s)")
#     fig.show()


# plot_q(sbc.Cance.mysmashmodel, optimized_gamma_model, g_node=380)
