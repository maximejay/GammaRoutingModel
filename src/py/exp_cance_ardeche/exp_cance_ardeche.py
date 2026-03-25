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

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mode", type=str)
parser.add_argument("-c", "--catchment", type=str)
parser.add_argument("-s", "--scenario_key", type=str)

args = parser.parse_args()

mode = args.mode
catchment = args.catchment
scenario_key = args.scenario_key

# mode = "smash_gamma"
# catchment = "Cance"
# scenario_key = "s1"

# Optimizer : SLSQP ?? au lieu de LBGSFB ?

scenarios_gamma = {
    "s1": {
        "control_vector": ["cp", "ct", "hydraulics_coefficient"],
        "ScaleGammaGradientsBySurface": False,
        "ScaleGradients": False,
        "gamma_settings": {
            "varying_spread": 0,
            "spreading_uniform": 1,
        },
    },
    "s2": {
        "control_vector": ["cp", "ct", "hydraulics_coefficient", "spreading"],
        "ScaleGammaGradientsBySurface": True,
        "ScaleGradients": True,
        "gamma_settings": {
            "varying_spread": 1,
            "spreading_uniform": 1,
        },
    },
    "s3": {
        "control_vector": ["cp", "ct", "hydraulics_coefficient", "spreading"],
        "ScaleGammaGradientsBySurface": True,
        "ScaleGradients": True,
        "gamma_settings": {
            "varying_spread": 1,
            "spreading_uniform": 0,
        },
    },
    "s4": {
        "control_vector": ["hydraulics_coefficient"],
        "ScaleGammaGradientsBySurface": True,
        "ScaleGradients": False,
        "gamma_settings": {
            "varying_spread": 0,
            "spreading_uniform": 1,
        },
    },
    "s5": {
        "control_vector": ["hydraulics_coefficient", "spreading"],
        "ScaleGammaGradientsBySurface": True,
        "ScaleGradients": False,
        "gamma_settings": {
            "varying_spread": 1,
            "spreading_uniform": 1,
        },
    },
    "s6": {
        "control_vector": ["hydraulics_coefficient", "spreading"],
        "ScaleGammaGradientsBySurface": True,
        "ScaleGradients": False,
        "scale"
        "gamma_settings": {
            "varying_spread": 1,
            "spreading_uniform": 0,
        },
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

working_path = "/home/maxime/DEV/GammaRoutingModel/src/py/exp_cance_ardeche/"
path_gamma_hdf5 = os.path.join(working_path, "Gamma_hdf5")
os.makedirs(path_gamma_hdf5, exist_ok=True)

start_time = "2008-01-01 00:00"
end_time = "2018-01-01 00:00"

# start_time = "2014-09-01 00:00"
# end_time = "2014-12-01 00:00"

# start_time = "2014-10-01 00:00"
# end_time = "2014-11-14 00:00"

val_start_time = "2018-01-01 00:00"
val_end_time = "2025-01-01 00:00"

# val_start_time = "2014-12-01 00:00"
# val_end_time = "2015-03-01 00:00"

gamma_dt = 900.0
warmup_time = 450

# pdt_start_optim = int(warmup_time * 24 * 3600 / gamma_dt)
# end_warmup = (
#     datetime.datetime.fromisoformat(start_time) + datetime.timedelta(days=warmup_time)
# ).strftime("%Y-%d-%m %H:%M")


gamma_settings = {
    "vmin": 0.1,
    "vmax": 10.0,
    "mode_discretization_step": 0.1,
    "spreading_discretization_step": 0.1,
    "ponderation_regul": 0.0,
    "velocity_computation": "qmm",
    "varying_spread": 1,
    "spreading_uniform": 1,
    "criteria": "nse",
    "ponderation_cost": 10000.0,
    "pdt_start_optim": int(warmup_time * 24 * 3600 / gamma_dt),
}

exu = {
    "Cance": ["V3524010", "V3515010", "V3517010"],
    "Ardeche": ["V5014010", "V5015210", "V5004030"],
}

sbc = sb.SmashBox()

sbc.myparam.list_param()
sbc.myparam.set_param("outletsID", exu[catchment])
sbc.myparam.set_param(
    "outlets_shapefile",
    "/home/maxime/DATA/BNBVlight/BNBV_light_bassins.shp",
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


def init_models(sbc, catchment, start_time, end_time, warmup_time):
    # create required base regional model
    sbc.newmodel(f"{catchment}_init")
    model = getattr(sbc, f"{catchment}_init")
    model.generate_mesh()
    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup({"start_time": start_time, "end_time": end_time})
    model.forward_run(warmup=warmup_time)

    model.save_model_container_hdf5(f"{working_path}/{catchment}_init.hdf5")
    model.save_model_container(f"{working_path}/{catchment}")

    # Gamma
    # for simul with gamma: create the required base model:
    sbc.copymodel(f"{catchment}_init", f"{catchment}_Gamma_init")
    model = getattr(sbc, f"{catchment}_Gamma_init")
    model.mysetup.update_setup({"routing_module": "zero", "return_opt_grad": "qe"})
    model.model()  # rebuild model because setup has changed but without reading the data
    sbc.copymodeldata(f"{catchment}_init", f"{catchment}_Gamma_init")  # copy the data
    model.forward_run()  # run anc copy warming states

    # configure the Gamma model
    model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel, dt=gamma_dt, **gamma_settings
    )

    # Set parameter
    model_gamma.routing_parameters_change(
        hydraulics_coefficient=1.0,
        spreading=2.0,
    )

    gamma.smashplug.RunCoupledModel(model.mysmashmodel, model_gamma)

    model.save_model_container_hdf5(f"{catchment}_Gamma_init.hdf5")
    model.save_model_container(f"./{catchment}")
    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_Gamma_init.hdf5", model_gamma
    )

    model.myplot.plot_mesh(fig_settings={"figname": "./figures/mesh_{catchment}.pdf"})


# sbc2 = sb.SmashBox()
# sbc2.load_containers(f"./{catchment}")
# sbc2.Cance.mymesh.load_catchment_polygon(param=sbc2.Cance._myparam.param)

# model.myplot.plot_hydrograph(columns=[0, 1, 2])
# model.myplot.multiplot_parameters()


def optimize_smash_model(
    sbc,
    catchment,
    s,
    scenario,
    bounds,
    warmup_time,
    start_time,
    end_time,
    val_start_time,
    val_end_time,
):
    # sbc.copymodel(catchment, f"{catchment}_{s}")
    sbc.newmodel(f"{catchment}_calibrated_{s}")
    model = getattr(sbc, f"{catchment}_calibrated_{s}")
    model.generate_mesh()
    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"

    model.mysetup.update_setup({"start_time": start_time, "end_time": end_time})
    model.model()

    # optimization including the warming period
    new_start_time = (
        datetime.datetime.fromisoformat(start_time) - datetime.timedelta(days=warmup_time)
    ).strftime("%Y-%d-%m %H:%M")

    new_bounds = {}
    for p in scenario[s]["control_vector"]:
        new_bounds.update({p: bounds[p]})

    optimize_options = {
        "parameters": scenario[s]["control_vector"],
        "bounds": new_bounds,
        "termination_crit": {
            "maxiter": 30,
            "factr": 1e6,
        },
    }
    cost_options = {
        "gauge": "dws",
        "end_warmup": start_time,
    }
    model.optimize(
        start_time=new_start_time,
        end_time=end_time,
        mapping="distributed",
        optimize_options=optimize_options,
        cost_options=cost_options,
    )

    # optimized_smash_parameters = model.optimize_model.rr_parameters.copy()

    # last forward run on the coalibrated period only
    model.forward_run(warmup=warmup_time)

    model.validate(start_time=val_start_time, end_time=val_end_time, warmup=warmup_time)

    model.save_model_container_hdf5(f"{working_path}/{catchment}_calibrated_{s}.hdf5")
    model.save_model_container(f"{working_path}/{catchment}")

    # validation temporelle => developp a validation model in smashbox with a function validation()
    # sbc.newmodel(f"{catchment}_calibrated_validation_{s}")
    # model = getattr(sbc, f"{catchment}_calibrated_validation_{s}")

    # model.generate_mesh()
    # model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    # model.mysetup.update_setup({"start_time": val_start_time, "end_time": val_end_time})
    # model.model()
    # model.mysmashmodel.rr_parameters = optimized_smash_parameters
    # model.forward_run(warmup=warmup_time)

    # model.save_model_container_hdf5(f"{working_path}/{catchment}.hdf5")
    # model.save_model_container(f"{working_path}/{catchment}")


def optimize_smash_gamma_model(
    sbc,
    catchment,
    s,
    scenario,
    bounds,
    gamma_settings,
    warmup_time,
    start_time,
    end_time,
    val_start_time,
    val_end_time,
):

    # update gamma settings
    gamma_settings.update(scenario[s]["gamma_settings"])

    sbc.newmodel(f"{catchment}_Gamma_optimize_{s}")
    model = getattr(sbc, f"{catchment}_Gamma_optimize_{s}")
    model.generate_mesh()
    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup({"routing_module": "zero", "return_opt_grad": "qe"})
    new_start_time = (
        datetime.datetime.fromisoformat(start_time) - datetime.timedelta(days=warmup_time)
    ).strftime("%Y-%m-%d %H:%M")
    model.mysetup.update_setup({"start_time": new_start_time, "end_time": end_time})
    model.model()

    new_bounds = {}
    for p in scenario[s]["control_vector"]:
        new_bounds.update({p: bounds[p]})

    # configure the Gamma model
    model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel, dt=gamma_dt, **gamma_settings
    )

    # Set parameter (initial guess)
    model_gamma.routing_parameters_change(
        hydraulics_coefficient=1.0,
        spreading=2.0,
    )

    model_gamma.routing_mesh.controlled_nodes = 0
    model_gamma.routing_mesh.controlled_nodes[0] = model_gamma.routing_mesh.gauge_nodes[
        np.argmax(
            model.mysmashmodel.mesh.area[model_gamma.routing_mesh.gauge_name_index - 1]
        )
    ]

    gamma.smashplug.RunCoupledModel(model.mysmashmodel, model_gamma)
    # get observed discharges:
    GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
        model.mysmashmodel.response_data.q,
        model_gamma,
        smash_dt=model.mysmashmodel.setup.dt,
    )

    # compute the initial cost
    model_gamma.cost_function(
        GammaGriddedObservation, model_gamma.routing_results.discharges
    )

    model_gamma.routing_results.costs

    BestControlVector, optimized_smash_model, optimized_gamma_model = (
        gamma.smashplug.OptimizeCoupledModel(
            model.mysmashmodel,
            model_gamma,
            control_parameters_list=scenario[s]["control_vector"],
            bounds=new_bounds,
            maxiter=20,
            maxfun=20,
            tol=1e-6,
            ScaleGradients=scenario[s]["ScaleGradients"],
            ScaleGammaGradientsBySurface=scenario[s]["ScaleGammaGradientsBySurface"],
            optim_type="local",
            local_optimizer="SLSQP",  # L-BFGS-B | trust-constr | SLSQP | TNC
        )
    )

    model.mysmashmodel = optimized_smash_model.copy()
    model.forward_run()  # no need to warmup (warmup period included)

    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_Gamma_optimize_{s}.hdf5", optimized_gamma_model
    )
    model.save_model_container_hdf5(f"{working_path}/{catchment}.hdf5")
    model.save_model_container(f"{working_path}/{catchment}")

    # run on the normal period + do a side warmup
    sbc.newmodel(f"{catchment}_Gamma_calibrated_{s}")
    model = getattr(sbc, f"{catchment}_Gamma_calibrated_{s}")
    model.generate_mesh()
    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup({"routing_module": "zero", "return_opt_grad": "qe"})
    model.mysetup.update_setup({"start_time": start_time, "end_time": end_time})
    model.model()
    model.optimize_model = optimized_smash_model.copy()
    model.forward_run(warmup=warmup_time)

    # configure the Gamma model warmup
    model_gamma_warmup = gamma.smashplug.ConfigureGammaWithSmash(
        model.warmup_model, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma_warmup.routing_parameters_change(
        hydraulics_coefficient=optimized_gamma_model.routing_parameters.hydraulics_coefficient,
        spreading=optimized_gamma_model.routing_parameters.spreading,
    )
    gamma.smashplug.RunCoupledModel(model.warmup_model, model_gamma_warmup)

    # configure Gamma model normal period
    model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma.routing_parameters_change(
        hydraulics_coefficient=optimized_gamma_model.routing_parameters.hydraulics_coefficient,
        spreading=optimized_gamma_model.routing_parameters.spreading,
    )

    model_gamma.routing_memory.remainder_init = (
        model_gamma_warmup.routing_memory.remainder
    )
    model_gamma.routing_memory.states_init = model_gamma_warmup.routing_memory.states
    model.mysmashmodel.rr_initial_states = model.warmup_model.rr_final_states.copy()

    gamma.smashplug.RunCoupledModel(model.mysmashmodel, model_gamma)

    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_Gamma_calibrated_{s}.hdf5", model_gamma
    )
    model.save_model_container_hdf5(f"{working_path}/{catchment}.hdf5")
    model.save_model_container(f"{working_path}/{catchment}")

    # validation temporelle
    sbc.newmodel(
        f"{catchment}_Gamma_calibrated_validation_{s}"
    )  # , "Cance_Gamma_validation")
    model = getattr(sbc, f"{catchment}_Gamma_calibrated_validation_{s}")
    model.generate_mesh()
    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup({"routing_module": "zero", "return_opt_grad": "qe"})
    model.mysetup.update_setup({"start_time": val_start_time, "end_time": val_end_time})
    model.model()
    model.mysmashmodel.rr_parameters = optimized_smash_model.rr_parameters.copy()
    model.forward_run(warmup=warmup_time)

    # configure the Gamma model (warmup)
    model_gamma_validation_warmup = gamma.smashplug.ConfigureGammaWithSmash(
        model.warmup_model, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma_validation_warmup.routing_parameters_change(
        hydraulics_coefficient=optimized_gamma_model.routing_parameters.hydraulics_coefficient,
        spreading=optimized_gamma_model.routing_parameters.spreading,
    )

    gamma.smashplug.RunCoupledModel(model.warmup_model, model_gamma_validation_warmup)

    # configure the Gamma model for validation only
    model_gamma_validation = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma_validation.routing_parameters_change(
        hydraulics_coefficient=optimized_gamma_model.routing_parameters.hydraulics_coefficient,
        spreading=optimized_gamma_model.routing_parameters.spreading,
    )
    model_gamma_validation.routing_memory.remainder_init = (
        model_gamma_validation_warmup.routing_memory.remainder
    )
    model_gamma_validation.routing_memory.states_init = (
        model_gamma_validation_warmup.routing_memory.states
    )
    model.mysmashmodel.rr_initial_states = model.warmup_model.rr_final_states.copy()

    gamma.smashplug.RunCoupledModel(model.mysmashmodel, model_gamma_validation)

    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_gamma_calibrated_validation_{s}.hdf5",
        model_gamma_validation,
    )
    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_gamma_calibrated_validation_warmup_{s}.hdf5",
        model_gamma_validation_warmup,
    )
    model.save_model_container_hdf5(f"{working_path}/{catchment}.hdf5")
    model.save_model_container(f"{working_path}/{catchment}")


# figures
def plot(catchment, s):
    sbc = sb.SmashBox()
    sbc.load_containers(f"{working_path}/{catchment}")

    model_smash = getattr(sbc, f"{catchment}_init")

    model_smash.myplot.multiplot_parameters(
        fig_settings={
            "figname": f"{working_path}/figures/{catchment}_init_parameters_{s}.pdf"
        }
    )

    model_smash = getattr(sbc, f"{catchment}_calibrated_{s}")

    model_smash.myplot.multiplot_parameters(
        fig_settings={
            "figname": f"{working_path}/figures/{catchment}_calibrated_parameters_{s}.pdf"
        }
    )

    model_gamma = pyhdf5_handler.read_hdf5file_as_dict(
        f"{path_gamma_hdf5}/{catchment}_Gamma_optimize_{s}.hdf5"
    )

    model_smash = getattr(sbc, f"{catchment}_Gamma_calibrated_{s}")

    fig, ax = gamma.smashplug.plot.multiplot_parameters(model_smash, model_gamma)
    sb.plot.plot.save_figure(
        fig,
        f"{working_path}/figures/{catchment}_Gamma_parameters_{s}.pdf",
        xsize=10,
        ysize=14,
    )

    model_smash_init = getattr(sbc, f"{catchment}_init")
    model_gamma_init = getattr(sbc, f"{catchment}_Gamma_init")
    model_gamma_calibrated = getattr(sbc, f"{catchment}_Gamma_calibrated_{s}")
    model_smash_calibrated = getattr(sbc, f"{catchment}_calibrated_{s}")
    fig, ax = gamma.smashplug.plot.plot_hydrograph(
        (
            {
                "model": model_smash_init.mysmashmodel,
                "fig_settings": {"color": "blue", "label": "Smash Regional"},
            },
            {
                "model": model_gamma_init.mysmashmodel,
                "fig_settings": {"color": "green", "label": "Smash Regional + Gamma"},
            },
            {
                "model": model_gamma_calibrated.mysmashmodel,
                "fig_settings": {
                    "color": "red",
                    "label": "Smash + Gamma-Calibrated",
                },
            },
            {
                "model": model_smash_calibrated.mysmashmodel,
                "fig_settings": {
                    "color": "orange",
                    "label": "Smash Calibrated",
                },
            },
        )
    )
    sb.plot.plot.save_figure(
        fig, f"{working_path}/figures/{catchment}_hydrograph_{s}.pdf", xsize=12, ysize=10
    )

    model_smash_validation = getattr(sbc, f"{catchment}_calibrated_{s}")
    model_gamma_validation = getattr(sbc, f"{catchment}_Gamma_calibrated_validation_{s}")
    fig, ax = gamma.smashplug.plot.plot_hydrograph(
        (
            {
                "model": model_smash_validation.validation_model,
                "fig_settings": {"color": "blue", "label": "Smash validation"},
                "outlets_name": ["V3524010"],
            },
            {
                "model": model_gamma_validation.mysmashmodel,
                "fig_settings": {
                    "color": "green",
                    "label": "Smash + Gamma validation",
                },
                "outlets_name": ["V3524010"],
            },
        )
    )
    sb.plot.plot.save_figure(
        fig,
        f"{working_path}/figures/{catchment}_hydrograph_validation_{s}.pdf",
        xsize=12,
        ysize=10,
    )


if mode == "init":
    init_models(sbc, catchment, start_time, end_time, warmup_time)

# Smash
# optimize smash model
if mode == "smash":
    optimize_smash_model(
        sbc,
        catchment,
        scenario_key,
        scenarios_smash,
        bounds,
        warmup_time,
        start_time,
        end_time,
        val_start_time,
        val_end_time,
    )


# Gamma
if mode == "smash_gamma":
    optimize_smash_gamma_model(
        sbc,
        catchment,
        scenario_key,
        scenarios_gamma,
        bounds,
        gamma_settings,
        warmup_time,
        start_time,
        end_time,
        val_start_time,
        val_end_time,
    )

if mode == "plot":
    plot(catchment, scenario_key)


# save figure


# configure the Gamma model for warmup (not used here)
# model_gamma_validation_warm = gamma.smashplug.ConfigureGammaWithSmash(
#     sbc.Cance_Gamma.warmup_model, dt=gamma_dt, **gamma_settings
# )
# model_gamma_validation_warm.routing_parameters.hydraulics_coefficient = (
#     optimized_gamma_model.routing_parameters.hydraulics_coefficient
# )
# model_gamma_validation_warm.routing_parameters.spreading = (
#     optimized_gamma_model.routing_parameters.spreading
# )

# gamma.smashplug.RunCoupledModel(sbc.Cance_Gamma.warmup_model, model_gamma_validation_warm)


# configure the Gamma model
# model_gamma_validation = gamma.smashplug.ConfigureGammaWithSmash(
#     sbc.Cance_Gamma.mysmashmodel, dt=gamma_dt, **gamma_settings
# )
# model_gamma_validation.routing_parameters.hydraulics_coefficient = (
#     optimized_gamma_model.routing_parameters.hydraulics_coefficient
# )
# model_gamma_validation.routing_parameters.spreading = (
#     optimized_gamma_model.routing_parameters.spreading
# )
# # Set initial_states : not useful, since the calibration involve a jump in discharge at the first time steps. may be better to start at 0.
# model_gamma_validation.routing_memory.remainder_init = (
#     model_gamma_validation_warm.routing_memory.remainder
# )
# model_gamma_validation.routing_memory.states_init = (
#     model_gamma_validation_warm.routing_memory.states
# )

# gamma.smashplug.RunCoupledModel(sbc.Cance_Gamma.warmup_model, model_gamma_validation)

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
