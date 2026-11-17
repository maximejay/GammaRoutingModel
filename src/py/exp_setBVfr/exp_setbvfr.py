#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:30:07 2026

@author: maxime
"""

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

scenario_key = "s1"


# Optimizer : SLSQP ?? au lieu de LBGSFB ?

scenarios_gamma = {
    "s1": {
        "control_vector": ["cp", "ct", "hc", "sc"],
        "mapping": "uniform",
        "gamma_settings": {
            "varying_spread": 1,
            "spreading_uniform": 1,
            "hydraulics_coef_uniform": 1,
        },
    },
    "s2": {
        "control_vector": ["cp", "ct", "hc", "sc"],
        "mapping": "distributed",
        "gamma_settings": {
            "varying_spread": 1,
            "spreading_uniform": 1,
        },
    },
}

scenarios_smash = {
    "s1": {
        "control_vector": ["cp", "ct", "llr"],
        "mapping": "uniform",
    },
    "s2": {
        "control_vector": ["cp", "ct", "llr"],
        "mapping": "distributed",
    },
}


bounds = {
    "cp": [1, 2000],
    "ct": [1, 1000],
    "llr": [1, 200],
    "hc": [0.3, 5],
    "sc": [0.5, 5],
}

working_path = "/home/maxime/DEV/GammaRoutingModel/src/py/exp_setBVfr/"
path_gamma_hdf5 = os.path.join(working_path, "Gamma_hdf5")
os.makedirs(path_gamma_hdf5, exist_ok=True)

start_time = "2012-01-01 00:00"
end_time = "2018-01-01 00:00"

start_time = "2014-09-01 00:00"
end_time = "2014-12-01 00:00"

start_time = "2014-10-01 00:00"
end_time = "2014-11-14 00:00"

val_start_time = "2018-01-01 00:00"
val_end_time = "2025-01-01 00:00"

val_start_time = "2014-12-01 00:00"
val_end_time = "2015-03-01 00:00"

gamma_dt = 900.0
warmup_time = 30


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
    "ponderation_cost": 10.0,
    "pdt_start_optim": int(warmup_time * 24 * 3600 / gamma_dt),
}

list_mesh = os.listdir(
    "/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh/"
)
mesh = list_mesh[0]
catchment = mesh[5:13]

sbc = sb.SmashBox()

sbc.myparam.set_param(
    "outlets_shapefile",
    "/home/maxime/DATA/BNBVlight/BNBV_light_bassins.shp",
)

sbc.myparam.set_param("smash_parameters_dt", 900)

os.makedirs(os.path.join(working_path, "figures"), exist_ok=True)


def init_models(
    sbc,
    catchment,
    start_time,
    end_time,
    warmup_time,
    val_start_time,
    val_end_time,
):
    # create required base regional model
    sbc.newmodel(f"{catchment}_init")
    model = getattr(sbc, f"{catchment}_init")

    model.mymesh.mesh = pyhdf5_handler.read_hdf5file_as_dict(
        f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh/mesh_{catchment}.hdf5"
    )

    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup(
        {"start_time": start_time, "end_time": end_time}
    )

    model.forward_run(
        warmup=warmup_time, set_parameters={"cp": 200, "ct": 100, "llr": 30}
    )

    model.validate(
        start_time=val_start_time, end_time=val_end_time, warmup=warmup_time
    )

    model.save_model_container(f"{working_path}/{catchment}")

    # Gamma
    # for simul with gamma: create the required base model:
    sbc.copymodel(f"{catchment}_init", f"{catchment}_Gamma_init")
    model = getattr(sbc, f"{catchment}_Gamma_init")
    model.mysetup.update_setup({"routing_module": "zero"})
    model.model(
        set_parameters={"cp": 200, "ct": 100}
    )  # rebuild model because setup has changed but without reading the data
    sbc.copymodeldata(
        f"{catchment}_init", f"{catchment}_Gamma_init"
    )  # copy the data
    model.forward_run()  # run anc copy warming states

    # configure the Gamma model
    model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel.smash, dt=gamma_dt, **gamma_settings
    )

    # Set parameter
    model_gamma.routing_parameters_change(
        hc=1.0,
        sc=2.0,
    )

    gamma.smashplug.RunCoupledModel(model.mysmashmodel.smash, model_gamma)

    model.save_model_container(f"{working_path}/{catchment}")
    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_Gamma_init.hdf5", model_gamma
    )

    model.myplot.plot_mesh(
        fig_settings={
            "figname": f"{working_path}/figures/mesh_{catchment}.pdf"
        }
    )


# sbc2 = sb.SmashBox()
# sbc2.load_containers(f"./{catchment}")
# sbc2.Cance.mymesh.load_catchment_polygon(param=sbc2.Cance._myparam.param)

# model.myplot.plot_hydrograph(columns=[0, 1, 2])
# model.myplot.multiplot_parameters()


def optimize_smash_model(
    sbc,
    catchment,
    scenario_key,
    scenario,
    bounds,
    warmup_time,
    start_time,
    end_time,
    val_start_time,
    val_end_time,
):
    # sbc.copymodel(catchment, f"{catchment}_{s}")
    sbc.newmodel(f"{catchment}_calibrated_{scenario_key}")
    model = getattr(sbc, f"{catchment}_calibrated_{scenario_key}")

    model.mymesh.mesh = pyhdf5_handler.read_hdf5file_as_dict(
        f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh/mesh_{catchment}.hdf5"
    )

    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"

    model.mysetup.update_setup(
        {"start_time": start_time, "end_time": end_time}
    )
    model.model(set_parameters={"cp": 200, "ct": 100, "llr": 30})

    # optimization including the warming period
    new_start_time = (
        datetime.datetime.fromisoformat(start_time)
        - datetime.timedelta(days=warmup_time)
    ).strftime("%Y-%d-%m %H:%M")

    new_bounds = {}
    for p in scenario[scenario_key]["control_vector"]:
        new_bounds.update({p: bounds[p]})

    optimize_options = {
        "parameters": scenario[scenario_key]["control_vector"],
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
        mapping=scenario[scenario_key]["mapping"],
        optimize_options=optimize_options,
        cost_options=cost_options,
    )

    # optimized_smash_parameters = model.optimize_model.rr_parameters.copy()

    # last forward run on the calibrated period only
    model.forward_run(warmup=warmup_time)

    model.validate(
        start_time=val_start_time, end_time=val_end_time, warmup=warmup_time
    )

    # model.save_model_container_hdf5(
    # f"{working_path}/{catchment}_calibrated_{s}.hdf5"
    # )
    model.save_model_container(f"{working_path}/{catchment}")

    # validation temporelle => developp a validation model in smashbox with a function validation()
    # sbc.newmodel(f"{catchment}_calibrated_validation_{s}")
    # model = getattr(sbc, f"{catchment}_calibrated_validation_{s}")

    # model.generate_mesh()
    # model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    # model.mysetup.update_setup({"start_time": val_start_time, "end_time": val_end_time})
    # model.model()
    # model.mysmashmodel.smash.rr_parameters = optimized_smash_parameters
    # model.forward_run(warmup=warmup_time)

    # model.save_model_container_hdf5(f"{working_path}/{catchment}.hdf5")
    # model.save_model_container(f"{working_path}/{catchment}")


def optimize_smash_gamma_model(
    sbc,
    catchment,
    scenario_key,
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
    gamma_settings.update(scenario[scenario_key]["gamma_settings"])

    sbc.newmodel(f"{catchment}_Gamma_optimize_{scenario_key}")
    model = getattr(sbc, f"{catchment}_Gamma_optimize_{scenario_key}")

    model.mymesh.mesh = pyhdf5_handler.read_hdf5file_as_dict(
        f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh/mesh_{catchment}.hdf5"
    )
    # model.generate_mesh()

    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup({"routing_module": "zero"})
    new_start_time = (
        datetime.datetime.fromisoformat(start_time)
        - datetime.timedelta(days=warmup_time)
    ).strftime("%Y-%m-%d %H:%M")
    model.mysetup.update_setup(
        {"start_time": new_start_time, "end_time": end_time}
    )
    model.model(set_parameters={"cp": 200, "ct": 100})

    new_bounds = {}
    for p in scenario[scenario_key]["control_vector"]:
        new_bounds.update({p: bounds[p]})

    # configure the Gamma model
    model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel.smash, dt=gamma_dt, **gamma_settings
    )

    # Set parameter (initial guess)
    model_gamma.routing_parameters_change(
        hc=1.0,
        sc=2.0,
    )

    model_gamma.routing_mesh.controlled_nodes = 0
    model_gamma.routing_mesh.controlled_nodes[0] = (
        model_gamma.routing_mesh.gauge_nodes[
            np.argmax(
                model.mysmashmodel.smash.mesh.area[
                    model_gamma.routing_mesh.gauge_name_index - 1
                ]
            )
        ]
    )

    gamma.smashplug.RunCoupledModel(model.mysmashmodel.smash, model_gamma)
    # get observed discharges:
    GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
        model.mysmashmodel.smash.response_data.q,
        model_gamma,
        smash_dt=model.mysmashmodel.smash.setup.dt,
    )

    # compute the initial cost
    model_gamma.cost_function(
        GammaGriddedObservation, model_gamma.routing_results.discharges
    )

    model_gamma.routing_results.costs

    BestControlVector, optimized_smash_model, optimized_gamma_model = (
        gamma.smashplug.OptimizeCoupledModel(
            model.mysmashmodel.smash,
            model_gamma,
            control_parameters_list=scenario[scenario_key]["control_vector"],
            bounds=new_bounds,
            maxiter=30,
            optim_type="uniform",
        )
    )

    model.mysmashmodel.smash = optimized_smash_model.copy()
    model.forward_run()  # no need to warmup (warmup period included)

    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_Gamma_optimize_{scenario_key}.hdf5",
        optimized_gamma_model,
    )

    model.save_model_container(f"{working_path}/{catchment}")

    # run on the normal period + do a side warmup
    sbc.newmodel(f"{catchment}_Gamma_calibrated_{scenario_key}")
    model = getattr(sbc, f"{catchment}_Gamma_calibrated_{scenario_key}")

    model.mymesh.mesh = pyhdf5_handler.read_hdf5file_as_dict(
        f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh/mesh_{catchment}.hdf5"
    )

    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup({"routing_module": "zero"})
    model.mysetup.update_setup(
        {"start_time": start_time, "end_time": end_time}
    )

    model.optimize_model.smash = optimized_smash_model.copy()
    model.model()

    model.forward_run(warmup=warmup_time)

    # configure the Gamma model warmup
    model_gamma_warmup = gamma.smashplug.ConfigureGammaWithSmash(
        model.warmup_model.smash, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma_warmup.routing_parameters_change(
        hc=optimized_gamma_model.routing_parameters.hc,
        sc=optimized_gamma_model.routing_parameters.sc,
    )
    gamma.smashplug.RunCoupledModel(
        model.warmup_model.smash, model_gamma_warmup
    )

    # configure Gamma model normal period
    model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel.smash, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma.routing_parameters_change(
        hc=optimized_gamma_model.routing_parameters.hc,
        sc=optimized_gamma_model.routing_parameters.sc,
    )

    model_gamma.routing_memory.remainder_init = (
        model_gamma_warmup.routing_memory.remainder
    )
    model_gamma.routing_memory.states_init = (
        model_gamma_warmup.routing_memory.states
    )
    model.mysmashmodel.smash.rr_initial_states = (
        model.warmup_model.smash.rr_final_states.copy()
    )

    gamma.smashplug.RunCoupledModel(model.mysmashmodel.smash, model_gamma)

    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_Gamma_calibrated_{scenario_key}.hdf5",
        model_gamma,
    )
    # model.save_model_container_hdf5(
    #     f"{working_path}/{catchment}_calibrated_{s}.hdf5"
    # )
    model.save_model_container(f"{working_path}/{catchment}")

    # validation temporelle
    sbc.newmodel(
        f"{catchment}_Gamma_calibrated_validation_{scenario_key}"
    )  # , "Cance_Gamma_validation")
    model = getattr(
        sbc, f"{catchment}_Gamma_calibrated_validation_{scenario_key}"
    )

    model.mymesh.mesh = pyhdf5_handler.read_hdf5file_as_dict(
        f"/home/maxime/DEV/GammaRoutingModel/src/py/Generation_mesh_bv_france/mesh/mesh_{catchment}.hdf5"
    )

    model.mysetup.setup_file = "setup_local_gr4_dt3600.yaml"
    model.mysetup.update_setup({"routing_module": "zero"})
    model.mysetup.update_setup(
        {"start_time": val_start_time, "end_time": val_end_time}
    )

    model.optimize_model.smash = optimized_smash_model.copy()
    model.model()

    model.forward_run(warmup=warmup_time)

    # configure the Gamma model (warmup)
    model_gamma_validation_warmup = gamma.smashplug.ConfigureGammaWithSmash(
        model.warmup_model.smash, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma_validation_warmup.routing_parameters_change(
        hc=optimized_gamma_model.routing_parameters.hc,
        sc=optimized_gamma_model.routing_parameters.sc,
    )

    gamma.smashplug.RunCoupledModel(
        model.warmup_model.smash, model_gamma_validation_warmup
    )

    # configure the Gamma model for validation only
    model_gamma_validation = gamma.smashplug.ConfigureGammaWithSmash(
        model.mysmashmodel.smash, dt=gamma_dt, **gamma_settings
    )
    # Set parameter
    model_gamma_validation.routing_parameters_change(
        hc=optimized_gamma_model.routing_parameters.hc,
        sc=optimized_gamma_model.routing_parameters.sc,
    )
    model_gamma_validation.routing_memory.remainder_init = (
        model_gamma_validation_warmup.routing_memory.remainder
    )
    model_gamma_validation.routing_memory.states_init = (
        model_gamma_validation_warmup.routing_memory.states
    )
    model.mysmashmodel.smash.rr_initial_states = (
        model.warmup_model.smash.rr_final_states.copy()
    )

    gamma.smashplug.RunCoupledModel(
        model.mysmashmodel.smash, model_gamma_validation
    )

    pyhdf5_handler.save_object_to_hdf5file(
        f"{path_gamma_hdf5}/{catchment}_gamma_calibrated_validation_{scenario_key}.hdf5",
        model_gamma_validation,
    )
    # pyhdf5_handler.save_object_to_hdf5file(
    #     f"{path_gamma_hdf5}/{catchment}_gamma_calibrated_validation_warmup_{scenario_key}.hdf5",
    #     model_gamma_validation_warmup,
    # )

    model.save_model_container(f"{working_path}/{catchment}")


# figures
def plot(catchment, s):

    list_dir = os.listdir(f"{working_path}/{catchment}")
    l = list(filter(lambda x: x.endswith("init"), list_dir))
    l.extend(list(filter(lambda x: x.endswith(f"_{s}"), list_dir)))

    sbc = sb.SmashBox()
    sbc.load_containers(
        f"{working_path}/{catchment}", load_catchment_polygon=True, only=l
    )

    model_smash = getattr(sbc, f"{catchment}_init")

    model_smash.myplot.multiplot_parameters(
        fig_settings={
            "figname": f"{working_path}/figures/{catchment}_init_parameters_{s}.pdf"
        }
    )

    if hasattr(sbc, f"{catchment}_calibrated_{s}"):
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

    fig, ax = gamma.smashplug.plot.multiplot_parameters(
        model_smash, model_gamma
    )
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
    tuple_plot = (
        {
            "model": model_smash_init.mysmashmodel.smash,
            "fig_settings": {"color": "blue", "label": "Smash Regional"},
        },
        {
            "model": model_gamma_init.mysmashmodel.smash,
            "fig_settings": {
                "color": "green",
                "label": "Smash Regional + Gamma",
            },
        },
        {
            "model": model_gamma_calibrated.mysmashmodel.smash,
            "fig_settings": {
                "color": "red",
                "label": "Smash + Gamma-Calibrated",
            },
        },
    )
    if hasattr(sbc, f"{catchment}_calibrated_{s}"):
        model_smash_calibrated = getattr(sbc, f"{catchment}_calibrated_{s}")
        tuple_plot = tuple_plot + (
            {
                "model": model_smash_calibrated.mysmashmodel.smash,
                "fig_settings": {
                    "color": "orange",
                    "label": "Smash Calibrated",
                },
            },
        )

    fig, ax = gamma.smashplug.plot.plot_hydrograph(tuple_plot)
    sb.plot.plot.save_figure(
        fig,
        f"{working_path}/figures/{catchment}_hydrograph_{s}.pdf",
        xsize=12,
        ysize=10,
    )

    model_smash_validation = getattr(sbc, f"{catchment}_calibrated_{s}")
    model_gamma_validation = getattr(
        sbc, f"{catchment}_Gamma_calibrated_validation_{s}"
    )

    tuple_plot = (
        {
            "model": model_smash_validation.validation_model.smash,
            "fig_settings": {"color": "blue", "label": "Smash validation"},
            "outlets_name": ["V3524010"],
        },
    )
    if hasattr(sbc, f"{catchment}_calibrated_{s}"):
        model_smash_validation = getattr(sbc, f"{catchment}_calibrated_{s}")
        tuple_plot = tuple_plot + (
            {
                "model": model_gamma_validation.mysmashmodel.smash,
                "fig_settings": {
                    "color": "green",
                    "label": "Smash + Gamma validation",
                },
                "outlets_name": ["V3524010"],
            },
        )

    model_gamma_validation = getattr(
        sbc, f"{catchment}_Gamma_calibrated_validation_{s}"
    )
    fig, ax = gamma.smashplug.plot.plot_hydrograph(tuple_plot)
    sb.plot.plot.save_figure(
        fig,
        f"{working_path}/figures/{catchment}_hydrograph_validation_{s}.pdf",
        xsize=12,
        ysize=10,
    )


def compute_stats(catchment, sbc=None, list_dir=None):
    if list_dir is None:
        list_dir = os.listdir(f"{working_path}/{catchment}")

    if sbc is None:
        sbc = sb.SmashBox()
        sbc.load_containers(
            f"{working_path}/{catchment}", load_catchment_polygon=False
        )

    for subdir in list_dir:
        res = getattr(sbc, subdir)
        res.mysmashmodel.mystats.foutlets_stats()
        res.mysmashmodel.mystats.fmisfit_stats()

        structure = (
            pyhdf5_handler.src.object_handler.generate_object_structure(
                res.mysmashmodel.mystats, include_method=False
            )
        )
        pyhdf5_handler.save_object_to_hdf5file(
            path_to_hdf5=os.path.join(
                working_path, catchment, subdir, "smashbox.hdf5"
            ),
            instance=res.mysmashmodel.mystats,
            keys_data=structure,
            location=f"./mysmashmodel/mystats",
        )


def plot_stats(catchment, sbc=None, list_dir=None):
    if list_dir is None:
        list_dir = os.listdir(f"{working_path}/{catchment}")

    if sbc is None:
        sbc = sb.SmashBox()
        sbc.load_containers(
            f"{working_path}/{catchment}", load_catchment_polygon=False
        )

    for subdir in list_dir:
        res = getattr(sbc, subdir)
        res.myplot.plot_misfit(
            fig_settings={
                "figname": f"{working_path}/figures/misfit_{subdir}.pdf"
            }
        )
        res.myplot.plot_misfit_map(
            fig_settings={
                "figname": f"{working_path}/figures/misfit_map_{subdir}.pdf"
            }
        )
        res.myplot.plot_outlet_stats(
            fig_settings={
                "figname": f"{working_path}/figures/outlets_stats_{subdir}.pdf"
            }
        )


if mode == "init":
    init_models(
        sbc,
        catchment,
        start_time,
        end_time,
        warmup_time,
        val_start_time,
        val_end_time,
    )

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

if mode == "compute_stats":
    compute_stats(catchment)

if mode == "plot_stats":
    plot_stats(catchment)

if mode == "stats":
    del sbc
    list_dir = os.listdir(f"{working_path}/{catchment}")
    for subdir in list_dir:
        print(f"</> run stats for {working_path}/{catchment}/{subdir}")
        sbc = sb.SmashBox()
        sbc.load_containers(
            f"{working_path}/{catchment}",
            load_catchment_polygon=False,
            only=[subdir],
        )
        compute_stats(catchment, sbc, list_dir=[subdir])
        plot_stats(catchment, sbc, list_dir=[subdir])
        del sbc


# save figure


# configure the Gamma model for warmup (not used here)
# model_gamma_validation_warm = gamma.smashplug.ConfigureGammaWithSmash(
#     sbc.Cance_Gamma.warmup_model, dt=gamma_dt, **gamma_settings
# )
# model_gamma_validation_warm.routing_parameters.hc = (
#     optimized_gamma_model.routing_parameters.hc
# )
# model_gamma_validation_warm.routing_parameters.sc = (
#     optimized_gamma_model.routing_parameters.sc
# )

# gamma.smashplug.RunCoupledModel(sbc.Cance_Gamma.warmup_model, model_gamma_validation_warm)


# configure the Gamma model
# model_gamma_validation = gamma.smashplug.ConfigureGammaWithSmash(
#     sbc.Cance_Gamma.mysmashmodel.smash, dt=gamma_dt, **gamma_settings
# )
# model_gamma_validation.routing_parameters.hc = (
#     optimized_gamma_model.routing_parameters.hc
# )
# model_gamma_validation.routing_parameters.sc = (
#     optimized_gamma_model.routing_parameters.sc
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
#     optimized_gamma_model.routing_parameters.hc,
#     optimized_smash_model,
#     optimized_gamma_model,
# )
# sc = gamma.smashplug.GammaVectorsToSmashGrid(
#     optimized_gamma_model.routing_parameters.sc,
#     optimized_smash_model,
#     optimized_gamma_model,
# )


# fig, ax = sb.plot.plot.plot_image(
#     sc,
#     mask=optimized_smash_model.mesh.active_cell,
#     ax_settings={"title": "sc coefficients", "title_fontsize": 14},
#     fig_settings={
#         "figname": f"ExpCalValRoutage_sc_VS{model_gamma.routing_setup.varying_spread}_SU{model_gamma.routing_setup.sc_uniform}_CAL{ctrl}.pdf"
#     },
# )
# fig.show()

# fig, ax = sb.plot.plot.plot_image(
#     Grid_Hydraulic_Coef,
#     mask=optimized_smash_model.mesh.active_cell,
#     ax_settings={"title": "Hydraulic coefficients", "title_fontsize": 14},
#     fig_settings={
#         "figname": f"CoefHydraulics_VS{model_gamma.routing_setup.varying_spread}_SU{model_gamma.routing_setup.sc_uniform}_CAL{ctrl}.pdf"
#     },
# )
# fig.show()


# fig, ax = sb.plot.plot.plot_hydrograph(
#     sbc.Cance.mysmashmodel.smash,
#     plot_rainfall=True,
#     plot_settings_sim={"color": "blue", "label": "Smash Regional"},
# )
# fig, ax = sb.plot.plot.plot_hydrograph(
#     sbc.Cance_Gamma.mysmashmodel.smash,
#     figure=(fig, *ax),
#     plot_rainfall=False,
#     plot_qobs=False,
#     plot_settings_sim={"color": "green", "label": "Smash Regional + Gamma"},
# )
# fig, ax = sb.plot.plot.plot_hydrograph(
#     sbc.Cance_Gamma_calibrated.mysmashmodel.smash,
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


# plot_q(sbc.Cance.mysmashmodel.smash, optimized_gamma_model, g_node=380)
