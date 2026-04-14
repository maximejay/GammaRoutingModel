#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:10:01 2026

@author: maxime
"""
import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma

import smash
import functions_smash_plot
import functions_smash_stats
import os
import rasterio
import re

dir_results = os.path.join("src", "py", "test", "figures")
os.makedirs(dir_results, exist_ok=True)


def read_ascii_data(filename):
    ascii_data = {}
    f = open(filename, "r")

    for line in f:

        strip_line = line.strip()
        split_line = strip_line.split()
        res = re.match("[a-zA-Z]*", split_line[0])
        if len(res.group()) > 0:
            ascii_data.update({split_line[0]: float(split_line[1])})
        else:
            break

    f.close()

    matrix = np.loadtxt(
        filename,
        dtype=float,
        comments="#",
        delimiter=None,
        converters=None,
        skiprows=6,
        usecols=None,
        unpack=False,
        ndmin=0,
        encoding=None,
        max_rows=None,
        quotechar=None,
        like=None,
    )

    ascii_data.update({"data": matrix})

    return ascii_data


def rasterio_write_tiff(filename, matrix, metadata):

    with rasterio.Env():
        with rasterio.open(filename, "w", compress="lzw", **metadata) as dst:
            dst.write(matrix, 1)


metadata = {
    "driver": "GTiff",
    "dtype": "float64",
    "nodata": None,
    "width": 10,
    "height": 2,
    "count": 1,
    "crs": None,
    "transform": rasterio.Affine(1000.0, 0.0, 1000.0, 0.0, -1000.0, 1000.0),
}

matrix = read_ascii_data("src/py/test/flwdir.ascii")
rasterio_write_tiff("src/py/test/flwdir.tif", np.atleast_2d(matrix["data"]), metadata)


mesh = smash.factory.generate_mesh(
    flwdir_path="src/py/test/flwdir.tif",
    x=[10000.0],
    y=[0.0],
    area=1e7,
    code=["gauge"],
    epsg="2154",
)

setup = {}
setup["hydrological_module"] = "gr4"
setup["routing_module"] = "lr"
setup["qobs_directory"] = "/home/maxime/DATA/QOBS_60min"
setup["prcp_directory"] = "/home/maxime/DATA/PLUIE"
setup["pet_directory"] = "/home/maxime/DATA/ETP-SFR-FRA-INTERA_L93"
setup["prcp_conversion_factor"] = 0.1
setup["read_prcp"] = False
setup["read_pet"] = False
setup["read_qobs"] = False
setup["daily_interannual_pet"] = False
setup["dt"] = 3600
setup["start_time"] = "2014-01-01 00:00"
setup["end_time"] = "2014-01-02 00:00"

model_smash_lr = smash.Model(setup, mesh)

setup["routing_module"] = "zero"

model_smash_true = smash.Model(setup, mesh)
model_smash_bck = smash.Model(setup, mesh)

prcp = np.zeros(shape=(1, 9, model_smash_true.setup.ntime_step))
prcp[0, 0:2, 0:5] = 50.0
model_smash_true.atmos_data.prcp = prcp
model_smash_lr.atmos_data.prcp = prcp
model_smash_bck.atmos_data.prcp = prcp
model_smash_true.atmos_data.pet = 0.0
model_smash_lr.atmos_data.pet = 0.0
model_smash_bck.atmos_data.pet = 0.0

model_smash_true.set_rr_parameters("cp", 200.0)
model_smash_true.set_rr_parameters("ct", 100.0)
model_smash_lr.set_rr_parameters("cp", 200.0)
model_smash_lr.set_rr_parameters("ct", 100.0)

model_smash_true.rr_initial_states.values = 0.5
model_smash_lr.rr_initial_states.values = 0.5
model_smash_bck.rr_initial_states.values = 0.5

model_smash_lr.rr_parameters.values[:, :, 4] = 35
model_smash_lr.forward_run()

return_var_true = model_smash_true.forward_run(
    return_options={"q_domain": True, "q_domain_kind_qt": True}
)

model_smash_bck.set_rr_parameters("cp", 100.0)
model_smash_bck.set_rr_parameters("ct", 50.0)

scenario = "1_sum"

return_var_bck = model_smash_bck.forward_run(
    return_options={"q_domain": True, "q_domain_kind_qt": True}
)


def smash_time_grid_to_time_vector_active_cells(grid, smash_model):

    vector = np.zeros(shape=(smash_model.mesh.nac, smash_model.setup.ntime_step))

    k = 0
    for col in range(smash_model.mesh.ncol):

        for row in range(smash_model.mesh.nrow):

            if smash_model.mesh.active_cell[row, col] > 0:

                vector[k, :] = grid[row, col, :]
                k = k + 1

    return vector


qt_true = smash_time_grid_to_time_vector_active_cells(
    return_var_true.q_domain, model_smash_true
)
qt_bck = smash_time_grid_to_time_vector_active_cells(
    return_var_bck.q_domain, model_smash_bck
)  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:10:01 2026

@author: maxime
"""
import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma

import smash
import functions_smash_plot
import functions_smash_stats
import os
import rasterio
import re

dir_results = os.path.join("src", "py", "test", "figures")
os.makedirs(dir_results, exist_ok=True)


def read_ascii_data(filename):
    ascii_data = {}
    f = open(filename, "r")

    for line in f:

        strip_line = line.strip()
        split_line = strip_line.split()
        res = re.match("[a-zA-Z]*", split_line[0])
        if len(res.group()) > 0:
            ascii_data.update({split_line[0]: float(split_line[1])})
        else:
            break

    f.close()

    matrix = np.loadtxt(
        filename,
        dtype=float,
        comments="#",
        delimiter=None,
        converters=None,
        skiprows=6,
        usecols=None,
        unpack=False,
        ndmin=0,
        encoding=None,
        max_rows=None,
        quotechar=None,
        like=None,
    )

    ascii_data.update({"data": matrix})

    return ascii_data


def rasterio_write_tiff(filename, matrix, metadata):

    with rasterio.Env():
        with rasterio.open(filename, "w", compress="lzw", **metadata) as dst:
            dst.write(matrix, 1)


metadata = {
    "driver": "GTiff",
    "dtype": "float64",
    "nodata": None,
    "width": 10,
    "height": 2,
    "count": 1,
    "crs": None,
    "transform": rasterio.Affine(1000.0, 0.0, 1000.0, 0.0, -1000.0, 1000.0),
}

matrix = read_ascii_data("src/py/test/flwdir.ascii")
rasterio_write_tiff("src/py/test/flwdir.tif", np.atleast_2d(matrix["data"]), metadata)


mesh = smash.factory.generate_mesh(
    flwdir_path="src/py/test/flwdir.tif",
    x=[10000.0],
    y=[0.0],
    area=1e7,
    code=["gauge"],
    epsg="2154",
)

setup = {}
setup["hydrological_module"] = "gr4"
setup["routing_module"] = "lr"
setup["qobs_directory"] = "/home/maxime/DATA/QOBS_60min"
setup["prcp_directory"] = "/home/maxime/DATA/PLUIE"
setup["pet_directory"] = "/home/maxime/DATA/ETP-SFR-FRA-INTERA_L93"
setup["prcp_conversion_factor"] = 0.1
setup["read_prcp"] = False
setup["read_pet"] = False
setup["read_qobs"] = False
setup["daily_interannual_pet"] = False
setup["dt"] = 3600
setup["start_time"] = "2014-01-01 00:00"
setup["end_time"] = "2014-01-02 00:00"

model_smash_lr = smash.Model(setup, mesh)

setup["routing_module"] = "zero"

model_smash_true = smash.Model(setup, mesh)
model_smash_bck = smash.Model(setup, mesh)

prcp = np.zeros(shape=(1, 9, model_smash_true.setup.ntime_step))
prcp[0, 0:2, 0:5] = 50.0
model_smash_true.atmos_data.prcp = prcp
model_smash_lr.atmos_data.prcp = prcp
model_smash_bck.atmos_data.prcp = prcp
model_smash_true.atmos_data.pet = 0.0
model_smash_lr.atmos_data.pet = 0.0
model_smash_bck.atmos_data.pet = 0.0

model_smash_true.set_rr_parameters("cp", 200.0)
model_smash_true.set_rr_parameters("ct", 100.0)
model_smash_lr.set_rr_parameters("cp", 200.0)
model_smash_lr.set_rr_parameters("ct", 100.0)

model_smash_true.rr_initial_states.values = 0.5
model_smash_lr.rr_initial_states.values = 0.5
model_smash_bck.rr_initial_states.values = 0.5

model_smash_lr.rr_parameters.values[:, :, 4] = 35
model_smash_lr.forward_run()

return_var_true = model_smash_true.forward_run(
    return_options={"q_domain": True, "q_domain_kind_qt": True}
)

model_smash_bck.set_rr_parameters("cp", 100.0)
model_smash_bck.set_rr_parameters("ct", 50.0)

scenario = "1_sum"

return_var_bck = model_smash_bck.forward_run(
    return_options={"q_domain": True, "q_domain_kind_qt": True}
)


def smash_time_grid_to_time_vector_active_cells(grid, smash_model):

    vector = np.zeros(shape=(smash_model.mesh.nac, smash_model.setup.ntime_step))

    k = 0
    for col in range(smash_model.mesh.ncol):

        for row in range(smash_model.mesh.nrow):

            if smash_model.mesh.active_cell[row, col] > 0:

                vector[k, :] = grid[row, col, :]
                k = k + 1

    return vector


qt_true = smash_time_grid_to_time_vector_active_cells(
    return_var_true.q_domain, model_smash_true
)
qt_bck = smash_time_grid_to_time_vector_active_cells(
    return_var_bck.q_domain, model_smash_bck
)


model_gamma_true = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash_true,
    dt=3600.0,
    velocity_computation="qmm",
    varying_spread=1,
    ponderation_regul=1000.0,
)

model_gamma_true.routing_parameters_change(hc=1.0, sc=3.0)
model_gamma_true.routing_memory.states_init = 1.0
model_gamma_true.run(qt_true.transpose())

model_gamma_background = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash_bck,
    dt=3600.0,
    velocity_computation="qmm",
    varying_spread=1,
    ponderation_regul=1000.0,
)

model_gamma_background.routing_parameters_change(hc=1.0, sc=3.0)
model_gamma_background.routing_memory.states_init = 1.0
model_gamma_background.run(qt_bck.transpose())

model_gamma_optimize = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash_bck,
    dt=3600.0,
    velocity_computation="qmm",
    varying_spread=0,
    ponderation_regul=1000.0,
    ponderation_cost=100.0,
)

model_gamma_optimize.routing_parameters_change(hc=1.0, sc=4.0)
model_gamma_optimize.routing_memory.states_init = 1.0
model_gamma_optimize.run(qt_bck.transpose())

model_gamma_optimize.routing_mesh_set_control(9)

model_smash_bck.response_data.q[0, :] = model_gamma_true.routing_results.discharges[:, 8]


# compute the initial cost
model_gamma_optimize.cost_function(
    model_gamma_true.routing_results.discharges,
    model_gamma_optimize.routing_results.discharges,
)

model_gamma_optimize.routing_results.costs


model_gamma_true = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash_true,
    dt=3600.0,
    velocity_computation="qmm",
    varying_spread=1,
    ponderation_regul=1000.0,
)

model_gamma_true.routing_parameters_change(hc=1.0, sc=3.6)
model_gamma_true.routing_memory.states_init = 0.1
model_gamma_true.run(qt_true.transpose())

model_gamma_background = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash_bck,
    dt=3600.0,
    velocity_computation="qmm",
    varying_spread=1,
    ponderation_regul=1000.0,
)

model_gamma_background.routing_parameters_change(hc=1.0, sc=3.6)
model_gamma_background.routing_memory.states_init = 0.1
model_gamma_background.run(qt_bck.transpose())

model_gamma_optimize = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash_bck,
    dt=3600.0,
    velocity_computation="qmm",
    varying_spread=1,
    ponderation_regul=1.0,
    ponderation_cost=1.0,
)

model_gamma_optimize.routing_parameters_change(hc=1.0, sc=3.6)
model_gamma_optimize.routing_memory.states_init = 0.1
model_gamma_optimize.run(qt_bck.transpose())

model_gamma_optimize.routing_mesh_set_control(9)

model_smash_bck.response_data.q[0, :] = model_gamma_true.routing_results.discharges[:, 8]


# compute the initial cost
model_gamma_optimize.cost_function(
    model_gamma_true.routing_results.discharges,
    model_gamma_optimize.routing_results.discharges,
)

model_gamma_optimize.routing_results.costs

res, model_smash_calibrated, model_gamma_calibrated = (
    gamma.smashplug.calibration_parameters_smash(
        model_smash_bck,
        model_gamma_optimize,
        bounds={
            "cp": [
                0.1,
                500.0,
            ],
            "ct": [0.1, 500.0],
        },
        observations=model_gamma_true.routing_results.discharges,
    )
)

model_smash_calibrated.rr_parameters.values


# BestControlVector, optimized_smash_model, model_gamma_optimize = (
#     gamma.smashplug.OptimizeCoupledModel(
#         model_smash_bck,
#         model_gamma_optimize,
#         model_gamma_true.routing_results.discharges,
#         control_parameters_list=["cp", "ct"],
#         bounds={
#             "cp": [
#                 0.1,
#                 500.0,
#             ],
#             "ct": [0.1, 500.0],
#         },
#         maxiter=30,
#         maxfun=100,
#         tol=None,
#         ScaleGradients=False,
#         ScaleGammaGradientsBySurface=False,
#         optim_type="local",
#         local_optimizer="L-BFGS-B",  # L-BFGS-B | trust-constr | SLSQP | TNC
#     )
# )


# fig = plt.figure(layout="constrained")
# fig, axes = plt.subplots(ncols=4, nrows=1, figsize=(15, 5))
# ax0, ax1, ax2, ax3 = axes[0], axes[1], axes[2], axes[3]

from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="constrained", figsize=(20, 10))
gs = GridSpec(
    3, 4, figure=fig, height_ratios=[1, 1, 1], width_ratios=[1, 1, 1, 1], hspace=0.2
)
ax0 = fig.add_subplot(gs[0, 2])  # Matrice 1
ax1 = fig.add_subplot(gs[0, 3])  # Matrice 1
ax8 = fig.add_subplot(gs[1:3, 2:4])  # Matrice 2

ax2 = fig.add_subplot(gs[0, 0])  # Matrice 2
ax3 = fig.add_subplot(gs[1, 0])  # Courbes (s'étend sur les 3 colonnes)
ax4 = fig.add_subplot(gs[2, 0])  # Courbes (s'étend sur les 3 colonnes)

ax5 = fig.add_subplot(gs[0, 1])  # Matrice 2
ax6 = fig.add_subplot(gs[1, 1])  # Courbes (s'étend sur les 3 colonnes)
ax7 = fig.add_subplot(gs[2, 1])  # Courbes (s'étend sur les 3 colonnes)

cax0 = ax0.matshow(
    prcp.transpose(),
    cmap="viridis",
    vmin=np.min(prcp),
    vmax=np.max(prcp),
    aspect="auto",
)
fig.colorbar(cax0, ax=ax0, fraction=0.046, pad=0.04, label="Prcp (mm)")

cax1 = ax1.matshow(
    qt_true.transpose(),
    cmap="viridis",
    vmin=np.min(qt_true),
    aspect="auto",
)
fig.colorbar(cax1, ax=ax1, fraction=0.046, pad=0.04, label="Inflows (True) (m3/s)")

vmin = min(
    np.min(model_smash_true.get_rr_parameters("cp")),
    np.min(model_smash_bck.get_rr_parameters("cp")),
    np.min(model_smash_calibrated.get_rr_parameters("cp")),
)
vmax = max(
    np.max(model_smash_true.get_rr_parameters("cp")),
    np.max(model_smash_bck.get_rr_parameters("cp")),
    np.max(model_smash_calibrated.get_rr_parameters("cp")),
)

cax2 = ax2.matshow(
    model_smash_true.get_rr_parameters("cp"),
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
)
fig.colorbar(cax2, ax=ax2, fraction=0.046, pad=0.04, label="Cp (True)")

cax3 = ax3.matshow(
    model_smash_bck.get_rr_parameters("cp"),
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
)
fig.colorbar(cax3, ax=ax3, fraction=0.046, pad=0.04, label="Cp (background)")

cax4 = ax4.matshow(
    model_smash_calibrated.get_rr_parameters("cp"),
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
)
fig.colorbar(cax4, ax=ax4, fraction=0.046, pad=0.04, label="Cp (calibrated)")

vmin = min(
    np.min(model_smash_true.get_rr_parameters("ct")),
    np.min(model_smash_bck.get_rr_parameters("ct")),
    np.min(model_smash_calibrated.get_rr_parameters("ct")),
)
vmax = max(
    np.max(model_smash_true.get_rr_parameters("ct")),
    np.max(model_smash_bck.get_rr_parameters("ct")),
    np.max(model_smash_calibrated.get_rr_parameters("ct")),
)

cax5 = ax5.matshow(
    model_smash_true.get_rr_parameters("ct"),
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
)
fig.colorbar(cax5, ax=ax5, fraction=0.046, pad=0.04, label="Ct (True)")

cax6 = ax6.matshow(
    model_smash_bck.get_rr_parameters("ct"),
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
)
fig.colorbar(cax6, ax=ax6, fraction=0.046, pad=0.04, label="Ct (background)")

cax7 = ax7.matshow(
    model_smash_calibrated.get_rr_parameters("ct"),
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
)
fig.colorbar(cax7, ax=ax7, fraction=0.046, pad=0.04, label="Ct (calibrated)")

ax2.set_title("Smash Cp (true - bck - optim)")
ax5.set_title("Smash Ct parameters (true - bck - optim)")

ax2.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax3.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax4.set_xlabel("X (upstream->downstream cells)", fontsize=10)

ax5.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax6.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax7.set_xlabel("X (upstream->downstream cells)", fontsize=10)


ax0.set_title("Prcp (mm)")  # TODO: Remplacer par votre titre
ax0.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax0.set_ylabel("Time-step (3600s)", fontsize=10)

ax1.set_title("Qt - Gamma Inflows True")  # TODO: Remplacer par votre titre
ax1.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax1.set_ylabel("Time-step (3600s)", fontsize=10)

time = np.arange(model_smash_lr.response.q.shape[1])
nodes = np.arange(model_gamma_true.routing_mesh.nb_nodes)

ax8.plot(
    time,
    model_gamma_true.routing_results.discharges[:, 8],
    label="Gamma  discharges (True)",
    lw=3,
    marker="o",
)
ax8.plot(
    time,
    model_gamma_background.routing_results.discharges[:, 8],
    label="Gamma  discharges (Background)",
)
ax8.plot(
    time,
    model_gamma_calibrated.routing_results.discharges[:, 8],
    label="Gamma  discharges (optimized)",
)
ax8.plot(time, model_smash_lr.response.q[0, :], label="True Smash discharges with llr")

ax8.set_title("Downstream discharges")
ax8.set_xlabel("Time-step (3600s)")
ax8.set_ylabel("Discharges (m3/s)")
ax8.legend()
ax8.grid(True)

plt.subplots_adjust(wspace=0.7)

# plt.show()

# fig.savefig(
#     os.path.join(dir_results, f"Exp_coupling_calibration_smash_{scenario}.pdf"),
#     bbox_inches="tight",
# )
