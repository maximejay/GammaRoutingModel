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
model_smash = smash.Model(setup, mesh)

prcp = np.zeros(shape=(1, 9, model_smash.setup.ntime_step))
prcp[0, 0:5, 0:10] = 100.0
model_smash.atmos_data.prcp = prcp
model_smash_lr.atmos_data.prcp = prcp
model_smash.atmos_data.pet = 0.0
model_smash_lr.atmos_data.pet = 0.0


model_smash.rr_initial_states.values = 0.5
model_smash_lr.rr_initial_states.values = 0.5

model = gamma.Model()

# parametre du model
model.routing_setup_init(
    npdt=24,
    dt=3600.0,
    vmin=0.1,
    vmax=10.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.2,
    ponderation_regul=10000.0,
    velocity_computation="qm3",
    varying_spread=1,
)

model.routing_mesh_init(nb_nodes=9)

model.routing_mesh_set_control(9)

model.routing_mesh.dx = 1000.0
model.routing_mesh_update()

model.routing_parameters_init(hc=1.0, sc=2.0)

model_smash_lr.rr_parameters.values[:, :, 4] = 35
model_smash_lr.forward_run()

return_var = model_smash.forward_run(
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


qt = smash_time_grid_to_time_vector_active_cells(return_var.q_domain, model_smash)

model.run(qt.transpose())


model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
    model_smash, dt=3600.0, velocity_computation="qm3"
)
model_gamma.routing_parameters_change(hc=1.0, sc=2.0)
gamma.smashplug.RunCoupledModel(model_smash, model_gamma)
# model_gamma.run(qt.transpose())

plt.plot(model.routing_results.discharges[:, 8])
plt.plot(model_gamma.routing_results.discharges[:, 8])
plt.plot(model_smash_lr.response.q[0, :])

fig = plt.figure(layout="constrained")
fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(10, 5))
ax0, ax1, ax2 = axes[0], axes[1], axes[2]

cax0 = ax0.matshow(prcp.transpose(), cmap="viridis", vmin=np.min(prcp), vmax=np.max(prcp))
fig.colorbar(cax0, ax=ax0, fraction=0.046, pad=0.04, label="Prcp (mm)")

cax1 = ax1.matshow(qt.transpose(), cmap="viridis", vmin=np.min(qt), vmax=np.max(qt))
fig.colorbar(cax1, ax=ax1, fraction=0.046, pad=0.04, label="Inflows (m3/s)")

ax0.set_title("Prcp (mm)")  # TODO: Remplacer par votre titre
ax0.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax0.set_ylabel("Time-step (3600s)", fontsize=10)

ax1.set_title("Qt - Gamma Inflows")  # TODO: Remplacer par votre titre
ax1.set_xlabel("X (upstream->downstream cells)", fontsize=10)
ax1.set_ylabel("Time-step (3600s)", fontsize=10)

time = np.arange(model_smash_lr.response.q.shape[1])
ax2.plot(
    time,
    model.routing_results.discharges[:, 8],
    label="Gamma discharges (manual interface)",
    marker="o",
    lw=3,
)
ax2.plot(
    time,
    model_gamma.routing_results.discharges[:, 8],
    label="Gamma  discharges (Smashplug module)",
)
ax2.plot(time, model_smash_lr.response.q[0, :], label="Smash discharges with llr")
ax2.set_title("Downstream discharges")
ax2.set_xlabel("Time-step (3600s)")
ax2.set_ylabel("Discharges")
ax2.legend()
ax2.grid(True)
plt.subplots_adjust(wspace=0.5)

plt.show()

fig.savefig(os.path.join(dir_results, "Exp_coupling_base.pdf"), bbox_inches="tight")
