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
setup["routing_module"] = "zero"
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

model_smash = smash.Model(setup, mesh)

prcp = np.zeros(shape=(1, 9, model_smash.setup.ntime_step))
prcp[0:5, 0:10] = 1000.0
model_smash.atmos_data.prcp = prcp
model_smash.rr_initial_states.values = 0.5

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

model.routing_parameters_init(hydraulics_coefficient=1.0, spreading=2.0)

return_var = model_smash.forward_run(return_options={"q_domain_kind": "qt"})

qt = model_smash.response.qac.copy()


model.run(qt.transpose())

plt.plot(model.routing_results.discharges[:, 8])
