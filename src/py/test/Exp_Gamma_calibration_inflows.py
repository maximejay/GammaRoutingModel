#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:55:30 2026

@author: maxime
"""

import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma


model = gamma.Model()

# parametre du model
model.routing_setup_init(
    npdt=100,
    dt=900.0,
    vmin=0.1,
    vmax=10.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.2,
    ponderation_regul=10000.0,
    velocity_computation="qm3",
    varying_spread=0,
)

model.routing_mesh_init(nb_nodes=10)

model.routing_mesh_set_control(10)

model.routing_mesh.dx = 1000.0
model.routing_mesh_update()

model.routing_parameters_init(hydraulics_coefficient=1.0, spreading=2.0)

inflows = np.zeros(shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes))

inflows[0:2, 0:10] = 4.0
# inflows[1:2, 0:10] = 2.0

model.run(inflows)

qobs = model.routing_results.discharges.copy()

inflows[0:2, 0:10] = 2.0

model.run(inflows)

model.cost_function(qobs, model.routing_results.discharges)

model.routing_results.costs

cost, gradient = model.run_backward_djdq(inflows.flatten(), inflows, qobs)

res = model.calibration_inflows(inflows, qobs)

k = 0
for i in range(inflows.shape[0]):
    inflows[i, :] = res.x[k : k + inflows.shape[1]]
    k = k + inflows.shape[1]

plt.imshow(inflows[0:10, :])
