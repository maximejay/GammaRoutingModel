import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma

import smash
import functions_smash_plot
import functions_smash_stats
import os

if "model" in locals():
    del model


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
    varying_spread=1,
)

# 4 nodes and 2 maximum upstream nodes
model.routing_mesh_init(nb_nodes=10)

# The controlloed nodes where obeervations are available : only one downstream control : node 4
model.routing_mesh_set_control(10)

# Create the mesh
model.routing_mesh.dx = 1000.0
model.routing_mesh_update()


# Simple twin experiment:
# creating an array of inflow
inflows = np.zeros(shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes))
inflows[0:2, 0] = 4.0

# Store the simulated discharge in an array: consider that is the true discharges, i.e our observation vector

# Initilaise the parameters
model.routing_parameters_init(hydraulics_coefficient=1.0, spreading=2.0)

# run the model (direct run)
model.run(inflows)

observations = np.zeros(shape=model.routing_results.discharges.shape)
observations[:, :] = model.routing_results.discharges[:, :].copy()
true_hc = model.routing_parameters.hydraulics_coefficient.copy()
true_sp = model.routing_parameters.spreading.copy()


# changing parameters
model.routing_parameters_change(hydraulics_coefficient=0.7, spreading=1.5)
back_hc = model.routing_parameters.hydraulics_coefficient.copy()
back_sp = model.routing_parameters.spreading.copy()

# run the model
model.run(inflows, states_init=0)

background_discharges = model.routing_results.discharges.copy()
# compute the cost function between the "observation" and the new simulated discharges
model.cost_function(observations, model.routing_results.discharges)

cost_initial = model.routing_results.costs.copy()

# calibrate the parameters to fit the "observed" discharges
model.calibration(inflows, observations, states_init=0)

optimal_discharges = model.routing_results.discharges.copy()
cost_final = model.routing_results.costs.copy()
optimal_hc = model.routing_parameters.hydraulics_coefficient.copy()
optimal_sp = model.routing_parameters.spreading.copy()

dir_results = os.path.join("src", "py", "test", "figures")
os.makedirs(dir_results, exist_ok=True)

from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="constrained", figsize=(10, 7))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.2)
ax1 = fig.add_subplot(gs[0, 0])  # Matrice 1
ax2 = fig.add_subplot(gs[0, 1])  # Matrice 2
ax3 = fig.add_subplot(gs[1, 0])  # Matrice 2
ax4 = fig.add_subplot(gs[1, 1])  # Matrice 2


time = np.arange(model.routing_setup.npdt)

ax1.plot(time, inflows[:, 0], label="Inflows", marker="", lw=2)
ax1.set_xlabel("Time-step (900s)")
ax1.set_ylabel("Discharges")
ax1.set_title("Upstream inflows")
ax1.legend()
ax1.grid(True)

ax2.plot(time, observations[:, 9], label="Downsrtream true discharges", marker="o", lw=3)
ax2.plot(time, background_discharges[:, 9], label="Downstream discharges (background)")
ax2.plot(time, optimal_discharges[:, 9], label="Downstream discharges (optimal)")
ax2.set_title("Downstream discharges")
ax2.set_xlabel("Time-step (900s)")
ax2.set_ylabel("Discharges")
ax2.legend()
ax2.grid(True)

cells = np.arange(model.routing_mesh.nb_nodes)
ax3.plot(cells, true_hc[:], label="True Hc")
ax3.plot(cells, back_hc[:], label="Background Hc")
ax3.plot(cells, optimal_hc[:], label="Optimal Hc")
ax3.set_xlabel("Cells (1000m)")
ax3.set_ylabel("Hc")
ax3.set_title("Coefficient Hc")
ax3.legend()
ax3.grid(True)

ax4.plot(cells, true_sp[:], label="True sp")
ax4.plot(cells, back_sp[:], label="Background sp")
ax4.plot(cells, optimal_sp[:], label="Optimal sp")
ax4.set_xlabel("Cells (1000m)")
ax4.set_ylabel("Sp")
ax4.set_title("Coefficient Sp")
ax4.legend()
ax4.grid(True)


plt.subplots_adjust(wspace=0.5)
# 7. Afficher le graphique
plt.show()

fig.savefig(os.path.join(dir_results, "Exp_Gamma_calibration.pdf"), bbox_inches="tight")
