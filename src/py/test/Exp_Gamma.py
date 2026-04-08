import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma

import smash
import functions_smash_plot
import functions_smash_stats
import os

if "model" in locals():
    del model

# Simple jonction test case for the Gamma routing model
# Here is the meshing o1-o4 are the nodes

# ~ !Confluence test case
# ~ !       o2
# ~ !  o1  /
# ~ !  \  /
# ~ !   o3
# ~ !   |
# ~ !   o4
# ~ !

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

# creation du maillage (Manual): Gamma have no meshing capabilities, you can use the interface to smash to build a mesh from a regular grid

# 4 nodes and 2 maximum upstream nodes
model.routing_mesh_init(nb_nodes=4, nb_upstream_nodes=2)
# Nodes computation order (node number)
model.routing_mesh.upstream_to_downstream_nodes = np.array([1, 2, 3, 4])
# Linker between one node with the upstream ones
model.routing_mesh.nodes_linker = 0
model.routing_mesh.nodes_linker[:, 2] = np.array(
    [1, 2]
)  # the third node (index=2) has upstream node 1 and 2 (index are different than node number ! This is due to the communication between fortran and Python)
model.routing_mesh.nodes_linker[1, 3] = 3  # the downstream node o4 targets o3

model.routing_mesh.surface = np.array(
    [1, 1, 1, 1]
)  # the drained surface in km2 à set to 1km (no influence). It is only important for the velocity computation vs Qspe
model.routing_mesh.dx = np.array(
    [1000.0, 7000.0, 2000.0, 1000.0]
)  # Distance vers le noeud aval. A long bief slows the model (Tabulation of the Gamma routing coefficients is long because the delay is big...). For better performances, we should split the bief.

# The controlloed nodes where obeervations are available : only one downstream control : node 4
model.routing_mesh_set_control(4)

# Create the mesh
model.routing_mesh_update()

# Testing combinaison of differents routing parameters
Xi = [0.5, 1.5]  # hydraulics coef
S = [1.0, 2.0]  # sc coef

dir_results = os.path.join("src", "py", "test", "figures")
os.makedirs(dir_results, exist_ok=True)

for i in range(2):

    for j in range(2):

        param_xi = Xi[i]
        param_S = S[j]
        # Initilaise the parameters
        model.routing_parameters_init(hc=param_xi, sc=param_S)

        # creating an array of inflow
        inflows = np.zeros(shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes))
        inflows[0:2, 0] = 2.0
        inflows[0:2, 1] = 4.0

        # run the model (direct run)
        model.run(inflows)

        fig, ax = plt.subplots()
        plt.rc("font", size=18)
        ax.plot(
            model.routing_results.discharges[:, 0], label="Gamma at node 01", lw=2
        )  # nodes are in Fortran index: in python array we need -1: how to fix that ?
        ax.plot(
            model.routing_results.discharges[:, 1], label="Gamma at node 02", lw=2
        )  # nodes are in Fortran index: in python array we need -1: how to fix that ?
        ax.plot(
            model.routing_results.discharges[:, 2], label="Gamma at node 03", lw=2
        )  # nodes are in Fortran index: in python array we need -1: how to fix that ?
        ax.plot(
            model.routing_results.discharges[
                :, model.routing_mesh.controlled_nodes[0] - 1
            ],
            label="Gamma at node 04",
            lw=2,
        )  # nodes are in Fortran index: in python array we need -1: how to fix that ?

        ax.legend(loc="upper left")
        ax.axes.grid(True, alpha=0.7, ls="--")
        ax.axes.set_xlabel("Time-Step (hours)")
        ax.axes.set_ylabel("Discharges (m^3/s)")
        ax.set_title(f"Hydraulics Coef={Xi[i]}, sc={S[j]}")
        fig.show()
        plot = (fig, ax)
        functions_smash_plot.save_figure(
            plot,
            figname=os.path.join(dir_results, f"Exp_Gamma_confluence_{i}{j}.pdf"),
            xlim=[0, 40],
            xsize=12,
            ysize=10,
        )


# Simple twin experiment:
# creating an array of inflow
inflows = np.zeros(shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes))
inflows[0:2, 0] = 2.0
inflows[0:2, 1] = 4.0
# Store the simulated discharge in an array: consider that is the true discharges, i.e our observation vector
# Initilaise the parameters
model.routing_parameters_init(hc=1.0, sc=2.0)
# run the model (direct run)
model.run(inflows)

observations = np.zeros(shape=model.routing_results.discharges.shape)
observations[:, :] = model.routing_results.discharges[:, :]
true_hc = model.routing_parameters.hc.copy()
true_sp = model.routing_parameters.sc.copy()


# changing parameters
model.routing_parameters_change(hc=0.7, sc=1.4)
back_hc = model.routing_parameters.hc.copy()
back_sp = model.routing_parameters.sc.copy()

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
optimal_hc = model.routing_parameters.hc.copy()
optimal_sp = model.routing_parameters.sc.copy()


from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="constrained", figsize=(10, 7))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.2)
ax1 = fig.add_subplot(gs[0, 0])  # Matrice 1
ax2 = fig.add_subplot(gs[0, 1])  # Matrice 2
ax3 = fig.add_subplot(gs[1, 0])  # Matrice 2
ax4 = fig.add_subplot(gs[1, 1])  # Matrice 2


time = np.arange(model.routing_setup.npdt)

ax1.plot(time, inflows[:, 0], label="Inflows node 0", marker="", lw=2)
ax1.plot(time, inflows[:, 1], label="Inflows node 1", marker="", lw=2)
ax1.set_xlabel("Time-step (900s)")
ax1.set_ylabel("Discharges")
ax1.set_title("Upstream inflows")
ax1.legend()
ax1.grid(True)

ax2.plot(time, observations[:, 3], label="Downsrtream true discharges", marker="o", lw=3)
ax2.plot(time, background_discharges[:, 3], label="Downstream discharges (background)")
ax2.plot(time, optimal_discharges[:, 3], label="Downstream discharges (optimal)")
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

fig.savefig(
    os.path.join(dir_results, "Exp_Gamma_confluence_calibration.pdf"), bbox_inches="tight"
)
