import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma

import smash
import functions_smash_plot
import functions_smash_stats
import os

if "model" in locals():
    del model

nb_nodes = 10
dt = 900.0
dx = 1000.0

model = gamma.Model()

# parametre du model
model.routing_setup_init(
    npdt=100,
    dt=dt,
    vmin=0.01,
    vmax=10.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.1,
    ponderation_regul=1.0,
    velocity_computation="qm3",
    varying_spread=0,
    sc_bnd=[0.5, 20.0],
)

# 4 nodes and 2 maximum upstream nodes
model.routing_mesh_init(nb_nodes=nb_nodes)

# The controlloed nodes where obeervations are available : only one downstream control : node 4
model.routing_mesh_set_control(nb_nodes)

# Create the mesh
model.routing_mesh.dx = dx
model.routing_mesh_update()


# Simple twin experiment:
# creating an array of inflow
inflows = np.zeros(shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes))
inflows[0:2, 0] = 4.0

# Store the simulated discharge in an array: consider that is the true discharges, i.e our observation vector

# Initilaise the parameters
model.routing_parameters_init(hc=1.0, sc=0.6)

# run the model (direct run)
model.run(inflows)


# Initilaise the parameters
model.routing_parameters_init(hc=1.0, sc=20.0)
# run the model (direct run)
model.run(inflows)

model.routing_states.tabulated_routing_coef


volume = {}
volume.update({"inflows": np.round(np.sum(inflows[:, :]), 3)})

hc_cst = 1.0
discharges = {}
for sc in range(6, 10, 1):
    model.routing_parameters_change(hc=hc_cst, sc=float(sc / 10))
    model.run(inflows, states_init=False)
    discharges.update({str(sc / 10): model.routing_results.discharges.copy()})
    volume.update(
        {
            str(sc / 10): np.round(
                np.sum(model.routing_results.discharges[:, nb_nodes - 1]), 3
            )
        }
    )

for sc in range(1, 21, 1):
    model.routing_parameters_change(hc=hc_cst, sc=float(sc))
    model.run(inflows, states_init=False)
    discharges.update({sc: model.routing_results.discharges.copy()})
    volume.update(
        {str(sc): np.round(np.sum(model.routing_results.discharges[:, nb_nodes - 1]), 3)}
    )


volume_hc = {}
volume_hc.update({"inflows": np.round(np.sum(inflows[:, :]), 3)})

sc_cst = 2.0
discharges_hc = {}
for hc in range(2, 10, 2):
    model.routing_parameters_change(hc=float(hc / 10), sc=sc_cst)
    model.run(inflows, states_init=False)
    discharges_hc.update({str(hc / 10): model.routing_results.discharges.copy()})
    volume_hc.update(
        {
            str(hc / 10): np.round(
                np.sum(model.routing_results.discharges[:, nb_nodes - 1]), 3
            )
        }
    )

for hc in range(1, 10, 1):
    model.routing_parameters_change(hc=hc, sc=sc_cst)
    model.run(inflows, states_init=False)
    discharges_hc.update({hc: model.routing_results.discharges.copy()})
    volume_hc.update(
        {str(hc): np.round(np.sum(model.routing_results.discharges[:, nb_nodes - 1]), 3)}
    )


dir_results = os.path.join("src", "py", "test", "figures")
os.makedirs(dir_results, exist_ok=True)

from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="constrained", figsize=(15, 10))
gs = GridSpec(
    2,
    3,
    figure=fig,
    height_ratios=[0.7, 0.3],
    width_ratios=[0.3, 1.0, 1.0],
    hspace=0.2,
)
ax1 = fig.add_subplot(gs[0, 0])  # Matrice 1
ax2 = fig.add_subplot(gs[0, 1])  # Matrice 2
ax3 = fig.add_subplot(gs[1, 0:2])  # Matrice 2
ax5 = fig.add_subplot(gs[1, 2])  # Matrice 2
ax4 = fig.add_subplot(gs[0, 2])  # Matrice 2

time = np.arange(model.routing_setup.npdt)

ax1.plot(time, inflows[:, 0], label="Inflows", marker="", lw=2)
ax1.set_xlabel("Time-step (900s)")
ax1.set_ylabel("Discharges")
ax1.set_title("Upstream inflows")
ax1.legend()
ax1.grid(True)

for key, value in discharges.items():
    ax2.plot(time, value[:, nb_nodes - 1], label=f"Dws discharges sc={key}")
ax2.set_title(
    f"Downstream discharges at node {nb_nodes} with dx={dx}\n for various sc values and hc={hc_cst}"
)
ax2.set_xlabel("Time-step (900s)")
ax2.set_ylabel("Discharges")
ax2.legend()
ax2.grid(True)

labels = list(volume.keys())  # ['A', 'B', 'C']
values = list(volume.values())  # [10, 20, 15]
ax3.bar(labels, values)

for i, v in enumerate(values):
    ax3.text(i, np.round(v - 1.0, 3), str(v), ha="center", va="bottom")

ax3.set_xlabel("")
ax3.set_ylabel("Volume en m3")
ax3.set_title(
    "Volume d'eau injecter (inflows) et router à l'exutoire aval\n pour chaque valeur de sc"
)


for key, value in discharges_hc.items():
    ax4.plot(time, value[:, nb_nodes - 1], label=f"Dws discharges hc={key}")

ax4.set_title(
    f"Downstream discharges at node {nb_nodes} with dx={dx}\n for various hc values and sc={sc_cst}"
)
ax4.set_xlabel("Time-step (900s)")
ax4.set_ylabel("Discharges")
ax4.legend()
ax4.grid(True)


labels = list(volume_hc.keys())  # ['A', 'B', 'C']
values = list(volume_hc.values())  # [10, 20, 15]
ax5.bar(labels, values)

for i, v in enumerate(values):
    ax5.text(i, np.round(v - 1.0, 3), str(v), ha="center", va="bottom")

ax5.set_xlabel("")
ax5.set_ylabel("Volume en m3")
ax5.set_title(
    "Volume d'eau injecter (inflows) et router à l'exutoire aval\n pour chaque valeur de hc"
)

plt.subplots_adjust(wspace=0.5, hspace=0.5)
# 7. Afficher le graphique
# plt.show()

# fig.savefig(os.path.join(dir_results, "Exp_Gamma_conservation.pdf"), bbox_inches="tight")
