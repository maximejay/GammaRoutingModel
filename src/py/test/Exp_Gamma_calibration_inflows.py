#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:55:30 2026

@author: maxime
"""

import numpy as np
import matplotlib.pyplot as plt
import gamma_routing_model as gamma
import os


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

inflows_origin = np.zeros(shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes))

inflows_origin[0:2, 0:10] = 4.0
# inflows[1:2, 0:10] = 2.0

model.run(inflows_origin)

discharges_origin = model.routing_results.discharges.copy()

qobs = model.routing_results.discharges.copy()

inflows_background = np.zeros(
    shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes)
)
inflows_background[0:2, 0:10] = 2.0

model.run(inflows_background)

discharges_background = model.routing_results.discharges.copy()

model.cost_function(qobs, model.routing_results.discharges)

model.routing_results.costs

# cost, gradient = model.run_backward_djdq(
#     inflows_background.flatten(), inflows_background, qobs
# )

res = model.calibration_inflows(inflows_background, qobs)

discharges_optimal = model.routing_results.discharges.copy()


inflows_optimal = np.zeros(shape=(model.routing_setup.npdt, model.routing_mesh.nb_nodes))
k = 0
for i in range(inflows_optimal.shape[0]):
    inflows_optimal[i, :] = res.x[k : k + inflows_optimal.shape[1]]
    k = k + inflows_optimal.shape[1]

# plt.imshow(inflows[0:10, :])

# 1. Créer une matrice numpy (exemple)
# data = np.random.rand(10, 10)  # TODO: Remplacer par votre matrice

dir_results = os.path.join("src", "py", "test", "figures")
os.makedirs(dir_results, exist_ok=True)

minval = min(inflows_origin.min(), inflows_background.min(), inflows_optimal.min())
maxval = max(inflows_origin.max(), inflows_background.max(), inflows_optimal.max())

# 2. Créer la figure et l'axe
# fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(10, 5))
# ax1, ax2, ax3 = axes[0, 0], axes[0, 1], axes[0, 2]
# ax4 = axes[1, 1]
# fig.delaxes(axes[1, 0])  # Supprimer le 4ème et 6eme sous-graphique
# fig.delaxes(axes[1, 2])  # Supprimer le 4ème et 6eme sous-graphique

from matplotlib.gridspec import GridSpec

fig = plt.figure(layout="constrained")
gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 1], width_ratios=[1, 1, 1], hspace=0.2)
ax1 = fig.add_subplot(gs[0, 0])  # Matrice 1
ax2 = fig.add_subplot(gs[0, 1])  # Matrice 2
ax3 = fig.add_subplot(gs[0, 2])  # Matrice 2
ax4 = fig.add_subplot(gs[1, :])  # Courbes (s'étend sur les 3 colonnes)


# 3. Tracer la matrice avec une échelle de couleur
cax1 = ax1.matshow(
    inflows_origin[0:10, :], cmap="viridis", vmin=minval, vmax=maxval
)  # TODO: Choisir une autre colormap si besoin
cax2 = ax2.matshow(
    inflows_background[0:10, :], cmap="viridis", vmin=minval, vmax=maxval
)  # TODO: Choisir une autre colormap si besoin
cax3 = ax3.matshow(
    inflows_optimal[0:10, :], cmap="viridis", vmin=minval, vmax=maxval
)  # TODO: Choisir une autre colormap si besoin


# 4. Ajouter une barre de couleur (échelle)
fig.colorbar(
    cax1, ax=ax1, fraction=0.046, pad=0.04, label="Inflows true (m3/s)"
)  # TODO: Changer le label si nécessaire
# 4. Ajouter une barre de couleur (échelle)
fig.colorbar(
    cax2, ax=ax2, fraction=0.046, pad=0.04, label="Inflows background (m3/s)"
)  # TODO: Changer le label si nécessaire
# 4. Ajouter une barre de couleur (échelle)
fig.colorbar(
    cax3, ax=ax3, fraction=0.046, pad=0.04, label="Inflows optimal (m3/s)"
)  # TODO: Changer le label si nécessaire


# 5. Ajouter un titre
ax1.set_title("True inflows")  # TODO: Remplacer par votre titre
ax2.set_title("Background inflows")  # TODO: Remplacer par votre titre
ax3.set_title("Optimal inflows")  # TODO: Remplacer par votre titre

# 6. Ajouter des légendes pour les axes
ax1.set_xlabel(
    "X (upstream->downstream cells)", fontsize=10
)  # TODO: Remplacer par votre légende
ax1.set_ylabel("Time-step (900q)", fontsize=10)  # TODO: Remplacer par votre légende
ax2.set_xlabel(
    "X (upstream->downstream cells)", fontsize=10
)  # TODO: Remplacer par votre légende
ax2.set_ylabel("Time-step (900q)", fontsize=10)  # TODO: Remplacer par votre légende
ax3.set_xlabel(
    "X (upstream->downstream cells)", fontsize=10
)  # TODO: Remplacer par votre légende
ax3.set_ylabel("Time-step (900q)", fontsize=10)  # TODO: Remplacer par votre légende

time = np.arange(discharges_origin.shape[0])
# 4. Tracer les 3 courbes temporelles sur ax3
ax4.plot(time, discharges_origin[:, 9], label="Origin discharges", marker="o", lw=3)
ax4.plot(time, discharges_background[:, 9], label="background discharges")
ax4.plot(time, discharges_optimal[:, 9], label="Optimal discharges")
ax4.set_title("Downstream discharges")
ax4.set_xlabel("Time-step (900s)")
ax4.set_ylabel("Discharges")
ax4.legend()
ax4.grid(True)

plt.subplots_adjust(wspace=0.5)
# 7. Afficher le graphique
plt.show()

fig.savefig(
    os.path.join(dir_results, "Exp_Gamma_calibration_inflows.pdf"), bbox_inches="tight"
)
