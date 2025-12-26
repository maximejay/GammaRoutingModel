import numpy as np
import matplotlib.pyplot as plt
import smash

import gamma_routing_model as gamma
import scipy


# import functions_smash_plot
# import functions_smash_stats

# Smash reference model structure gr-a
setup_cance, mesh_cance = smash.factory.load_dataset("Cance")
setup_cance["routing_module"] = "rm_zero"
smash_model = smash.Model(setup_cance, mesh_cance)
smash_model.forward_run()

# smash_model.response.qt.shape

# smash_model.optimize(mapping="distributed")

# smash compilé avec les options de debuggage
# At line 11781 of file ../smash/fcore/forward/forward_openmp_db.f90
# Fortran runtime error: Assignment of scalar to unallocated array
# parameters_b%control%x = 0.0_4
# Tips : Ou est ajouté le call de la la fonciton cout ? il y a option_b en plus ... le crash intervient dans la fonction CLASSICAL_COMPUTE_JREG_B


# configure the Gamma model
model_gamma = gamma.smashplug.ConfigureGammaWithSmash(
    smash_model,
    dt=900,
    vmin=0.1,
    vmax=5.0,
    mode_discretization_step=0.1,
    spreading_discretization_step=0.1,
    ponderation_regul=0.0,
    velocity_computation="qm3",
    varying_spread=1,
    spreading_uniform=1,
    criteria="nse",
    ponderation_cost=10000.0,
    pdt_start_optim=1600,
)

# Set parameter
model_gamma.routing_parameters.hydraulics_coefficient = 0.5
model_gamma.routing_parameters.spreading = 1.0
model_gamma.routing_mesh.controlled_nodes[1:3] = 0

# direct run of the coupled model
gamma.smashplug.RunCoupledModel(smash_model, model_gamma)

# get observed discharges:
GammaGriddedObservation = gamma.smashplug.SmashDataVectorsToGrid(
    smash_model.response_data.q, model_gamma, smash_dt=smash_model.setup.dt
)

# compute the initial cost
initial_cost = model_gamma.cost_function(
    GammaGriddedObservation, model_gamma.routing_results.discharges
)

# Get the control vector
ControlVector = gamma.smashplug.VectorizeModelParameters(
    smash_model,
    model_gamma,
    control_parameters_list=["cp", "ct", "hydraulics_coefficient", "spreading"],
)

ControlVector.update(
    {
        "bounds": {
            "cp": [0.1, 1000.0],
            "ct": [0.1, 1000.0],
            "hydraulics_coefficient": [0.3, 5.0],
            "spreading": [0.5, 3.0],
        }
    }
)

# get the inflows
GammaInflows = gamma.smashplug.GetGammaInflowFromSmash(smash_model, dt=900)


cost, grad = gamma.smashplug.ComputeCostAndGradients(
    ControlVector["X"],
    ControlVector,
    smash_model,
    model_gamma,
    GammaGriddedObservation,
    True,
    True,
)

# Compute the gradient
gradient = gamma.smashplug.ComputeModelGradients(
    ControlVector, smash_model, model_gamma, GammaInflows, GammaGriddedObservation
)

# optimize parameters
# Here smash parameter are not normalized. That cause slow convergence. We try some gradient normalisation technics : ScaleGradientsByBounds and ScaleGammaGradients
# the spreading coefficient will uniformly calibrated => the gradient of this variable is scaled by the nb of nodes !
BestControlVector, optimized_smash_model, optimized_gamma_model = (
    gamma.smashplug.OptimizeCoupledModel(
        smash_model,
        model_gamma,
        GammaGriddedObservation,
        control_parameters_list=["cp", "ct", "hydraulics_coefficient", "spreading"],
        bounds={
            "cp": [0.1, 1000.0],
            "ct": [0.1, 1000.0],
            "hydraulics_coefficient": [0.3, 5.0],
            "spreading": [0.5, 3.0],
        },
        maxiter=15,
        tol=0.00001,
        ScaleGradientsByBounds=True,
        ScaleGammaGradients=True,
    )
)
