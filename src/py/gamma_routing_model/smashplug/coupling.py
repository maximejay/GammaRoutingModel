import smash
import gamma_routing_model as gamma
import numpy as np
import scipy

# ~ GammaRouting is a conceptual flow propagation model
# ~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

# ~ This file is part of GammaRouting.

# ~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# ~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# ~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

# ~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


def MatrixCoordsToVectorIndexes(row, col, ncol):
    # Convert Matrix coordinates (row,col) to a vector of indexes
    # if row, col are scalars, it return a scalar
    # if row, col are vector, it return a vector

    return row * ncol + col


def VectorIndexesOfActiveSmashCells(smash_model):
    # Calcul des indexes des pixels actifs du modèle Smash (grille/matrix) et stockage dans un vecteur dont la dimension correspond aux nombres de pixels actifs de Smash
    # Input: un objet modelè smash
    # output: un vecteur contenant les indexes des pixels actifs

    ncol = smash_model.mesh.ncol
    nrow = smash_model.mesh.nrow
    nb_active_cell = np.sum(smash_model.mesh.active_cell)

    index_vector_matrix = np.zeros(nb_active_cell)

    # Index des pixels actifs des matrices Smash stockés dans un vecteur
    k = 0
    for row in range(nrow):
        for col in range(ncol):

            if smash_model.mesh.active_cell[row, col] > 0:

                index_vector_matrix[k] = row * ncol + col
                k = k + 1

    return index_vector_matrix


def VectorIndexesOfSmashGauges(smash_model):
    # Compute the indexes of each smash gauges. Smash gauges coordinates are defined as (row,col) coordinates system
    # Return a vector of the indexes of each gauge computed as row*ncol+col

    ncol = smash_model.mesh.ncol
    nb_active_cell = np.sum(smash_model.mesh.active_cell)
    # index_gauge_pos=smash_model.mesh.gauge_pos[:,0]*(ncol)+smash_model.mesh.gauge_pos[:,1]

    index_gauge_pos = MatrixCoordsToVectorIndexes(
        smash_model.mesh.gauge_pos[:, 0], smash_model.mesh.gauge_pos[:, 1], ncol
    )
    index_vector_matrix = VectorIndexesOfActiveSmashCells(smash_model)

    gauge_pos = np.zeros(nb_active_cell) - 1

    for i in range(len(index_gauge_pos)):

        pos = np.where(index_gauge_pos[i] == index_vector_matrix)[0][0]
        gauge_pos[i] = pos

    return gauge_pos


def SmashMesh2DToVector(smash_model):
    # Convert some matrixes of the smash mesh into vectors suitable for the gamma meshing
    # It convert the following matrixes:
    # flwacc : accumulation flow
    # active_cell: the active cells
    # It also build a vector of the indexes of the cells of the smash matrixes: index_vector_matrix[k]=row*ncol+col
    # This function return a dictionary with these 3 vectors

    ncol = smash_model.mesh.ncol
    nrow = smash_model.mesh.nrow
    nb_active_cell = np.sum(smash_model.mesh.active_cell)

    index_vector_matrix = np.zeros(nb_active_cell)
    flw_acc_lin = np.zeros(nb_active_cell)
    active_cell_lin = np.zeros(nb_active_cell)
    dx = np.zeros(nb_active_cell)

    k = 0
    for row in range(nrow):

        for col in range(ncol):

            if smash_model.mesh.active_cell[row, col] > 0:

                index_vector_matrix[k] = row * ncol + col
                flw_acc_lin[k] = smash_model.mesh.flwacc[row, col]
                active_cell_lin[k] = smash_model.mesh.active_cell[row, col]
                dx[k] = smash_model.mesh.dx[row, col]
                k = k + 1

    VectorMesh1D = {
        "MatrixIndexes": index_vector_matrix,
        "FlowAcc": flw_acc_lin,
        "ActivesCells": active_cell_lin,
        "dx": dx,
    }

    return VectorMesh1D


def ComputeNodeLinker(smash_model):
    # This function compute the nodes linker path vector. It return a 2D vector with the corresponding upstream nodes for each node.
    # This mesh is build from a grid mesh where the water flow in 8 possibles directions
    # Return the node linker vector

    nb_upstream_nodes = 8  # number of possible flow directions
    dcol = [0, -1, -1, -1, 0, 1, 1, 1]  # moving directions in column
    drow = [1, 1, 0, -1, -1, -1, 0, 1]  # moving direction in rows
    dkind = [1, 2, 3, 4, 5, 6, 7, 8]  # coding directions

    ncol = smash_model.mesh.ncol
    nrow = smash_model.mesh.nrow
    nb_active_cell = np.sum(smash_model.mesh.active_cell)

    index_vector_matrix = VectorIndexesOfActiveSmashCells(smash_model)

    nodes_linker = np.zeros(shape=(nb_upstream_nodes, nb_active_cell))

    k = 0
    for row in range(nrow):

        for col in range(ncol):

            if smash_model.mesh.active_cell[row, col] > 0:

                l = 0
                for up_nd in range(nb_upstream_nodes):

                    row_imd = int(row + drow[up_nd])
                    col_imd = int(col + dcol[up_nd])

                    if (
                        (col_imd >= 0)
                        and (col_imd < ncol)
                        and (row_imd >= 0)
                        and (row_imd < nrow)
                    ):

                        if smash_model.mesh.flwdir[row_imd, col_imd] == dkind[up_nd]:

                            pos_grid = row_imd * (ncol) + col_imd
                            pos_vector = np.where(pos_grid == index_vector_matrix)
                            nodes_linker[l, k] = pos_vector[0][0] + 1
                            l = l + 1

                k = k + 1

    return nodes_linker


def SmashMeshToGammaMesh(smash_model, model_gamma):
    # Build the Gamma mesh from the one of the Smash model.
    # The smash model use a regular grid for the mesh. It is build with python/fortran matrixes
    # Gamma use a custom mesh defined with nodes and reach. It is build with a python/fortran vectors
    # We can convert regular matrixes mesh to vector with this function. We get nodes and reaches with a uniform spatial discretization step.
    # Input:
    # smash_model: the smash model object
    # model_gamma: the  gamma object model to be updated with a new mesh

    nb_nodes = np.sum(smash_model.mesh.active_cell)
    nb_upstream_nodes = 8  # properties of gridded mesh
    model_gamma.routing_mesh_init(
        nb_nodes=nb_nodes, nb_upstream_nodes=nb_upstream_nodes
    )  # initialise the gamma mesh

    VectorMesh1D = SmashMesh2DToVector(smash_model)
    nodes_linker = ComputeNodeLinker(smash_model)
    gauge_pos = VectorIndexesOfSmashGauges(smash_model)

    # Filling the Gamma mesh properties
    model_gamma.routing_mesh.upstream_to_downstream_nodes = (
        np.argsort(VectorMesh1D["FlowAcc"]) + 1
    )
    model_gamma.routing_mesh.nodes_linker = nodes_linker
    model_gamma.routing_mesh.controlled_nodes = gauge_pos + 1
    model_gamma.routing_mesh.surface = VectorMesh1D["dx"] ** 2.0 / 1000.0**2.0  #!km²
    model_gamma.routing_mesh.dx = VectorMesh1D["dx"]

    # Update and finalise the gamma mesh for the model_gamma object
    model_gamma.routing_mesh_update()


def ConfigureGammaWithSmash(smash_model, dt=None, **kwargs):
    # Configure Gamma from the smash model structure
    # Input:
    # smash_model: smash model object
    # dt: time step of the gamma model
    # **kwargs: any arguments for gamma.routing_setup

    if smash_model.setup.routing_module != "rm_zero":
        print("Warnings: setting smash_model.setup.routing_module to rm_zero")
        smash_model.setup.structure = "rm_zero"

    # initialise the gamma model object
    model_gamma = gamma.Model()

    if "dt" in kwargs:
        kwargs.pop("dt")

    if "npdt" in kwargs:
        kwargs.pop("npdt")
        print("Warnings: npdt cannot set by the user")

    if dt is None:
        dt = smash_model.setup.dt
    else:
        if dt > smash_model.setup.dt:
            print("Warning: Gamma dt cannot be greater than Smash dt")
            dt = smash_model.setup.dt

    npdt = int(smash_model.setup.ntime_step * smash_model.setup.dt / dt)

    # Initialise the gamma model setup
    model_gamma.routing_setup_init(npdt=npdt, dt=dt, **kwargs)

    # Build the gamma mesh from the smash mesh
    SmashMeshToGammaMesh(smash_model, model_gamma)

    # initialise the parameter of the gamma model
    model_gamma.routing_parameters_init()

    # intialise the states of the gamma model
    model_gamma.routing_states_init()

    # return the complete gamma model object, ready for run
    return model_gamma


def GetGammaInflowFromSmash(smash_model, dt=None):
    # Get the current inflow from Smash outflow
    # Interpolate if dt gamma is lower
    # return the the matrix 'interpolated_inflows' ready to run the gamma model with gamma.run(interpolated_inflows)

    if dt is None:
        dt = smash_model.setup.dt

    npdt = int(smash_model.setup.ntime_step * smash_model.setup.dt / dt)
    nb_nodes = np.sum(smash_model.mesh.active_cell)
    ncol = smash_model.mesh.ncol
    nrow = smash_model.mesh.nrow

    # interface inflows
    # inflows_smash = smash_model.output.qsim_domain.copy()
    inflows_smash = smash_model.response.qt.copy()
    # inflows_domain = np.zeros(shape=(smash_model.setup.ntime_step, nb_nodes))
    inflows_domain = inflows_smash.transpose()

    # interpolation si pdt gamma inférieur
    if dt < smash_model.setup.dt:

        frac_dt = int(smash_model.setup.dt / dt)
        interpolated_inflows = np.zeros(shape=(npdt, inflows_domain.shape[1]))

        # Le pdt entre 0 à 1 pour gamma correspond au premier pdt de smash
        interpolated_inflows[0:frac_dt, :] = inflows_domain[0, :]

        itg = frac_dt
        for ts in range(smash_model.setup.ntime_step - 1):

            tg = ts + 1.0 / frac_dt

            for itfrac in range(frac_dt):

                interpolated_inflows[itg, :] = (
                    inflows_domain[ts + 1, :] - inflows_domain[ts, :]
                ) * (tg - ts) / (ts + 1 - ts) + inflows_domain[ts, :]

                itg = itg + 1
                tg = tg + 1.0 / frac_dt

    else:
        interpolated_inflows = inflows_domain[:, :].copy()

    return interpolated_inflows


def InterpolatedObservationsAtGauge(smash_model, dt):
    # Interpolate the observations of the smash model at every gauge for a new time-step

    npdt = int(smash_model.setup._ntime_step * smash_model.setup.dt / dt)
    # interface observations
    observations = smash_model.input_data.qobs.copy()

    # interpolation si pdt gamma inférieur
    if dt < smash_model.setup.dt:

        frac_dt = int(smash_model.setup.dt / dt)
        interpolated_observations = np.zeros(shape=(npdt, observations.shape[0]))
        itg = frac_dt

        # cdt initial
        interpolated_observations[0:itg, :] = observations[:, 0]

        for ts in range(smash_model.setup._ntime_step - 1):

            tg = ts + 1.0 / frac_dt

            for itfrac in range(frac_dt):

                interpolated_observations[itg, :] = np.transpose(
                    (observations[:, ts + 1] - observations[:, ts])
                    * (tg - ts)
                    / (ts + 1 - ts)
                    + observations[:, ts]
                )

                itg = itg + 1
                tg = tg + 1.0 / frac_dt

    else:
        interpolated_observations = observations[:, :].copy()

    return interpolated_observations


def FitDataToNewDt(data, dt, new_dt):
    # interpolate a data vector to a new time step

    if new_dt < dt:

        npdt = int(data.shape[0] * dt / new_dt)
        new_data = np.zeros(npdt)
        frac_dt = int(dt / new_dt)

        i_nt = frac_dt

        # cdt initial
        new_data[0:i_nt] = data[0]

        for ts in range(data.shape[0] - 1):

            i_t = ts + 1.0 / frac_dt

            for itfrac in range(frac_dt):

                new_data[i_nt] = (data[ts + 1] - data[ts]) * (i_t - ts) / (
                    ts + 1 - ts
                ) + data[ts]

                i_nt = i_nt + 1
                i_t = i_t + 1.0 / frac_dt

    elif new_dt > dt:

        npdt = int(data.shape[0] * dt / new_dt)
        new_data = np.zeros(npdt)
        frac_dt = int(new_dt / dt)

        n_ts = 0

        for ts in range(npdt):

            new_data[ts] = np.mean(data[n_ts : n_ts + frac_dt])
            n_ts = n_ts + frac_dt

    else:

        new_data = data.copy()

    return new_data


def GammaVectorsToSmashGrid(vector, smash_model, model_gamma):
    # Convert vector from gamma model to smash grid
    # vector[nb_nodes]

    grid = np.zeros(shape=(smash_model.mesh.nrow, smash_model.mesh.ncol))

    k = 0
    for row in range(smash_model.mesh.nrow):

        for col in range(smash_model.mesh.ncol):

            if smash_model.mesh.active_cell[row, col] > 0:

                grid[row, col] = vector[k]
                k = k + 1

    return grid


def SmashGridToGammaVectors(grid, smash_model, model_gamma):
    # Convert smash gridded data to vectors for gamma model
    # grid[nbx,nby]

    nb_gauge = len(np.where(model_gamma.routing_mesh.controlled_nodes > 0)[0])

    vector = np.zeros(shape=(model_gamma.routing_mesh.nb_nodes))

    k = 0
    for row in range(smash_model.mesh.nrow):

        for col in range(smash_model.mesh.ncol):

            if smash_model.mesh.active_cell[row, col] > 0:

                vector[k] = grid[row, col]
                k = k + 1

    return vector


def SmashDataVectorsToGrid(observations, model_gamma, smash_dt):
    # observations[nb_gauge,npdt]
    # gridded_observation[npdt,nb_cells]

    gridded_observation = np.zeros(
        shape=(model_gamma.routing_setup.npdt, model_gamma.routing_mesh.nb_nodes)
    )

    pos_observation = 0
    for i in range(model_gamma.routing_mesh.nb_nodes):

        k = model_gamma.routing_mesh.controlled_nodes[i] - 1

        if (k > 0) and (k <= model_gamma.routing_mesh.nb_nodes):
            # print(f"fitting at obs {pos_observation} for grid {k} from {smash_dt} to {model_gamma.routing_setup.dt}")
            gridded_observation[:, k] = FitDataToNewDt(
                observations[pos_observation, :], smash_dt, model_gamma.routing_setup.dt
            )
            pos_observation = pos_observation + 1

    return gridded_observation


def GriddedDataToSmashVectors(gridded_observation, model_gamma, smash_dt):
    # gridded_observation[npdt,nb_cells]
    # observations[nb_gauge,npdt]

    nb_gauge = len(np.where(model_gamma.routing_mesh.controlled_nodes > 0)[0])

    observations = np.zeros(
        shape=(
            nb_gauge,
            int(gridded_observation.shape[0] * model_gamma.routing_setup.dt / smash_dt),
        )
    )

    pos_observation = 0
    for i in range(model_gamma.routing_mesh.nb_nodes):

        k = model_gamma.routing_mesh.controlled_nodes[i] - 1

        if (k > 0) and (k <= model_gamma.routing_mesh.nb_nodes):

            observations[pos_observation, :] = FitDataToNewDt(
                gridded_observation[:, k], model_gamma.routing_setup.dt, smash_dt
            )
            pos_observation = pos_observation + 1

    return observations


def GetCorrespondingTimeStepRange(smash_model, model_gamma):
    # Get a vector of time-steps of gamma and smash (compare to smash)
    # This is useful for quick ploting
    # Time-Step of Smash is counting from 1 to N
    # Time-Step of Gamma is counting from 1/dt_smash to N+1/dt_smash

    frac_dt = int(smash_model.setup.dt / model_gamma.routing_setup.dt)

    xsmash = np.arange(1, smash_model.setup._ntime_step + 1)
    xgamma = np.arange(
        1 / frac_dt, smash_model.setup._ntime_step + 1 / frac_dt, 1 / frac_dt
    )

    return {"ts_smash": xsmash, "ts_gamma": xgamma}


def VectorizeModelParameters(smash_model, model_gamma, control_parameters_list=list()):
    # Get all controlled parameters from Smash and Gamma and vectorize them
    # Return a the control vector dictionary with component:
    #'X' : the control vector
    #'ParamList' : the list of controled parameter
    #'NbNodes' : the number of active cells/nodes for each controlloed parameter
    # X has the following order:
    # - Gradient of Smash in the order of the control_vector
    # - Gradient of Gamma : 1st hydraulics_coefficient, 2nd spreading

    if len(control_parameters_list) == 0:

        control_parameters_list = list()
        control_parameters_list.append(
            smash._constant.STRUCTURE_RR_PARAMETERS[smash_model.setup.structure]
        )

        control_parameters_list.append("hydraulics_coefficient")
        if model_gamma.setup.varying_spread == True:
            control_parameters_list.append("spreading")

    else:

        for ctrl_var in control_parameters_list:

            if not (
                ctrl_var in smash_model.rr_parameters.keys
                or ctrl_var in smash_model.rr_initial_states.keys
                or (hasattr(model_gamma.routing_parameters, ctrl_var))
            ):

                print("Error: Wrong parameter/states in control vector: " + ctrl_var)
                return

    SizeOfControl = len(control_parameters_list)
    LinearizedParameters = np.zeros(SizeOfControl * model_gamma.routing_mesh.nb_nodes)

    k = 0  # indexe of the current cell

    for ctrl_var in control_parameters_list:

        if ctrl_var in smash_model.rr_parameters.keys:

            index_param = list(smash_model.rr_parameters.keys).index(ctrl_var)
            param_smash = smash_model.rr_parameters.values[:, :, index_param]

            for row in range(smash_model.mesh.nrow):

                for col in range(smash_model.mesh.ncol):

                    if smash_model.mesh.active_cell[row, col] > 0:

                        LinearizedParameters[k] = param_smash[row, col]
                        k = k + 1

        if ctrl_var in smash_model.rr_initial_states.keys:

            index_state = list(smash_model.rr_initial_states).keys.index(ctrl_var)
            states_smash = smash_model.rr_initial_states.values[:, :, index_state]

            for row in range(smash_model.mesh.nrow):

                for col in range(smash_model.mesh.ncol):

                    if smash_model.mesh.active_cell[row, col] > 0:

                        LinearizedParameters[k] = states_smash[row, col]
                        k = k + 1

    nb_nodes = model_gamma.routing_mesh.nb_nodes

    for ctrl_var in control_parameters_list:

        if hasattr(model_gamma.routing_parameters, ctrl_var):

            if ctrl_var == "hydraulics_coefficient":

                LinearizedParameters[k : k + nb_nodes] = (
                    model_gamma.routing_parameters.hydraulics_coefficient[:]
                )
                k = k + nb_nodes

            if ctrl_var == "spreading":

                LinearizedParameters[k : k + nb_nodes] = (
                    model_gamma.routing_parameters.spreading[:]
                )
                k = k + nb_nodes

    control_vector = {
        "X": LinearizedParameters,
        "ParamList": control_parameters_list,
        "NbNodes": nb_nodes,
    }

    return control_vector


def SetVectorizedModelParameters(control_vector, smash_model, model_gamma):
    # Update the smash model object and the gamma model object from the control_vector dictionary
    # control_vector['X'] has the following order:
    # - Gradient of Smash in the order of the control_vector["ParamList"]
    # - Gradient of Gamma : 1st hydraulics_coefficient, 2nd spreading

    if not (
        ("X" in control_vector)
        and ("ParamList" in control_vector)
        and ("NbNodes" in control_vector)
    ):
        print(
            "Error: the control vector is an incomplete dictionary. You must provide a control vector build with the function VectorizeModelParameters"
        )
        raise ValueError(
            "Error: the control vector is an incomplete dictionary. You must provide a control vector build with the function VectorizeModelParameters"
        )

    k = 0

    # corresponding smash function : smash/core/simulation/_optimize/_x_to_parameters_states
    for ctrl_var in control_vector["ParamList"]:

        if hasattr(smash_model.parameters, ctrl_var):

            MatrixParameters = getattr(smash_model.parameters, ctrl_var)

            for sub_row in range(smash_model.mesh.nrow):

                for sub_col in range(smash_model.mesh.ncol):

                    if smash_model.mesh.active_cell[sub_row, sub_col] > 0:

                        MatrixParameters[sub_row, sub_col] = control_vector["X"][k]

                        k = k + 1

            setattr(smash_model.parameters, ctrl_var, MatrixParameters)

        if hasattr(smash_model.states, ctrl_var):

            MatrixStates = getattr(smash_model.states, ctrl_var)

            for sub_row in range(smash_model.mesh.nrow):

                for sub_col in range(smash_model.mesh.ncol):

                    if smash_model.mesh.active_cell[sub_row, sub_col] > 0:

                        MatrixStates[sub_row, sub_col] = control_vector["X"][k]

                        k = k + 1

            setattr(smash_model.states, ctrl_var, MatrixStates)

    # Upgrade for Gamma only
    nb_nodes = control_vector["NbNodes"]

    for ctrl_var in control_vector["ParamList"]:

        if hasattr(model_gamma.routing_parameters, ctrl_var):

            if ctrl_var == "hydraulics_coefficient":

                if model_gamma.routing_setup.hydraulics_coef_uniform == 1:
                    model_gamma.routing_parameters.hydraulics_coefficient[:] = (
                        control_vector["X"][k]
                    )
                else:
                    model_gamma.routing_parameters.hydraulics_coefficient[:] = (
                        control_vector["X"][k : k + nb_nodes]
                    )

                k = k + nb_nodes

            if ctrl_var == "spreading":

                if model_gamma.routing_setup.spreading_uniform == 1:
                    model_gamma.routing_parameters.spreading[:] = control_vector["X"][k]
                else:
                    model_gamma.routing_parameters.spreading[:] = control_vector["X"][
                        k : k + nb_nodes
                    ]

                k = k + nb_nodes


def ComputeModelGradients(
    control_vector,
    smash_model,
    model_gamma,
    interpolated_inflows,
    observations,
    set_param=False,
    ScaleGammaGradients=True,
    ScaleGradientsByBounds=True,
):
    # Compute the Gradients dCost/dX for model M(X)=M(Xgamma,Xsmash)
    # That computation required:
    # 1)dCost/dXgamma
    # 2)dCost/dYsmash
    # 3)dYsmash/dXsmash
    # This function return a single vector with all aggregated gradient in the following order:
    # - Gradient of Smash in the order of the control_vector
    # - Gradient of Gamma : 1st hydraulics_coefficient, 2nd spreading

    # Update model with the control vector if required
    # ~ if set_param:
    # ~ SetVectorizedModelParameters(control_vector,smash_model,model_gamma)

    # Gamma Side gradients computation :

    # Here we should normalize parameters
    # (p-lb) / (ub-lb)
    # set gamma_model.routing_setup.denormalize_parameter=1
    # implement in forward this condition to denormalize parameter

    # Gradients of COST with respect to routing_parameters
    Grad_dCOST_dROUTINGPARAMETERS = np.zeros(
        shape=(2, model_gamma.routing_mesh.nb_nodes), order="F", dtype="float32"
    )

    print(
        "call gamma.routing_gamma_forward_adjoint_b: Gradients of COST with respect to varying inputs ROUTING_PARAMETERS"
    )
    gamma.libfgamma.Mod_Gamma_Interface.routing_gamma_forward_adjoint_b(
        model_gamma.routing_setup,
        model_gamma.routing_mesh,
        model_gamma.routing_parameters,
        interpolated_inflows,
        observations,
        model_gamma.routing_states,
        model_gamma.routing_results,
        Grad_dCOST_dROUTINGPARAMETERS,
    )

    # Gradients of COST with respect to inflows (interpolated_inflows)
    Grad_dCOST_dINFLOWS = np.zeros(
        shape=(model_gamma.routing_setup.npdt, model_gamma.routing_mesh.nb_nodes),
        order="F",
        dtype="float32",
    )

    print(
        "call gamma.routing_gamma_forward_adjoint_b0: Gradients of COST with respect to varying input INFLOWS (interpolated_inflows)"
    )
    gamma.libfgamma.Mod_Gamma_Interface.routing_gamma_forward_adjoint_b0(
        model_gamma.routing_setup,
        model_gamma.routing_mesh,
        model_gamma.routing_parameters,
        interpolated_inflows,
        observations,
        model_gamma.routing_states,
        model_gamma.routing_results,
        Grad_dCOST_dINFLOWS,
    )

    # Intégratioon over the time for Gradients of COST with respect to inflows (interpolated_inflows): Cost and Inflows have not the same dimension. Inflows vary in time but not Cost as it is alrady integrated. Thus Tapenade compte the gradients for all time step. We need to integrate it manually.
    Int_Grad_dCOST_dINFLOWS = np.zeros(Grad_dCOST_dINFLOWS.shape[1])
    for i in range(Grad_dCOST_dINFLOWS.shape[1]):
        Int_Grad_dCOST_dINFLOWS[i] = np.mean(Grad_dCOST_dINFLOWS[:, i])

    # SMASH side gradients computation:

    # Gradients of OUTPUT.QSIM_DOMAIN with respect to PARAMETERS and STATES

    # Here we should normalize parameters
    # (p-lb) / (ub-lb)
    # set smash_model.setup._optimize.denormalize_parameter=1
    # set smash_model.setup._optimize.lb_parameters=
    # set smash_model.setup._optimize.ub_parameters=
    # set smash_model.setup._optimize.optim_parameters=
    # set for states also

    # initialisation des types dérivés fortran
    Grad_dQSIMDOMAIN_dPARAMETERS = smash.fcore._mwd_parameters.ParametersDT(
        smash_model.setup, smash_model.mesh
    )
    Grad_dQSIMDOMAIN_dOUPUT = smash.fcore._mwd_output.OutputDT(
        smash_model.setup, smash_model.mesh
    )
    Grad_dQSIMDOMAIN_dOUPUT.response.qt = np.float32(1)

    # Options = smash.default_optimize_options(smash_model, mapping="distributed")
    cost_options = (
        smash.core.simulation._standardize._standardize_simulation_cost_options(
            smash_model, "optimize", None
        )
    )
    smash.core.simulation._standardize._standardize_simulation_cost_options_finalize(
        smash_model, "optimize", cost_options
    )

    return_options = (
        smash.core.simulation._standardize._standardize_simulation_return_options(
            smash_model, "optimize", None
        )
    )
    smash.core.simulation._standardize._standardize_simulation_return_options_finalize(
        smash_model, return_options
    )

    print(
        "call smash.solver._mw_forward.wrapped_forward_b0: Gradients of OUTPUT.QSIM_DOMAIN with respect to varying inputs PARAMETERS and STATES"
    )
    smash.fcore._mw_forward.forward_run_b0(
        smash_model.setup,
        smash_model.mesh,
        smash_model._input_data,
        smash_model._parameters,
        Grad_dQSIMDOMAIN_dPARAMETERS,
        smash_model._output,
        Grad_dQSIMDOMAIN_dOUPUT,
        smash.fcore._mwd_options.OptionsDT(
            smash_model.setup,
            smash_model.mesh,
            cost_options["njoc"],
            cost_options["njrc"],
        ),
        smash.fcore._mwd_returns.ReturnsDT(
            smash_model.setup,
            smash_model.mesh,
            return_options["nmts"],
            return_options["fkeys"],
        ),
    )

    # Smash gradients dSimDomain/dParameter are already integrated in time.But we need to vectorize these gradient and keep only those on the active cells. MoreOver we need to compute the full gradient dCost/dParameters (dCost/dXsmash) given by the product : dQSIMDOMAIN/dPARAMETERS x dCost/dInflows
    # linear_grandient_cp=np.zeros(Int_Grad_dCOST_dINFLOWS.shape)

    # ~ grad_param_smash_cp=LinearizedGradientsForSmash[:]*Int_Grad_dCOST_dINFLOWS[:]

    list_param_smash = control_vector["ParamList"].copy()
    if "hydraulics_coefficient" in list_param_smash:
        list_param_smash.remove("hydraulics_coefficient")

    if "spreading" in list_param_smash:
        list_param_smash.remove("spreading")

    nb_param_smash = len(list_param_smash)
    LinearizedGradientsForSmash = np.zeros(
        model_gamma.routing_mesh.nb_nodes * nb_param_smash
    )

    # print(control_vector['ParamList'],list_param_smash,len(LinearizedGradientsForSmash))

    k = 0
    for ctrl_var in control_vector["ParamList"]:

        # check if parameters type
        if hasattr(smash_model.parameters, ctrl_var):

            MatrixGradients = getattr(Grad_dQSIMDOMAIN_dPARAMETERS, ctrl_var)
            i = 0

            for sub_row in range(smash_model.mesh.nrow):

                for sub_col in range(smash_model.mesh.ncol):

                    if smash_model.mesh.active_cell[sub_row, sub_col] > 0:

                        LinearizedGradientsForSmash[k] = (
                            MatrixGradients[sub_row, sub_col] * Int_Grad_dCOST_dINFLOWS[i]
                        )
                        k = k + 1
                        i = i + 1

        # check if states type
        # if hasattr(smash_model.states, ctrl_var):

        #     MatrixGradients = getattr(Grad_dQSIMDOMAIN_dSTATES, ctrl_var)
        #     i = 0

        #     for sub_row in range(smash_model.mesh.nrow):

        #         for sub_col in range(smash_model.mesh.ncol):

        #             if smash_model.mesh.active_cell[sub_row, sub_col] > 0:

        #                 LinearizedGradientsForSmash[k] = (
        #                     MatrixGradients[sub_row, sub_col] * Int_Grad_dCOST_dINFLOWS[i]
        #                 )
        #                 k = k + 1
        #                 i = i + 1

    # At this step we own all gradients given by Tapenade.
    # Here we must normalise the gradients, so that the gradients of Gamma and Smash have the same magnitude (same physical meaning)
    # The Smash gradients are given for 1 cell and 1 input, thus for 1 km²
    # The Gamma gradients are given along the network with increasing surface from upstream to downstream.
    # We choose to normalise the Gamma gradients by the cumulated surface.
    if ScaleGammaGradients and model_gamma.routing_setup.hydraulics_coef_uniform == 0:
        Grad_dCOST_dROUTINGPARAMETERS[0, :] = (
            Grad_dCOST_dROUTINGPARAMETERS[0, :]
            / model_gamma.routing_mesh.cumulated_surface[:]
        )

    if ScaleGammaGradients and model_gamma.routing_setup.spreading_uniform == 0:
        Grad_dCOST_dROUTINGPARAMETERS[1, :] = (
            Grad_dCOST_dROUTINGPARAMETERS[1, :]
            / model_gamma.routing_mesh.cumulated_surface[:]
        )

    # Here is an attempt to normalize by the mean discharge
    # ~ mean_gamma_discharges=np.zeros(model_gamma.routing_mesh.nb_nodes)
    # ~ for i in range(model_gamma.routing_mesh.nb_nodes):
    # ~ mean_gamma_discharges[i]=np.mean(model_gamma.routing_results.discharges[:,i])
    # normalisation par mean discharges
    # ~ Grad_dCOST_dROUTINGPARAMETERS[0,:]=Grad_dCOST_dROUTINGPARAMETERS[0,:]/mean_gamma_discharges[:]
    # ~ Grad_dCOST_dROUTINGPARAMETERS[1,:]=Grad_dCOST_dROUTINGPARAMETERS[1,:]/mean_gamma_discharges[:]

    # Parameters/states are not normalized... We should normalize the gradient ?

    # Finally we aggregate all computed gradients in a single vector given the control_vector
    NbControledGammaParameter = 0
    if "hydraulics_coefficient" in control_vector["ParamList"]:
        NbControledGammaParameter = NbControledGammaParameter + 1

    if "spreading" in control_vector["ParamList"]:
        NbControledGammaParameter = NbControledGammaParameter + 1

    OutputsGradients = np.zeros(
        len(LinearizedGradientsForSmash)
        + NbControledGammaParameter * control_vector["NbNodes"]
    )

    OutputsGradients[0 : len(LinearizedGradientsForSmash)] = LinearizedGradientsForSmash

    position = len(LinearizedGradientsForSmash)

    for ctrl_var in control_vector["ParamList"]:

        if hasattr(model_gamma.routing_parameters, ctrl_var):

            if ctrl_var == "hydraulics_coefficient":

                if model_gamma.routing_setup.hydraulics_coef_uniform == 1:
                    OutputsGradients[position : position + control_vector["NbNodes"]] = (
                        Grad_dCOST_dROUTINGPARAMETERS[0, :]
                        / float(control_vector["NbNodes"])
                    )
                else:
                    OutputsGradients[position : position + control_vector["NbNodes"]] = (
                        Grad_dCOST_dROUTINGPARAMETERS[0, :]
                    )

                position = position + control_vector["NbNodes"]

            if ctrl_var == "spreading":

                if model_gamma.routing_setup.spreading_uniform == 1:
                    OutputsGradients[position : position + control_vector["NbNodes"]] = (
                        Grad_dCOST_dROUTINGPARAMETERS[1, :]
                        / float(control_vector["NbNodes"])
                    )
                else:
                    OutputsGradients[position : position + control_vector["NbNodes"]] = (
                        Grad_dCOST_dROUTINGPARAMETERS[1, :]
                    )

                position = position + control_vector["NbNodes"]

    # Normalize the output gradients with respect to the boundaries
    if ("bounds" in control_vector) and (ScaleGradientsByBounds == True):

        position = 0
        max_amplitude = 0.0
        # sum_amplitude=0.
        for param in control_vector["ParamList"]:

            amplitude = (
                control_vector["bounds"][param][1] - control_vector["bounds"][param][0]
            )

            if amplitude > max_amplitude:
                max_amplitude = amplitude
                # sum_amplitude=sum_amplitude+amplitude

            if amplitude > 0:
                OutputsGradients[position : position + control_vector["NbNodes"]] = (
                    OutputsGradients[position : position + control_vector["NbNodes"]]
                    / amplitude
                )

            position = position + control_vector["NbNodes"]

        OutputsGradients = OutputsGradients * max_amplitude
        # OutputsGradients=OutputsGradients* sum_amplitude / len(control_vector["ParamList"])

    return OutputsGradients


def RunCoupledModel(smash_model, model_gamma):
    # Run the direct problem of the coupled model Smash and Gamma

    print("Run of the Smash model")
    smash_model.forward_run()

    print(
        f"Getting the interpolated Smash outflow at time-step {model_gamma.routing_setup.dt}"
    )
    interpolated_inflows = GetGammaInflowFromSmash(
        smash_model, dt=model_gamma.routing_setup.dt
    )

    # run the model
    print("Run the Gamma Routing model")
    model_gamma.run(interpolated_inflows, states_init=0)


def ComputeCostAndGradients(
    X,
    control_vector,
    smash_model,
    model_gamma,
    observations,
    ScaleGammaGradients=True,
    ScaleGradientsByBounds=True,
):
    # Compute the gamma cost function and the gradients of dCost/(dXsmash,dXgamma)

    # Get the new control vector
    print("Updating the control vector")
    control_vector["X"] = X
    # Set the new control vector to the Smash and Gamma Model
    SetVectorizedModelParameters(control_vector, smash_model, model_gamma)

    print("smash_model_run")
    smash_model.run(inplace=True)

    print("interolate inflows")
    interpolated_inflows = GetGammaInflowFromSmash(
        smash_model, dt=model_gamma.routing_setup.dt
    )

    # run the model
    print("model_gamma_run")
    model_gamma.run(interpolated_inflows, states_init=0)

    print("compute cost")
    model_gamma.cost_function(observations, model_gamma.routing_results.discharges)
    cost = model_gamma.routing_results.costs[0]

    print("compute gradients")
    gradient = ComputeModelGradients(
        control_vector,
        smash_model,
        model_gamma,
        interpolated_inflows,
        observations,
        ScaleGammaGradients=ScaleGammaGradients,
        ScaleGradientsByBounds=ScaleGradientsByBounds,
    )

    return cost, gradient


def OptimizeCoupledModel(
    smash_model,
    model_gamma,
    observations,
    control_parameters_list=["cp", "hydraulics_coefficient", "spreading"],
    bounds={
        "cp": [0.1, 1000.0],
        "hydraulics_coefficient": [0.5, 10],
        "spreading": [1.0, 3.0],
    },
    maxiter=10,
    tol=0.00001,
    ScaleGammaGradients=True,
    ScaleGradientsByBounds=True,
    inplace=False,
):

    # Optimize the distributed parameters of the coupled model using the lbfgsb controler

    optimized_gamma = model_gamma.copy()
    optimized_smash = smash_model.copy()

    if observations.shape != optimized_gamma.routing_results.discharges.shape:
        raise ValueError(
            f"Observations vectors has not the same shape ({observations.shape}) than the outputs discharges of Gamma ({optimized_gamma.routing_results.discharges.shape})"
        )

    # define the control vector
    ControlVector = gamma.smashplug.VectorizeModelParameters(
        optimized_smash, optimized_gamma, control_parameters_list=control_parameters_list
    )

    boundaries = np.zeros(shape=(len(ControlVector["X"]), 2))

    # define the boundaris
    if len(control_parameters_list) != len(bounds):
        raise ValueError(
            "Error: The number of boundaries is different than the number of controls"
        )

    position = 0
    for control in control_parameters_list:

        # set only smash boundaries first
        if hasattr(optimized_smash.parameters, control):

            if not control in bounds:
                raise ValueError(
                    f"Error: The controlled variable {control} is not in the boundaries definition"
                )

            boundaries[position : position + ControlVector["NbNodes"], 0] = bounds[
                control
            ][0]
            boundaries[position : position + ControlVector["NbNodes"], 1] = bounds[
                control
            ][1]
            position = position + ControlVector["NbNodes"]

        if hasattr(optimized_smash.states, control):

            if not control in bounds:
                raise ValueError(
                    f"Error: The controlled variable {control} is not in the boundaries definition"
                )

            boundaries[position : position + ControlVector["NbNodes"], 0] = bounds[
                control
            ][0]
            boundaries[position : position + ControlVector["NbNodes"], 1] = bounds[
                control
            ][1]
            position = position + ControlVector["NbNodes"]

    for control in control_parameters_list:

        # set only gamma boundaries
        if hasattr(optimized_gamma.routing_parameters, control):

            if control == "hydraulics_coefficient":

                if not "hydraulics_coefficient" in bounds:
                    raise ValueError(
                        f"Error: The controlled variable hydraulics_coefficient is not in the boundaries definition"
                    )

                boundaries[position : position + ControlVector["NbNodes"], 0] = bounds[
                    "hydraulics_coefficient"
                ][0]
                boundaries[position : position + ControlVector["NbNodes"], 1] = bounds[
                    "hydraulics_coefficient"
                ][1]
                position = position + ControlVector["NbNodes"]

            if control == "spreading":

                if not "spreading" in bounds:
                    raise ValueError(
                        f"Error: The controlled variable spreading is not in the boundaries definition"
                    )

                boundaries[position : position + ControlVector["NbNodes"], 0] = bounds[
                    "spreading"
                ][0]
                boundaries[position : position + ControlVector["NbNodes"], 1] = bounds[
                    "spreading"
                ][1]
                position = position + ControlVector["NbNodes"]

    # update the control vector with bounds
    if ScaleGradientsByBounds:
        ControlVector.update({"bounds": bounds})

    res = scipy.optimize.minimize(
        gamma.smashplug.ComputeCostAndGradients,
        ControlVector["X"],
        args=(
            ControlVector,
            optimized_smash,
            optimized_gamma,
            observations,
            ScaleGammaGradients,
            ScaleGradientsByBounds,
        ),
        method="L-BFGS-B",
        jac=True,
        bounds=boundaries,
        tol=tol,
        options={"disp": True, "iprint": 101, "maxiter": maxiter},
    )

    # get results
    new_parameters = res.x
    # set the optimal parameter to the Smash and Gamma Model
    ControlVector["X"] = new_parameters
    gamma.smashplug.SetVectorizedModelParameters(
        ControlVector, optimized_smash, optimized_gamma
    )

    # Run the direct coupled model Gamma o Smash-
    gamma.smashplug.RunCoupledModel(optimized_smash, optimized_gamma)

    # convert Vector discharges from Smash to Gamma gridded discharges
    # ~ GammaGriddedObservation=gamma.smashplug.SmashDataVectorsToGrid(optimized_smash.input_data.qobs,optimized_gamma,smash_dt=optimized_smash.setup.dt)

    # compute the final cost
    # ~ final_cost=model_gamma.cost_function(GammaGriddedObservation,optimized_gamma.routing_results.discharges)

    if inplace == True:
        model_gamma = optimized_gamma.copy()
        smash_model = optimized_smash.copy()
        return ControlVector

    if inplace == False:
        return (ControlVector, optimized_smash, optimized_gamma)
