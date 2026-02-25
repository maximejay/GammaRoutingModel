# ~ GammaRouting is a conceptual flow propagation model
# ~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

# ~ This file is part of GammaRouting.

# ~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# ~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# ~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

# ~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


from __future__ import annotations

import traceback
import os
import numpy as np

from gamma_routing_model.libfgamma import (
    Mod_Gamma_Routing_Setup,
    Mod_Gamma_Routing_Mesh,
    Mod_Gamma_Routing_Parameters,
    Mod_Gamma_Routing_States,
    Mod_Gamma_Routing_Memory,
    Mod_Gamma_Routing_Results,
    Mod_Gamma_Interface,
)


__all__ = ["Model"]


class Model(object):
    """
    Primary data structure of the hydrological model `GammaRouting`.

    Additional help can be found here for configuration file options : (user guide)

    Parameters
    ----------
    **kwargs:
            Arbitrary keyword arguments where keys must be attributes of the
            derived type `model.routing_setup`.

    Examples
    --------
    >>> import gamma
    >>> model = model=gamma.Model()

    See Also
    --------
    gamma.routing_setup : Fortan wrapped derived type.
    """

    def __init__(self, **kwargs):

        print("")
        print(
            "GammaRouting  Copyright (C) 2022,2023  Hydris-hydrologie, Maxime Jay-Allemand (contact: maxime.jay.allemand@hydris-hydrologie.fr)"
        )
        print(
            "This program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions."
        )
        print("")

        # __init__
        # self.kwargs = kwargs

        # check all arg and set default value if not set
        # then pass it to Type_Routing_Setup

        # ~ try:
        self._routing_setup = Mod_Gamma_Routing_Setup.type_routing_setup()
        self._routing_mesh = Mod_Gamma_Routing_Mesh.type_routing_mesh()
        self._routing_parameters = Mod_Gamma_Routing_Parameters.type_routing_parameter()
        self._routing_states = Mod_Gamma_Routing_States.type_routing_states()
        self._routing_memory = Mod_Gamma_Routing_Memory.type_routing_memory()
        self._routing_results = Mod_Gamma_Routing_Results.type_routing_results()

        # ~ except:
        # ~ print("failed")

        return

    def __repr__(self):

        ret = [""]

        ret = "".join(ret)

        return ret

    @property
    def routing_setup(self):

        return self._routing_setup

    @routing_setup.setter
    def routing_setup(self, value):
        self._routing_setup = value

    @property
    def routing_mesh(self):

        return self._routing_mesh

    @routing_mesh.setter
    def routing_mesh(self, value):
        self._routing_mesh = value

    @property
    def routing_states(self):

        return self._routing_states

    @routing_states.setter
    def routing_states(self, value):
        self._routing_states = value

    @property
    def routing_memory(self):

        return self._routing_memory

    @routing_memory.setter
    def routing_memory(self, value):
        self._routing_memory = value

    @property
    def routing_parameters(self):
        return self._routing_parameters

    @routing_parameters.setter
    def routing_parameters(self, value):
        self._routing_parameters = value

    def routing_setup_init(self, **kwargs):

        self.routing_setup = Mod_Gamma_Routing_Setup.type_routing_setup()

        options = kwargs

        if options is None:
            options = {}

        for key, value in options.items():
            setattr(self.routing_setup, key, value)

        # bounds calculation : default bounds
        self.routing_setup.hydrau_coef_boundaries[0] = 0.1
        self.routing_setup.hydrau_coef_boundaries[1] = 10.0
        self.routing_setup.spreading_boundaries[0] = 0.5
        self.routing_setup.spreading_boundaries[1] = 5.0
        # ~ self.routing_setup.spreading_boundaries[1]=self.routing_setup.dt+300.*(self.routing_setup.dt**0.5)

    def routing_mesh_init(self, **kwargs):

        self.routing_mesh = Mod_Gamma_Routing_Mesh.type_routing_mesh()

        options = kwargs

        if options is None:
            options = {}

        for key, value in options.items():
            setattr(self.routing_mesh, key, value)

        Mod_Gamma_Routing_Mesh.routing_mesh_self_initialisation(
            self.routing_mesh,
            self.routing_mesh.nb_nodes,
            self.routing_mesh.nb_upstream_nodes,
        )

        # read mesh here... and then...
        # call mesh_update(routing_mesh)

    def routing_mesh_update(self):

        Mod_Gamma_Routing_Mesh.mesh_update(self.routing_mesh)

    def routing_parameters_init(self, **kwargs):

        self.routing_parameters = Mod_Gamma_Routing_Parameters.type_routing_parameter()

        Mod_Gamma_Routing_Parameters.routing_parameter_self_initialisation(
            self.routing_parameters, self.routing_setup, self.routing_mesh
        )

        options = kwargs

        if options is None:
            options = {}

        for key, value in options.items():
            setattr(self.routing_parameters, key, value)

    def routing_parameters_change(self, **kwargs):

        options = kwargs

        if options is None:
            options = {}

        for key, value in options.items():
            if key == "hydraulics_coefficient":
                self.routing_parameters.hydraulics_coefficient = value

            if key == "spreading":
                self.routing_parameters.spreading = value

        self.routing_states_init()

    def routing_states_init(self):

        self.routing_states = Mod_Gamma_Routing_States.type_routing_states()
        self.routing_memory = Mod_Gamma_Routing_Memory.type_routing_memory()

        Mod_Gamma_Routing_States.routing_state_self_initialisation(
            self.routing_setup,
            self.routing_mesh,
            self.routing_parameters,
            self.routing_states,
        )

        Mod_Gamma_Interface.routing_gamma_precomputation(
            self.routing_setup,
            self.routing_mesh,
            self.routing_states,
            self.routing_memory,
        )

    def routing_memory_reset(self):

        Mod_Gamma_Routing_Memory.routing_memory_reset(self.routing_memory)

    def routing_results_init(self):

        self.routing_results = Mod_Gamma_Routing_Results.type_routing_results()

        Mod_Gamma_Routing_Results.routing_results_self_initialisation(
            self.routing_setup, self.routing_mesh, self.routing_results
        )

    def run(self, inflows, states_init=True, memory_reset=True):

        # first check if we need to force states_init
        if self.routing_setup.varying_spread == 1:

            max_spreading = self.routing_setup.spreading_boundaries[1]

            nb_spreads = (
                int(max_spreading / self.routing_setup.spreading_discretization_step) + 1
            )

        else:
            max_spreading = np.max(self.routing_parameters.spreading)
            nb_spreads = 1

        if max_spreading != self.routing_states.max_spreading:
            print(
                "Force states_init, reason change of max_spreading {max_spreading} != {self.routing_states.max_spreading}"
            )
            states_init = True

        if nb_spreads != self.routing_states.nb_spreads:
            print(
                f"Force states_init, reason change of nb_spreads {nb_spreads} != {self.routing_states.nb_spreads}"
            )
            states_init = True

        # first initialise states
        # initialise states
        if states_init:
            self.routing_states_init()

        # reset discharges states to zeros
        if memory_reset:
            self.routing_memory_reset()

        # initialise routing_results
        self.routing_results_init()

        inflows = np.array(inflows, order="F", dtype="float32")

        Mod_Gamma_Interface.routing_gamma_run(
            self.routing_setup,
            self.routing_mesh,
            self.routing_parameters,
            inflows,
            self.routing_states,
            self.routing_memory,
            self.routing_results,
        )

    def calibration(self, inflows, observations, states_init=True, memory_reset=True):

        # first initialise states
        # initialise states
        if states_init:
            self.routing_states_init()

        # reset discharges states to zeros
        if memory_reset:
            self.routing_memory_reset()

        # initialise routing_results
        self.routing_results_init()

        inflows = np.array(inflows, order="F", dtype="float32")
        observations = np.array(observations, order="F", dtype="float32")

        Mod_Gamma_Interface.routing_gamma_control(
            self.routing_setup,
            self.routing_mesh,
            self.routing_parameters,
            inflows,
            observations,
            self.routing_states,
            self.routing_memory,
            self.routing_results,
        )

    def cost_function(self, observations, qnetwork):

        # tab_cost=np.zeros(shape=(3),order='F',dtype="float32")

        # print(np.where(observations-qnetwork!=0))

        Mod_Gamma_Interface.routing_gamma_cost_function(
            self.routing_setup,
            self.routing_mesh,
            self.routing_parameters,
            observations,
            qnetwork,
            self.routing_results,
        )

    def copy(self):
        """
        Make a deepcopy of the Model.

        Returns
        -------
        Model
            A copy of Model.
        """

        copy = Model()
        copy.routing_setup = Mod_Gamma_Routing_Setup.routing_setup_copy(
            self.routing_setup
        )
        copy.routing_mesh = Mod_Gamma_Routing_Mesh.routing_mesh_copy(self.routing_mesh)
        copy.routing_parameters = Mod_Gamma_Routing_Parameters.routing_parameter_copy(
            self.routing_parameters
        )
        copy.routing_states = Mod_Gamma_Routing_States.routing_states_copy(
            self.routing_states
        )
        copy.routing_memory = Mod_Gamma_Routing_Memory.routing_memory_copy(
            self.routing_memory
        )
        copy.routing_results = Mod_Gamma_Routing_Results.routing_results_copy(
            self.routing_results
        )

        return copy
