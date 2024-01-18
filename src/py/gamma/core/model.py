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

from gamma.wrapping import (
    mod_gamma_routing_setup,
    mod_gamma_routing_mesh,
    mod_gamma_routing_parameters,
    mod_gamma_routing_states,
    mod_gamma_routing_results,
    mod_gamma_interface
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
        print("GammaRouting  Copyright (C) 2022,2023  Hydris-hydrologie, Maxime Jay-Allemand (contact: maxime.jay.allemand@hydris-hydrologie.fr)")
        print("This program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions.")
        print("")
        
        # __init__
        #self.kwargs = kwargs
        
        
        # check all arg and set default value if not set
        #then pass it to Type_Routing_Setup
        
        # ~ try:
        self._routing_setup = mod_gamma_routing_setup.type_routing_setup()
        self._routing_mesh = mod_gamma_routing_mesh.type_routing_mesh()
        self._routing_parameters = mod_gamma_routing_parameters.type_routing_parameter()
        self._routing_states = mod_gamma_routing_states.type_routing_states()
        self._routing_results = mod_gamma_routing_results.type_routing_results()
        
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
    def routing_parameters(self):
        
        return self._routing_parameters
    
    @routing_parameters.setter
    def routing_parameters(self, value):
        self._routing_parameters = value
    
    
    def routing_setup_init(self, **kwargs):
        
        self.routing_setup = mod_gamma_routing_setup.type_routing_setup()
        
        options=kwargs
        
        if options is None:
            options = {}
            
        for key, value in options.items():
            setattr(self.routing_setup, key, value)
        
        #bounds calculation : default bounds
        self.routing_setup.hydrau_coef_boundaries[0]=0.1
        self.routing_setup.hydrau_coef_boundaries[1]=10.
        self.routing_setup.spreading_boundaries[0]=0.5
        self.routing_setup.spreading_boundaries[1]=5.
        # ~ self.routing_setup.spreading_boundaries[1]=self.routing_setup.dt+300.*(self.routing_setup.dt**0.5)
        
    
    def routing_mesh_init(self,**kwargs):
        
        self.routing_mesh = mod_gamma_routing_mesh.type_routing_mesh()
        
        options=kwargs
        
        if options is None:
            options = {}
            
        for key, value in options.items():
            setattr(self.routing_mesh, key, value)
        
        mod_gamma_routing_mesh.routing_mesh_self_initialisation(self.routing_mesh,self.routing_mesh.nb_nodes,self.routing_mesh.nb_upstream_nodes)
        
        #read mesh here... and then...
        #call mesh_update(routing_mesh)
    
    def routing_mesh_update(self):
        
        mod_gamma_routing_mesh.mesh_update(self.routing_mesh)
        
    
    def routing_parameters_init(self,states_init=True,**kwargs):
        
        self.routing_parameters = mod_gamma_routing_parameters.type_routing_parameter()
        
        mod_gamma_routing_parameters.routing_parameter_self_initialisation(self.routing_parameters,self.routing_setup,self.routing_mesh)
        
        options=kwargs
        
        if options is None:
            options = {}
            
        for key, value in options.items():
            setattr(self.routing_parameters, key, value)
        
        
    
    def routing_parameters_change(self,**kwargs):
        
        options=kwargs
        
        if options is None:
            options = {}
            
        for key, value in options.items():
            if (key=="hydraulics_coefficient") :
                self.routing_parameters.hydraulics_coefficient=value
            
            if (key=="spreading") :
                self.routing_parameters.spreading=value
            
        self.routing_states_init()
    
    
    def routing_states_init(self):
        
        self.routing_states = mod_gamma_routing_states.type_routing_states()
        
        mod_gamma_routing_states.routing_state_self_initialisation(self.routing_setup,self.routing_mesh,self.routing_parameters,self.routing_states)
        
        mod_gamma_interface.routing_gamma_precomputation(self.routing_setup,self.routing_mesh,self.routing_states)
    
    
    def routing_states_reset(self):
        
        mod_gamma_routing_states.routing_states_reset(self.routing_states)
    
    
    def routing_results_init(self):
        
        self.routing_results = mod_gamma_routing_results.type_routing_results()
        
        mod_gamma_routing_results.routing_results_self_initialisation(self.routing_setup,self.routing_mesh,self.routing_results)
    
    
    def run(self,inflows,states_init=True,states_reset=True):
        
        #first initialise states
        #initialise states
        if (states_init) :
            self.routing_states_init()
        
        #reset discharges states to zeros
        if (states_reset) :
            self.routing_states_reset()
        
        #initialise routing_results
        self.routing_results_init()
        
        inflows=np.array(inflows,order='F',dtype="float32")
        
        mod_gamma_interface.routing_gamma_run(self.routing_setup,self.routing_mesh,self.routing_parameters,
        inflows,self.routing_states,self.routing_results)
        
    
    
    def calibration(self,inflows,observations,states_init=True,states_reset=True):
        
        #first initialise states
        #initialise states
        if (states_init) :
            self.routing_states_init()
        
        #reset discharges states to zeros
        if (states_reset) :
            self.routing_states_reset()
        
        #initialise routing_results
        self.routing_results_init()
        
        inflows=np.array(inflows,order='F',dtype="float32")
        observations=np.array(observations,order='F',dtype="float32")
        
        mod_gamma_interface.routing_gamma_control(self.routing_setup,self.routing_mesh,self.routing_parameters,inflows,observations,self.routing_states,self.routing_results)
        
    
    
    def cost_function(self,observations,qnetwork):
        
        #tab_cost=np.zeros(shape=(3),order='F',dtype="float32")
        
        #print(np.where(observations-qnetwork!=0))
        
        mod_gamma_interface.routing_gamma_cost_function(self.routing_setup,self.routing_mesh,self.routing_parameters,observations,qnetwork,self.routing_results)
        
        
    
    
    def copy(self):
        """
        Make a deepcopy of the Model.
        
        Returns
        -------
        Model
            A copy of Model.
        """
        
        copy = Model()
        copy.routing_setup = mod_gamma_routing_setup.routing_setup_copy(self.routing_setup)
        copy.routing_mesh = mod_gamma_routing_mesh.routing_mesh_copy(self.routing_mesh)
        copy.routing_parameters = mod_gamma_routing_parameters.routing_parameter_copy(self.routing_parameters)
        copy.routing_states = mod_gamma_routing_states.routing_states_copy(self.routing_states)
        copy.routing_results = mod_gamma_routing_results.routing_results_copy(self.routing_results)
        
        return copy

