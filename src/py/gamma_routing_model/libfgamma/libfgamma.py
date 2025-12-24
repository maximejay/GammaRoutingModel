from __future__ import print_function, absolute_import, division
from . import _libfgamma
import f90wrap.runtime
import logging
import numpy
import warnings
import weakref

class Mod_Gamma_Routing_Setup(f90wrap.runtime.FortranModule):
    """
    Module mod_gamma_routing_setup
    Defined at mod_routing_setup.f90 lines 8-156
    """
    @f90wrap.runtime.register_class("libfgamma.type_routing_setup")
    class type_routing_setup(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=type_routing_setup)
        Defined at mod_routing_setup.f90 lines 10-30
        """
        def __init__(self, handle=None):
            """
            Automatically generated constructor for type_routing_setup
            
            self = Type_Routing_Setup()
            Defined at mod_routing_setup.f90 lines 10-30
            
            Returns
            -------
            this : Type_Routing_Setup
                Object to be constructed
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            if isinstance(handle, numpy.ndarray) and handle.ndim == 1 and handle.dtype.num \
                == 5:
                self._handle = handle
                self._alloc = True
            else:
                result = \
                    _libfgamma.f90wrap_mod_gamma_routing_setup__type_routing_setup_initialise()
                self._handle = result[0] if isinstance(result, tuple) else result
                self._alloc = True
            self._setup_finalizer()
        
        def _setup_finalizer(self):
            """Set up weak reference destructor to prevent Fortran memory leaks."""
            if self._alloc:
                destructor = getattr(_libfgamma, \
                    "f90wrap_mod_gamma_routing_setup__type_routing_setup_finalise")
                self._finalizer = weakref.finalize(self, destructor, self._handle)
        
        @property
        def npdt(self):
            """
            Element npdt ftype=integer  pytype=int32
            Defined at mod_routing_setup.f90 line 11
            """
            return _libfgamma.f90wrap_type_routing_setup__get__npdt(self._handle)
        
        @npdt.setter
        def npdt(self, npdt):
            _libfgamma.f90wrap_type_routing_setup__set__npdt(self._handle, npdt)
        
        @property
        def dt(self):
            """
            Element dt ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 12
            """
            return _libfgamma.f90wrap_type_routing_setup__get__dt(self._handle)
        
        @dt.setter
        def dt(self, dt):
            _libfgamma.f90wrap_type_routing_setup__set__dt(self._handle, dt)
        
        @property
        def vmin(self):
            """
            Element vmin ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 13
            """
            return _libfgamma.f90wrap_type_routing_setup__get__vmin(self._handle)
        
        @vmin.setter
        def vmin(self, vmin):
            _libfgamma.f90wrap_type_routing_setup__set__vmin(self._handle, vmin)
        
        @property
        def vmax(self):
            """
            Element vmax ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 14
            """
            return _libfgamma.f90wrap_type_routing_setup__get__vmax(self._handle)
        
        @vmax.setter
        def vmax(self, vmax):
            _libfgamma.f90wrap_type_routing_setup__set__vmax(self._handle, vmax)
        
        @property
        def spreading_gamma_threshold(self):
            """
            Element spreading_gamma_threshold ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 15
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__spreading_gamma_threshold(self._handle)
        
        @spreading_gamma_threshold.setter
        def spreading_gamma_threshold(self, spreading_gamma_threshold):
            _libfgamma.f90wrap_type_routing_setup__set__spreading_gamma_threshold(self._handle, \
                spreading_gamma_threshold)
        
        @property
        def elongation_factor(self):
            """
            Element elongation_factor ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 16
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__elongation_factor(self._handle)
        
        @elongation_factor.setter
        def elongation_factor(self, elongation_factor):
            _libfgamma.f90wrap_type_routing_setup__set__elongation_factor(self._handle, \
                elongation_factor)
        
        @property
        def mode_discretization_step(self):
            """
            Element mode_discretization_step ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 17
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__mode_discretization_step(self._handle)
        
        @mode_discretization_step.setter
        def mode_discretization_step(self, mode_discretization_step):
            _libfgamma.f90wrap_type_routing_setup__set__mode_discretization_step(self._handle, \
                mode_discretization_step)
        
        @property
        def spreading_discretization_step(self):
            """
            Element spreading_discretization_step ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 18
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__spreading_discretization_step(self._handle)
        
        @spreading_discretization_step.setter
        def spreading_discretization_step(self, spreading_discretization_step):
            _libfgamma.f90wrap_type_routing_setup__set__spreading_discretization_step(self._handle, \
                spreading_discretization_step)
        
        @property
        def hydrau_coef_boundaries(self):
            """
            Element hydrau_coef_boundaries ftype=real pytype=float array
            Defined at mod_routing_setup.f90 line 19
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_setup__array__hydrau_coef_boundaries(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            hydrau_coef_boundaries = self._arrays.get(array_hash)
            if hydrau_coef_boundaries is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if hydrau_coef_boundaries.ctypes.data != array_handle:
                    hydrau_coef_boundaries = None
            if hydrau_coef_boundaries is None:
                try:
                    hydrau_coef_boundaries = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_setup__array__hydrau_coef_boundaries)
                except TypeError:
                    hydrau_coef_boundaries = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = hydrau_coef_boundaries
            return hydrau_coef_boundaries
        
        @hydrau_coef_boundaries.setter
        def hydrau_coef_boundaries(self, hydrau_coef_boundaries):
            self.hydrau_coef_boundaries[...] = hydrau_coef_boundaries
        
        @property
        def spreading_boundaries(self):
            """
            Element spreading_boundaries ftype=real pytype=float array
            Defined at mod_routing_setup.f90 line 20
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_setup__array__spreading_boundaries(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            spreading_boundaries = self._arrays.get(array_hash)
            if spreading_boundaries is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if spreading_boundaries.ctypes.data != array_handle:
                    spreading_boundaries = None
            if spreading_boundaries is None:
                try:
                    spreading_boundaries = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_setup__array__spreading_boundaries)
                except TypeError:
                    spreading_boundaries = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = spreading_boundaries
            return spreading_boundaries
        
        @spreading_boundaries.setter
        def spreading_boundaries(self, spreading_boundaries):
            self.spreading_boundaries[...] = spreading_boundaries
        
        @property
        def ponderation_regul(self):
            """
            Element ponderation_regul ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 21
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__ponderation_regul(self._handle)
        
        @ponderation_regul.setter
        def ponderation_regul(self, ponderation_regul):
            _libfgamma.f90wrap_type_routing_setup__set__ponderation_regul(self._handle, \
                ponderation_regul)
        
        @property
        def ponderation_cost(self):
            """
            Element ponderation_cost ftype=real  pytype=float32
            Defined at mod_routing_setup.f90 line 22
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__ponderation_cost(self._handle)
        
        @ponderation_cost.setter
        def ponderation_cost(self, ponderation_cost):
            _libfgamma.f90wrap_type_routing_setup__set__ponderation_cost(self._handle, \
                ponderation_cost)
        
        @property
        def pdt_start_optim(self):
            """
            Element pdt_start_optim ftype=integer  pytype=int32
            Defined at mod_routing_setup.f90 line 23
            """
            return _libfgamma.f90wrap_type_routing_setup__get__pdt_start_optim(self._handle)
        
        @pdt_start_optim.setter
        def pdt_start_optim(self, pdt_start_optim):
            _libfgamma.f90wrap_type_routing_setup__set__pdt_start_optim(self._handle, \
                pdt_start_optim)
        
        @property
        def auto_reg(self):
            """
            Element auto_reg ftype=integer  pytype=int32
            Defined at mod_routing_setup.f90 line 24
            """
            return _libfgamma.f90wrap_type_routing_setup__get__auto_reg(self._handle)
        
        @auto_reg.setter
        def auto_reg(self, auto_reg):
            _libfgamma.f90wrap_type_routing_setup__set__auto_reg(self._handle, auto_reg)
        
        @property
        def hydraulics_coef_uniform(self):
            """
            Element hydraulics_coef_uniform ftype=integer  pytype=int32
            Defined at mod_routing_setup.f90 line 25
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__hydraulics_coef_uniform(self._handle)
        
        @hydraulics_coef_uniform.setter
        def hydraulics_coef_uniform(self, hydraulics_coef_uniform):
            _libfgamma.f90wrap_type_routing_setup__set__hydraulics_coef_uniform(self._handle, \
                hydraulics_coef_uniform)
        
        @property
        def spreading_uniform(self):
            """
            Element spreading_uniform ftype=integer  pytype=int32
            Defined at mod_routing_setup.f90 line 26
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__spreading_uniform(self._handle)
        
        @spreading_uniform.setter
        def spreading_uniform(self, spreading_uniform):
            _libfgamma.f90wrap_type_routing_setup__set__spreading_uniform(self._handle, \
                spreading_uniform)
        
        @property
        def iter_max(self):
            """
            Element iter_max ftype=integer  pytype=int32
            Defined at mod_routing_setup.f90 line 27
            """
            return _libfgamma.f90wrap_type_routing_setup__get__iter_max(self._handle)
        
        @iter_max.setter
        def iter_max(self, iter_max):
            _libfgamma.f90wrap_type_routing_setup__set__iter_max(self._handle, iter_max)
        
        @property
        def velocity_computation(self):
            """
            Element velocity_computation ftype=character(3) pytype=str
            Defined at mod_routing_setup.f90 line 28
            """
            return \
                _libfgamma.f90wrap_type_routing_setup__get__velocity_computation(self._handle)
        
        @velocity_computation.setter
        def velocity_computation(self, velocity_computation):
            _libfgamma.f90wrap_type_routing_setup__set__velocity_computation(self._handle, \
                velocity_computation)
        
        @property
        def criteria(self):
            """
            Element criteria ftype=character(4) pytype=str
            Defined at mod_routing_setup.f90 line 29
            """
            return _libfgamma.f90wrap_type_routing_setup__get__criteria(self._handle)
        
        @criteria.setter
        def criteria(self, criteria):
            _libfgamma.f90wrap_type_routing_setup__set__criteria(self._handle, criteria)
        
        @property
        def varying_spread(self):
            """
            Element varying_spread ftype=logical pytype=bool
            Defined at mod_routing_setup.f90 line 30
            """
            return _libfgamma.f90wrap_type_routing_setup__get__varying_spread(self._handle)
        
        @varying_spread.setter
        def varying_spread(self, varying_spread):
            _libfgamma.f90wrap_type_routing_setup__set__varying_spread(self._handle, \
                varying_spread)
        
        def __str__(self):
            ret = ['<type_routing_setup>{\n']
            ret.append('    npdt : ')
            ret.append(repr(self.npdt))
            ret.append(',\n    dt : ')
            ret.append(repr(self.dt))
            ret.append(',\n    vmin : ')
            ret.append(repr(self.vmin))
            ret.append(',\n    vmax : ')
            ret.append(repr(self.vmax))
            ret.append(',\n    spreading_gamma_threshold : ')
            ret.append(repr(self.spreading_gamma_threshold))
            ret.append(',\n    elongation_factor : ')
            ret.append(repr(self.elongation_factor))
            ret.append(',\n    mode_discretization_step : ')
            ret.append(repr(self.mode_discretization_step))
            ret.append(',\n    spreading_discretization_step : ')
            ret.append(repr(self.spreading_discretization_step))
            ret.append(',\n    hydrau_coef_boundaries : ')
            ret.append(repr(self.hydrau_coef_boundaries))
            ret.append(',\n    spreading_boundaries : ')
            ret.append(repr(self.spreading_boundaries))
            ret.append(',\n    ponderation_regul : ')
            ret.append(repr(self.ponderation_regul))
            ret.append(',\n    ponderation_cost : ')
            ret.append(repr(self.ponderation_cost))
            ret.append(',\n    pdt_start_optim : ')
            ret.append(repr(self.pdt_start_optim))
            ret.append(',\n    auto_reg : ')
            ret.append(repr(self.auto_reg))
            ret.append(',\n    hydraulics_coef_uniform : ')
            ret.append(repr(self.hydraulics_coef_uniform))
            ret.append(',\n    spreading_uniform : ')
            ret.append(repr(self.spreading_uniform))
            ret.append(',\n    iter_max : ')
            ret.append(repr(self.iter_max))
            ret.append(',\n    velocity_computation : ')
            ret.append(repr(self.velocity_computation))
            ret.append(',\n    criteria : ')
            ret.append(repr(self.criteria))
            ret.append(',\n    varying_spread : ')
            ret.append(repr(self.varying_spread))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def routing_setup_self_initialisation(self, npdt=None, dt=None, vmin=None, \
        vmax=None, elongation_factor=None, mode_discretization_step=None, \
        spreading_discretization_step=None, ponderation_regul=None, \
        ponderation_cost=None, auto_reg=None, hydraulics_coef_uniform=None, \
        spreading_uniform=None, iter_max=None, velocity_computation=None, \
        criteria=None, varying_spread=None, interface_call=False):
        """
        routing_setup_self_initialisation(self[, npdt, dt, vmin, vmax, \
            elongation_factor, mode_discretization_step, spreading_discretization_step, \
            ponderation_regul, ponderation_cost, auto_reg, hydraulics_coef_uniform, \
            spreading_uniform, iter_max, velocity_computation, criteria, \
            varying_spread])
        Defined at mod_routing_setup.f90 lines 36-135
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        npdt : int32
        dt : float32
        vmin : float32
        vmax : float32
        elongation_factor : float32
        mode_discretization_step : float32
        spreading_discretization_step : float32
        ponderation_regul : float32
        ponderation_cost : float32
        auto_reg : int32
        hydraulics_coef_uniform : int32
        spreading_uniform : int32
        iter_max : int32
        velocity_computation : str
        criteria : str
        varying_spread : bool
        """
        _libfgamma.f90wrap_mod_gamma_routing_setup__routing_setup_self_initiala6ab(routing_setup=self._handle, \
            npdt=npdt, dt=dt, vmin=vmin, vmax=vmax, elongation_factor=elongation_factor, \
            mode_discretization_step=mode_discretization_step, \
            spreading_discretization_step=spreading_discretization_step, \
            ponderation_regul=ponderation_regul, ponderation_cost=ponderation_cost, \
            auto_reg=auto_reg, hydraulics_coef_uniform=hydraulics_coef_uniform, \
            spreading_uniform=spreading_uniform, iter_max=iter_max, \
            velocity_computation=velocity_computation, criteria=criteria, \
            varying_spread=varying_spread)
    
    @staticmethod
    def routing_setup_clear(self, interface_call=False):
        """
        routing_setup_clear(self)
        Defined at mod_routing_setup.f90 lines 137-151
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        """
        _libfgamma.f90wrap_mod_gamma_routing_setup__routing_setup_clear(routing_setup=self._handle)
    
    @staticmethod
    def routing_setup_copy(self, interface_call=False):
        """
        object_copy = routing_setup_copy(self)
        Defined at mod_routing_setup.f90 lines 153-156
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        
        Returns
        -------
        object_copy : Type_Routing_Setup
        """
        object_copy = \
            _libfgamma.f90wrap_mod_gamma_routing_setup__routing_setup_copy(routing_setup=self._handle)
        object_copy = \
            f90wrap.runtime.lookup_class("libfgamma.type_routing_setup").from_handle(object_copy, \
            alloc=True)
        object_copy._setup_finalizer()
        return object_copy
    
    _dt_array_initialisers = []
    
    

mod_gamma_routing_setup = Mod_Gamma_Routing_Setup()

class Mod_Gamma_Routing_Mesh(f90wrap.runtime.FortranModule):
    """
    Module mod_gamma_routing_mesh
    Defined at mod_routing_mesh.f90 lines 8-220
    """
    @f90wrap.runtime.register_class("libfgamma.type_routing_mesh")
    class type_routing_mesh(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=type_routing_mesh)
        Defined at mod_routing_mesh.f90 lines 10-23
        """
        def __init__(self, handle=None):
            """
            Automatically generated constructor for type_routing_mesh
            
            self = Type_Routing_Mesh()
            Defined at mod_routing_mesh.f90 lines 10-23
            
            Returns
            -------
            this : Type_Routing_Mesh
                Object to be constructed
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            if isinstance(handle, numpy.ndarray) and handle.ndim == 1 and handle.dtype.num \
                == 5:
                self._handle = handle
                self._alloc = True
            else:
                result = \
                    _libfgamma.f90wrap_mod_gamma_routing_mesh__type_routing_mesh_initialise()
                self._handle = result[0] if isinstance(result, tuple) else result
                self._alloc = True
            self._setup_finalizer()
        
        def _setup_finalizer(self):
            """Set up weak reference destructor to prevent Fortran memory leaks."""
            if self._alloc:
                destructor = getattr(_libfgamma, \
                    "f90wrap_mod_gamma_routing_mesh__type_routing_mesh_finalise")
                self._finalizer = weakref.finalize(self, destructor, self._handle)
        
        @property
        def nb_nodes(self):
            """
            Element nb_nodes ftype=integer  pytype=int32
            Defined at mod_routing_mesh.f90 line 11
            """
            return _libfgamma.f90wrap_type_routing_mesh__get__nb_nodes(self._handle)
        
        @nb_nodes.setter
        def nb_nodes(self, nb_nodes):
            _libfgamma.f90wrap_type_routing_mesh__set__nb_nodes(self._handle, nb_nodes)
        
        @property
        def nb_upstream_nodes(self):
            """
            Element nb_upstream_nodes ftype=integer  pytype=int32
            Defined at mod_routing_mesh.f90 line 12
            """
            return \
                _libfgamma.f90wrap_type_routing_mesh__get__nb_upstream_nodes(self._handle)
        
        @nb_upstream_nodes.setter
        def nb_upstream_nodes(self, nb_upstream_nodes):
            _libfgamma.f90wrap_type_routing_mesh__set__nb_upstream_nodes(self._handle, \
                nb_upstream_nodes)
        
        @property
        def dx(self):
            """
            Element dx ftype=real pytype=float array
            Defined at mod_routing_mesh.f90 line 13
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__dx(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            dx = self._arrays.get(array_hash)
            if dx is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if dx.ctypes.data != array_handle:
                    dx = None
            if dx is None:
                try:
                    dx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__dx)
                except TypeError:
                    dx = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = dx
            return dx
        
        @dx.setter
        def dx(self, dx):
            self.dx[...] = dx
        
        @property
        def index_varying_dx(self):
            """
            Element index_varying_dx ftype=integer pytype=int array
            Defined at mod_routing_mesh.f90 line 14
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__index_varying_dx(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            index_varying_dx = self._arrays.get(array_hash)
            if index_varying_dx is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if index_varying_dx.ctypes.data != array_handle:
                    index_varying_dx = None
            if index_varying_dx is None:
                try:
                    index_varying_dx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__index_varying_dx)
                except TypeError:
                    index_varying_dx = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = index_varying_dx
            return index_varying_dx
        
        @index_varying_dx.setter
        def index_varying_dx(self, index_varying_dx):
            self.index_varying_dx[...] = index_varying_dx
        
        @property
        def varying_dx(self):
            """
            Element varying_dx ftype=real pytype=float array
            Defined at mod_routing_mesh.f90 line 15
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__varying_dx(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            varying_dx = self._arrays.get(array_hash)
            if varying_dx is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if varying_dx.ctypes.data != array_handle:
                    varying_dx = None
            if varying_dx is None:
                try:
                    varying_dx = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__varying_dx)
                except TypeError:
                    varying_dx = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = varying_dx
            return varying_dx
        
        @varying_dx.setter
        def varying_dx(self, varying_dx):
            self.varying_dx[...] = varying_dx
        
        @property
        def nodes_indexes(self):
            """
            Element nodes_indexes ftype=integer pytype=int array
            Defined at mod_routing_mesh.f90 line 16
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__nodes_indexes(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            nodes_indexes = self._arrays.get(array_hash)
            if nodes_indexes is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if nodes_indexes.ctypes.data != array_handle:
                    nodes_indexes = None
            if nodes_indexes is None:
                try:
                    nodes_indexes = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__nodes_indexes)
                except TypeError:
                    nodes_indexes = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = nodes_indexes
            return nodes_indexes
        
        @nodes_indexes.setter
        def nodes_indexes(self, nodes_indexes):
            self.nodes_indexes[...] = nodes_indexes
        
        @property
        def nodes_names(self):
            """
            Element nodes_names ftype=character(100) pytype=str array
            Defined at mod_routing_mesh.f90 line 17
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__nodes_names(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            nodes_names = self._arrays.get(array_hash)
            if nodes_names is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if nodes_names.ctypes.data != array_handle:
                    nodes_names = None
            if nodes_names is None:
                try:
                    nodes_names = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__nodes_names)
                except TypeError:
                    nodes_names = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = nodes_names
            return nodes_names
        
        @nodes_names.setter
        def nodes_names(self, nodes_names):
            self.nodes_names[...] = nodes_names
        
        @property
        def surface(self):
            """
            Element surface ftype=real pytype=float array
            Defined at mod_routing_mesh.f90 line 18
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__surface(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            surface = self._arrays.get(array_hash)
            if surface is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if surface.ctypes.data != array_handle:
                    surface = None
            if surface is None:
                try:
                    surface = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__surface)
                except TypeError:
                    surface = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = surface
            return surface
        
        @surface.setter
        def surface(self, surface):
            self.surface[...] = surface
        
        @property
        def cumulated_surface(self):
            """
            Element cumulated_surface ftype=real pytype=float array
            Defined at mod_routing_mesh.f90 line 19
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__cumulated_surface(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            cumulated_surface = self._arrays.get(array_hash)
            if cumulated_surface is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if cumulated_surface.ctypes.data != array_handle:
                    cumulated_surface = None
            if cumulated_surface is None:
                try:
                    cumulated_surface = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__cumulated_surface)
                except TypeError:
                    cumulated_surface = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = cumulated_surface
            return cumulated_surface
        
        @cumulated_surface.setter
        def cumulated_surface(self, cumulated_surface):
            self.cumulated_surface[...] = cumulated_surface
        
        @property
        def nodes_linker(self):
            """
            Element nodes_linker ftype=integer pytype=int array
            Defined at mod_routing_mesh.f90 line 20
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__nodes_linker(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            nodes_linker = self._arrays.get(array_hash)
            if nodes_linker is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if nodes_linker.ctypes.data != array_handle:
                    nodes_linker = None
            if nodes_linker is None:
                try:
                    nodes_linker = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__nodes_linker)
                except TypeError:
                    nodes_linker = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = nodes_linker
            return nodes_linker
        
        @nodes_linker.setter
        def nodes_linker(self, nodes_linker):
            self.nodes_linker[...] = nodes_linker
        
        @property
        def upstream_to_downstream_nodes(self):
            """
            Element upstream_to_downstream_nodes ftype=integer pytype=int array
            Defined at mod_routing_mesh.f90 line 21
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__upstream_to_downstream_nodes(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            upstream_to_downstream_nodes = self._arrays.get(array_hash)
            if upstream_to_downstream_nodes is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if upstream_to_downstream_nodes.ctypes.data != array_handle:
                    upstream_to_downstream_nodes = None
            if upstream_to_downstream_nodes is None:
                try:
                    upstream_to_downstream_nodes = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__upstream_to_downstream_nodes)
                except TypeError:
                    upstream_to_downstream_nodes = f90wrap.runtime.direct_c_array(array_type, \
                        array_shape, array_handle)
                self._arrays[array_hash] = upstream_to_downstream_nodes
            return upstream_to_downstream_nodes
        
        @upstream_to_downstream_nodes.setter
        def upstream_to_downstream_nodes(self, upstream_to_downstream_nodes):
            self.upstream_to_downstream_nodes[...] = upstream_to_downstream_nodes
        
        @property
        def cum_node_index(self):
            """
            Element cum_node_index ftype=integer pytype=int array
            Defined at mod_routing_mesh.f90 line 22
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__cum_node_index(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            cum_node_index = self._arrays.get(array_hash)
            if cum_node_index is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if cum_node_index.ctypes.data != array_handle:
                    cum_node_index = None
            if cum_node_index is None:
                try:
                    cum_node_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__cum_node_index)
                except TypeError:
                    cum_node_index = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = cum_node_index
            return cum_node_index
        
        @cum_node_index.setter
        def cum_node_index(self, cum_node_index):
            self.cum_node_index[...] = cum_node_index
        
        @property
        def controlled_nodes(self):
            """
            Element controlled_nodes ftype=integer pytype=int array
            Defined at mod_routing_mesh.f90 line 23
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_mesh__array__controlled_nodes(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            controlled_nodes = self._arrays.get(array_hash)
            if controlled_nodes is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if controlled_nodes.ctypes.data != array_handle:
                    controlled_nodes = None
            if controlled_nodes is None:
                try:
                    controlled_nodes = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_mesh__array__controlled_nodes)
                except TypeError:
                    controlled_nodes = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = controlled_nodes
            return controlled_nodes
        
        @controlled_nodes.setter
        def controlled_nodes(self, controlled_nodes):
            self.controlled_nodes[...] = controlled_nodes
        
        def __str__(self):
            ret = ['<type_routing_mesh>{\n']
            ret.append('    nb_nodes : ')
            ret.append(repr(self.nb_nodes))
            ret.append(',\n    nb_upstream_nodes : ')
            ret.append(repr(self.nb_upstream_nodes))
            ret.append(',\n    dx : ')
            ret.append(repr(self.dx))
            ret.append(',\n    index_varying_dx : ')
            ret.append(repr(self.index_varying_dx))
            ret.append(',\n    varying_dx : ')
            ret.append(repr(self.varying_dx))
            ret.append(',\n    nodes_indexes : ')
            ret.append(repr(self.nodes_indexes))
            ret.append(',\n    nodes_names : ')
            ret.append(repr(self.nodes_names))
            ret.append(',\n    surface : ')
            ret.append(repr(self.surface))
            ret.append(',\n    cumulated_surface : ')
            ret.append(repr(self.cumulated_surface))
            ret.append(',\n    nodes_linker : ')
            ret.append(repr(self.nodes_linker))
            ret.append(',\n    upstream_to_downstream_nodes : ')
            ret.append(repr(self.upstream_to_downstream_nodes))
            ret.append(',\n    cum_node_index : ')
            ret.append(repr(self.cum_node_index))
            ret.append(',\n    controlled_nodes : ')
            ret.append(repr(self.controlled_nodes))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def routing_mesh_self_initialisation(self, nb_nodes=None, \
        nb_upstream_nodes=None, dx=None, interface_call=False):
        """
        routing_mesh_self_initialisation(self[, nb_nodes, nb_upstream_nodes, dx])
        Defined at mod_routing_mesh.f90 lines 26-78
        
        Parameters
        ----------
        routing_mesh : Type_Routing_Mesh
        nb_nodes : int32
        nb_upstream_nodes : int32
        dx : float32
        """
        _libfgamma.f90wrap_mod_gamma_routing_mesh__routing_mesh_self_initialis7f71(routing_mesh=self._handle, \
            nb_nodes=nb_nodes, nb_upstream_nodes=nb_upstream_nodes, dx=dx)
    
    @staticmethod
    def routing_mesh_clear(self, interface_call=False):
        """
        routing_mesh_clear(self)
        Defined at mod_routing_mesh.f90 lines 80-94
        
        Parameters
        ----------
        routing_mesh : Type_Routing_Mesh
        """
        _libfgamma.f90wrap_mod_gamma_routing_mesh__routing_mesh_clear(routing_mesh=self._handle)
    
    @staticmethod
    def mesh_update(self, interface_call=False):
        """
        mesh_update(self)
        Defined at mod_routing_mesh.f90 lines 96-112
        
        Parameters
        ----------
        routing_mesh : Type_Routing_Mesh
        """
        _libfgamma.f90wrap_mod_gamma_routing_mesh__mesh_update(routing_mesh=self._handle)
    
    @staticmethod
    def mesh_compute_cumulated_surface(self, interface_call=False):
        """
        mesh_compute_cumulated_surface(self)
        Defined at mod_routing_mesh.f90 lines 114-142
        
        Parameters
        ----------
        routing_mesh : Type_Routing_Mesh
        """
        _libfgamma.f90wrap_mod_gamma_routing_mesh__mesh_compute_cumulated_surface(routing_mesh=self._handle)
    
    @staticmethod
    def mesh_compute_cumulated_node_index(self, interface_call=False):
        """
        mesh_compute_cumulated_node_index(self)
        Defined at mod_routing_mesh.f90 lines 144-172
        
        Parameters
        ----------
        routing_mesh : Type_Routing_Mesh
        """
        _libfgamma.f90wrap_mod_gamma_routing_mesh__mesh_compute_cumulated_node3f82(routing_mesh=self._handle)
    
    @staticmethod
    def mesh_uniq_dx(self, interface_call=False):
        """
        mesh_uniq_dx(self)
        Defined at mod_routing_mesh.f90 lines 174-215
        
        Parameters
        ----------
        routing_mesh : Type_Routing_Mesh
        """
        _libfgamma.f90wrap_mod_gamma_routing_mesh__mesh_uniq_dx(routing_mesh=self._handle)
    
    @staticmethod
    def routing_mesh_copy(self, interface_call=False):
        """
        object_copy = routing_mesh_copy(self)
        Defined at mod_routing_mesh.f90 lines 217-220
        
        Parameters
        ----------
        routing_mesh : Type_Routing_Mesh
        
        Returns
        -------
        object_copy : Type_Routing_Mesh
        """
        object_copy = \
            _libfgamma.f90wrap_mod_gamma_routing_mesh__routing_mesh_copy(routing_mesh=self._handle)
        object_copy = \
            f90wrap.runtime.lookup_class("libfgamma.type_routing_mesh").from_handle(object_copy, \
            alloc=True)
        object_copy._setup_finalizer()
        return object_copy
    
    _dt_array_initialisers = []
    
    

mod_gamma_routing_mesh = Mod_Gamma_Routing_Mesh()

class Mod_Gamma_Routing_Parameters(f90wrap.runtime.FortranModule):
    """
    Module mod_gamma_routing_parameters
    Defined at mod_routing_parameters.f90 lines 8-78
    """
    @f90wrap.runtime.register_class("libfgamma.type_routing_parameter")
    class type_routing_parameter(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=type_routing_parameter)
        Defined at mod_routing_parameters.f90 lines 10-12
        """
        def __init__(self, handle=None):
            """
            Automatically generated constructor for type_routing_parameter
            
            self = Type_Routing_Parameter()
            Defined at mod_routing_parameters.f90 lines 10-12
            
            Returns
            -------
            this : Type_Routing_Parameter
                Object to be constructed
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            if isinstance(handle, numpy.ndarray) and handle.ndim == 1 and handle.dtype.num \
                == 5:
                self._handle = handle
                self._alloc = True
            else:
                result = \
                    _libfgamma.f90wrap_mod_gamma_routing_parameters__type_routing_paramete1ecf()
                self._handle = result[0] if isinstance(result, tuple) else result
                self._alloc = True
            self._setup_finalizer()
        
        def _setup_finalizer(self):
            """Set up weak reference destructor to prevent Fortran memory leaks."""
            if self._alloc:
                destructor = getattr(_libfgamma, \
                    "f90wrap_mod_gamma_routing_parameters__type_routing_paramete5302")
                self._finalizer = weakref.finalize(self, destructor, self._handle)
        
        @property
        def hydraulics_coefficient(self):
            """
            Element hydraulics_coefficient ftype=real pytype=float array
            Defined at mod_routing_parameters.f90 line 11
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_parameter__array__hydraulics_coefficient(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            hydraulics_coefficient = self._arrays.get(array_hash)
            if hydraulics_coefficient is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if hydraulics_coefficient.ctypes.data != array_handle:
                    hydraulics_coefficient = None
            if hydraulics_coefficient is None:
                try:
                    hydraulics_coefficient = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_parameter__array__hydraulics_coefficient)
                except TypeError:
                    hydraulics_coefficient = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = hydraulics_coefficient
            return hydraulics_coefficient
        
        @hydraulics_coefficient.setter
        def hydraulics_coefficient(self, hydraulics_coefficient):
            self.hydraulics_coefficient[...] = hydraulics_coefficient
        
        @property
        def spreading(self):
            """
            Element spreading ftype=real pytype=float array
            Defined at mod_routing_parameters.f90 line 12
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_parameter__array__spreading(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            spreading = self._arrays.get(array_hash)
            if spreading is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if spreading.ctypes.data != array_handle:
                    spreading = None
            if spreading is None:
                try:
                    spreading = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_parameter__array__spreading)
                except TypeError:
                    spreading = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = spreading
            return spreading
        
        @spreading.setter
        def spreading(self, spreading):
            self.spreading[...] = spreading
        
        def __str__(self):
            ret = ['<type_routing_parameter>{\n']
            ret.append('    hydraulics_coefficient : ')
            ret.append(repr(self.hydraulics_coefficient))
            ret.append(',\n    spreading : ')
            ret.append(repr(self.spreading))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def routing_parameter_self_initialisation(self, routing_setup, routing_mesh, \
        hydraulics_coefficient=None, spreading=None, interface_call=False):
        """
        routing_parameter_self_initialisation(self, routing_setup, routing_mesh[, \
            hydraulics_coefficient, spreading])
        Defined at mod_routing_parameters.f90 lines 16-57
        
        Parameters
        ----------
        routing_parameter : Type_Routing_Parameter
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        hydraulics_coefficient : float32
        spreading : float32
        """
        _libfgamma.f90wrap_mod_gamma_routing_parameters__routing_parameter_selcac8(routing_parameter=self._handle, \
            routing_setup=routing_setup._handle, routing_mesh=routing_mesh._handle, \
            hydraulics_coefficient=hydraulics_coefficient, spreading=spreading)
    
    @staticmethod
    def routing_parameter_clear(self, interface_call=False):
        """
        routing_parameter_clear(self)
        Defined at mod_routing_parameters.f90 lines 59-73
        
        Parameters
        ----------
        routing_parameter : Type_Routing_Parameter
        """
        _libfgamma.f90wrap_mod_gamma_routing_parameters__routing_parameter_clear(routing_parameter=self._handle)
    
    @staticmethod
    def routing_parameter_copy(self, interface_call=False):
        """
        object_copy = routing_parameter_copy(self)
        Defined at mod_routing_parameters.f90 lines 75-78
        
        Parameters
        ----------
        routing_parameter : Type_Routing_Parameter
        
        Returns
        -------
        object_copy : Type_Routing_Parameter
        """
        object_copy = \
            _libfgamma.f90wrap_mod_gamma_routing_parameters__routing_parameter_copy(routing_parameter=self._handle)
        object_copy = \
            f90wrap.runtime.lookup_class("libfgamma.type_routing_parameter").from_handle(object_copy, \
            alloc=True)
        object_copy._setup_finalizer()
        return object_copy
    
    _dt_array_initialisers = []
    
    

mod_gamma_routing_parameters = Mod_Gamma_Routing_Parameters()

class Mod_Gamma_Routing_States(f90wrap.runtime.FortranModule):
    """
    Module mod_gamma_routing_states
    Defined at mod_routing_states.f90 lines 8-123
    """
    @f90wrap.runtime.register_class("libfgamma.type_routing_states")
    class type_routing_states(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=type_routing_states)
        Defined at mod_routing_states.f90 lines 10-27
        """
        def __init__(self, handle=None):
            """
            Automatically generated constructor for type_routing_states
            
            self = Type_Routing_States()
            Defined at mod_routing_states.f90 lines 10-27
            
            Returns
            -------
            this : Type_Routing_States
                Object to be constructed
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            if isinstance(handle, numpy.ndarray) and handle.ndim == 1 and handle.dtype.num \
                == 5:
                self._handle = handle
                self._alloc = True
            else:
                result = \
                    _libfgamma.f90wrap_mod_gamma_routing_states__type_routing_states_initi9a74()
                self._handle = result[0] if isinstance(result, tuple) else result
                self._alloc = True
            self._setup_finalizer()
        
        def _setup_finalizer(self):
            """Set up weak reference destructor to prevent Fortran memory leaks."""
            if self._alloc:
                destructor = getattr(_libfgamma, \
                    "f90wrap_mod_gamma_routing_states__type_routing_states_finalise")
                self._finalizer = weakref.finalize(self, destructor, self._handle)
        
        @property
        def window_length(self):
            """
            Element window_length ftype=integer pytype=int array
            Defined at mod_routing_states.f90 line 11
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__window_length(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            window_length = self._arrays.get(array_hash)
            if window_length is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if window_length.ctypes.data != array_handle:
                    window_length = None
            if window_length is None:
                try:
                    window_length = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__window_length)
                except TypeError:
                    window_length = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = window_length
            return window_length
        
        @window_length.setter
        def window_length(self, window_length):
            self.window_length[...] = window_length
        
        @property
        def nb_mode(self):
            """
            Element nb_mode ftype=integer  pytype=int32
            Defined at mod_routing_states.f90 line 12
            """
            return _libfgamma.f90wrap_type_routing_states__get__nb_mode(self._handle)
        
        @nb_mode.setter
        def nb_mode(self, nb_mode):
            _libfgamma.f90wrap_type_routing_states__set__nb_mode(self._handle, nb_mode)
        
        @property
        def nb_spreads(self):
            """
            Element nb_spreads ftype=integer  pytype=int32
            Defined at mod_routing_states.f90 line 13
            """
            return _libfgamma.f90wrap_type_routing_states__get__nb_spreads(self._handle)
        
        @nb_spreads.setter
        def nb_spreads(self, nb_spreads):
            _libfgamma.f90wrap_type_routing_states__set__nb_spreads(self._handle, \
                nb_spreads)
        
        @property
        def max_mode(self):
            """
            Element max_mode ftype=real  pytype=float32
            Defined at mod_routing_states.f90 line 14
            """
            return _libfgamma.f90wrap_type_routing_states__get__max_mode(self._handle)
        
        @max_mode.setter
        def max_mode(self, max_mode):
            _libfgamma.f90wrap_type_routing_states__set__max_mode(self._handle, max_mode)
        
        @property
        def min_mode(self):
            """
            Element min_mode ftype=real  pytype=float32
            Defined at mod_routing_states.f90 line 15
            """
            return _libfgamma.f90wrap_type_routing_states__get__min_mode(self._handle)
        
        @min_mode.setter
        def min_mode(self, min_mode):
            _libfgamma.f90wrap_type_routing_states__set__min_mode(self._handle, min_mode)
        
        @property
        def max_spreading(self):
            """
            Element max_spreading ftype=real  pytype=float32
            Defined at mod_routing_states.f90 line 16
            """
            return _libfgamma.f90wrap_type_routing_states__get__max_spreading(self._handle)
        
        @max_spreading.setter
        def max_spreading(self, max_spreading):
            _libfgamma.f90wrap_type_routing_states__set__max_spreading(self._handle, \
                max_spreading)
        
        @property
        def window_shift(self):
            """
            Element window_shift ftype=real  pytype=float32
            Defined at mod_routing_states.f90 line 17
            """
            return _libfgamma.f90wrap_type_routing_states__get__window_shift(self._handle)
        
        @window_shift.setter
        def window_shift(self, window_shift):
            _libfgamma.f90wrap_type_routing_states__set__window_shift(self._handle, \
                window_shift)
        
        @property
        def scale_coef(self):
            """
            Element scale_coef ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 19
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__scale_coef(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            scale_coef = self._arrays.get(array_hash)
            if scale_coef is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if scale_coef.ctypes.data != array_handle:
                    scale_coef = None
            if scale_coef is None:
                try:
                    scale_coef = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__scale_coef)
                except TypeError:
                    scale_coef = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = scale_coef
            return scale_coef
        
        @scale_coef.setter
        def scale_coef(self, scale_coef):
            self.scale_coef[...] = scale_coef
        
        @property
        def param_normalisation(self):
            """
            Element param_normalisation ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 20
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__param_normalisation(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            param_normalisation = self._arrays.get(array_hash)
            if param_normalisation is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if param_normalisation.ctypes.data != array_handle:
                    param_normalisation = None
            if param_normalisation is None:
                try:
                    param_normalisation = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__param_normalisation)
                except TypeError:
                    param_normalisation = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = param_normalisation
            return param_normalisation
        
        @param_normalisation.setter
        def param_normalisation(self, param_normalisation):
            self.param_normalisation[...] = param_normalisation
        
        @property
        def quantile(self):
            """
            Element quantile ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 21
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__quantile(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            quantile = self._arrays.get(array_hash)
            if quantile is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if quantile.ctypes.data != array_handle:
                    quantile = None
            if quantile is None:
                try:
                    quantile = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__quantile)
                except TypeError:
                    quantile = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = quantile
            return quantile
        
        @quantile.setter
        def quantile(self, quantile):
            self.quantile[...] = quantile
        
        @property
        def tabulated_delay(self):
            """
            Element tabulated_delay ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 22
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__tabulated_delay(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            tabulated_delay = self._arrays.get(array_hash)
            if tabulated_delay is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if tabulated_delay.ctypes.data != array_handle:
                    tabulated_delay = None
            if tabulated_delay is None:
                try:
                    tabulated_delay = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__tabulated_delay)
                except TypeError:
                    tabulated_delay = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = tabulated_delay
            return tabulated_delay
        
        @tabulated_delay.setter
        def tabulated_delay(self, tabulated_delay):
            self.tabulated_delay[...] = tabulated_delay
        
        @property
        def tabulated_spreading(self):
            """
            Element tabulated_spreading ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 23
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__tabulated_spreading(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            tabulated_spreading = self._arrays.get(array_hash)
            if tabulated_spreading is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if tabulated_spreading.ctypes.data != array_handle:
                    tabulated_spreading = None
            if tabulated_spreading is None:
                try:
                    tabulated_spreading = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__tabulated_spreading)
                except TypeError:
                    tabulated_spreading = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = tabulated_spreading
            return tabulated_spreading
        
        @tabulated_spreading.setter
        def tabulated_spreading(self, tabulated_spreading):
            self.tabulated_spreading[...] = tabulated_spreading
        
        @property
        def tabulated_routing_coef(self):
            """
            Element tabulated_routing_coef ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 25
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__tabulated_routing_coef(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            tabulated_routing_coef = self._arrays.get(array_hash)
            if tabulated_routing_coef is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if tabulated_routing_coef.ctypes.data != array_handle:
                    tabulated_routing_coef = None
            if tabulated_routing_coef is None:
                try:
                    tabulated_routing_coef = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__tabulated_routing_coef)
                except TypeError:
                    tabulated_routing_coef = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = tabulated_routing_coef
            return tabulated_routing_coef
        
        @tabulated_routing_coef.setter
        def tabulated_routing_coef(self, tabulated_routing_coef):
            self.tabulated_routing_coef[...] = tabulated_routing_coef
        
        @property
        def states(self):
            """
            Element states ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 26
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__states(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            states = self._arrays.get(array_hash)
            if states is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if states.ctypes.data != array_handle:
                    states = None
            if states is None:
                try:
                    states = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__states)
                except TypeError:
                    states = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = states
            return states
        
        @states.setter
        def states(self, states):
            self.states[...] = states
        
        @property
        def remainder(self):
            """
            Element remainder ftype=real pytype=float array
            Defined at mod_routing_states.f90 line 27
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_states__array__remainder(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            remainder = self._arrays.get(array_hash)
            if remainder is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if remainder.ctypes.data != array_handle:
                    remainder = None
            if remainder is None:
                try:
                    remainder = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_states__array__remainder)
                except TypeError:
                    remainder = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = remainder
            return remainder
        
        @remainder.setter
        def remainder(self, remainder):
            self.remainder[...] = remainder
        
        def __str__(self):
            ret = ['<type_routing_states>{\n']
            ret.append('    window_length : ')
            ret.append(repr(self.window_length))
            ret.append(',\n    nb_mode : ')
            ret.append(repr(self.nb_mode))
            ret.append(',\n    nb_spreads : ')
            ret.append(repr(self.nb_spreads))
            ret.append(',\n    max_mode : ')
            ret.append(repr(self.max_mode))
            ret.append(',\n    min_mode : ')
            ret.append(repr(self.min_mode))
            ret.append(',\n    max_spreading : ')
            ret.append(repr(self.max_spreading))
            ret.append(',\n    window_shift : ')
            ret.append(repr(self.window_shift))
            ret.append(',\n    scale_coef : ')
            ret.append(repr(self.scale_coef))
            ret.append(',\n    param_normalisation : ')
            ret.append(repr(self.param_normalisation))
            ret.append(',\n    quantile : ')
            ret.append(repr(self.quantile))
            ret.append(',\n    tabulated_delay : ')
            ret.append(repr(self.tabulated_delay))
            ret.append(',\n    tabulated_spreading : ')
            ret.append(repr(self.tabulated_spreading))
            ret.append(',\n    tabulated_routing_coef : ')
            ret.append(repr(self.tabulated_routing_coef))
            ret.append(',\n    states : ')
            ret.append(repr(self.states))
            ret.append(',\n    remainder : ')
            ret.append(repr(self.remainder))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def routing_state_self_initialisation(self, routing_mesh, routing_parameter, \
        routing_states, interface_call=False):
        """
        routing_state_self_initialisation(self, routing_mesh, routing_parameter, \
            routing_states)
        Defined at mod_routing_states.f90 lines 30-83
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_parameter : Type_Routing_Parameter
        routing_states : Type_Routing_States
        """
        _libfgamma.f90wrap_mod_gamma_routing_states__routing_state_self_initiaf922(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, \
            routing_parameter=routing_parameter._handle, \
            routing_states=routing_states._handle)
    
    @staticmethod
    def routing_states_reset(self, interface_call=False):
        """
        routing_states_reset(self)
        Defined at mod_routing_states.f90 lines 85-101
        
        Parameters
        ----------
        routing_states : Type_Routing_States
        """
        _libfgamma.f90wrap_mod_gamma_routing_states__routing_states_reset(routing_states=self._handle)
    
    @staticmethod
    def routing_states_clear(self, interface_call=False):
        """
        routing_states_clear(self)
        Defined at mod_routing_states.f90 lines 103-118
        
        Parameters
        ----------
        routing_states : Type_Routing_States
        """
        _libfgamma.f90wrap_mod_gamma_routing_states__routing_states_clear(routing_states=self._handle)
    
    @staticmethod
    def routing_states_copy(self, interface_call=False):
        """
        object_copy = routing_states_copy(self)
        Defined at mod_routing_states.f90 lines 120-123
        
        Parameters
        ----------
        routing_states : Type_Routing_States
        
        Returns
        -------
        object_copy : Type_Routing_States
        """
        object_copy = \
            _libfgamma.f90wrap_mod_gamma_routing_states__routing_states_copy(routing_states=self._handle)
        object_copy = \
            f90wrap.runtime.lookup_class("libfgamma.type_routing_states").from_handle(object_copy, \
            alloc=True)
        object_copy._setup_finalizer()
        return object_copy
    
    _dt_array_initialisers = []
    
    

mod_gamma_routing_states = Mod_Gamma_Routing_States()

class Mod_Gamma_Routing_Results(f90wrap.runtime.FortranModule):
    """
    Module mod_gamma_routing_results
    Defined at mod_routing_results.f90 lines 8-120
    """
    @f90wrap.runtime.register_class("libfgamma.type_routing_results")
    class type_routing_results(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=type_routing_results)
        Defined at mod_routing_results.f90 lines 10-21
        """
        def __init__(self, handle=None):
            """
            Automatically generated constructor for type_routing_results
            
            self = Type_Routing_Results()
            Defined at mod_routing_results.f90 lines 10-21
            
            Returns
            -------
            this : Type_Routing_Results
                Object to be constructed
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            if isinstance(handle, numpy.ndarray) and handle.ndim == 1 and handle.dtype.num \
                == 5:
                self._handle = handle
                self._alloc = True
            else:
                result = \
                    _libfgamma.f90wrap_mod_gamma_routing_results__type_routing_results_inieb68()
                self._handle = result[0] if isinstance(result, tuple) else result
                self._alloc = True
            self._setup_finalizer()
        
        def _setup_finalizer(self):
            """Set up weak reference destructor to prevent Fortran memory leaks."""
            if self._alloc:
                destructor = getattr(_libfgamma, \
                    "f90wrap_mod_gamma_routing_results__type_routing_results_fin4f66")
                self._finalizer = weakref.finalize(self, destructor, self._handle)
        
        @property
        def discharges(self):
            """
            Element discharges ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 11
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__discharges(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            discharges = self._arrays.get(array_hash)
            if discharges is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if discharges.ctypes.data != array_handle:
                    discharges = None
            if discharges is None:
                try:
                    discharges = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__discharges)
                except TypeError:
                    discharges = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = discharges
            return discharges
        
        @discharges.setter
        def discharges(self, discharges):
            self.discharges[...] = discharges
        
        @property
        def velocities(self):
            """
            Element velocities ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 12
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__velocities(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            velocities = self._arrays.get(array_hash)
            if velocities is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if velocities.ctypes.data != array_handle:
                    velocities = None
            if velocities is None:
                try:
                    velocities = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__velocities)
                except TypeError:
                    velocities = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = velocities
            return velocities
        
        @velocities.setter
        def velocities(self, velocities):
            self.velocities[...] = velocities
        
        @property
        def costs(self):
            """
            Element costs ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 13
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__costs(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            costs = self._arrays.get(array_hash)
            if costs is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if costs.ctypes.data != array_handle:
                    costs = None
            if costs is None:
                try:
                    costs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__costs)
                except TypeError:
                    costs = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = costs
            return costs
        
        @costs.setter
        def costs(self, costs):
            self.costs[...] = costs
        
        @property
        def allgrads_hydraulic_coef(self):
            """
            Element allgrads_hydraulic_coef ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 14
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__allgrads_hydraulic_coef(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            allgrads_hydraulic_coef = self._arrays.get(array_hash)
            if allgrads_hydraulic_coef is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if allgrads_hydraulic_coef.ctypes.data != array_handle:
                    allgrads_hydraulic_coef = None
            if allgrads_hydraulic_coef is None:
                try:
                    allgrads_hydraulic_coef = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__allgrads_hydraulic_coef)
                except TypeError:
                    allgrads_hydraulic_coef = f90wrap.runtime.direct_c_array(array_type, \
                        array_shape, array_handle)
                self._arrays[array_hash] = allgrads_hydraulic_coef
            return allgrads_hydraulic_coef
        
        @allgrads_hydraulic_coef.setter
        def allgrads_hydraulic_coef(self, allgrads_hydraulic_coef):
            self.allgrads_hydraulic_coef[...] = allgrads_hydraulic_coef
        
        @property
        def allgrads_spreading(self):
            """
            Element allgrads_spreading ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 15
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__allgrads_spreading(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            allgrads_spreading = self._arrays.get(array_hash)
            if allgrads_spreading is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if allgrads_spreading.ctypes.data != array_handle:
                    allgrads_spreading = None
            if allgrads_spreading is None:
                try:
                    allgrads_spreading = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__allgrads_spreading)
                except TypeError:
                    allgrads_spreading = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = allgrads_spreading
            return allgrads_spreading
        
        @allgrads_spreading.setter
        def allgrads_spreading(self, allgrads_spreading):
            self.allgrads_spreading[...] = allgrads_spreading
        
        @property
        def gradients_hydraulic_coef(self):
            """
            Element gradients_hydraulic_coef ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 16
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__gradients_hydraulic_coef(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            gradients_hydraulic_coef = self._arrays.get(array_hash)
            if gradients_hydraulic_coef is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if gradients_hydraulic_coef.ctypes.data != array_handle:
                    gradients_hydraulic_coef = None
            if gradients_hydraulic_coef is None:
                try:
                    gradients_hydraulic_coef = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__gradients_hydraulic_coef)
                except TypeError:
                    gradients_hydraulic_coef = f90wrap.runtime.direct_c_array(array_type, \
                        array_shape, array_handle)
                self._arrays[array_hash] = gradients_hydraulic_coef
            return gradients_hydraulic_coef
        
        @gradients_hydraulic_coef.setter
        def gradients_hydraulic_coef(self, gradients_hydraulic_coef):
            self.gradients_hydraulic_coef[...] = gradients_hydraulic_coef
        
        @property
        def gradients_spreading(self):
            """
            Element gradients_spreading ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 17
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__gradients_spreading(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            gradients_spreading = self._arrays.get(array_hash)
            if gradients_spreading is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if gradients_spreading.ctypes.data != array_handle:
                    gradients_spreading = None
            if gradients_spreading is None:
                try:
                    gradients_spreading = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__gradients_spreading)
                except TypeError:
                    gradients_spreading = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = gradients_spreading
            return gradients_spreading
        
        @gradients_spreading.setter
        def gradients_spreading(self, gradients_spreading):
            self.gradients_spreading[...] = gradients_spreading
        
        @property
        def initial_hydraulic_coef(self):
            """
            Element initial_hydraulic_coef ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 18
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__initial_hydraulic_coef(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            initial_hydraulic_coef = self._arrays.get(array_hash)
            if initial_hydraulic_coef is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if initial_hydraulic_coef.ctypes.data != array_handle:
                    initial_hydraulic_coef = None
            if initial_hydraulic_coef is None:
                try:
                    initial_hydraulic_coef = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__initial_hydraulic_coef)
                except TypeError:
                    initial_hydraulic_coef = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = initial_hydraulic_coef
            return initial_hydraulic_coef
        
        @initial_hydraulic_coef.setter
        def initial_hydraulic_coef(self, initial_hydraulic_coef):
            self.initial_hydraulic_coef[...] = initial_hydraulic_coef
        
        @property
        def initial_spreading(self):
            """
            Element initial_spreading ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 19
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__initial_spreading(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            initial_spreading = self._arrays.get(array_hash)
            if initial_spreading is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if initial_spreading.ctypes.data != array_handle:
                    initial_spreading = None
            if initial_spreading is None:
                try:
                    initial_spreading = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__initial_spreading)
                except TypeError:
                    initial_spreading = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = initial_spreading
            return initial_spreading
        
        @initial_spreading.setter
        def initial_spreading(self, initial_spreading):
            self.initial_spreading[...] = initial_spreading
        
        @property
        def final_hydraulic_coef(self):
            """
            Element final_hydraulic_coef ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 20
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__final_hydraulic_coef(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            final_hydraulic_coef = self._arrays.get(array_hash)
            if final_hydraulic_coef is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if final_hydraulic_coef.ctypes.data != array_handle:
                    final_hydraulic_coef = None
            if final_hydraulic_coef is None:
                try:
                    final_hydraulic_coef = \
                        f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__final_hydraulic_coef)
                except TypeError:
                    final_hydraulic_coef = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = final_hydraulic_coef
            return final_hydraulic_coef
        
        @final_hydraulic_coef.setter
        def final_hydraulic_coef(self, final_hydraulic_coef):
            self.final_hydraulic_coef[...] = final_hydraulic_coef
        
        @property
        def final_spreading(self):
            """
            Element final_spreading ftype=real pytype=float array
            Defined at mod_routing_results.f90 line 21
            """
            array_ndim, array_type, array_shape, array_handle = \
                _libfgamma.f90wrap_type_routing_results__array__final_spreading(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            final_spreading = self._arrays.get(array_hash)
            if final_spreading is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if final_spreading.ctypes.data != array_handle:
                    final_spreading = None
            if final_spreading is None:
                try:
                    final_spreading = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _libfgamma.f90wrap_type_routing_results__array__final_spreading)
                except TypeError:
                    final_spreading = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = final_spreading
            return final_spreading
        
        @final_spreading.setter
        def final_spreading(self, final_spreading):
            self.final_spreading[...] = final_spreading
        
        def __str__(self):
            ret = ['<type_routing_results>{\n']
            ret.append('    discharges : ')
            ret.append(repr(self.discharges))
            ret.append(',\n    velocities : ')
            ret.append(repr(self.velocities))
            ret.append(',\n    costs : ')
            ret.append(repr(self.costs))
            ret.append(',\n    allgrads_hydraulic_coef : ')
            ret.append(repr(self.allgrads_hydraulic_coef))
            ret.append(',\n    allgrads_spreading : ')
            ret.append(repr(self.allgrads_spreading))
            ret.append(',\n    gradients_hydraulic_coef : ')
            ret.append(repr(self.gradients_hydraulic_coef))
            ret.append(',\n    gradients_spreading : ')
            ret.append(repr(self.gradients_spreading))
            ret.append(',\n    initial_hydraulic_coef : ')
            ret.append(repr(self.initial_hydraulic_coef))
            ret.append(',\n    initial_spreading : ')
            ret.append(repr(self.initial_spreading))
            ret.append(',\n    final_hydraulic_coef : ')
            ret.append(repr(self.final_hydraulic_coef))
            ret.append(',\n    final_spreading : ')
            ret.append(repr(self.final_spreading))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def routing_results_self_initialisation(self, routing_mesh, routing_results, \
        interface_call=False):
        """
        routing_results_self_initialisation(self, routing_mesh, routing_results)
        Defined at mod_routing_results.f90 lines 24-74
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_results : Type_Routing_Results
        """
        _libfgamma.f90wrap_mod_gamma_routing_results__routing_results_self_iniec9c(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, routing_results=routing_results._handle)
    
    @staticmethod
    def routing_results_clear(self, interface_call=False):
        """
        routing_results_clear(self)
        Defined at mod_routing_results.f90 lines 76-90
        
        Parameters
        ----------
        routing_results : Type_Routing_Results
        """
        _libfgamma.f90wrap_mod_gamma_routing_results__routing_results_clear(routing_results=self._handle)
    
    @staticmethod
    def routing_results_reset(self, interface_call=False):
        """
        routing_results_reset(self)
        Defined at mod_routing_results.f90 lines 92-115
        
        Parameters
        ----------
        routing_results : Type_Routing_Results
        """
        _libfgamma.f90wrap_mod_gamma_routing_results__routing_results_reset(routing_results=self._handle)
    
    @staticmethod
    def routing_results_copy(self, interface_call=False):
        """
        object_copy = routing_results_copy(self)
        Defined at mod_routing_results.f90 lines 117-120
        
        Parameters
        ----------
        routing_results : Type_Routing_Results
        
        Returns
        -------
        object_copy : Type_Routing_Results
        """
        object_copy = \
            _libfgamma.f90wrap_mod_gamma_routing_results__routing_results_copy(routing_results=self._handle)
        object_copy = \
            f90wrap.runtime.lookup_class("libfgamma.type_routing_results").from_handle(object_copy, \
            alloc=True)
        object_copy._setup_finalizer()
        return object_copy
    
    _dt_array_initialisers = []
    
    

mod_gamma_routing_results = Mod_Gamma_Routing_Results()

class Mod_Gamma_Interface(f90wrap.runtime.FortranModule):
    """
    Module mod_gamma_interface
    Defined at mod_gamma_interface.f90 lines 8-355
    """
    @staticmethod
    def auto_compute_boundaries(self, routing_mesh, observed_discharges, \
        interface_call=False):
        """
        auto_compute_boundaries(self, routing_mesh, observed_discharges)
        Defined at mod_gamma_interface.f90 lines 18-87
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        observed_discharges : float array
        """
        _libfgamma.f90wrap_mod_gamma_interface__auto_compute_boundaries(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, observed_discharges=observed_discharges)
    
    @staticmethod
    def routing_gamma_change_parameters(self, routing_states, routing_setup, \
        routing_mesh, hydraulics_coefficient=None, spreading=None, \
        interface_call=False):
        """
        routing_gamma_change_parameters(self, routing_states, routing_setup, \
            routing_mesh[, hydraulics_coefficient, spreading])
        Defined at mod_gamma_interface.f90 lines 89-127
        
        Parameters
        ----------
        routing_parameter : Type_Routing_Parameter
        routing_states : Type_Routing_States
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        hydraulics_coefficient : float32
        spreading : float32
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_change_parameters(routing_parameter=self._handle, \
            routing_states=routing_states._handle, routing_setup=routing_setup._handle, \
            routing_mesh=routing_mesh._handle, \
            hydraulics_coefficient=hydraulics_coefficient, spreading=spreading)
    
    @staticmethod
    def routing_gamma_precomputation(self, routing_mesh, routing_states, \
        interface_call=False):
        """
        routing_gamma_precomputation(self, routing_mesh, routing_states)
        Defined at mod_gamma_interface.f90 lines 129-148
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_states : Type_Routing_States
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_precomputation(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, routing_states=routing_states._handle)
    
    @staticmethod
    def routing_states_update(self, routing_setup, routing_mesh, routing_states, \
        interface_call=False):
        """
        routing_states_update(self, routing_setup, routing_mesh, routing_states)
        Defined at mod_gamma_interface.f90 lines 150-172
        
        Parameters
        ----------
        routing_parameter : Type_Routing_Parameter
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_states : Type_Routing_States
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_states_update(routing_parameter=self._handle, \
            routing_setup=routing_setup._handle, routing_mesh=routing_mesh._handle, \
            routing_states=routing_states._handle)
    
    @staticmethod
    def routing_gamma_run(self, routing_mesh, routing_parameter, inflows, \
        routing_states, routing_results, interface_call=False):
        """
        routing_gamma_run(self, routing_mesh, routing_parameter, inflows, \
            routing_states, routing_results)
        Defined at mod_gamma_interface.f90 lines 174-202
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_parameter : Type_Routing_Parameter
        inflows : float array
        routing_states : Type_Routing_States
        routing_results : Type_Routing_Results
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_run(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, \
            routing_parameter=routing_parameter._handle, inflows=inflows, \
            routing_states=routing_states._handle, \
            routing_results=routing_results._handle)
    
    @staticmethod
    def routing_gamma_control(self, routing_mesh, routing_parameter, inflows, \
        observations, routing_states, routing_results, interface_call=False):
        """
        routing_gamma_control(self, routing_mesh, routing_parameter, inflows, \
            observations, routing_states, routing_results)
        Defined at mod_gamma_interface.f90 lines 204-238
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_parameter : Type_Routing_Parameter
        inflows : float array
        observations : float array
        routing_states : Type_Routing_States
        routing_results : Type_Routing_Results
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_control(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, \
            routing_parameter=routing_parameter._handle, inflows=inflows, \
            observations=observations, routing_states=routing_states._handle, \
            routing_results=routing_results._handle)
    
    @staticmethod
    def routing_gamma_forward_adjoint_b(self, routing_mesh, routing_parameter, \
        inflows, observations, routing_states, routing_results, gradients, \
        interface_call=False):
        """
        routing_gamma_forward_adjoint_b(self, routing_mesh, routing_parameter, inflows, \
            observations, routing_states, routing_results, gradients)
        Defined at mod_gamma_interface.f90 lines 240-265
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_parameter : Type_Routing_Parameter
        inflows : float array
        observations : float array
        routing_states : Type_Routing_States
        routing_results : Type_Routing_Results
        gradients : float array
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_forward_adjoint_b(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, \
            routing_parameter=routing_parameter._handle, inflows=inflows, \
            observations=observations, routing_states=routing_states._handle, \
            routing_results=routing_results._handle, gradients=gradients)
    
    @staticmethod
    def routing_gamma_forward_adjoint_b0(self, routing_mesh, routing_parameter, \
        inflows, observations, routing_states, routing_results, gradients, \
        interface_call=False):
        """
        routing_gamma_forward_adjoint_b0(self, routing_mesh, routing_parameter, inflows, \
            observations, routing_states, routing_results, gradients)
        Defined at mod_gamma_interface.f90 lines 267-291
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_parameter : Type_Routing_Parameter
        inflows : float array
        observations : float array
        routing_states : Type_Routing_States
        routing_results : Type_Routing_Results
        gradients : float array
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_forward_adjoint_b0(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, \
            routing_parameter=routing_parameter._handle, inflows=inflows, \
            observations=observations, routing_states=routing_states._handle, \
            routing_results=routing_results._handle, gradients=gradients)
    
    @staticmethod
    def routing_gamma_cost_function(self, routing_mesh, routing_parameter, \
        observations, qnetwork, routing_results, interface_call=False):
        """
        routing_gamma_cost_function(self, routing_mesh, routing_parameter, observations, \
            qnetwork, routing_results)
        Defined at mod_gamma_interface.f90 lines 293-321
        
        Parameters
        ----------
        routing_setup : Type_Routing_Setup
        routing_mesh : Type_Routing_Mesh
        routing_parameter : Type_Routing_Parameter
        observations : float array
        qnetwork : float array
        routing_results : Type_Routing_Results
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_cost_function(routing_setup=self._handle, \
            routing_mesh=routing_mesh._handle, \
            routing_parameter=routing_parameter._handle, observations=observations, \
            qnetwork=qnetwork, routing_results=routing_results._handle)
    
    @staticmethod
    def routing_gamma_linear_interpolation(delay, index_varying_dx, routing_states, \
        gamma_coefficient, interface_call=False):
        """
        routing_gamma_linear_interpolation(delay, index_varying_dx, routing_states, \
            gamma_coefficient)
        Defined at mod_gamma_interface.f90 lines 324-330
        
        Parameters
        ----------
        delay : float32
        index_varying_dx : int32
        routing_states : Type_Routing_States
        gamma_coefficient : float array
        """
        _libfgamma.f90wrap_mod_gamma_interface__routing_gamma_linear_interpolation(delay=delay, \
            index_varying_dx=index_varying_dx, routing_states=routing_states._handle, \
            gamma_coefficient=gamma_coefficient)
    
    @staticmethod
    def interface_compute_gamma_coefficient(scale, mode, quantile, window_shift, \
        density_function, gamma_coefficient, interface_call=False):
        """
        interface_compute_gamma_coefficient(scale, mode, quantile, window_shift, \
            density_function, gamma_coefficient)
        Defined at mod_gamma_interface.f90 lines 334-342
        
        Parameters
        ----------
        scale : float32
        mode : float32
        quantile : float array
        window_shift : float32
        density_function : str
        gamma_coefficient : float array
        """
        _libfgamma.f90wrap_mod_gamma_interface__interface_compute_gamma_coeffiac06(scale=scale, \
            mode=mode, quantile=quantile, window_shift=window_shift, \
            density_function=density_function, gamma_coefficient=gamma_coefficient)
    
    @staticmethod
    def interface_compute_routing_coefficient(scale, mode, quantile, window_shift, \
        density_function, gamma_coefficient, interface_call=False):
        """
        interface_compute_routing_coefficient(scale, mode, quantile, window_shift, \
            density_function, gamma_coefficient)
        Defined at mod_gamma_interface.f90 lines 346-354
        
        Parameters
        ----------
        scale : float32
        mode : float32
        quantile : float array
        window_shift : float32
        density_function : str
        gamma_coefficient : float array
        """
        _libfgamma.f90wrap_mod_gamma_interface__interface_compute_routing_coefbd47(scale=scale, \
            mode=mode, quantile=quantile, window_shift=window_shift, \
            density_function=density_function, gamma_coefficient=gamma_coefficient)
    
    _dt_array_initialisers = []
    

mod_gamma_interface = Mod_Gamma_Interface()

