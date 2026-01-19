!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


module mod_gamma_routing_states
    
    implicit none
    
    type type_routing_states
        integer,dimension(:),allocatable :: window_length ! length of the delay windows (model memory)
        integer :: nb_mode !number of discretization mode
        integer :: nb_spreads!number of discretization spreads
        real :: max_mode !the highest mode i.e for vmin and dx max
        real :: min_mode !the lowest mode i.e for vmax and dx min
        real :: max_spreading !the highest spreading values in s/m
        real :: window_shift !shift so that the peak of the PDF is located at x=dmin/vmax
        !real :: max_scale !uniform scale coefficient computed thanks to the max_spread coeff
        real,dimension(:),allocatable :: scale_coef !non uniform scale coefficient computed thanks to the max_spread coeff
        real, dimension(2) :: param_normalisation !Array of factor to normaize the model parameters (hydraulic_coeff, spreading)
        real, dimension(:), allocatable :: quantile !quantiles series to compute the Gamma pdf/cdf 
        real, dimension(:), allocatable :: tabulated_delay !tabulated delay (or mode) to locate the Gamma pdf
        real, dimension(:), allocatable :: tabulated_spreading !tabulated spreading to spread the Gamma pdf
        !real,dimension(:,:,:), allocatable :: tabulated_routing_coef !tabulated routing coefficient for the unit hydrogram
        real,dimension(:,:,:,:), allocatable :: tabulated_routing_coef !tabulated routing coefficient for the unit hydrogram
        real,dimension(:,:), allocatable :: states !state of the system at t0
        real,dimension(:,:),allocatable :: remainder!remainder for the routing scheme
    end type type_routing_states
    

    contains
    
    subroutine routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
        
        ! Notes
        ! -----
        ! **routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)** :
        !
        ! - Initialise the routing_states derived type, allocate all components and precompute some variables
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_parameter``                   routing_parameter Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! =============================           ===================================
        
        use mod_gamma_routing_setup
        use mod_gamma_routing_mesh
        use mod_gamma_routing_parameters
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        type(type_routing_states), intent(inout) :: routing_states
        
        integer :: i
        real :: max_spreading
        
        if (allocated(routing_states%window_length)) then
            deallocate(routing_states%window_length)
        end if
        allocate(routing_states%window_length(routing_mesh%nb_nodes))
        
        if (allocated(routing_states%scale_coef)) then
            deallocate(routing_states%scale_coef)
        end if
        allocate(routing_states%scale_coef(size(routing_mesh%varying_dx)))
        
        
        if (routing_setup%varying_spread>0) then
            routing_states%max_spreading=routing_setup%spreading_boundaries(2)
            routing_states%nb_spreads=int(routing_states%max_spreading/routing_setup%spreading_discretization_step)+1
        else
            routing_states%nb_spreads=1
!~             max_spreading=0.0
!~             max_spreading=maxval(routing_parameter%spreading)
!~             do i=1,routing_mesh%nb_nodes
!~                 max_spreading=max(max_spreading,routing_parameter%spreading(i))
!~             end do
!~             routing_states%max_spreading=max_spreading
            routing_states%max_spreading=maxval(routing_parameter%spreading)
!~             routing_states%max_spreading=routing_setup%spreading_boundaries(2)
        end if
        
        
        routing_states%max_mode=maxval(routing_mesh%dx)/routing_setup%vmin/routing_setup%dt
        routing_states%min_mode=minval(routing_mesh%dx)/routing_setup%vmax/routing_setup%dt
        routing_states%nb_mode=int((routing_states%max_mode+1.0)/routing_setup%mode_discretization_step)+1
        !routing_states%nb_mode=int((routing_states%max_mode-routing_states%min_mode+1.0)/routing_setup%mode_discretization_step)+1
        
        !Condition à respecter dx/dt >= 1 => courant pour v=1m/s
        
        if (routing_states%max_mode>routing_setup%npdt) then
            write(*,*) "Warning : npdt lower than than the maximum model delay:",routing_setup%npdt,"<"&
            &,routing_states%max_mode, "this could provoke a calibration failure or weird outputs."
        end if
        
        routing_states%param_normalisation=1.0
        
    end subroutine routing_state_self_initialisation
    
    
    subroutine routing_states_reset(routing_states)
        
        ! Notes
        ! -----
        ! **routing_states_reset(routing_states)** :
        !
        ! - Reset the derived type routing_states, set to zeros the states and the remainder components
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_states``                      routing_states Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_states), intent(inout) :: routing_states
        
        !default value
        routing_states%remainder=0.
        routing_states%states=0.
        
    end subroutine routing_states_reset
    
    subroutine routing_states_clear(routing_states)
        
        ! Notes
        ! -----
        ! **routing_states_clear(routing_states)** :
        !
        ! - Clear the derived type routing_states
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_states``                      routing_states Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_states), intent(inout) :: routing_states
        
        type(type_routing_states) :: routing_states_new
        
        routing_states=routing_states_new
    endsubroutine routing_states_clear
    
    
    subroutine routing_states_copy(routing_states, object_copy)
    
        type(type_routing_states), intent(in) :: routing_states
        type(type_routing_states), intent(out) :: object_copy
        
        object_copy=routing_states
    
    end subroutine routing_states_copy
    
end module mod_gamma_routing_states
