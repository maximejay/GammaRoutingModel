!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


module mod_gamma_routing_memory
    
    implicit none
    
    type type_routing_memory
        real,dimension(:,:), allocatable :: states !state of the system at t0
        real,dimension(:,:), allocatable :: remainder!remainder for the routing scheme
    end type type_routing_memory
    
    contains
    
    subroutine routing_memory_self_initialisation(routing_mesh,routing_states,routing_memory)
        
        ! Notes
        ! -----
        ! **routing_state_self_initialisation(routing_mesh,routing_states,routing_memory)** :
        !
        ! - Initialise the routing_memory derived type, allocate all components and precompute some variables
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_parameter``                   routing_parameter Derived Type (in)
        ! ``routing_states``                      routing_states Derived Type (inout)
        ! ``routing_memory``                      routing_memory Derived Type (inout)
        ! =============================           ===================================
        
        use mod_gamma_routing_mesh
        use mod_gamma_routing_states
        
        implicit none
        
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(in) :: routing_states
        type(type_routing_memory), intent(inout) :: routing_memory
        
        if (allocated(routing_memory%remainder)) then
            deallocate(routing_memory%remainder)
        end if
        if (allocated(routing_memory%states)) then
            deallocate(routing_memory%states)
        end if
        
        allocate(routing_memory%remainder(int(maxval(routing_states%window_length)),routing_mesh%nb_nodes))
        allocate(routing_memory%states(int(maxval(routing_states%window_length)),routing_mesh%nb_nodes))
        
        !default value
        routing_memory%remainder=0.
        routing_memory%states=0.
        
    end subroutine routing_memory_self_initialisation
    
    
    subroutine routing_memory_reset(routing_memory)
        
        ! Notes
        ! -----
        ! **routing_memory_reset(routing_memory)** :
        !
        ! - Reset the derived type routing_memory, set to zeros the states and the remainder components
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_memory``                      routing_memory Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_memory), intent(inout) :: routing_memory
        
        !default value
        routing_memory%remainder=0.
        routing_memory%states=0.
        
    end subroutine routing_memory_reset
    
    subroutine routing_memory_clear(routing_memory)
        
        ! Notes
        ! -----
        ! **routing_memory_clear(routing_memory)** :
        !
        ! - Clear the derived type routing_memory
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_memory``                      routing_memory Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_memory), intent(inout) :: routing_memory
        
        type(type_routing_memory) :: routing_memory_new
        
        routing_memory=routing_memory_new
    endsubroutine routing_memory_clear
    
    
    subroutine routing_memory_copy(routing_memory, object_copy)
    
        type(type_routing_memory), intent(in) :: routing_memory
        type(type_routing_memory), intent(out) :: object_copy
        
        object_copy=routing_memory
    
    end subroutine routing_memory_copy
    
end module mod_gamma_routing_memory
