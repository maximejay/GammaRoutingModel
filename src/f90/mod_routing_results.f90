!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


module mod_gamma_routing_results
    
    implicit none
    
    type type_routing_results
        real, dimension(:,:), allocatable :: discharges 
        real, dimension(:,:), allocatable :: velocities
        real, dimension(3) :: costs
        real, dimension(:,:), allocatable :: allgrads_hydraulic_coef
        real, dimension(:,:), allocatable :: allgrads_spreading
        real, dimension(:), allocatable :: gradients_hydraulic_coef
        real, dimension(:), allocatable :: gradients_spreading
        real, dimension(:), allocatable :: initial_hydraulic_coef
        real, dimension(:), allocatable :: initial_spreading
        real, dimension(:), allocatable :: final_hydraulic_coef
        real, dimension(:), allocatable :: final_spreading
    end type type_routing_results
    
    contains
    
    
    subroutine routing_results_self_initialisation(routing_setup,routing_mesh,routing_results)
        
        ! Notes
        ! -----
        ! **routing_results_self_initialisation(routing_setup,routing_mesh,routing_results)** :
        !
        ! - Initialise the routing_results derived type, allocate all components
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_results``                     routing_results Derived Type (in)
        ! =============================           ===================================
        
        
        use mod_gamma_routing_setup
        use mod_gamma_routing_mesh
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_results), intent(inout) :: routing_results
        
        if (.not.allocated(routing_results%discharges)) then
            allocate(routing_results%discharges(routing_setup%npdt,routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%velocities)) then
            allocate(routing_results%velocities(routing_setup%npdt,routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%allgrads_hydraulic_coef)) then
            allocate(routing_results%allgrads_hydraulic_coef(routing_setup%iter_max,routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%allgrads_spreading)) then
            allocate(routing_results%allgrads_spreading(routing_setup%iter_max,routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%gradients_hydraulic_coef)) then
            allocate(routing_results%gradients_hydraulic_coef(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%gradients_spreading)) then
            allocate(routing_results%gradients_spreading(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%initial_hydraulic_coef)) then
            allocate(routing_results%initial_hydraulic_coef(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%initial_spreading)) then
            allocate(routing_results%initial_spreading(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%final_hydraulic_coef)) then
            allocate(routing_results%final_hydraulic_coef(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_results%final_spreading)) then
            allocate(routing_results%final_spreading(routing_mesh%nb_nodes))
        end if
        
        call routing_results_reset(routing_results)
        
    end subroutine routing_results_self_initialisation
    
    
    
    subroutine routing_results_clear(routing_results)
        
        ! Notes
        ! -----
        ! **routing_results_clear(routing_results)** :
        !
        ! - Clear the derived type routing_results
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_results``                     routing_results Derived Type (inout)
        ! =============================           ===================================
        
        type(type_routing_results), intent(inout) :: routing_results
        
        type(type_routing_results) :: new_routing_results
        
        routing_results=new_routing_results
        
    endsubroutine routing_results_clear
    
    
    subroutine routing_results_reset(routing_results)
        
        ! Notes
        ! -----
        ! **routing_results_reset(routing_results)** :
        !
        ! - Reset the derived type routing_results
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_results``                     routing_results Derived Type (inout)
        ! =============================           ===================================
        
        type(type_routing_results), intent(inout) :: routing_results
        
        routing_results%discharges=0. 
        routing_results%velocities=0.
        routing_results%costs=0.
        routing_results%allgrads_hydraulic_coef=0.
        routing_results%allgrads_spreading=0.
        routing_results%gradients_hydraulic_coef=0.
        routing_results%gradients_spreading=0.
        routing_results%initial_hydraulic_coef=0.
        routing_results%initial_spreading=0.
        routing_results%final_hydraulic_coef=0.
        routing_results%final_spreading=0.
        
    endsubroutine routing_results_reset
    
    
    subroutine routing_results_copy(routing_results,object_copy)
    
        type(type_routing_results), intent(in) :: routing_results
        type(type_routing_results), intent(out) :: object_copy
        
        object_copy=routing_results
    
    end subroutine routing_results_copy
    
end module mod_gamma_routing_results
