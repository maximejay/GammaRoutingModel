!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


module mod_gamma_routing_parameters
    
    implicit none
    
    type type_routing_parameter
        real, dimension(:), allocatable :: hydraulics_coefficient 
        real, dimension(:), allocatable :: spreading ! damping coefficient in seconds (s/m): spreading of the Gamma law
    end type type_routing_parameter
    
    contains
    
    
    subroutine routing_parameter_self_initialisation(routing_parameter,routing_setup,routing_mesh,&
    &hydraulics_coefficient,spreading)
        
        ! Notes
        ! -----
        ! **routing_parameter_self_initialisation(routing_parameter,routing_setup,routing_mesh,hydraulics_coefficient,spreading)** :
        !
        ! - Initialise the routing_parameter derived type with user values, allocate all components and set user values for all nodes
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_parameter``                   routing_parameter Derived Type (inout)
        ! ``routing_setup``                       routing_setup Derived Type (inout)
        ! ``routing_mesh``                        Routing_mesh Derived Type (inout)
        ! ``hydraulics_coefficient=1.``           Value of the hydraulic coefficient (optional)
        ! ``spreading=dt./dx.``                   Value of the spreading coefficient, default is set to dt (in second/m)(optional)
        ! =============================           ===================================
        
        
        use mod_gamma_routing_setup
        use mod_gamma_routing_mesh
        
        implicit none
        
        type(type_routing_parameter), intent(inout) :: routing_parameter
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        real,optional,intent(in) :: hydraulics_coefficient
        real,optional,intent(in) :: spreading
        
        integer :: i
        
        if (.not.allocated(routing_parameter%hydraulics_coefficient)) then
            allocate(routing_parameter%hydraulics_coefficient(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_parameter%spreading)) then
            allocate(routing_parameter%spreading(routing_mesh%nb_nodes))
        end if
        
        if (present(hydraulics_coefficient) .AND. hydraulics_coefficient>0) then
            
            routing_parameter%hydraulics_coefficient=hydraulics_coefficient
        else
            routing_parameter%hydraulics_coefficient=1.0 !default value
        end if
        
        if (present(spreading) .AND. spreading>0) then
            
            routing_parameter%spreading=spreading ! given in s/m 
        else
            
            do i=1,routing_mesh%nb_nodes
                routing_parameter%spreading(i)=routing_setup%dt/routing_mesh%dx(i)  !default value
            end do
            
        end if
        
        !reading parameter
        
        !setting parameter
        
    end subroutine routing_parameter_self_initialisation
    
    
    
    subroutine routing_parameter_clear(routing_parameter)
        
        ! Notes
        ! -----
        ! **routing_parameter_clear(routing_parameter)** :
        !
        ! - Clear the derived type routing_parameter
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_parameter``                   routing_parameter Derived Type (inout)
        ! =============================           ===================================
        
        type(type_routing_parameter), intent(inout) :: routing_parameter
        
        type(type_routing_parameter) :: routing_parameter_new
        
        routing_parameter=routing_parameter_new
    endsubroutine routing_parameter_clear
    
    subroutine routing_parameter_copy(routing_parameter,object_copy)
    
        type(type_routing_parameter), intent(in) :: routing_parameter
        type(type_routing_parameter), intent(out) :: object_copy
        
        object_copy=routing_parameter
    
    end subroutine routing_parameter_copy
    
end module mod_gamma_routing_parameters
