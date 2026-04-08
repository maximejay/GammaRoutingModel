!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


module mod_gamma_routing_parameters
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    
    implicit none
    
    type type_routing_parameter
        real, dimension(:), allocatable :: hc 
        real, dimension(:), allocatable :: sc ! damping coefficient in seconds (s/m): spreading of the Gamma law
        real, dimension(:), allocatable :: hc_n 
        real, dimension(:), allocatable :: sc_n ! damping coefficient in seconds (s/m): spreading of the Gamma law
        integer :: normalized
    end type type_routing_parameter
    
    contains
    
    
    subroutine routing_parameter_self_initialisation(routing_parameter,routing_setup,routing_mesh,&
    &hc,sc)
        
        ! Notes
        ! -----
        ! **routing_parameter_self_initialisation(routing_parameter,routing_setup,routing_mesh,hc,sc)** :
        !
        ! - Initialise the routing_parameter derived type with user values, allocate all components and set user values for all nodes
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_parameter``                   routing_parameter Derived Type (inout)
        ! ``routing_setup``                       routing_setup Derived Type (inout)
        ! ``routing_mesh``                        Routing_mesh Derived Type (inout)
        ! ``hc=1.``                               Value of the hydraulic coefficient (optional)
        ! ``sc=dt./dx.``                          Value of the spreading coefficient, default is set to dt (in second/m)(optional)
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_parameter), intent(inout) :: routing_parameter
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        real,optional,intent(in) :: hc
        real,optional,intent(in) :: sc
        
        integer :: i
        real :: hc_copy, sc_copy
        
        if (.not.allocated(routing_parameter%hc)) then
            allocate(routing_parameter%hc(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_parameter%sc)) then
            allocate(routing_parameter%sc(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_parameter%hc_n)) then
            allocate(routing_parameter%hc_n(routing_mesh%nb_nodes))
        end if
        if (.not.allocated(routing_parameter%sc_n)) then
            allocate(routing_parameter%sc_n(routing_mesh%nb_nodes))
        end if
        
        hc_copy=-99.
        sc_copy=-99.
        if (present(hc)) then
            hc_copy=hc
            if (hc_copy<routing_setup%hydrau_coef_boundaries(1)) then
            hc_copy=routing_setup%hydrau_coef_boundaries(1)
            end if
            if (hc_copy>routing_setup%hydrau_coef_boundaries(2)) then
                hc_copy=routing_setup%hydrau_coef_boundaries(2)
            end if
        end if
        if (present(sc)) then
            sc_copy=sc
            if (sc_copy<routing_setup%spreading_boundaries(1)) then
            sc_copy=routing_setup%spreading_boundaries(1)
            end if
            if (sc_copy>routing_setup%spreading_boundaries(2)) then
                sc_copy=routing_setup%spreading_boundaries(2)
            end if
        end if
        
        if (hc_copy>0.0) then
            routing_parameter%hc=hc_copy
        else
            routing_parameter%hc=1.0 !default value
        end if
        
        if (sc_copy>0.0) then
            routing_parameter%sc=sc_copy ! given in s/m 
        else
            do i=1,routing_mesh%nb_nodes
                routing_parameter%sc(i)=routing_setup%dt/routing_mesh%dx(i)  !default value
            end do
        end if
        
        call normalize_routing_parameters(routing_parameter, routing_setup, routing_mesh)
        
    end subroutine routing_parameter_self_initialisation
    
    subroutine normalize_routing_parameters(routing_parameter, routing_setup, routing_mesh)
        
        type(type_routing_parameter), intent(inout) :: routing_parameter
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        
        routing_parameter%hc_n=(routing_parameter%hc -&
        & routing_setup%hydrau_coef_boundaries(1)) /&
        &(routing_setup%hydrau_coef_boundaries(2)-routing_setup%hydrau_coef_boundaries(1))
        
        routing_parameter%sc_n=(routing_parameter%sc -&
        & routing_setup%spreading_boundaries(1)) /&
        &(routing_setup%spreading_boundaries(2)-routing_setup%spreading_boundaries(1))

    end subroutine normalize_routing_parameters
    
    subroutine unnormalize_routing_parameters(routing_parameter, routing_setup, routing_mesh)
        
        type(type_routing_parameter), intent(inout) :: routing_parameter
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        
        routing_parameter%hc=routing_parameter%hc_n*&
        &(routing_setup%hydrau_coef_boundaries(2)-routing_setup%hydrau_coef_boundaries(1))&
        &+routing_setup%hydrau_coef_boundaries(1)
        
        routing_parameter%sc=routing_parameter%sc_n*&
        &(routing_setup%spreading_boundaries(2)-routing_setup%spreading_boundaries(1))&
        &+routing_setup%spreading_boundaries(1)
        
    end subroutine unnormalize_routing_parameters
    
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
