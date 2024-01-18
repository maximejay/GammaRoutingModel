!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr



module mod_gamma_routing_mesh
    implicit none
    
    type type_routing_mesh
        integer :: nb_nodes=2 !nbx*nby
        integer :: nb_upstream_nodes=1 !maximum number of upstrems contribution
        real, dimension(:), allocatable :: dx !discretization-space in meter (m)
        integer, dimension(:), allocatable :: index_varying_dx ! index of different dx values in varying_dx
        real, dimension(:), allocatable :: varying_dx !different dx values in the mesh in ascend order
        integer, dimension(:), allocatable :: nodes_indexes ! nodes indexes in the order of reading
        character(100), dimension(:), allocatable :: nodes_names ! nodes names in the same order than nodes indexes
        real, dimension(:), allocatable :: surface ! surface of each node, dimension(nb_nodes)
        real, dimension(:), allocatable :: cumulated_surface ! cumulated surface at each node, dimension(nb_nodes)
        integer, dimension(:,:), allocatable :: nodes_linker ! link between the current nodes and the upstreams nodes : dimension(nb_upstream_nodes,nb_nodes)
        integer, dimension(:), allocatable :: upstream_to_downstream_nodes ! flow path, dimension(nb_nodes)
        integer, dimension(:), allocatable :: cum_node_index ! for optimization only : cumulated node from upstream to downstreamm
        integer, dimension(:), allocatable :: controlled_nodes ! for optimization only : indexes of the controlled nodes
    end type type_routing_mesh
    
    
    contains
    
    subroutine routing_mesh_self_initialisation(routing_mesh,nb_nodes,nb_upstream_nodes,dx)
        
        ! Notes
        ! -----
        ! **routing_mesh_self_initialisation(routing_mesh,nb_nodes,nb_upstream_nodes)** :
        !
        ! - Initialise the routing_mesh derived type with user values, allocate all components of the mesh structure et set default values
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_mesh``                        Routing_mesh Derived Type (inout)
        ! ``nb_nodes=2``                          Number of computation nodes (or cell) (optional)
        ! ``nb_upstream_nodes=1.``                Maximum number of upstream conribution at each node (optional)
        ! ``dx=1000.``                            step-size (optional)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_mesh), intent(inout) :: routing_mesh
        integer, optional :: nb_nodes
        integer, optional :: nb_upstream_nodes
        real, optional ::dx
        
        integer :: i
        
        if (present(nb_nodes)) then 
            routing_mesh%nb_nodes=nb_nodes
        end if
        if (present(nb_upstream_nodes)) then 
            routing_mesh%nb_upstream_nodes=nb_upstream_nodes
        end if
        
        allocate(routing_mesh%dx(routing_mesh%nb_nodes))
        if (present(dx)) then 
            routing_mesh%dx=dx
        else
            routing_mesh%dx=1000.0
        end if
        
        allocate(routing_mesh%index_varying_dx(routing_mesh%nb_nodes))
        allocate(routing_mesh%nodes_linker(routing_mesh%nb_upstream_nodes,routing_mesh%nb_nodes))
        allocate(routing_mesh%upstream_to_downstream_nodes(routing_mesh%nb_nodes))
        allocate(routing_mesh%surface(routing_mesh%nb_nodes))
        allocate(routing_mesh%cumulated_surface(routing_mesh%nb_nodes))
        allocate(routing_mesh%cum_node_index(routing_mesh%nb_nodes))
        allocate(routing_mesh%nodes_indexes(routing_mesh%nb_nodes))
        allocate(routing_mesh%nodes_names(routing_mesh%nb_nodes))
        allocate(routing_mesh%controlled_nodes(routing_mesh%nb_nodes))
        
        !default value
        routing_mesh%surface=1. !km2
        routing_mesh%nodes_indexes=(/ (i, i=1,routing_mesh%nb_nodes) /)
        routing_mesh%nodes_names=""
        routing_mesh%upstream_to_downstream_nodes=(/ (i, i=1,routing_mesh%nb_nodes) /)
        routing_mesh%nodes_linker=0
        routing_mesh%nodes_linker(1,:)=(/ (i, i=0,routing_mesh%nb_nodes-1) /)
        routing_mesh%controlled_nodes=-1
        
        call mesh_uniq_dx(routing_mesh)
        call mesh_compute_cumulated_surface(routing_mesh)
        call mesh_compute_cumulated_node_index(routing_mesh)
        
    end subroutine routing_mesh_self_initialisation
    
    
    subroutine routing_mesh_clear(routing_mesh)
        
        ! Notes
        ! -----
        ! **routing_mesh_clear(routing_mesh)** :
        !
        ! - Clear the derived type routing_mesh
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_mesh``                        routing_mesh Derived Type (inout)
        ! =============================           ===================================
        
        type(type_routing_mesh), intent(inout) :: routing_mesh
        
        type(type_routing_mesh) :: routing_mesh_new
        
        routing_mesh=routing_mesh_new
    endsubroutine routing_mesh_clear
    
    
    subroutine mesh_update(routing_mesh)
        
        ! Notes
        ! -----
        ! **mesh_update(routing_mesh)** :
        !
        ! - update the mesh structure after initialisation and modification. The cumulated surface and the cumulative node index are computed
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_mesh``                        routing_mesh Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_mesh), intent(inout) :: routing_mesh
        
        call mesh_uniq_dx(routing_mesh)
        call mesh_compute_cumulated_surface(routing_mesh)
        call mesh_compute_cumulated_node_index(routing_mesh)
        
    end subroutine mesh_update
    
    
    subroutine mesh_compute_cumulated_surface(routing_mesh)
        
        ! Notes
        ! -----
        ! **mesh_compute_cumulated_surface(routing_mesh)** :
        !
        ! Compute the cumulated surface 
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_mesh``                        routing_mesh Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_mesh), intent(inout) :: routing_mesh
        
        integer :: i,j
        integer :: current_node, upstream_node
        real :: cum_upstream_surface
        
        routing_mesh%cumulated_surface=0.0 !km2
        do i=1,routing_mesh%nb_nodes
        
            current_node=routing_mesh%upstream_to_downstream_nodes(i)
            
            cum_upstream_surface=0.
            
            do j=1,routing_mesh%nb_upstream_nodes
                upstream_node=routing_mesh%nodes_linker(j,current_node)
                if (upstream_node>0) then
                    cum_upstream_surface=cum_upstream_surface+routing_mesh%cumulated_surface(upstream_node)
                end if
            end do
            
            routing_mesh%cumulated_surface(current_node)=routing_mesh%surface(current_node)+cum_upstream_surface
            
        end do
    
    end subroutine mesh_compute_cumulated_surface
    
    
    subroutine mesh_compute_cumulated_node_index(routing_mesh)
        
        ! Notes
        ! -----
        ! **mesh_compute_cumulated_surface(routing_mesh)** :
        !
        ! Compute the cumulative node index
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_mesh``                        routing_mesh Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_mesh), intent(inout) :: routing_mesh
        
        integer :: i,j
        integer :: current_node, upstream_node
        integer :: cum_node_index
        
        routing_mesh%cum_node_index=1 !km2
        do i=1,routing_mesh%nb_nodes
        
            current_node=routing_mesh%upstream_to_downstream_nodes(i)
            
            cum_node_index=0
            
            do j=1,routing_mesh%nb_upstream_nodes
            
                upstream_node=routing_mesh%nodes_linker(j,current_node)
                
                if (upstream_node>0) then
                    cum_node_index=cum_node_index+routing_mesh%cum_node_index(upstream_node)
                end if
                
            end do
            
            routing_mesh%cum_node_index(current_node)=routing_mesh%cum_node_index(current_node)+cum_node_index
            
        end do
    
    end subroutine mesh_compute_cumulated_node_index
    
    
    subroutine mesh_uniq_dx(routing_mesh)
    
        implicit none
        
        type(type_routing_mesh), intent(inout) :: routing_mesh
                
        real, dimension(routing_mesh%nb_nodes) ::unique
        real, dimension(routing_mesh%nb_nodes) ::dx_copy
        integer :: i,j
        real :: min_val, max_val
        
        unique=0.
        min_val = minval(routing_mesh%dx)-1
        max_val = maxval(routing_mesh%dx)
        
        dx_copy=routing_mesh%dx
        i=0
        
        !write(*,*) i, min_val, max_val,dx_copy
        
        do while (min_val<max_val)
            
            i = i+1
            min_val = minval(dx_copy)
            unique(i) = min_val
            
            do j=1,routing_mesh%nb_nodes
            
                if (dx_copy(j)<=min_val) then
                    dx_copy(j)=max_val+1.
                end if
                
            end do
            
            
            !write(*,*) i, min_val, max_val,dx_copy
            !pause
            
        end do
        !min_val = minval(routing_mesh%dx, mask=(routing_mesh%dx>min_val))
        
        !write(*,*) unique
        
        if (allocated(routing_mesh%varying_dx)) then
            deallocate(routing_mesh%varying_dx)
        end if
        
        allocate(routing_mesh%varying_dx(i))   !<-- Or, just use unique(1:i) 
        routing_mesh%varying_dx(:)=unique(1:i)
        
        !fill index_varying_dx
        do i=1,routing_mesh%nb_nodes
            
            do j=1,size(routing_mesh%varying_dx)
            
                if ( routing_mesh%dx(i) == routing_mesh%varying_dx(j) ) then
                    
                    routing_mesh%index_varying_dx(i)=j
                    exit
                    
                end if
            
            end do
            
        end do
        
        !write(*,*) routing_mesh%index_varying_dx
        
    end subroutine mesh_uniq_dx
    
    
    subroutine routing_mesh_copy(routing_mesh,object_copy)
    
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_mesh), intent(out) :: object_copy
        
        object_copy=routing_mesh
    
    end subroutine routing_mesh_copy
    
end module mod_gamma_routing_mesh
