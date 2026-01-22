!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr



module mod_gamma_routing
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    use mod_gamma_routing_parameters
    use mod_gamma_routing_states
    use mod_gamma_routing_results
    
    implicit none
    
    ! Creation of a local type useful for the routing model (memory) ! this is a trick for the differentiation of the model. This type need to be allocated before use: type(routing_memory), dimension(nb_nodes) : routingmem
    type routing_memory
        real :: states
        real :: remainder
    end type routing_memory
    
    contains
    
    function x_unn(flag,lb,ub,x)
        integer :: flag
        real :: x_unn
        real, intent(in) :: lb
        real, intent(in) :: ub
        real, intent(in) :: x
        
        if (flag == 1) then
            x_unn=x*(ub-lb)+lb
        else 
            x_unn=x
        endif
        
    end function x_unn
    
    subroutine routing_hydrogram(routing_setup,routing_mesh,routing_parameter,&
    &inflows,routing_states,routing_results)
        
        ! Notes
        ! -----
        ! **routing_hydrogram(routing_setup,routing_mesh,routing_parameter,inflows,routing_states,routing_results)** :
        !
        ! - Run the model an propagate the hydrogram thanks to the inflows
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_parameter``                   routing_parameter Derived Type (in)
        ! ``inflows``                             Inflows, array(npdt,nb_nodes) (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! ``routing_results``                     Routing_results Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: inflows
        type(type_routing_states), intent(inout) :: routing_states
        type(type_routing_results), intent(inout) :: routing_results
        
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: qnetwork
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: vnetwork
        
        integer :: i
        real,dimension(routing_mesh%nb_nodes) :: qmesh
        real,dimension(routing_mesh%nb_nodes) :: velocities
        real,dimension(routing_mesh%nb_nodes) :: inflow
        
        type(routing_memory),dimension(size(routing_states%quantile),routing_mesh%nb_nodes) :: routingmem
        
        routingmem(:,:)%states=routing_states%states
        routingmem(:,:)%remainder=routing_states%remainder
        
        do i=1,routing_setup%npdt
            velocities=0.
            qmesh=0.
            inflow=inflows(i,1:routing_mesh%nb_nodes)
            call routing_flow(routing_setup,routing_mesh,routing_parameter,inflow,routing_states,routingmem,qmesh,velocities)
            qnetwork(i,:)=qmesh
            vnetwork(i,:)=velocities
        end do
        
        routing_states%states=routingmem(:,:)%states
        routing_states%remainder=routingmem(:,:)%remainder
        
        !storing results
        routing_results%discharges=qnetwork
        routing_results%velocities=vnetwork
        
    end subroutine routing_hydrogram
    
    subroutine routing_flow(routing_setup,routing_mesh,routing_parameter,inflows,routing_states,routingmem,qnetwork,velocities)
        
        ! Notes
        ! -----
        ! **routing_flow(routing_setup,routing_mesh,routing_parameter,inflows,routing_states,qnetwork,vnetwork)** :
        !
        ! - Routing the flow for one time-step only
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_parameter``                   routing_parameter Derived Type (in)
        ! ``inflows``                             Inflows, array(npdt,nb_nodes) (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! ``qnetwork``                            Discharges in the network, array(npdt,nb_nodes) (inout)
        ! ``vnetwork``                            Velocities in the network, array(npdt,nb_nodes) (inout)
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        real, dimension(routing_mesh%nb_nodes), intent(in) :: inflows
        type(type_routing_states), intent(inout) :: routing_states
        type(routing_memory),dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: routingmem
        real, dimension(routing_mesh%nb_nodes), intent(out) :: qnetwork
        real, dimension(routing_mesh%nb_nodes), intent(out) :: velocities
        
        real, dimension(size(routing_states%quantile)) :: gamma_coefficient
        real :: velocity
        real :: qcell
        real :: mode
        real :: hydro_param_unn, spreading_unn
        integer :: i
        integer :: current_node
        integer :: index_varying_dx
        
        do i=1,routing_mesh%nb_nodes
            
            !order matters
            current_node=routing_mesh%upstream_to_downstream_nodes(i)
            
            !upstream routed discharges + inflow (m3/s)
            call get_discharges(routing_mesh,routing_states,routingmem,current_node,inflows(current_node),qcell)
            
            !qnetwork is always in m3 => output discharges
            qnetwork(current_node)=qcell
            
            !write(*,*) i,current_node,qcell
            hydro_param_unn=x_unn(routing_parameter%normalized,routing_setup%hydrau_coef_boundaries(1),&
            &routing_setup%hydrau_coef_boundaries(2),routing_parameter%hydraulics_coefficient(current_node))
            
            call compute_velocity(hydro_param_unn,&
            &routing_setup,routing_mesh,routing_states,current_node,qcell,velocity)
            
            velocities(current_node)=velocity
            mode=((routing_mesh%dx(current_node)/(velocity))/routing_setup%dt)
            index_varying_dx=routing_mesh%index_varying_dx(current_node)
            
!~             write(*,*) i,qcell,velocity, mode, routing_mesh%dx(current_node),routing_setup%dt
            
            if (routing_setup%varying_spread==1) then
                
                spreading_unn=x_unn(routing_parameter%normalized,routing_setup%spreading_boundaries(1),&
                &routing_setup%spreading_boundaries(2),routing_parameter%spreading(current_node))
!~                 spreading_unn=1.0
                call interpolated_routing_coefficients_bilinear(mode,spreading_unn,&
                &index_varying_dx,routing_states,gamma_coefficient)
                
            else
                
                call interpolated_routing_coefficients_linear(mode,index_varying_dx,routing_states,gamma_coefficient)
                
            end if
            
            !if (current_node==10) then 
            !    write(*,*) velocity,mode,gamma_coefficient 
            !endif
            
            call LocalMemStorage(routing_mesh,routing_states,gamma_coefficient,qcell,current_node,routingmem)
            
            call MemMassTransfert(routing_states,routing_setup,routing_mesh,current_node,routingmem)
            
        end do
        
        
    end subroutine routing_flow
    
    
    
    !get the discharge at current node
     subroutine get_discharges(routing_mesh,routing_states,routingmem,current_node,inflow,qcell)
        
        ! Notes
        ! -----
        ! **get_discharges(routing_mesh,routing_states,routingmem,current_node,inflow,qcell)** :
        !
        ! - Compute the discharge at the current node, get the discharges from upstream nodes
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (in)
        ! ``routingmem``                          routingmem Derived Type (in)
        ! ``current_node``                        Indexe of the current node (in)
        ! ``inflow``                              Input discharge / rainfall in m3/s
        ! ``qcell``                               Discharge at the current cell in m3/s
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(in) :: routing_states
        type(routing_memory), dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(in) :: routingmem
        integer, intent(in) :: current_node
        real, intent(in) :: inflow
        real, intent(out) :: qcell
        
        integer :: i
        integer :: upstream_node
        real :: qrout
        
        qcell = 0.
        qrout = 0.
        
        do i=1,routing_mesh%nb_upstream_nodes
            
            upstream_node=routing_mesh%nodes_linker(i,current_node)
            
            if (upstream_node>0) then
                qrout = qrout + routingmem(1,upstream_node)%states!ici on est en m3/s
            endif
            
        end do
        
        qcell = inflow + qrout ! en m3 
        
    endsubroutine get_discharges
    
    
    pure subroutine compute_velocity(hydraulics_coefficient,routing_setup,routing_mesh,routing_states,current_node,&
    &incoming_discharges,velocity)
        
        ! Notes
        ! -----
        ! **compute_velocity(hydraulics_coefficient,routing_setup,routing_mesh,routing_states,current_node,incoming_discharges,velocity)** :
        !
        ! - Compute the flow velocity for the current node. The velocity is computed with the specific discharge (mm) or with the discharge (m3/s)
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``hydraulics_coefficient``              Hydraulic coefficient (in)
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (in)
        ! ``current_node``                        Indexe of the current node (in)
        ! ``incoming_discharges``                 Input discharge / rainfall in m3/s
        ! ``velocity``                            velocity at the current cell in m/s
        ! =============================           ===================================
        
        implicit none
        real, intent(in) :: hydraulics_coefficient
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(in) :: routing_states
        integer, intent(in) :: current_node
        real, intent(inout) :: incoming_discharges
        real, intent(out) :: velocity
        
        real :: baseflow
        
        !baseflow is mandatory, even very small, to prevent NAN computation in the adjoint code due to very very low incomming discharges
        !It is not related to the velocity computation, but may be du to the spreading in memory ??
        baseflow=0.0000001
        if (incoming_discharges < baseflow ) then
            incoming_discharges=baseflow
        endif
        
        velocity=routing_setup%vmax
        
        if (routing_setup%velocity_computation.eq."qmm") then
            velocity =  hydraulics_coefficient * &
            &( incoming_discharges * routing_setup%dt * 1000. / &
            &(routing_mesh%cumulated_surface(current_node) * 1000.0**2.) )**0.4
        end if
        
        if (routing_setup%velocity_computation.eq."qm3") then
            velocity =  hydraulics_coefficient * incoming_discharges**0.4
        end if
        
        
        if (velocity < routing_setup%vmin ) then
            velocity=routing_setup%vmin
        endif
        
        if (velocity > routing_setup%vmax ) then
            velocity=routing_setup%vmax
        endif
        
    endsubroutine compute_velocity
    
    subroutine interpolated_routing_coefficients_linear(delay,index_dx,routing_states,gamma_coefficient)
        
        ! Notes
        ! -----
        ! **interpolated_routing_coefficients_linear(delay,routing_states,gamma_coefficient)** :
        !
        ! - Linear interpolation of the routing coefficient for a givent delay (i.e mode of the Gamma PDF)
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``delay``                               Hydraulic coefficient (in)
        ! ``routing_states``                      Routing_mesh Derived Type (in)
        ! ``gamma_coefficient``                   Interpolated gamma coefficient for the quantiles (inout)
        ! =============================           ===================================
        
        implicit none
        real, intent(in) :: delay
        integer, intent(in) :: index_dx
        type(type_routing_states),intent(in) :: routing_states
        real, dimension(size(routing_states%quantile)), intent(out) :: gamma_coefficient
        
        integer :: i,lim_sup_ind_mu,lim_inf_ind_mu
        real :: delay_inf,delay_sup
        real, dimension(size(routing_states%quantile)) :: diffcoef
        
        lim_inf_ind_mu=1
        lim_sup_ind_mu=2
        
        !search indices
        do i=1,routing_states%nb_mode
            
            if (delay<routing_states%tabulated_delay(i)) then
                lim_inf_ind_mu=i-1
                lim_sup_ind_mu=i
                exit
            endif
        end do
        
        delay_inf=routing_states%tabulated_delay(lim_inf_ind_mu)
        delay_sup=routing_states%tabulated_delay(lim_sup_ind_mu)
        
        !index_dx=routing_mesh%index_varying_dx(current-node)
        
        !trick for differentiation
        diffcoef=routing_states%tabulated_routing_coef(:,lim_sup_ind_mu,1,index_dx)-&
        &routing_states%tabulated_routing_coef(:,lim_inf_ind_mu,1,index_dx)
        
        gamma_coefficient=routing_states%tabulated_routing_coef(:,lim_inf_ind_mu,1,index_dx)+&
        &(delay-delay_inf)*(diffcoef) / (delay_sup-delay_inf)
        
        
    end subroutine interpolated_routing_coefficients_linear
    
    
    subroutine interpolated_routing_coefficients_bilinear(mode,spreading,index_dx,routing_states,gamma_coefficient)
        
        ! Notes
        ! -----
        ! **interpolated_routing_coefficients_bilinear(delay,spreading,index_dx,routing_states,gamma_coefficient)** :
        !
        ! - Bi-Linear interpolation of the routing coefficient for a givent delay (i.e mode of the Gamma PDF) and the spreading (i.e the spreading of the gamma pdf)
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``delay``                               Hydraulic coefficient (in)
        ! ``spreading``                           Spreading coefficient (in)
        ! ``index_dx``                            index of uniq dx (sorted in ascendant order) (in)
        ! ``routing_states``                      Routing_mesh Derived Type (in)
        ! ``gamma_coefficient``                   gamma coefficient for all quantile (inout)
        ! =============================           ===================================
        
        implicit none
        real, intent(in) :: mode
        real, intent(in) :: spreading
        integer, intent(in) :: index_dx
        type(type_routing_states),intent(in) :: routing_states
        real, dimension(size(routing_states%quantile)), intent(out) :: gamma_coefficient
        
        integer :: i
        integer :: ix1,iy1,ix2,iy2
        real :: x1,y1,x2,y2,fx,fy
!~         real, dimension(size(routing_states%quantile)) :: diffcoefXY1, diffcoefXY2, diffcoefYX
!~         real, dimension(size(routing_states%quantile)) :: coefXY2,coefXY1,coefYX
        !index_dx=routing_mesh%index_varying_dx(current-node)
        
        !Avoid chocs by increasing the spreading with the mode: we could calibrate the exponnent instead of the spreading... This works very well, less chocs when v is low (hydro_coef is low too)
        !true_spreading=spreading*routing_states%param_normalisation(2)+spreading*routing_states%param_normalisation(2)*mode**0.5
        
        !search indices for delay x
        ix1=1
        ix2=2
        do i=2,routing_states%nb_mode
            if (mode<routing_states%tabulated_delay(i)) then
                ix1=i-1
                ix2=i
                exit
            endif
        end do
        
        !search indices for spreading y
        iy1=1
        iy2=2
        do i=2,routing_states%nb_spreads
            if (spreading<routing_states%tabulated_spreading(i)) then
                iy1=i-1
                iy2=i
                exit
            endif
        end do
        
        x1=routing_states%tabulated_delay(ix1)
        x2=routing_states%tabulated_delay(ix2)
        y1=routing_states%tabulated_spreading(iy1)
        y2=routing_states%tabulated_spreading(iy2)
        
        
        
        !trick for differentiation
!~         diffcoefXY1=routing_states%tabulated_routing_coef(:,ix2,iy1,index_dx)-&
!~         &routing_states%tabulated_routing_coef(:,ix1,iy1,index_dx)
        
!~         diffcoefXY2=routing_states%tabulated_routing_coef(:,ix2,iy2,index_dx)-&
!~         &routing_states%tabulated_routing_coef(:,ix1,iy2,index_dx)
        
!~         coefXY1=routing_states%tabulated_routing_coef(:,ix1,iy1,index_dx)+&
!~         &(mode-x1)*(diffcoefXY1) / (x2-x1)
        
!~         coefXY2=routing_states%tabulated_routing_coef(:,ix1,iy2,index_dx)+&
!~         &(mode-x1)*(diffcoefXY2) / (x2-x1)
        
!~         diffcoefYX=coefXY2-coefXY1
        
!~         gamma_coefficient=coefXY2+&
!~         &(spreading-y1)*(diffcoefXY2) / (y2-y1)
        
!~         write(*,*) mode,spreading,gamma_coefficient
        
        gamma_coefficient=&
        & ((mode-x2)/(x1-x2)) * ((spreading-y2)/(y1-y2)) * &
        &routing_states%tabulated_routing_coef(:,ix1,iy1,index_dx)&
        &+((mode-x1)/(x2-x1)) * ((spreading-y2)/(y1-y2)) * &
        &routing_states%tabulated_routing_coef(:,ix2,iy1,index_dx)&
        &+((mode-x2)/(x1-x2)) * ((spreading-y1)/(y2-y1)) * &
        &routing_states%tabulated_routing_coef(:,ix1,iy2,index_dx)&
        &+((mode-x1)/(x2-x1)) * ((spreading-y1)/(y2-y1)) * &
        &routing_states%tabulated_routing_coef(:,ix2,iy2,index_dx)
        
!~         write(*,*) spreading,y1,y2, routing_states%tabulated_routing_coef(:,ix1,iy1,index_dx)
    
    end subroutine interpolated_routing_coefficients_bilinear
    
    
    
    
    !>Store and spread in memory the delayed discharge
     subroutine LocalMemStorage(routing_mesh,routing_states,gamma_coefficient,qcell,current_node,routingmem)
        
        ! Notes
        ! -----
        ! **LocalMemStorage(routing_mesh,routing_states,gamma_coefficient,qcell,current_node,routingmem)** :
        !
        ! - Spreading in memory the discharge thank to the interpolated Gamma coefficients
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_mesh``                        routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (in)
        ! ``gamma_coefficient``                   gamma coefficient for all quantile (in)
        ! ``qcell``                               Discharge at the current cell in m3/s (in)
        ! ``current_node``                        Current node index (in)
        ! ``routingmem``                          routingmem derived type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_mesh),intent(in) :: routing_mesh
        type(type_routing_states),intent(inout) :: routing_states
        real, dimension(size(routing_states%quantile)),intent(in) :: gamma_coefficient
        real, intent(in) :: qcell
        integer, intent(in) :: current_node
        type(routing_memory),dimension(size(routing_states%quantile),routing_mesh%nb_nodes),intent(inout) :: routingmem
        
        integer :: t
        
        do t=1,routing_states%window_length(current_node)
            routingmem(t,current_node)%states=routingmem(t,current_node)%states + gamma_coefficient(t)*qcell
        enddo
        
    endsubroutine LocalMemStorage
    
    
    
    !>Switch up in time the local memory_storage array at position ix,iy
    pure subroutine MemMassTransfert(routing_states,routing_setup,routing_mesh,current_node,routingmem)
        
        ! Notes
        ! -----
        ! **MemMassTransfert(routing_states,routing_setup,routing_mesh,current_node,routingmem)** :
        !
        ! - Switch up in time the memory storage array
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_states``                      Routing_mesh Derived Type (in)
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        routing_mesh Derived Type (in)
        ! ``gamma_coefficient``                   gamma coefficient for all quantile (in)
        ! ``current_node``                        Current node index (in)
        ! ``routingmem``                          routingmem derived type (inout)
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_states),intent(inout) :: routing_states
        type(type_routing_setup),intent(in) :: routing_setup
        type(type_routing_mesh),intent(in) :: routing_mesh
        integer, intent(in) :: current_node
        type(routing_memory),dimension(size(routing_states%quantile),routing_mesh%nb_nodes),intent(inout) :: routingmem
        
        integer :: nbmemcell
        integer :: upstream_node
        integer :: t
        integer :: i
        real :: dqfill
        real :: inv_elongation_factor
        
        if (routing_mesh%cum_node_index(current_node) .gt. 1) then
            
            inv_elongation_factor=1./routing_setup%elongation_factor
            
            do i=1,routing_mesh%nb_upstream_nodes
                
                upstream_node=routing_mesh%nodes_linker(i,current_node)
                
                if (upstream_node>0) then
                    
                    nbmemcell=routing_states%window_length(upstream_node)
                    
                    !loop over delay in memory
                    do t=1,nbmemcell-1
                    
                        dqfill=inv_elongation_factor * routingmem(t+1,upstream_node)%states
                        
                        routingmem(t,upstream_node)%states= dqfill + routingmem(t,upstream_node)%remainder
                        
                        routingmem(t,upstream_node)%remainder=routingmem(t+1,upstream_node)%states-dqfill
                        
                    enddo
                    
                    !last time step in memory
                    routingmem(nbmemcell,upstream_node)%states=0.
                    routingmem(nbmemcell,upstream_node)%remainder=0.
                
                end if
                
            end do
            
        end if
        
    endsubroutine MemMassTransfert
    
end module mod_gamma_routing
    
