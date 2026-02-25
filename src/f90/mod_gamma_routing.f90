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
    use mod_gamma_routing_memory
    use mod_gamma_routing_results
    
    implicit none
    
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
    &inflows,routing_states,routing_memory,routing_results)
        
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
        type(type_routing_memory), intent(inout) :: routing_memory
        type(type_routing_results), intent(inout) :: routing_results
        
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: qnetwork
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: vnetwork
        
        integer :: i
        real,dimension(routing_mesh%nb_nodes) :: qmesh
        real,dimension(routing_mesh%nb_nodes) :: velocities
        real,dimension(routing_mesh%nb_nodes) :: inflow
                
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes) :: remainder
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes) :: states
        
        remainder=routing_memory%remainder
        states=routing_memory%states
        
        do i=1,routing_setup%npdt
            velocities=0.
            qmesh=0.
            inflow=inflows(i,1:routing_mesh%nb_nodes)
            call routing_flow(routing_setup,routing_mesh,routing_parameter,inflow,routing_states,remainder,&
            &states,qmesh,velocities)
            qnetwork(i,:)=qmesh
            vnetwork(i,:)=velocities
        end do
        
        routing_memory%remainder=remainder
        routing_memory%states=states
        
        !storing results
        routing_results%discharges=qnetwork
        routing_results%velocities=vnetwork
        
    end subroutine routing_hydrogram
    
    subroutine routing_flow(routing_setup,routing_mesh,routing_parameter,inflows,routing_states,remainder,states,&
    &qnetwork,velocities)
        
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
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: remainder
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: states
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
            call get_discharges(routing_mesh,routing_states,remainder,states,current_node,&
            &inflows(current_node),qcell)
            
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
                call interpolated_routing_coefficients_bilinear(mode,spreading_unn,&
                &index_varying_dx,routing_states,gamma_coefficient)
!~                 call interpolated_routing_coefficients_bicubic(mode,spreading_unn,&
!~                 &index_varying_dx,routing_states,gamma_coefficient)
                
            else
                
                call interpolated_routing_coefficients_linear(mode,index_varying_dx,routing_states,gamma_coefficient)
                
            end if
            
            !if (current_node==10) then 
            !    write(*,*) velocity,mode,gamma_coefficient 
            !endif
            
            call LocalMemStorage(routing_mesh,routing_states,gamma_coefficient,qcell,current_node,&
            &remainder,states)
            
            call MemMassTransfert(routing_states,routing_setup,routing_mesh,current_node,remainder,states)
            
        end do
        
        
    end subroutine routing_flow
    
    
    
    !get the discharge at current node
     subroutine get_discharges(routing_mesh,routing_states,remainder,states,current_node,inflow,qcell)
        
        ! Notes
        ! -----
        ! **get_discharges(routing_mesh,routing_states,current_node,inflow,qcell)** :
        !
        ! - Compute the discharge at the current node, get the discharges from upstream nodes
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (in)
        ! ``routing_memory``                      routing_memory Derived Type (in)
        ! ``current_node``                        Indexe of the current node (in)
        ! ``inflow``                              Input discharge / rainfall in m3/s
        ! ``qcell``                               Discharge at the current cell in m3/s
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(in) :: routing_states
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: remainder
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: states
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
                qrout = qrout + states(1,upstream_node)!ici on est en m3/s
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
        real :: x1,y1,x2,y2
        real :: wx1, wx2,wy1, wy2
        real :: x,y
        !real :: spreading
        
        !Avoid chocs by increasing the spreading with the mode: we could calibrate the exponnent instead of the spreading... This works very well, less chocs when v is low (hydro_coef is low too)
        !spreading=spreading_origin+spreading_origin*mode**0.5
        
        x = max(routing_states%tabulated_delay(1), min(mode, routing_states%tabulated_delay(routing_states%nb_mode)))
        y = max(routing_states%tabulated_spreading(1), &
        &min(spreading, routing_states%tabulated_spreading(routing_states%nb_spreads)))
        
        !search indices for delay x
        ix1=1
        do i=1,routing_states%nb_mode-1
            if (x<routing_states%tabulated_delay(i+1)) then
                ix1=i
                exit
            endif
        end do
        ix1 = max(1, min(ix1, routing_states%nb_mode-1))
        ix2 = min(ix1 + 1, routing_states%nb_mode)
        
        !search indices for spreading y
        iy1=1
        do i=1,routing_states%nb_spreads-1
            if (y<routing_states%tabulated_spreading(i+1)) then
                iy1=i
                exit
            endif
        end do
        iy1 = max(1, min(iy1, routing_states%nb_spreads-1))
        iy2 = min(iy1 + 1, routing_states%nb_spreads)
        
        x1=routing_states%tabulated_delay(ix1)
        x2=routing_states%tabulated_delay(ix2)
        y1=routing_states%tabulated_spreading(iy1)
        y2=routing_states%tabulated_spreading(iy2)
        
!~         gamma_coefficient=&
!~         & ((mode-x2)/(x1-x2)) * ((spreading-y2)/(y1-y2)) * &
!~         &routing_states%tabulated_routing_coef(:,ix1,iy1,index_dx)&
!~         &+((mode-x1)/(x2-x1)) * ((spreading-y2)/(y1-y2)) * &
!~         &routing_states%tabulated_routing_coef(:,ix2,iy1,index_dx)&
!~         &+((mode-x2)/(x1-x2)) * ((spreading-y1)/(y2-y1)) * &
!~         &routing_states%tabulated_routing_coef(:,ix1,iy2,index_dx)&
!~         &+((mode-x1)/(x2-x1)) * ((spreading-y1)/(y2-y1)) * &
!~         &routing_states%tabulated_routing_coef(:,ix2,iy2,index_dx)
        
        wx1 = (x2 - x)/(x2 - x1)
        wx2 = (x - x1)/(x2 - x1)
        wy1 = (y2 - y)/(y2 - y1)
        wy2 = (y - y1)/(y2 - y1)
        
        gamma_coefficient = &
            wx1*wy1 * routing_states%tabulated_routing_coef(:,ix1,iy1,index_dx) &
          + wx2*wy1 * routing_states%tabulated_routing_coef(:,ix2,iy1,index_dx) &
          + wx1*wy2 * routing_states%tabulated_routing_coef(:,ix1,iy2,index_dx) &
          + wx2*wy2 * routing_states%tabulated_routing_coef(:,ix2,iy2,index_dx)
    
    end subroutine interpolated_routing_coefficients_bilinear
    
    
    subroutine interpolated_routing_coefficients_bicubic(mode,spreading,index_dx,routing_states,gamma_coefficient)
        implicit none
        real, intent(in) :: mode, spreading
        integer, intent(in) :: index_dx
        type(type_routing_states), intent(in) :: routing_states
        real, dimension(size(routing_states%quantile)), intent(out) :: gamma_coefficient

        integer :: i, j, ix, iy
        real :: x, y
        real :: x1, x2, x0, x3
        real :: y1, y2, y0, y3
        real :: t, u
        real, dimension(4) :: px, py
        real, dimension(4,4) :: f
        real :: cubic_interp

        !-------------------------------
        ! Clip mode and spreading to bounds
        !-------------------------------
        x = max(routing_states%tabulated_delay(1), min(mode, routing_states%tabulated_delay(routing_states%nb_mode)))
        y = max(routing_states%tabulated_spreading(1), &
        &min(spreading, routing_states%tabulated_spreading(routing_states%nb_spreads)))

        !-------------------------------
        ! Find lower index for mode (x)
        !-------------------------------
        ix = 2
        do i = 2,routing_states%nb_mode
            if (x < routing_states%tabulated_delay(i)) then
                ix = i - 1
                exit
            end if
        end do
        ix = max(2, min(ix, routing_states%nb_mode-2)) ! Ensure we have 4 points for bicubic
        x0 = routing_states%tabulated_delay(ix-1)
        x1 = routing_states%tabulated_delay(ix)
        x2 = routing_states%tabulated_delay(ix+1)
        x3 = routing_states%tabulated_delay(ix+2)

        t = (x - x1) / (x2 - x1)

        !-------------------------------
        ! Find lower index for spreading (y)
        !-------------------------------
        iy = 2
        do i = 2,routing_states%nb_spreads
            if (y < routing_states%tabulated_spreading(i)) then
                iy = i - 1
                exit
            end if
        end do
        iy = max(2, min(iy, routing_states%nb_spreads-2))
        y0 = routing_states%tabulated_spreading(iy-1)
        y1 = routing_states%tabulated_spreading(iy)
        y2 = routing_states%tabulated_spreading(iy+1)
        y3 = routing_states%tabulated_spreading(iy+2)

        u = (y - y1) / (y2 - y1)

        !-------------------------------
        ! Extract 4x4 neighborhood
        !-------------------------------
        do i = 1, size(gamma_coefficient)
            f(1,1) = routing_states%tabulated_routing_coef(i, ix-1, iy-1, index_dx)
            f(2,1) = routing_states%tabulated_routing_coef(i, ix,   iy-1, index_dx)
            f(3,1) = routing_states%tabulated_routing_coef(i, ix+1, iy-1, index_dx)
            f(4,1) = routing_states%tabulated_routing_coef(i, ix+2, iy-1, index_dx)

            f(1,2) = routing_states%tabulated_routing_coef(i, ix-1, iy, index_dx)
            f(2,2) = routing_states%tabulated_routing_coef(i, ix,   iy, index_dx)
            f(3,2) = routing_states%tabulated_routing_coef(i, ix+1, iy, index_dx)
            f(4,2) = routing_states%tabulated_routing_coef(i, ix+2, iy, index_dx)

            f(1,3) = routing_states%tabulated_routing_coef(i, ix-1, iy+1, index_dx)
            f(2,3) = routing_states%tabulated_routing_coef(i, ix,   iy+1, index_dx)
            f(3,3) = routing_states%tabulated_routing_coef(i, ix+1, iy+1, index_dx)
            f(4,3) = routing_states%tabulated_routing_coef(i, ix+2, iy+1, index_dx)

            f(1,4) = routing_states%tabulated_routing_coef(i, ix-1, iy+2, index_dx)
            f(2,4) = routing_states%tabulated_routing_coef(i, ix,   iy+2, index_dx)
            f(3,4) = routing_states%tabulated_routing_coef(i, ix+1, iy+2, index_dx)
            f(4,4) = routing_states%tabulated_routing_coef(i, ix+2, iy+2, index_dx)

            !-------------------------------
            ! Bicubic interpolation in x then y
            !-------------------------------
            do j = 1,4
                px(j) = cubic_hermite(f(1,j), f(2,j), f(3,j), f(4,j), t)
            end do

            gamma_coefficient(i) = cubic_hermite(px(1), px(2), px(3), px(4), u)
        end do

    end subroutine interpolated_routing_coefficients_bicubic
    
    !-------------------------------------------------
    ! Cubic Hermite interpolation (C1 continuous)
    !-------------------------------------------------
    real function cubic_hermite(f0,f1,f2,f3,s)
        implicit none
        real, intent(in) :: f0,f1,f2,f3,s
        real :: a0,a1,a2,a3
        a0 = f1
        a1 = 0.5*(f2 - f0)
        a2 = f0 - 2.5*f1 + 2.0*f2 - 0.5*f3
        a3 = -0.5*f0 + 1.5*f1 - 1.5*f2 + 0.5*f3
        cubic_hermite = a0 + a1*s + a2*s*s + a3*s*s*s
    end function cubic_hermite
    
    
    !>Store and spread in memory the delayed discharge
     subroutine LocalMemStorage(routing_mesh,routing_states,gamma_coefficient,qcell,current_node,remainder,states)
        
        ! Notes
        ! -----
        ! **LocalMemStorage(routing_mesh,routing_states,gamma_coefficient,qcell,current_node)** :
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
        ! ``routing_memory``                      routing_memory derived type (inout)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_mesh),intent(in) :: routing_mesh
        type(type_routing_states),intent(inout) :: routing_states
        real, dimension(size(routing_states%quantile)),intent(in) :: gamma_coefficient
        real, intent(in) :: qcell
        integer, intent(in) :: current_node
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: remainder
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: states
        
        integer :: t
        
        do t=1,routing_states%window_length(current_node)
            states(t,current_node)=states(t,current_node) + gamma_coefficient(t)*qcell
        enddo
        
    endsubroutine LocalMemStorage
    
    
    
    !>Switch up in time the local memory_storage array at position ix,iy
    pure subroutine MemMassTransfert(routing_states,routing_setup,routing_mesh,current_node,remainder,states)
        
        ! Notes
        ! -----
        ! **MemMassTransfert(routing_states,routing_setup,routing_mesh,current_node)** :
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
        ! ``routing_memory``                      routing_memory derived type (inout)
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_states),intent(inout) :: routing_states
        type(type_routing_setup),intent(in) :: routing_setup
        type(type_routing_mesh),intent(in) :: routing_mesh
        integer, intent(in) :: current_node
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: remainder
        real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes), intent(inout) :: states
        
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
                    
                        dqfill=inv_elongation_factor * states(t+1,upstream_node)
                        
                        states(t,upstream_node)= dqfill + remainder(t,upstream_node)
                        remainder(t,upstream_node)=states(t+1,upstream_node)-dqfill
                        
                    enddo
                    
                    !last time step in memory
                    states(nbmemcell,upstream_node)=0.
                    remainder(nbmemcell,upstream_node)=0.
                
                end if
                
            end do
            
        end if
        
    endsubroutine MemMassTransfert
    
end module mod_gamma_routing
    
