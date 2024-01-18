!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr



module mod_gamma_interface
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    use mod_gamma_routing_parameters
    use mod_gamma_routing_states
    use mod_gamma_routing_results
    use mod_gamma_routing
    use mod_gamma_function
    use mod_gamma_sorting
    
    contains
    
    
    subroutine auto_compute_boundaries(routing_setup,routing_mesh,observed_discharges)
        
        ! Notes
        ! -----
        ! **auto_compute_boundaries(routing_setup,routing_mesh,observed_discharges)** :
        !
        ! - Auto-compute the bounds for parameters. Experimental, do not use this function.
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``observed_discharges``                 Observed discharges, array(npdt,nb_nodes) (in)
        ! =============================           ===================================
        
        implicit none
        type(type_routing_setup), intent(inout) :: routing_setup
        type(type_routing_mesh), intent(inout) :: routing_mesh
        real,dimension(routing_setup%npdt,routing_mesh%nb_nodes),intent(in) :: observed_discharges
        
        
        real,dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: observed_discharges_sorted
        real,dimension(routing_setup%npdt) :: discharge_at_station
        integer :: i,j,k
        real :: qspe,qspe_min_obs,qspe_max_obs,quant_10,quant_90
        
        
        observed_discharges_sorted=observed_discharges
        
        k=0
        do i=1,size(observed_discharges,2)
        
            discharge_at_station=observed_discharges_sorted(:,i)
            call sort(discharge_at_station)
            observed_discharges_sorted(:,i)=discharge_at_station
            
            do j=1,size(observed_discharges,1)
                
                k=k+1
                
                if (observed_discharges(j,i)>0.1) then
                
                    if (routing_setup%velocity_computation.eq."qmm") then
                
                        qspe = observed_discharges(j,i) * routing_setup%dt * 1000. &
                        &/ (routing_mesh%cumulated_surface(i) * 1000.0**2.)
                        
                    end if
                    
                    if (routing_setup%velocity_computation.eq."qm3") then
                    
                        qspe = observed_discharges(j,i)
                        
                    end if
                    
                    if (k==1) then
                        qspe_min_obs=qspe
                        qspe_max_obs=qspe
                    else
                        qspe_min_obs = min(qspe_min_obs,qspe)
                        qspe_max_obs = max(qspe_max_obs,qspe)
                    end if
                
                end if
                
            end do
            
        end do
        
        
        
        quant_90=0.
        quant_10=100.
        do i=1,size(observed_discharges_sorted,2)
        
            if (routing_setup%velocity_computation.eq."qmm") then
                
                observed_discharges_sorted(:,i) = observed_discharges_sorted(:,i) * routing_setup%dt * 1000. &
                &/ (routing_mesh%cumulated_surface(i) * 1000.0**2.)
            
            end if
            !remove null q ?
            quant_10=min(observed_discharges_sorted(int(10.*size(observed_discharges_sorted,1)/100.),i),quant_10)
            quant_90=max(observed_discharges_sorted(int(90.*size(observed_discharges_sorted,1)/100.),i),quant_90)
            
        end do
        
        write(*,*) "Quantile 10/90",quant_10,quant_90
        
        if (qspe_min_obs<=0. .or. qspe_max_obs<=0.) then
            write(*,*) "Cannot compute boundaries... left with its default values."
            return
        end if
        
        !write(*,*) qspe_min_obs,qspe_max_obs,routing_setup%vmin,routing_setup%vmax
        routing_setup%hydrau_coef_boundaries(1)=routing_setup%vmin/qspe_min_obs**0.4
        routing_setup%hydrau_coef_boundaries(2)=routing_setup%vmax/qspe_max_obs**0.4
        
        routing_setup%spreading_boundaries(1)=routing_setup%dt
        routing_setup%spreading_boundaries(2)=routing_setup%dt+100.*sqrt(routing_setup%dt)
        
    end subroutine auto_compute_boundaries
    
    
    subroutine routing_gamma_change_parameters(routing_parameter,routing_states,routing_setup,routing_mesh,&
    &hydraulics_coefficient,spreading)
        
        ! Notes
        ! -----
        ! **routing_gamma_change_parameters(routing_parameter,routing_states,routing_setup,routing_mesh,hydraulics_coefficient,spreading)** :
        !
        ! - Change the parameters (routing_parameter) and update routing_states if necessary
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_parameter``                   routing_parameter Derived Type (inout)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``hydraulics_coefficient``              Hydraulic coefficient, real (in)
        ! ``spreading``                           Spreading coefficient (in)
        ! =============================           ===================================
        
        
        implicit none
        
        type(type_routing_parameter), intent(inout) :: routing_parameter
        type(type_routing_states), intent(inout) :: routing_states
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        real,optional,intent(in) :: hydraulics_coefficient
        real,optional,intent(in) :: spreading
        
        if (present(hydraulics_coefficient)) then
            routing_parameter%hydraulics_coefficient=hydraulics_coefficient
        end if
        
        if (present(spreading)) then
        
            routing_parameter%spreading=spreading
            
            if (routing_setup%varying_spread) then
                !routing_states%max_spreading=routing_setup%spreading_boundaries(2)
                !routing_states%nb_spreads=int(routing_states%max_spreading/routing_setup%spreading_discretization_step)+1
            else
                routing_states%nb_spreads=1
                routing_states%max_spreading=maxval(routing_parameter%spreading)
                !recompute and reallocate some variables in routing states
                call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
            end if
            
        end if
        
    end subroutine routing_gamma_change_parameters
    
    
    subroutine routing_gamma_precomputation(routing_setup,routing_mesh,routing_states)
        
        ! Notes
        ! -----
        ! **routing_gamma_precomputation(routing_setup,routing_mesh,routing_states)** :
        !
        ! - Compute all Gamma parameters and the Unit-Hydrogram coefficients (fill routing_states)
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! =============================           ==================================
        
        use mod_gamma_function
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        
        call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
        
    end subroutine routing_gamma_precomputation
    
    
    subroutine routing_states_update(routing_parameter,routing_setup,routing_mesh,routing_states)
        
        ! Notes
        ! -----
        ! **routing_states_update(routing_parameter,routing_setup,routing_mesh,routing_states)** :
        !
        ! - Re-initialise routing_states and re-compute all gamma parameters
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_parameter``                   routing_parameter Derived Type (in)
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_parameter), intent(in) :: routing_parameter
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        
        call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
        call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
        
    end subroutine routing_states_update
    
    
    subroutine routing_gamma_run(routing_setup,routing_mesh,routing_parameter,&
    &inflows,routing_states,routing_results)
        
        ! Notes
        ! -----
        ! **routing_gamma_run(routing_setup,routing_mesh,routing_parameter,inflows,routing_states,routing_results)** :
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
        !real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(inout) :: qnetwork
        !real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(inout) :: vnetwork
        
        call routing_hydrogram(routing_setup,routing_mesh,routing_parameter,&
        &inflows,routing_states,routing_results)
        
    end subroutine routing_gamma_run
    
    
    subroutine routing_gamma_control(routing_setup,routing_mesh,routing_parameter,&
    &inflows,observations,routing_states,routing_results)
        
        ! Notes
        ! -----
        ! **routing_gamma_control(routing_setup,routing_mesh,routing_parameter,inflows,observations,routing_states,qnetwork,vnetwork,cost)** :
        !
        ! - Estimate the parameters of the model with a variationnal algorithm
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_parameter``                   routing_parameter Derived Type (in)
        ! ``inflows``                             Inflows, array(npdt,nb_nodes) (in)
        ! ``observations``                        Discharges observations, array(npdt,nb_nodes) (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! ``qnetwork``                            Discharges in the network, array(npdt,nb_nodes) (inout)
        ! ``vnetwork``                            Velocities in the network, array(npdt,nb_nodes) (inout)
        ! ``cost``                                Cost, function evaluation, real (inout)
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: inflows
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: observations
        type(type_routing_states), intent(inout) :: routing_states
        type(type_routing_results), intent(inout) :: routing_results
        !real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(out) :: qnetwork
        !real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(out) :: vnetwork
        !real, dimension(3), intent(inout) :: tab_cost
        !real, intent(out) :: cost
        
        call control(routing_setup,routing_mesh,routing_parameter,&
            &inflows,observations,routing_states,routing_results)
        
    end subroutine routing_gamma_control
    
    
    subroutine routing_gamma_forward_adjoint_b(routing_setup,routing_mesh,routing_parameter,&
    &inflows,observations,routing_states,routing_results,gradients)
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: inflows
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: observations
        type(type_routing_states), intent(inout) :: routing_states
        type(type_routing_results), intent(inout) :: routing_results
        real,dimension(2,routing_mesh%nb_nodes),intent(inout) :: gradients
        
        type(type_routing_parameter) :: routing_parameterb
        real :: cost
        real :: costb
        
        call routing_parameter_self_initialisation(routing_parameter=routing_parameterb,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,hydraulics_coefficient=0.,spreading=0.)
        
        call routing_states_reset(routing_states)
        
        gradients=0.
        cost=0.
        costb=1.
        call routing_hydrogram_forward_b(routing_setup, &
                &   routing_mesh, routing_parameter, routing_parameterb, inflows, &
                &   observations, routing_states, routing_results, cost, &
                &   costb)
        
        gradients(1,:)=routing_parameterb%hydraulics_coefficient
        gradients(2,:)=routing_parameterb%spreading
        
    end subroutine routing_gamma_forward_adjoint_b
    
    
    
    subroutine routing_gamma_forward_adjoint_b0(routing_setup,routing_mesh,routing_parameter,&
    &inflows,observations,routing_states,routing_results,gradients)
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: inflows
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: observations
        type(type_routing_states), intent(inout) :: routing_states
        type(type_routing_results), intent(inout) :: routing_results
        real,dimension(routing_setup%npdt,routing_mesh%nb_nodes),intent(inout) :: gradients
        
        type(type_routing_parameter) :: routing_parameterb
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: inflowsb
        
        real :: cost
        real :: costb
        
        call routing_parameter_self_initialisation(routing_parameter=routing_parameterb,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,hydraulics_coefficient=0.,spreading=0.)
        
        call routing_states_reset(routing_states)
        
        gradients=0.
        cost=0.
        costb=1.
        call routing_hydrogram_forward_b0(routing_setup, &
                &   routing_mesh, routing_parameter,routing_parameterb, inflows, inflowsb,&
                &   observations, routing_states, routing_results, cost, &
                &   costb)
        
        gradients=inflowsb
        
    end subroutine routing_gamma_forward_adjoint_b0
    
    
    subroutine routing_gamma_cost_function(routing_setup,routing_mesh,routing_parameter,observations,qnetwork,routing_results)
        
        ! Notes
        ! -----
        ! **routing_gamma_cost_function(npdt,routing_mesh,observations,qnetwork,routing_results)** :
        !
        ! - Compute the cost function and return the cost (roots mean square)
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``npdt``                                nimber of time-step, integer (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``observations``                        Discharges observations, array(npdt,nb_nodes) (in)
        ! ``qnetwork``                            Discharges in the network, array(npdt,nb_nodes) (in)
        ! ``routing_results``                     routing_results derived type (inout)
        ! =============================           ===================================
        
        implicit none
    
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: observations
        real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: qnetwork
        type(type_routing_results), intent(inout) :: routing_results
        
        real, dimension(3) :: tab_cost
        real :: cost
        
        
        cost=0.
        tab_cost=0.
        call cost_function(routing_setup,routing_mesh,routing_parameter,observations,qnetwork,tab_cost,cost)
        
        routing_results%costs=tab_cost
        
    end subroutine routing_gamma_cost_function
    
    !interface to interpolated_routing_coefficients_linear
    subroutine routing_gamma_linear_interpolation(delay,index_varying_dx,routing_states,gamma_coefficient)
        implicit none
        real, intent(in) :: delay
        integer, intent(in) :: index_varying_dx
        type(type_routing_states),intent(in) :: routing_states
        real, dimension(size(routing_states%quantile)), intent(out) :: gamma_coefficient
        
        call interpolated_routing_coefficients_linear(delay,index_varying_dx,routing_states,gamma_coefficient)
        
    end subroutine routing_gamma_linear_interpolation
    
    
    !interface to compute_gamma_coefficient
    subroutine interface_compute_gamma_coefficient(scale,mode,quantile,window_shift,&
    &density_function,gamma_coefficient)
        implicit none
        
        real,intent(in) :: scale
        real,intent(in) :: mode
        real,dimension(:),intent(in) :: quantile
        real,intent(in) :: window_shift
        character(3),intent(in) :: density_function
        real,dimension(size(quantile)),intent(out) :: gamma_coefficient
        
        call compute_gamma_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
        
    end subroutine interface_compute_gamma_coefficient
    
    
    !interface to compute_gamma_routing_coefficient
    subroutine interface_compute_routing_coefficient(scale,mode,quantile,window_shift,&
    &density_function,gamma_coefficient)
        implicit none
        
        real,intent(in) :: scale
        real,intent(in) :: mode
        real,dimension(:),intent(in) :: quantile
        real,intent(in) :: window_shift
        character(3),intent(in) :: density_function
        real,dimension(size(quantile)),intent(out) :: gamma_coefficient
        
        call compute_gamma_routing_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
        
    end subroutine interface_compute_routing_coefficient
    
end module mod_gamma_interface

