!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr

program routing
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    use mod_gamma_routing_parameters
    use mod_gamma_routing_states
    use mod_gamma_routing_results
    use mod_gamma_function
    use mod_gamma_interface
    use mod_gamma_routing
    
    implicit none
    
    !config
    real :: dt,dx,vmin,vmax,spread,epsilon
    real :: mode
    
    !variable
    real :: scale=0.5
    integer :: window_length=10
    real :: window_shift
    real :: velocity
    real :: delay
    real :: elongation_coefficient
    
    real, dimension(:,:), allocatable :: inflows,qnetwork,vnetwork,observations,random_serie
    real, dimension(:), allocatable :: inflow,qmesh,velocities
    real :: qcell
    real :: cost
    integer :: i, current_node
    integer :: npdt
    
    real,dimension(:),allocatable :: gamma_values,gamma_values_cdf,gamma_coefficient,gamma_coefficient_cdf,&
    &interpolated_gamma_coefficient
    real, dimension(:,:,:), allocatable :: tabulated_gamma_coefficient
    real, dimension(:,:,:), allocatable :: tabulated_gamma_coefficient_3D
    real, dimension(:), allocatable :: tabulated_delay
    real, dimension(:), allocatable :: tabulated_spreading
    real,dimension(:),allocatable :: quantile
    real,dimension(:),allocatable :: adjusted_quantile
    character(300) :: filename
    logical :: old_test
    
    type(type_routing_setup) :: routing_setup
    type(type_routing_mesh) :: routing_mesh
    type(type_routing_states) :: routing_states
    type(type_routing_states) :: routing_states_3D
    type(type_routing_results) :: routing_results
    type(type_routing_parameter) :: routing_parameter
    type(type_routing_parameter) :: routing_parameter_3D
    
    real :: alpha,beta
    
    write(*,*) "hello Gamma-routing"
    write(*,*) ""
    
    call system("mkdir out/")
    
!~     write(*,*) "manual_gradient_test..."
!~     write(*,*) ""
!~     call manual_gradient_test()
    
!~     pause
    
    !simple test case
    
    write(*,*) "routing_setup_self_initialisation..."
    write(*,*) ""
    call routing_setup_self_initialisation(routing_setup,npdt=100,dt=900.,vmin=0.1,vmax=10.,&
    &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=10000.0&
    &,velocity_computation="qm3",varying_spread=1)
    
    write(*,*) "routing_mesh_self_initialisation..."
    write(*,*) ""
    call routing_mesh_self_initialisation(routing_mesh,nb_nodes=10,nb_upstream_nodes=1)
    routing_mesh%controlled_nodes(1)=10
    call mesh_update(routing_mesh)
        
    write(*,*) "routing_parameter_self_initialisation..."
    write(*,*) ""
    call routing_parameter_self_initialisation(routing_parameter=routing_parameter,routing_setup=routing_setup,&
    &routing_mesh=routing_mesh,&
    &hydraulics_coefficient=0.5,spreading=1.)
    
    
    write(*,*) "routing_state_self_initialisation..."
    write(*,*) ""
    call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
    
    
    write(*,*) "compute_gamma_parameters..."
    write(*,*) ""
    call compute_gamma_parameters(routing_setup,routing_mesh,routing_states) !#todo introduire cela dans routing_gamma_run ou dans routing_state_self_initialisation ?
    
    write(*,*) "routing_states:"
    write(*,*) "routing_states%max_mode=",routing_states%max_mode
    write(*,*) "routing_states%min_mode=",routing_states%min_mode
    write(*,*) "routing_states%nb_mode=",routing_states%nb_mode
    write(*,*) "routing_states%scale_coef=",routing_states%scale_coef
    write(*,*) "routing_states%window_length=",routing_states%window_length    
    write(*,*) "routing_states%nb_spreads=",routing_states%nb_spreads    
    write(*,*) "routing_states%max_spreading=",routing_states%max_spreading
    write(*,*) "routing_states%window_shift=",routing_states%window_shift 
    write(*,*) "routing_states%param_normalisation=",routing_states%param_normalisation
    write(*,*) "routing_states%quantile=",routing_states%quantile
    write(*,*) "routing_states%tabulated_delay=",routing_states%tabulated_delay
    write(*,*) "routing_states%tabulated_routing_coef=",routing_states%tabulated_routing_coef(1,1,:,1)
    write(*,*) "routing_states%states=",routing_states%states
    write(*,*) "routing_states%remainder=",routing_states%remainder
    write(*,*) ""
    
    filename="out/tabulated_gamma_coefficient_3D_pdf.txt"
    call write_tabulated_coefficients_3D(routing_states%tabulated_routing_coef(:,:,:,1),routing_states%quantile,filename)
    pause
    
    write(*,*) "routing_results_self_initialisation..."
    call routing_results_self_initialisation(routing_setup,routing_mesh,routing_results)
    
    allocate(inflows(routing_setup%npdt,routing_mesh%nb_nodes))
    !inputs : initialisation
    inflows=0.0
    inflows(1,1)=10.
    
    write(*,*) ""
    write(*,*) "routing_parameter%hydraulics_coefficient=",routing_parameter%hydraulics_coefficient
    write(*,*) "routing_parameter%spreading=",routing_parameter%spreading
    write(*,*) ""
    
    write(*,*) "routing_gamma_run..."
    call routing_gamma_run(routing_setup,routing_mesh,routing_parameter,inflows,&
    &routing_states,routing_results)
    
    write(*,*) "routing_gamma_cost_function..."
    allocate(observations(routing_setup%npdt,routing_mesh%nb_nodes))
    allocate(qnetwork(routing_setup%npdt,routing_mesh%nb_nodes))
    allocate(random_serie(routing_setup%npdt,routing_mesh%nb_nodes))
    call random_init(.true., .true.)
    call random_number(random_serie)
    
    observations=routing_results%discharges*random_serie
    qnetwork=routing_results%discharges
    call routing_gamma_cost_function(routing_setup,routing_mesh,routing_parameter,observations,qnetwork,routing_results)
    
    write(*,*) routing_results%costs
    
    call routing_setup_clear(routing_setup)
    call routing_states_clear(routing_states)
    call routing_mesh_clear(routing_mesh)
    call routing_parameter_clear(routing_parameter)
    call routing_results_clear(routing_results)
    deallocate(inflows)
    deallocate(qnetwork)
    deallocate(observations)
    
    pause
    !dt=3600
    
    write(*,*) "routing_setup_self_initialisation..."
    write(*,*) ""
    call routing_setup_self_initialisation(routing_setup,npdt=100,dt=3600.,vmin=0.01,vmax=10.,&
    &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=10000.0&
    &,velocity_computation="qm3",varying_spread=1)
    
    write(*,*) "routing_mesh_self_initialisation..."
    write(*,*) ""
    call routing_mesh_self_initialisation(routing_mesh,nb_nodes=10,nb_upstream_nodes=1)
    
    write(*,*) "routing_parameter_self_initialisation..."
    write(*,*) ""
    call routing_parameter_self_initialisation(routing_parameter=routing_parameter,routing_setup=routing_setup,&
    &routing_mesh=routing_mesh,&
    &hydraulics_coefficient=0.5,spreading=1.)
    
    write(*,*) "routing_state_self_initialisation..."
    write(*,*) ""
    call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
    
    write(*,*) "compute_gamma_parameters..."
    write(*,*) ""
    call compute_gamma_parameters(routing_setup,routing_mesh,routing_states) 
    
    filename="out/tabulated_gamma_coefficient_3600s_pdf.txt"
    call write_tabulated_coefficients_3D(routing_states%tabulated_routing_coef(:,:,:,1),routing_states%quantile,filename)
    
    call routing_setup_clear(routing_setup)
    call routing_states_clear(routing_states)
    call routing_mesh_clear(routing_mesh)
    call routing_parameter_clear(routing_parameter)
    
    
    !dx=250m
    
    write(*,*) "routing_setup_self_initialisation..."
    write(*,*) ""
    call routing_setup_self_initialisation(routing_setup,npdt=100,dt=900.,vmin=0.01,vmax=10.,&
    &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=10000.0&
    &,velocity_computation="qm3",varying_spread=1)
    
    write(*,*) "routing_mesh_self_initialisation..."
    write(*,*) ""
    call routing_mesh_self_initialisation(routing_mesh,nb_nodes=10,nb_upstream_nodes=1,dx=250.)
    
    write(*,*) "routing_parameter_self_initialisation..."
    write(*,*) ""
    call routing_parameter_self_initialisation(routing_parameter=routing_parameter,routing_setup=routing_setup,&
    &routing_mesh=routing_mesh,&
    &hydraulics_coefficient=0.5,spreading=1.)
    
    write(*,*) "routing_state_self_initialisation..."
    write(*,*) ""
    call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
    
    write(*,*) "compute_gamma_parameters..."
    write(*,*) ""
    call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
    
    filename="out/tabulated_gamma_coefficient_250m_pdf.txt"
    call write_tabulated_coefficients_3D(routing_states%tabulated_routing_coef(:,:,:,1),routing_states%quantile,filename)
    
    call routing_setup_clear(routing_setup)
    call routing_states_clear(routing_states)
    call routing_mesh_clear(routing_mesh)
    call routing_parameter_clear(routing_parameter)
    
    
    !dx=4000m
    
    write(*,*) "routing_setup_self_initialisation..."
    write(*,*) ""
    call routing_setup_self_initialisation(routing_setup,npdt=100,dt=3600.,vmin=0.005,vmax=10.,&
    &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=0.1,ponderation_regul=10000.0&
    &,velocity_computation="qm3",varying_spread=1)
    
    write(*,*) "routing_mesh_self_initialisation..."
    write(*,*) ""
    call routing_mesh_self_initialisation(routing_mesh,nb_nodes=10,nb_upstream_nodes=1,dx=250.)
    
    write(*,*) "routing_parameter_self_initialisation..."
    write(*,*) ""
    call routing_parameter_self_initialisation(routing_parameter=routing_parameter,routing_setup=routing_setup,&
    &routing_mesh=routing_mesh,&
    &hydraulics_coefficient=0.5,spreading=1.)
    
    write(*,*) "routing_state_self_initialisation..."
    write(*,*) ""
    call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
    
    write(*,*) "compute_gamma_parameters..."
    write(*,*) ""
    call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
    
    filename="out/tabulated_gamma_coefficient_3600s250m_pdf.txt"
    call write_tabulated_coefficients_3D(routing_states%tabulated_routing_coef(:,:,:,1),routing_states%quantile,filename)
    
    call routing_setup_clear(routing_setup)
    call routing_states_clear(routing_states)
    call routing_mesh_clear(routing_mesh)
    call routing_parameter_clear(routing_parameter)
    
    
    pause
    
    !complexe test case
    !      o
    !  o  /
    !  \ /
    !   o
    !   |
    !   o
    !
    write(*,*) "routing_setup_self_initialisation..."
    write(*,*) ""
    call routing_setup_self_initialisation(routing_setup,npdt=100,dt=900.,vmin=0.1,vmax=10.,&
    &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=0.2,ponderation_regul=10000.0&
    &,velocity_computation="qm3",varying_spread=1)
    
    
    
    write(*,*) "routing_mesh_self_initialisation..."
    write(*,*) ""
    call routing_mesh_self_initialisation(routing_mesh,nb_nodes=4,nb_upstream_nodes=2)
    
    !meshing : we should import from runoff or text_file ! 2 functions !
    routing_mesh%upstream_to_downstream_nodes=(/1,2,3,4/)
    
    routing_mesh%nodes_linker=0
    routing_mesh%nodes_linker(:,3)=(/1,2/)
    routing_mesh%nodes_linker(1,4)=3
    
    routing_mesh%surface=(/1,1,1,1/)
    routing_mesh%dx=(/1000.0,7000.0,2000.0,1000./)
    
    routing_mesh%controlled_nodes(1)=4
    
    !update and auto compute mesh...
    write(*,*) "mesh_update_new..."
    write(*,*) ""
    call mesh_update(routing_mesh)
    
    
    allocate(inflows(routing_setup%npdt,routing_mesh%nb_nodes))
    allocate(inflow(routing_mesh%nb_nodes))
    allocate(qmesh(routing_mesh%nb_nodes))
    allocate(velocities(routing_mesh%nb_nodes))
    allocate(qnetwork(routing_setup%npdt,routing_mesh%nb_nodes))
    allocate(vnetwork(routing_setup%npdt,routing_mesh%nb_nodes))
    
    !inputs : initialisation
    inflows=0.0
    inflows(1,1)=2.
    inflows(1,2)=5.
    
    qmesh=0.
    velocities=0.0
    qnetwork=0.
    vnetwork=0.
    
    
    write(*,*) "routing_parameter_self_initialisation..."
    write(*,*) ""
    call routing_parameter_self_initialisation(routing_parameter=routing_parameter,routing_setup=routing_setup,&
    &routing_mesh=routing_mesh,&
    &hydraulics_coefficient=0.5,spreading=1.)
    
    
    write(*,*) "routing_state_self_initialisation..."
    write(*,*) ""
    call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
    
    
    write(*,*) "compute_gamma_parameters..."
    write(*,*) ""
    call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
    
    
    write(*,*) "routing_states:"
    !write(*,*) "routing_states%max_scale=",routing_states%max_scale
    write(*,*) "routing_states%max_scale=",routing_states%scale_coef
    write(*,*) "routing_states%window_length=",routing_states%window_length    
    write(*,*) "routing_states%max_mode=",routing_states%max_mode    
    write(*,*) "routing_states%nb_mode=",routing_states%nb_mode    
    write(*,*) "routing_states%nb_spreads=",routing_states%nb_spreads    
    write(*,*) "routing_states%max_spreading=",routing_states%max_spreading
    write(*,*) "routing_states%window_shift=",routing_states%window_shift 
    write(*,*) ""
    
    
    !!!!!!!!!!!!!!CONTROL !!!!!!!!!!!!!!!!!!!!
    
    write(*,*) "routing_results_self_initialisation..."
    call routing_results_self_initialisation(routing_setup,routing_mesh,routing_results)
    
    
    write(*,*) "routing_hydrogram_forward..."
    write(*,*) ""
    allocate(observations(routing_setup%npdt,routing_mesh%nb_nodes))
    observations=0.0
    call routing_states_reset(routing_states)
    call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,inflows,observations,&
    &routing_states,routing_results,cost)
    
    
    observations=routing_results%discharges
    
    filename="out/qnetwork_jonction.txt"
    call write_matrix(routing_setup%npdt,routing_mesh%nb_nodes,routing_results%discharges,filename)
    
    filename="out/vnetwork_jonction.txt"
    call write_matrix(routing_setup%npdt,routing_mesh%nb_nodes,routing_results%velocities,filename)
    
    
    
    write(*,*) "routing_gamma_change_parameters..."
    write(*,*) ""
    call routing_gamma_change_parameters(routing_parameter,routing_states,routing_setup,&
    &routing_mesh,hydraulics_coefficient=0.7,spreading=1.4)
    
    
    write(*,*) "routing_hydrogram_forward..."
    write(*,*) ""
    call routing_states_reset(routing_states)
    call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,inflows,observations,&
    &routing_states,routing_results,cost)
    
    filename="out/qnetwork_initial.txt"
    call write_matrix(routing_setup%npdt,routing_mesh%nb_nodes,routing_results%discharges,filename)
    
    filename="out/vnetwork_initial.txt"
    call write_matrix(routing_setup%npdt,routing_mesh%nb_nodes,routing_results%velocities,filename)
    
    
    write(*,*) "control..."
    write(*,*) ""
    call routing_states_reset(routing_states)
    
    call control(routing_setup,routing_mesh,routing_parameter,&
    &inflows,observations,routing_states,routing_results)
    
    
    filename="out/qnetwork_evaluation.txt"
    call write_matrix(routing_setup%npdt,routing_mesh%nb_nodes,routing_results%discharges,filename)
    
    filename="out/vnetwork_evaluation.txt"
    call write_matrix(routing_setup%npdt,routing_mesh%nb_nodes,routing_results%velocities,filename)
    
    
    deallocate(inflows,qnetwork,vnetwork,inflow,qmesh,velocities)
    call routing_setup_clear(routing_setup)
    call routing_mesh_clear(routing_mesh)
    call routing_states_clear(routing_states)
    call routing_parameter_clear(routing_parameter)
    call routing_results_clear(routing_results)
    
    
    
    
    
    
    
    
!~     !call routing_setup_self_initialisation(routing_setup,dt=900.,vmin=0.1,vmax=10.,&
!~     !&elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=300.)
    
!~     call routing_setup_self_initialisation(routing_setup,dt=900.,vmin=0.1,vmax=10.,&
!~     &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=300.,&
!~     &velocity_computation="qmm")
    
!~     write(*,*) ""
!~     write(*,*) "Routing_setup:"
!~     write(*,*) "routing_setup%dx=",routing_setup%dx
!~     write(*,*) "routing_setup%dt=",routing_setup%dt
!~     write(*,*) "routing_setup%mode_discretization_step=",routing_setup%mode_discretization_step
!~     write(*,*) "routing_setup%upstream_to_downstream_nodes",routing_setup%upstream_to_downstream_nodes
!~     write(*,*) "routing_setup%nodes_linker",routing_setup%nodes_linker
!~     write(*,*) "routing_setup%surface",routing_setup%surface
!~     write(*,*) "routing_setup%cumulated_surface",routing_setup%cumulated_surface
!~     write(*,*) ""
    
!~     call routing_parameter_self_initialisation(routing_parameter,routing_setup,routing_mesh)
    
!~     call routing_state_self_initialisation(routing_setup,routing_mesh,routing_states)
    
!~     write(*,*) "routing_states:"
!~     write(*,*) "routing_states%max_scale=",routing_states%max_scale
!~     write(*,*) "routing_states%window_length=",routing_states%window_length    
!~     write(*,*) "routing_states%max_mode=",routing_states%max_mode    
!~     write(*,*) "routing_states%nb_mode=",routing_states%nb_mode    
!~     write(*,*) "routing_states%nb_spreads=",routing_states%nb_spreads    
!~     write(*,*) "routing_states%max_spreading=",routing_states%max_spreading
!~     write(*,*) "routing_states%window_shift=",routing_states%window_shift 
!~     write(*,*) ""
      
!~     write(*,*) "routing_states%quantile=",routing_states%quantile
!~     write(*,*) ""
!~     write(*,*) "routing_states%tabulated_delay=",routing_states%tabulated_delay
!~     write(*,*) ""
!~     write(*,*) "routing_states%tabulated_routing_coef(:,1)=",routing_states%tabulated_routing_coef(:,1,1) 
    
!~     filename="out/new_tabulated_routing_coefficient_pdf.txt"
!~     call write_tabulated_coefficients(routing_states%tabulated_routing_coef,routing_states%quantile,filename)
    
    
    
!~     write(*,*) ""
!~     write(*,*) "------------ Routing the flow... ---------------"
!~     write(*,*) ""
    
!~     npdt=30
    
!~     call routing_setup_clear(routing_setup)
!~     call routing_setup_self_initialisation(routing_setup,dt=900.,vmin=0.1,vmax=10.,&
!~     &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=1800.,&
!~     &velocity_computation="qmm")
    
!~     !routing_mesh%cumulated_surface=1000.**2./(900.*1000.0)
!~     write(*,*) "routing_mesh%surface",routing_setup%surface
!~     write(*,*) "routing_mesh%cumulated_surface",routing_setup%cumulated_surface
    
!~     call routing_parameter_clear(routing_parameter)
!~     call routing_parameter_self_initialisation(routing_parameter=routing_parameter,routing_setup=routing_setup,&
!~     &routing_mesh=routing_mesh,hydraulics_coefficient=1,spreading=1800.)
    
    
!~     call routing_states_clear(routing_states)
!~     call routing_state_self_initialisation(routing_setup,routing_mesh,routing_states)
    
!~     allocate(inflow(routing_mesh%nb_nodes))
!~     allocate(qmesh(routing_mesh%nb_nodes))
!~     allocate(velocities(routing_mesh%nb_nodes))
!~     allocate(gamma_coefficient(size(routing_states%quantile)))
    
!~     allocate(inflows(routing_setup%npdt,routing_mesh%nb_nodes))
!~     allocate(qnetwork(routing_setup%npdt,routing_mesh%nb_nodes))
!~     allocate(vnetwork(routing_setup%npdt,routing_mesh%nb_nodes))
    
!~     inflows=0.
!~     inflows(1,1)=1.
!~     qnetwork=0.
!~     vnetwork=0.
    
!~     inflow=0.
!~     inflow(1)=1.
!~     qmesh=0.
!~     velocities=0.0
    
    
!~     call routing_hydrogram(npdt,routing_setup,routing_parameter,inflows,routing_states,qnetwork,vnetwork)
    
!~     filename="out/qnetwork.txt"
!~     call write_matrix(npdt,routing_mesh%nb_nodes,qnetwork,filename)
    
!~     filename="out/vnetwork.txt"
!~     call write_matrix(npdt,routing_mesh%nb_nodes,vnetwork,filename)
    
    
!~     write(*,*) ""
!~     write(*,*) "------------------"
!~     write(*,*) "qnetwork=",qnetwork(:,routing_mesh%upstream_to_downstream_nodes(routing_mesh%nb_nodes))
!~     write(*,*) "vnetwork=",vnetwork(:,routing_mesh%upstream_to_downstream_nodes(routing_mesh%nb_nodes))
!~     write(*,*) "------------------"
!~     write(*,*) ""
    
!~     call routing_states_clear(routing_states)
!~     call routing_state_self_initialisation(routing_setup,routing_mesh,routing_states)
!~     call routing_flow(routing_setup,routing_parameter,inflow,routing_states,qmesh,velocities)
!~     inflow=0.0
!~     call routing_flow(routing_setup,routing_parameter,inflow,routing_states,qmesh,velocities)
!~     call routing_flow(routing_setup,routing_parameter,inflow,routing_states,qmesh,velocities)
!~     call routing_flow(routing_setup,routing_parameter,inflow,routing_states,qmesh,velocities)
!~     call routing_flow(routing_setup,routing_parameter,inflow,routing_states,qmesh,velocities)
    
    
!~     write(*,*) ""
!~     write(*,*) "------------------"
!~     write(*,*) "qnetwork=",qmesh(:)
!~     write(*,*) "velocities=",velocities(:)
!~     write(*,*) "------------------"
!~     write(*,*) ""
    
    
!~     call routing_setup_clear(routing_setup)
!~     call routing_mesh_clear(routing_mesh)
!~     call routing_states_clear(routing_states)
!~     call routing_parameter_clear(routing_parameter)
    
!~     scale=0.1169
!~     mode=100.0
!~     alpha=dble(1.0+(mode+1.0)/scale) !position parameter ot the pdf, here the pdf is centered on its mode. The mode is shifted +1
!~     beta=dble(1./scale) 
!~     write(*,*) alpha,beta
!~     write(*,*) GammaCDF(dble(100.5), dble(alpha), dble(beta))-GammaCDF(dble(99.5), dble(alpha), dble(beta))
!~     write(*,*) GammaPDF(dble(40.0), dble(alpha), dble(beta))
    
    !pause
    
!~     write(*,*) ""
!~     write(*,*) "---------- 3D tests cases ------------"
!~     write(*,*) ""
    
!~     call routing_setup_clear(routing_setup)
!~     call routing_setup_self_initialisation(routing_setup,dt=900.,vmin=0.1,vmax=10.,&
!~     &max_spreading=9000.0,elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=300.&
!~     &,velocity_computation="qmm",varying_spread=1)
    
!~     routing_mesh%dx=1000.0
    
!~     call routing_parameter_self_initialisation(routing_setup,routing_mesh,routing_parameter_3D)
    
!~     call routing_state_self_initialisation(routing_setup,routing_mesh,routing_states_3D)
    
!~     write(*,*) "routing_states:"
!~     write(*,*) "routing_states%max_scale=",routing_states_3D%max_scale
!~     write(*,*) "routing_states%window_length=",routing_states_3D%window_length    
!~     write(*,*) "routing_states%max_mode=",routing_states_3D%max_mode    
!~     write(*,*) "routing_states%nb_mode=",routing_states_3D%nb_mode    
!~     write(*,*) "routing_states%nb_spreads=",routing_states_3D%nb_spreads    
!~     write(*,*) "routing_states%max_spreading=",routing_states_3D%max_spreading
!~     write(*,*) "routing_states%window_shift=",routing_states_3D%window_shift 
!~     write(*,*) ""
      
!~     write(*,*) "routing_states%quantile=",routing_states_3D%quantile
!~     write(*,*) ""
!~     write(*,*) "routing_states%tabulated_delay=",routing_states_3D%tabulated_delay
!~     write(*,*) ""
!~     write(*,*) "routing_states%tabulated_spreading=",routing_states_3D%tabulated_spreading
!~     write(*,*) ""
!~     write(*,*) "routing_states%tabulated_routing_coef(:,1,1)=",routing_states_3D%tabulated_routing_coef(:,1,1) 
    
!~     filename="out/new_tabulated_routing_coefficient_3D_pdf.txt"
!~     call write_tabulated_coefficients_3D(routing_states_3D%tabulated_routing_coef,routing_states_3D%quantile,filename)
    
!~     velocity=0.5
!~     mode=routing_mesh%dx(1)/(velocity*routing_setup%dt)
!~     routing_parameter_3D%spreading=900.
!~     elongation_coefficient=1.
    
!~     allocate(interpolated_gamma_coefficient(maxval(routing_states_3D%window_length)))
    
!~     !call interpolated_routing_coefficients_3D(mode,routing_parameter_3D%spreading(1),maxval(routing_states_3D%window_length)&
!~     !&,size(routing_states_3D%tabulated_delay),size(routing_states_3D%tabulated_spreading),&
!~     !&routing_states_3D%tabulated_delay,routing_states_3D%tabulated_spreading,&
!~     !&routing_states_3D%tabulated_routing_coef,interpolated_gamma_coefficient)
    
!~     call new_interpolated_routing_coefficients_linear(mode,routing_states_3D,interpolated_gamma_coefficient)
!~     !call new_interpolated_routing_coefficients_3D(mode,routing_parameter_3D%spreading(1),routing_states_3D,&
!~     !&interpolated_gamma_coefficient)
    
!~     filename="out/new_interpolated_routing_coefficients_bilinear.txt"
!~     call write_coefficients(interpolated_gamma_coefficient,routing_states_3D%quantile,filename)
!~     write(*,*) ""
!~     write(*,*) "New_Interpolated_gamma_coefficient_3D=",interpolated_gamma_coefficient
!~     write(*,*) ""
    
!~     deallocate(interpolated_gamma_coefficient)
    
    
    
    
    
    
    old_test=.false.
    if (old_test) then
        write(*,*) ""
        write(*,*) "---------- OLD TESTS ------------"
        write(*,*) ""
        !config
        dx=1000.0
        dt=900.0
        vmin=0.1
        vmax=10.0
        spread=1800.0
        epsilon=0.1
        velocity=0.8
        elongation_coefficient=1.0  !work for 1<elongation_coefficient<1.15 above value change the shape of the gamma function => for higher delay the spreading decrease. Its an artefact...
        
        !Test functions  compute_gamma_scale and compute_gamma_window
        
        call generic_compute_gamma_scale(dx=dx,dt=dt,vmax=vmax,spread=spread,epsilon=epsilon,scale=scale)
        
        write(*,*) ""
        write(*,*) "scale=",scale
        write(*,*) ""
        
        call generic_compute_gamma_window(dx=dx,dt=dt,vmin=vmin,vmax=vmax,scale=scale,spread=spread,window_length=window_length,&
        &window_shift=window_shift,quantile=quantile)
        
        write(*,*) "window_length=",window_length
        write(*,*) "window_shift=",window_shift
        write(*,*) "quantile=",quantile
        write(*,*) ""
        
        call generic_gamma_elongation_cells(window_length,elongation_coefficient,window_length,&
            &adjusted_quantile)
        
        write(*,*) "New window_length (new window_length with elongated time)=",window_length
        write(*,*) "adjusted_quantile=",adjusted_quantile
        write(*,*) ""
        
        
        ! test functions compute_gamma_coefficient and compute_gamma_routing_coefficient
        deallocate(gamma_coefficient)
        allocate(gamma_values(window_length))
        allocate(gamma_values_cdf(window_length))
        allocate(gamma_coefficient(window_length))
        allocate(gamma_coefficient_cdf(window_length))
            
        mode=dx/(vmin*dt)
        call compute_gamma_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_values)
        call compute_gamma_routing_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_coefficient)
        write(*,*) "gamma_values_vmin=",gamma_values
        write(*,*) "gamma_coefficient_vmin=",gamma_coefficient
        write(*,*) ""
        
        filename="out/gamma_coefficient_vmin.txt"
        call write_coefficients(gamma_coefficient,adjusted_quantile,filename)
        
        mode=dx/(0.5*dt)
        call compute_gamma_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_values)
        call compute_gamma_routing_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_coefficient)
        write(*,*) "gamma_values_vmoy=",gamma_values
        write(*,*) "gamma_coefficient_vmoy=",gamma_coefficient
        write(*,*) ""
        
        filename="out/gamma_coefficient_vmoy.txt"
        call write_coefficients(gamma_coefficient,adjusted_quantile,filename)
        
        mode=dx/(5.0*dt)
        call compute_gamma_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_values)
        call compute_gamma_routing_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_coefficient)
        write(*,*) "gamma_values_vmax=",gamma_values
        write(*,*) "gamma_coefficient_vmax=",gamma_coefficient
        write(*,*) ""
        
        filename="out/gamma_coefficient_vmax.txt"
        call write_coefficients(gamma_coefficient,adjusted_quantile,filename)
        
        
        mode=dx/(0.5*dt)
        call compute_gamma_coefficient(scale,mode,adjusted_quantile,window_shift,"cdf",gamma_values_cdf)
        call compute_gamma_routing_coefficient(scale,mode,adjusted_quantile,window_shift,"cdf",gamma_coefficient_cdf)
        call compute_gamma_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_values)
        call compute_gamma_routing_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_coefficient)
        write(*,*) "gamma_values_cdf=",gamma_values_cdf
        write(*,*) "gamma_coefficient_cdf=",gamma_coefficient_cdf
        write(*,*) "Delta_gamma_values_cdf-pdf=",gamma_values_cdf-gamma_values
        write(*,*) "Delta_gamma_coefficient_cdf-pdf=",gamma_coefficient_cdf-gamma_coefficient
        write(*,*) ""
        
        filename="out/gamma_coefficient_cdf.txt"
        call write_coefficients(gamma_coefficient_cdf,adjusted_quantile,filename)
        
        
        
        !Test functions tabulated_routing_coefficients_2D 
        
        call generic_tabulated_routing_coefficients_2D(scale,adjusted_quantile,window_shift,(dx/vmin/dt),0.1,"pdf",&
        &tabulated_gamma_coefficient)
        filename="out/tabulated_routing_coefficient_pdf.txt"
        call write_tabulated_coefficients(tabulated_gamma_coefficient,adjusted_quantile,filename)
        
        call generic_tabulated_routing_coefficients_2D(scale,adjusted_quantile,window_shift,(dx/vmin/dt),0.1,"cdf",&
        &tabulated_gamma_coefficient)
        filename="out/tabulated_routing_coefficient_cdf.txt"
        call write_tabulated_coefficients(tabulated_gamma_coefficient,adjusted_quantile,filename)
        
        !test tabulated_delay_for_gamma computation
        
        call tabulated_delay_for_gamma(int((dx/vmin/dt)/0.1),0.1,tabulated_delay)
        write(*,*) "tabulated_delay=",tabulated_delay
        write(*,*) ""
        
        
        !Test comparaison gamma_table and tabulated_gamma_function_2D cdf and pdf
        
        !call gamma_table(window_length,ceiling(((dx/vmin/dt)+1)/0.1),0.1,scale,&
        !&elongation_coefficient,tabulated_delay,tabulated_gamma_coefficient(:,:,1))
        !filename="out/tabulated_gamma_coefficient_Igor.txt"
        !call write_tabulated_coefficients(tabulated_gamma_coefficient,adjusted_quantile,filename)
        !write(*,*) tabulated_delay
        
        call generic_tabulated_gamma_function_2D(scale,adjusted_quantile,window_shift,(dx/vmin/dt),0.1,"cdf",&
        &tabulated_gamma_coefficient)
        filename="out/tabulated_gamma_coefficient_cdf.txt"
        call write_tabulated_coefficients(tabulated_gamma_coefficient,adjusted_quantile,filename)
        
        call generic_tabulated_gamma_function_2D(scale,adjusted_quantile,window_shift,(dx/vmin/dt),0.1,"pdf",&
        &tabulated_gamma_coefficient)
        filename="out/tabulated_gamma_coefficient_pdf.txt"
        call write_tabulated_coefficients(tabulated_gamma_coefficient,adjusted_quantile,filename)
        
        
        
        !Test function interpoalted_gamma_coefficients
        
        allocate(interpolated_gamma_coefficient(window_length))
        
        call generic_tabulated_routing_coefficients_2D(scale,adjusted_quantile,window_shift,(dx/vmin/dt),0.1,"pdf",&
        &tabulated_gamma_coefficient)
        
        call tabulated_delay_for_gamma(int((dx/vmin/dt)/0.1),0.1,tabulated_delay)
        
        velocity=0.5
        mode=(dx/velocity)/dt
        
        call compute_gamma_routing_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_coefficient)
        filename="out/gamma_coefficient_true.txt"
        call write_coefficients(gamma_coefficient,adjusted_quantile,filename)
        
        write(*,*) ""
        write(*,*) "Interpolated gamma coefficient for delay=",mode
        write(*,*) ""
        
        call generic_interpolated_routing_coefficients_linear(mode,1,window_length,size(tabulated_delay),1,tabulated_delay,&
            &tabulated_gamma_coefficient,interpolated_gamma_coefficient)
        filename="out/interpolated_gamma_coefficient.txt"
        call write_coefficients(interpolated_gamma_coefficient,adjusted_quantile,filename)
        write(*,*) ""
        write(*,*) "Interpolated_gamma_coefficient=",interpolated_gamma_coefficient
        write(*,*) ""
        write(*,*) ""
        write(*,*) "True_gamma_coefficient=",gamma_coefficient
        write(*,*) ""
        
        
        !test function tabulated_routing_coefficients_3D
        
        velocity=0.5
        mode=dx/(velocity*dt)
        spread=2100.
        elongation_coefficient=1.
        
        deallocate(quantile,adjusted_quantile)
        
        call generic_compute_gamma_scale(dx=dx,dt=dt,vmax=vmax,spread=10.*dt,epsilon=epsilon,scale=scale)
        call generic_compute_gamma_window(dx=dx,dt=dt,vmin=vmin,vmax=vmax,scale=scale,spread=10.*dt,window_length=window_length,&
        &window_shift=window_shift,quantile=quantile)
        call generic_gamma_elongation_cells(window_length,elongation_coefficient,window_length,adjusted_quantile)
            
        call generic_tabulated_routing_coefficients_3D(dx,dt,vmax,epsilon,adjusted_quantile,window_shift, dx/vmin/dt, 10.*dt, &
            & 0.1,300., "pdf", tabulated_gamma_coefficient_3D)
            
        call tabulated_spreading_for_gamma(int(10.*dt/300),300.,tabulated_spreading)
        write(*,*) "tabulated_spreading=",tabulated_spreading
        write(*,*) ""
        
        filename="out/tabulated_gamma_coefficient_3D_pdf.txt"
        call write_tabulated_coefficients_3D(tabulated_gamma_coefficient_3D,adjusted_quantile,filename)
        
        
        
        quantile=0.
        adjusted_quantile=0.
        deallocate(quantile,adjusted_quantile,gamma_coefficient,interpolated_gamma_coefficient)
        
        call generic_compute_gamma_scale(dx=dx,dt=dt,vmax=vmax,spread=10.*dt,epsilon=epsilon,scale=scale)
        call generic_compute_gamma_window(dx=dx,dt=dt,vmin=vmin,vmax=vmax,scale=scale,spread=10.*dt,window_length=window_length,&
        &window_shift=window_shift,quantile=quantile)
        call generic_gamma_elongation_cells(window_length,elongation_coefficient,window_length,adjusted_quantile)
            
        allocate(gamma_coefficient(window_length))
        
        call generic_compute_gamma_scale(dx=dx,dt=dt,vmax=vmax,spread=spread,epsilon=epsilon,scale=scale)
        call compute_gamma_routing_coefficient(scale,mode,adjusted_quantile,window_shift,"pdf",gamma_coefficient)
        filename="out/gamma_coefficient_3D_true.txt"
        call write_coefficients(gamma_coefficient,adjusted_quantile,filename)
        
        write(*,*) ""
        write(*,*) "Gamma_coefficient_3D_true=",gamma_coefficient !for mode 2.2222222
        write(*,*) ""
        
        write(*,*) ""
        write(*,*) "tabulated_gamma_coefficient_3D_true=",tabulated_gamma_coefficient_3D(:,23,7) !for mode 2.2
        write(*,*) ""
        
        write(*,*) ""
        write(*,*) "Interpolated gamma coefficient 3D for delay=",mode !interpolation
        write(*,*) ""
        
        allocate(interpolated_gamma_coefficient(window_length))
        call generic_interpolated_routing_coefficients_bilinear(mode,spread,1,window_length,size(tabulated_delay),&
        &size(tabulated_spreading),1,tabulated_delay,tabulated_spreading,tabulated_gamma_coefficient_3D,&
        &interpolated_gamma_coefficient)
        filename="out/interpolated_routing_coefficients_bilinear.txt"
        call write_coefficients(interpolated_gamma_coefficient,adjusted_quantile,filename)
        write(*,*) ""
        write(*,*) "Interpolated_gamma_coefficient_3D=",interpolated_gamma_coefficient
        write(*,*) ""
    
    
    endif
    
    !plot validation graphics
    write(*,*) "plotting results..."
    call system("gnuplot ../src/gnuplot/plot_coefficients.py")
    
end program routing
