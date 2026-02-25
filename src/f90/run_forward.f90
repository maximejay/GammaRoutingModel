!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


subroutine routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,&
&inflows,observations,routing_states,routing_memory,routing_results,cost)
    
    ! Notes
    ! -----
    ! **routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,inflows,observations,routing_states,routing_memory,routing_results,cost)** :
    !
    ! - Run the model an propagate the hydrogram thanks to the inflows and compute the cost function. This subroutine is differentiable
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
    ! ``routing_memory``                      routing_memory Derived Type (inout)
    ! ``routing_results``                     Routing_results Derived Type (inout)
    ! ``cost``                                Cost function evaluation (inout)
    ! =============================           ===================================
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    use mod_gamma_routing_parameters
    use mod_gamma_routing_states
    use mod_gamma_routing_memory
    use mod_gamma_routing_results
    use mod_gamma_routing
    
    implicit none
    
    type(type_routing_setup), intent(in) :: routing_setup
    type(type_routing_mesh), intent(in) :: routing_mesh
    type(type_routing_parameter), intent(inout) :: routing_parameter
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: inflows
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: observations
    type(type_routing_states), intent(inout) :: routing_states
    type(type_routing_memory), intent(inout) :: routing_memory
    type(type_routing_results), intent(inout) :: routing_results
    real, intent(out) :: cost
    
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: qnetwork
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes) :: vnetwork
    real, dimension(3) :: tab_cost
    integer :: i
    real,dimension(routing_mesh%nb_nodes) :: qmesh
    real,dimension(routing_mesh%nb_nodes) :: velocities
    real,dimension(routing_mesh%nb_nodes) :: inflow
        
    real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes) :: remainder
    real, dimension(size(routing_states%quantile),routing_mesh%nb_nodes) :: states
    
    if (routing_setup%hydraulics_coef_uniform==1) then
        do i=1,routing_mesh%nb_nodes
            routing_parameter%hydraulics_coefficient(i)=routing_parameter%hydraulics_coefficient(1)
        end do
    end if
    if (routing_setup%spreading_uniform==1) then
        do i=1,routing_mesh%nb_nodes
            routing_parameter%spreading(i)=routing_parameter%spreading(1)
        end do
    end if
    
    !Here we should denormalise parameter if nedeed:
    !if routing_setup_normalized=True:
    !   call denormalized(routing_parameter)
    
    remainder=routing_memory%remainder
    states=routing_memory%states
    
    do i=1,routing_setup%npdt
        velocities=0.
        qmesh=0.
        inflow=inflows(i,1:routing_mesh%nb_nodes)
        call routing_flow(routing_setup,routing_mesh,routing_parameter,inflow,routing_states,&
        &remainder,states,qmesh,velocities)
        
        qnetwork(i,:)=qmesh
        vnetwork(i,:)=velocities
        
    end do
    
    call cost_function(routing_setup,routing_mesh,routing_parameter,observations,qnetwork,tab_cost,cost)
    
    !storing results
    routing_results%costs=tab_cost
    routing_results%discharges=qnetwork
    routing_results%velocities=vnetwork
    
    routing_memory%remainder=remainder
    routing_memory%states=states

end subroutine routing_hydrogram_forward

