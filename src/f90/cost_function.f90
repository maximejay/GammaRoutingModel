!~ This file is part of GammaRouting.

!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


!>costfun : Fonction cout
subroutine cost_function(routing_setup,routing_mesh,routing_parameter,observations,qnetwork,tab_cost,cost_final)
    
    ! Notes
    ! -----
    ! **cost_function(npdt,routing_mesh,observations,qnetwork,cost)
    !
    ! - Compute the cost function and return the cost (roots mean square)
    !        
    ! =============================           ===================================
    ! Parameters                              Description
    ! =============================           ===================================
    ! ``npdt``                                nimber of time-step, integer (in)
    ! ``routing_mesh``                        Routing_mesh Derived Type (in)
    ! ``observations``                        Discharges observations, array(npdt,nb_nodes) (in)
    ! ``qnetwork``                            Discharges in the network, array(npdt,nb_nodes) (inout)
    ! ``cost``                                Cost, function evaluation, real (inout)
    ! =============================           ===================================
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    use mod_gamma_routing_parameters
    
    implicit none
    
    type(type_routing_setup), intent(in) :: routing_setup
    type(type_routing_mesh), intent(in) :: routing_mesh
    type(type_routing_parameter), intent(in) :: routing_parameter
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: observations
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(in) :: qnetwork
    real, dimension(3), intent(inout) :: tab_cost
    real, intent(inout) :: cost_final
    
    real :: penalty,cost
    real :: num,den,meanobs,sumobs
    integer :: j,i,k,numobs
    
    cost=0.
    tab_cost=0.
    cost_final=0.
    
    
    !write(*,*) any((observations-qnetwork)**2.>0)
    
    
    
    
    select case(trim(routing_setup%criteria))
    case("rmse")
        
        cost=0.
        do j=1,routing_mesh%nb_nodes
            
            k=routing_mesh%controlled_nodes(j)
            
            if ( (k>0) .and. (k<=routing_mesh%nb_nodes) ) then
                
                do i=routing_setup%pdt_start_optim,routing_setup%npdt
                    
                    !write(*,*) i,k,(observations(i,k)-qnetwork(i,k))**2.,observations(i,k),qnetwork(i,k)
!~                     if (observations(i,k).ne.qnetwork(i,k))then
                    
!~                         write(*,*) i,k,(observations(i,k)-qnetwork(i,k))**2.,observations(i,k),qnetwork(i,k)
                        
!~                     end if
                    
                    cost=cost+(observations(i,k)-qnetwork(i,k))**2.
                
                end do
                
            endif
            
        enddo
        
    case("nse")
        
        cost=0.
        do j=1,routing_mesh%nb_nodes
            
            k=routing_mesh%controlled_nodes(j)
            
            if ( (k>0) .and. (k<=routing_mesh%nb_nodes) ) then
                
                sumobs=0.
                num=0.
                den=0.
                numobs=0
                
                do i=routing_setup%pdt_start_optim,routing_setup%npdt
                    
                    sumobs=sumobs+observations(i,k)
                    numobs=numobs+1
                    
                end do
                
                if (numobs>0) then
                    meanobs=sumobs/numobs
                endif
                
                do i=routing_setup%pdt_start_optim,routing_setup%npdt
                    
                    num=num+(observations(i,k)-qnetwork(i,k))**2.
                    den=den+(observations(i,k)-meanobs)**2.
                    
                end do
                
                if (den>0.) then
                    cost=cost+num/den
                end if
                
            endif
            
        end do
        
    end select
    
    penalty=0.
    call regularization(routing_mesh,routing_parameter,penalty)
    
    cost_final=routing_setup%ponderation_cost*cost + routing_setup%ponderation_regul*penalty
    
    tab_cost(1)=cost_final
    tab_cost(2)=cost
    tab_cost(3)=penalty
    
!~     write(*,*) ""
!~     write(*,*) "-----------------------------------------------------"
!~     write(*,*) "cost=",cost_final," ; j0=",cost," ; penalty=",routing_setup%ponderation_regul*penalty
!~     write(*,*) "-----------------------------------------------------"
!~     write(*,*) ""
    
    
endsubroutine cost_function
    

subroutine regularization(routing_mesh,routing_parameter,penalty)
    use mod_gamma_routing_parameters
    use mod_gamma_routing_mesh
    
    implicit none
    
    type(type_routing_mesh), intent(in) :: routing_mesh
    type(type_routing_parameter), intent(in) :: routing_parameter
    real,intent(inout) :: penalty
    
    integer :: i, current_node, next_node, previous_node
    
    penalty=0.
    do i=1,routing_mesh%nb_nodes
        !order matters
        current_node=routing_mesh%upstream_to_downstream_nodes(i)
        if (i<routing_mesh%nb_nodes) then
            next_node=routing_mesh%upstream_to_downstream_nodes(i+1)
        else
            next_node=current_node
        endif
        
        if (i>1) then
            previous_node=routing_mesh%upstream_to_downstream_nodes(i-1)
        else
            previous_node=current_node
        endif
        
        penalty=penalty+&
        &0.5*((routing_parameter%hydraulics_coefficient(current_node)-&
        &routing_parameter%hydraulics_coefficient(previous_node))**2.&
        &-(routing_parameter%hydraulics_coefficient(current_node)-&
        &routing_parameter%hydraulics_coefficient(next_node))**2.)**2.
        
        penalty=penalty+&
        &0.5*((routing_parameter%spreading(current_node)-routing_parameter%spreading(previous_node))**2.&
        &-(routing_parameter%spreading(current_node)-routing_parameter%spreading(next_node))**2.)**2.
        
    end do
    
endsubroutine regularization

