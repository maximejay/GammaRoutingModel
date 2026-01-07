!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022 Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr



subroutine manual_gradient_test()
        
        use mod_gamma_routing_setup
        use mod_gamma_routing_mesh
        use mod_gamma_routing_parameters
        use mod_gamma_routing_states
        use mod_gamma_function
        use mod_gamma_routing
        
        !use mod_gamma_routing_diff
        !use mod_gamma_routing_diff_d
        !use MOD_GAMMA_ROUTING_DIFF
        !use MOD_GAMMA_ROUTING_DIFF_D
        
        implicit none
        
        type(type_routing_setup) :: routing_setup
        type(type_routing_mesh) :: routing_mesh
        type(type_routing_parameter) :: routing_parameter
        type(type_routing_states) :: routing_states
        real :: cost
        
        real, dimension(:,:), allocatable :: inflows,qnetwork,vnetwork,observations
        
        type(type_routing_parameter) :: routing_parameterd
        real, dimension(:,:), allocatable :: qnetworkd
        real :: costd
        
        type(type_routing_parameter) :: routing_parameterb
        real, dimension(:,:), allocatable :: qnetworkb
        real :: costb
        
       
        
        
        
        
        CHARACTER(200) :: FILENAME
        CHARACTER(2) :: int_to_str

        INTEGER :: ni,iter
        INTEGER :: VALIDATION ! Statut de la validation du test du gradient 1 OK, 0 Failed

        REAL :: pdx
        REAL , DIMENSION(3) :: COST_F ! Valeur du cout : ecart entre simule et observe : Tableau de dimension NITER MAX =100 (devra etre un parametre)
        REAL :: THRESHOLD_GRAD_VALIDATION !Valeur du seuil de validation du test du gradient

        REAL, dimension(4) :: dx_tab !Tableau Perturbation de l element
        REAL, dimension(4) :: dx_tab1 !Tableau Perturbation de l element
        REAL, dimension(4) :: dx_tab2 !Tableau Perturbation de l element

        REAL, dimension(2) :: grad !Tableau des gradients de l'adjoint
        
        integer :: var_spread

        DATA dx_tab1 /0.001,0.002,0.005,0.01/ 
        DATA dx_tab2 /0.1,1.,10.,100./ 
        
        WRITE(*,*) "Test du gradient : Les Resultats devraient satisfaire : derivative = costd = Cost_adjoint"
        
        ! Test du gradient
        ni=1 ! Numero du paramètre perturbe

        if (ni.eq.1) then
            dx_tab=dx_tab1
            var_spread=0
        elseif(ni.eq.2)then
            dx_tab=dx_tab2
            var_spread=1
        endif

        write(*,*) dx_tab
        
        write(*,*) "routing_setup_self_initialisation..."
        write(*,*) ""
        call routing_setup_self_initialisation(routing_setup,npdt=30,dt=900.,vmin=0.1,vmax=10.,&
        &elongation_factor=1.0,mode_discretization_step=0.1,spreading_discretization_step=100.&
        &,velocity_computation="qm3",varying_spread=var_spread)
        
        !routing_setup%hydrau_coef_boundaries=5.0
        !routing_setup%spreading_boundaries=7200.0
        
        write(*,*) "routing_mesh_self_initialisation..."
        write(*,*) ""
        call routing_mesh_self_initialisation(routing_mesh,nb_nodes=5,nb_upstream_nodes=1)
        
        
        
        allocate(inflows(routing_setup%npdt,routing_mesh%nb_nodes))
        allocate(qnetwork(routing_setup%npdt,routing_mesh%nb_nodes))
        allocate(qnetworkd(routing_setup%npdt,routing_mesh%nb_nodes))
        allocate(qnetworkb(routing_setup%npdt,routing_mesh%nb_nodes))
        allocate(vnetwork(routing_setup%npdt,routing_mesh%nb_nodes))
        
        !inputs : initialisation
        inflows=0.0
        inflows(1,1)=10.
        qnetwork=0.
        vnetwork=0.
        
        write(*,*) "routing_parameter_self_initialisation..."
        write(*,*) ""
        
        if (ni.eq.1) then
            call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
            &routing_setup=routing_setup,routing_mesh=routing_mesh,&
            &hydraulics_coefficient=0.5,spreading=1800.)
        elseif(ni.eq.2)then
            call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
            &routing_setup=routing_setup,routing_mesh=routing_mesh,&
            &hydraulics_coefficient=1.0,spreading=2100.)
        endif
        
        
        
        write(*,*) "routing_state_self_initialisation..."
        write(*,*) ""
        call routing_state_self_initialisation(routing_setup,routing_mesh,routing_parameter,routing_states)
        
        write(*,*) "compute_gamma_parameters..."
        write(*,*) ""
        call compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
        
        write(*,*) "routing_hydrogram_forward..."
        write(*,*) ""
        allocate(observations(routing_setup%npdt,routing_mesh%nb_nodes))
        observations=0.0
        call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,inflows,observations,&
        &routing_states,qnetwork,vnetwork,cost)
        
        observations=qnetwork
        
        
        
        
        do iter=1,size(dx_tab)
        
            pdx=dx_tab(iter)

            WRITE(*,*) "Perturbation du parametre ni = ", ni
            WRITE(*,*) "Perturbation pdx = ", pdx
            WRITE(*,*) ""
            WRITE(*,*) "Call grd1d x3 :"
            
            
            call routing_states_reset(routing_states)
            call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
            &routing_setup=routing_setup,routing_mesh=routing_mesh,&
            &hydraulics_coefficient=1.0,spreading=1800.)
               
            ! run grd accoring spatial parameters  param and spatial time-varying states 
            cost=0.
            call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,&
            &inflows,observations,routing_states,qnetwork,vnetwork,cost)
            COST_F(1)=cost
            write(*,*) "At starting point, initial cost=",COST_F(1)
                        
            THRESHOLD_GRAD_VALIDATION=abs(0.001*COST_F(1)) ! Valeur du seuil du succès du test 0.1% du cout

            ! On perturbe l element +pdx
            ! On appelle ensuite le modele et on stocke le cout

            
            call routing_states_reset(routing_states)
            if (ni.eq.1) then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0+pdx,spreading=1800.)
            elseif(ni.eq.2)then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0,spreading=1800.+pdx)
            endif
           
            

           ! run grd accoring spatial parameters  param and spatial time-varying states 
            cost=0.
            call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,&
            &inflows,observations,routing_states,qnetwork,vnetwork,cost)
            COST_F(2)=cost
            write(*,*) "cost +pdx=",COST_F(2)
                        
            
            ! On perturbe l element +pdx
            ! On appelle ensuite le modele et on stocke le cout

            
    
            call routing_states_reset(routing_states)
            if (ni.eq.1) then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0-pdx,spreading=1800.)
            elseif(ni.eq.2)then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0,spreading=1800.-pdx)
            endif
           
           

           ! run grd accoring spatial parameters  param and spatial time-varying states 
            cost=0.
            call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,&
            &inflows,observations,routing_states,qnetwork,vnetwork,cost)
            COST_F(3)=cost
            write(*,*) "cost -pdx=",COST_F(3)
            
                        
            
            !Calcul de la variation du cout:
            WRITE(*,*) "Direct Model, derivative = ", (COST_F(2)-COST_F(3))/2.0/pdx
            WRITE(*,*) ""



            ! Calcul du Lineaire tangent
            ! Param
            WRITE(*,*) "Call TANGENT LINEAR MODEL SIMU_FORWARD_D :"
            
            
            
            call routing_states_reset(routing_states)
            if (ni.eq.1) then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0,spreading=1800.)
                
                call routing_parameter_self_initialisation(routing_parameter=routing_parameterd,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=pdx,spreading=0.0)
                
            elseif(ni.eq.2)then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0,spreading=1800.)
                
                call routing_parameter_self_initialisation(routing_parameter=routing_parameterd,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=0.0,spreading=pdx)
            endif
            
            ! run grd accoring spatial parameters  param and spatial time-varying states 
            cost=0.
            !call tangeant
            costd=1.
            call ROUTING_HYDROGRAM_FORWARD_D( routing_setup, &
            &   routing_mesh, routing_parameter, routing_parameterd, inflows, &
            &   observations, routing_states, qnetwork, qnetworkd, vnetwork, cost, &
            &   costd)
            write(*,*) "New Cost = ",cost
            
            WRITE(*,*) "Lineaire Tangent Cost_D = ", costd/pdx ! Ecriture de la variation du cout par le lineaire tangent :


            WRITE(*,*) ""
            
            
            WRITE(*,*) "CALL ADJOINT MODEL ... :"
            
            call routing_states_reset(routing_states)
            if (ni.eq.1) then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0,spreading=1800.)
                
                call routing_parameter_self_initialisation(routing_parameter=routing_parameterb,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=0.,spreading=0.)
                
            elseif(ni.eq.2)then
                call routing_parameter_self_initialisation(routing_parameter=routing_parameter,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=1.0,spreading=1800.)
                
                call routing_parameter_self_initialisation(routing_parameter=routing_parameterb,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,&
                &hydraulics_coefficient=0.,spreading=0.)
            endif
            
            !call adjoint
            cost=0.
            costb=1.
            
            call ROUTING_HYDROGRAM_FORWARD_B( routing_setup, &
            &   routing_mesh, routing_parameter, routing_parameterb, inflows, &
            &   observations, routing_states, qnetwork, qnetworkb, vnetwork, cost, &
            &   costb)
            
            write(*,*) "New Cost = ",cost
            grad(1)=sum(routing_parameterb%hydraulics_coefficient)
            grad(2)=sum(routing_parameterb%spreading)
            write(*,*)"Adjoint Gradient = ",grad


            WRITE(*,*) ""

            if (abs((COST_F(2)-COST_F(3))/2.0/pdx-costd/pdx).le.THRESHOLD_GRAD_VALIDATION .AND. &
            & abs(costd/pdx-grad(ni)).le.THRESHOLD_GRAD_VALIDATION .AND. &
            & abs((COST_F(2)-COST_F(3))/2.0/pdx-grad(ni)).le.THRESHOLD_GRAD_VALIDATION) then
                write(*,*) "Gradient test OK"
                VALIDATION=1
            else
                write(*,*) "Gradient test failed"
                write(*,*) "Threshold = ",THRESHOLD_GRAD_VALIDATION
                VALIDATION=0
                if (abs((COST_F(2)-COST_F(3))/2.0/pdx-costd/pdx).gt.THRESHOLD_GRAD_VALIDATION) then
                    write(*,*) "Delta COST - Cost_d = ", abs((COST_F(2)-COST_F(3))/2.0/pdx-costd/pdx)
                endif
                 if (abs(costd/pdx-grad(ni)).gt.THRESHOLD_GRAD_VALIDATION) then
                    write(*,*) "Cost_d - Adjoint_Gradient = ", abs(costd/pdx-grad(ni))
                endif
                 if (abs((COST_F(2)-COST_F(3))/2.0/pdx-grad(ni)).gt.THRESHOLD_GRAD_VALIDATION) then
                    write(*,*) "Delta COST - Adjoint_Gradient = ", abs((COST_F(2)-COST_F(3))/2.0/pdx-grad(ni))
                endif
            endif

            !Ecriture dans un fichier des resultats
            WRITE(*,*) ""
            write(int_to_str,'(I2)') ni
            filename='Test_gradient_param_'//trim(int_to_str)//'.txt'
            if (iter.eq.1) then
                open(10,FILE=trim(filename),STATUS='REPLACE')
                write(10,*) "iter   pdx  Direct  TLM    Adjoint    Validation"
            else
                open(10,FILE=trim(filename),access='append')
            endif
            write(10,*) iter,pdx,(COST_F(2)-COST_F(3))/2.0/pdx,costd/pdx,grad(ni),VALIDATION
            close(10)

        enddo
            
    endsubroutine manual_gradient_test

