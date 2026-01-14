!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr



!subroutine controleur
!>@states : spatial and time-dependent variables
!>@param : spatial parameters
subroutine control(routing_setup,routing_mesh,routing_parameter,&
&inflows,observations,routing_states,routing_results)
    
    ! Notes
    ! -----
    ! **control(routing_setup,routing_mesh,routing_parameter,inflows,observations,routing_states,qnetwork,vnetwork,cost)** :
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
    ! ``routing_results``                     Routing_results Derived Type (inout)
    ! =============================           ===================================
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    use mod_gamma_routing_parameters
    use mod_gamma_routing_states
    use mod_gamma_routing_results
    use mod_gamma_routing
    use mod_gamma_interface
    
    implicit none
    
    type(type_routing_setup), intent(inout) :: routing_setup
    type(type_routing_mesh), intent(inout) :: routing_mesh
    type(type_routing_parameter), intent(inout) :: routing_parameter
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(inout) :: inflows
    real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(inout) :: observations
    type(type_routing_states), intent(inout) :: routing_states
    type(type_routing_results), intent(inout) :: routing_results
    !real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(inout) :: qnetwork
    !real, dimension(routing_setup%npdt,routing_mesh%nb_nodes), intent(inout) :: vnetwork
    !real, dimension(3) :: tab_cost
    !real, intent(inout) :: cost
    
    !type(type_routing_setup) :: routing_setup
    !type(type_routing_mesh) :: routing_mesh
    type(type_routing_parameter) :: routing_parameterb
    type(type_routing_parameter) :: routing_parameter_initial
    real :: cost
    real, dimension(:,:), allocatable :: qnetworkb
    real :: costb
    real, dimension(2) :: hydrau_coef_bounds
    real, dimension(2) :: spreading_bounds
    
    integer :: loop,nbloop,iter
    
    character(400) :: costfilename
    logical :: qui
    
    integer    :: nn , mm
    integer :: iprint 
    integer,  parameter    :: dp = kind(1.0d0)
    real(dp)    :: factr , pgtol 
    character(len=60)      :: task, csave
    logical                :: lsave(4)
    integer                :: isave(44)
    real(dp)               :: f
    real(dp)               :: dsave(29)
    integer,dimension(:),  allocatable  :: nbd, iwa
    real(dp),dimension(:), allocatable  :: xx, l, u, g, wa
    integer                :: i
    
    nn = 2*routing_mesh%nb_nodes
    mm = 10

    allocate( nbd(nn), xx(nn), l(nn), u(nn), g(nn) )
    allocate ( iwa(3*nn) )
    allocate ( wa(2*mm*nn + 5*nn + 11*mm*mm + 8*mm) )

!   Declare a few additional variables for the sample problem.
!   General use arrays/parameters
!   MODEL VARIABLES ---
    
    !store initial parameters
    routing_results%initial_hydraulic_coef=routing_parameter%hydraulics_coefficient
    routing_results%initial_spreading=routing_parameter%spreading
    
    
    !--------------- Normalisation des paramètres et des bornes -----------------------
!~     routing_states%param_normalisation(1)=routing_setup%hydrau_coef_boundaries(2)
!~     routing_states%param_normalisation(2)=routing_setup%spreading_boundaries(2)
    
!~     routing_parameter%hydraulics_coefficient=routing_parameter%hydraulics_coefficient/routing_states%param_normalisation(1)
!~     routing_parameter%spreading=routing_parameter%spreading/routing_states%param_normalisation(2)
    
!~     spreading_bounds=routing_setup%spreading_boundaries/routing_setup%spreading_boundaries(2)
!~     hydrau_coef_bounds=routing_setup%hydrau_coef_boundaries/routing_setup%hydrau_coef_boundaries(2)
    ! --------------------------------------------------------------------------------
    
    call normalize_routing_parameters(routing_setup, routing_mesh, routing_parameter)
    hydrau_coef_bounds(1)=0.
    hydrau_coef_bounds(2)=1.
    spreading_bounds(1)=0.
    spreading_bounds(2)=1.
    
    call routing_parameter_self_initialisation(routing_parameter=routing_parameter_initial,routing_setup=routing_setup,&
    &routing_mesh=routing_mesh,hydraulics_coefficient=0.,spreading=0.)
    
    call routing_parameter_self_initialisation(routing_parameter=routing_parameterb,routing_setup=routing_setup,&
    &routing_mesh=routing_mesh,hydraulics_coefficient=0.,spreading=0.)
    
    allocate(qnetworkb(routing_setup%npdt,routing_mesh%nb_nodes))
    
    !initialisation des autres variables
    factr=1.e+1
    pgtol=1.e-14 
    iprint=1!01
    
    costfilename="out/cost.txt"
    inquire(file = trim("out/"),exist=qui)
    if(.not.qui) then
        call system("mkdir -p "//"out/")
    endif
    open(101,file=trim(costfilename))
    
    
    write(*,*) "nparam=",nn
    
    !write(*,*)routing_setup%hydrau_coef_boundaries
    !write(*,*)routing_setup%spreading_boundaries
    
    !pause
    
    !save the parameters
    routing_parameter_initial=routing_parameter
    
    if (routing_setup%auto_reg == 1) then 
        nbloop=2
        routing_setup%ponderation_regul=0.
    else
        nbloop=1
    end if
    
    do loop=1,nbloop
        
        routing_parameter=routing_parameter_initial
        
        if (routing_setup%auto_reg == 1) then
            if (loop.eq.2) then
                if (routing_results%costs(3)>0.) then
                    routing_setup%ponderation_regul=&
                    &routing_setup%ponderation_cost*routing_results%costs(2)/routing_results%costs(3)  ! ici c'est par rapport à un minimum local trouvé à l'aide d'un premier cylce avec weight_reg=0 (pareil qu' avec l'iterative regularisation'
                else
                    exit
                endif
            endif
        end if
        
        call affect_bound(hydrau_coef_bounds,spreading_bounds,routing_mesh,nbd,l,u)
        
        
        cost=0.
        
        call routing_states_reset(routing_states)
        call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,inflows,observations,&
        &routing_states,routing_results,cost)
        
        
        write(*,*) "At starting point, initial cost=",cost
        
        write(101,*) "cost"," J0"," BetaJreg"," Jreg","  NbCurIter","  TotalNbFuncEval","  NbFunEvalCurIter","  ProjGrad"
        write(101,*) cost,0,0,0,0
        
        xx=0.
        call linearise_control_vector(nn,routing_mesh,routing_parameter,xx)
        
        task = 'START'
        
        iter=1
        do while((task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
                &task.eq.'START')) 
             
            !call optimiseur
            call setulb(nn,mm,xx,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)
            
            !New x
            call unlinearise_control_vector(nn,routing_mesh,xx,routing_parameter)
            
            
             if (task(1:2) .eq. 'FG') then
                
                call routing_parameter_self_initialisation(routing_parameter=routing_parameterb,&
                &routing_setup=routing_setup,routing_mesh=routing_mesh,hydraulics_coefficient=0.,spreading=0.)
                
                cost=0.
                costb=1.
                call routing_states_reset(routing_states)
                
                call routing_hydrogram_forward_b(routing_setup, &
                &   routing_mesh, routing_parameter, routing_parameterb, inflows, &
                &   observations, routing_states, routing_results, cost, &
                &   costb)
                
                !store all gradients
                routing_results%allgrads_hydraulic_coef(iter,:)=routing_parameterb%hydraulics_coefficient
                routing_results%allgrads_spreading(iter,:)=routing_parameterb%spreading
                
                routing_results%gradients_hydraulic_coef=routing_parameterb%hydraulics_coefficient
                routing_results%gradients_spreading=routing_parameterb%spreading
                
                write(*,*) "New Cost = ",cost
                
                !Direct cost function evaluation
                f=real(cost,kind(f))
                
                !Gradient
                call make_gradient_vector(nn,routing_mesh,routing_parameterb,g)
                
            endif
            

            if (task(1:5) .eq. 'NEW_X') then   

                write(101,*) cost,&
                &isave(30),isave(34),isave(36),dsave(13)
                
                iter=iter+1
    !    
    !       the minimization routine has returned with a new iterate.
    !       At this point have the opportunity of stopping the iteration 
    !       or observing the values of certain parameters
    !
    !       First are two examples of stopping tests.

    !       Note: task(1:4) must be assigned the value 'STOP' to terminate  
    !         the iteration and ensure that the final results are
    !         printed in the default format. The rest of the character
    !         string TASK may be used to store other information.

    !       1) Terminate if the total number of f and g evaluations
    !            exceeds 99.

                 if (isave(34) .ge. 1000)  &
                    task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'

    !       2) Terminate if  |proj g|/(1+|f|) < 1.0d-10, where 
    !          "proj g" denoted the projected gradient

                 if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f))) &
                   task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
    
    !       3)Terminate if iter >= iter_max
                 if (iter>=routing_setup%iter_max) then
                    task='STOP: TOTAL ITERATIONS EXCEEDS LIMIT'
                 end if
                 
    !       We now wish to print the following information at each
    !       iteration:
    !       
    !         1) the current iteration number, isave(30),
    !         2) the total number of f and g evaluations, isave(34),
    !         3) the value of the objective function f,
    !         4) the norm of the projected gradient,  dsve(13)
    !
    !       See the comments at the end of driver1 for a description
    !       of the variables isave and dsave.
             
                write (6,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate' &
                   , isave(30),'nfg =',isave(34),'f =',f,'|proj g| =',dsave(13)

    !        If the run is to be terminated, we print also the information
    !        contained in task as well as the final value of x.

                if (task(1:4) .eq. 'STOP') then
                   write (6,*) task  
                   write (6,*) 'Final X='
                   write (6,'((1x,1p, 6(1x,d11.4)))') (xx(i),i = 1,2*routing_mesh%nb_nodes)
                end if
            
            end if
            
        enddo
        
    end do
    !New x
    call unlinearise_control_vector(nn,routing_mesh,xx,routing_parameter)
    
    
    call routing_states_reset(routing_states)
    call routing_hydrogram_forward(routing_setup,routing_mesh,routing_parameter,inflows,observations,&
    &routing_states,routing_results,cost)
    
    
    write(101,*) cost,&
    &isave(30),isave(34),isave(36),dsave(13)
    close(101)
    
    ! ---------------- Un-Normalisation des paramètres -------------------------------
!~     routing_parameter%hydraulics_coefficient=routing_parameter%hydraulics_coefficient*routing_states%param_normalisation(1)
!~     routing_parameter%spreading=routing_parameter%spreading*routing_states%param_normalisation(2)
!~     routing_states%param_normalisation(1)=1.
!~     routing_states%param_normalisation(2)=1.
    ! --------------------------------------------------------------------------------
    call unnormalize_routing_parameters(routing_setup, routing_mesh, routing_parameter)
    
    !store final parameters
    routing_results%final_hydraulic_coef=routing_parameter%hydraulics_coefficient
    routing_results%final_spreading=routing_parameter%spreading
    
    write(*,*) "--------- Final estimation ----------"
    write(*,*) ""
    write(*,*) "Hydraulics parameters:"
    write (*,'((1x,1p, 6(1x,d11.4)))') (routing_parameter%hydraulics_coefficient(i),i = 1,routing_mesh%nb_nodes)
    write(*,*) "spreading coefficients:"
    write (*,'((1x,1p, 6(1x,d11.4)))') (routing_parameter%spreading(i),i = 1,routing_mesh%nb_nodes)
    write(*,*) ""
    write(*,*) "-------------------------------------"
    
    
    contains
    
    !subroutine starting_point
    !linearise param_dist vers la variable X en focntion des paramètres à optimiser
    !input/output : nbd,l,u
    ! input z=nb params or states
    subroutine affect_bound(hydrau_coef_bounds,spreading_bounds,routing_mesh,nbd,l,u)
        implicit none
        real,dimension(2), intent(in) :: hydrau_coef_bounds
        real,dimension(2), intent(in) :: spreading_bounds
        type(type_routing_mesh), intent(in) :: routing_mesh
        integer,dimension(routing_mesh%nb_nodes*2) :: nbd
        real*8,dimension(routing_mesh%nb_nodes*2) :: l,u
        
        integer :: p,i!,k,i,j,z
        nbd=0
        l=0.
        u=1.
        
        p=1
        do i=1,routing_mesh%nb_nodes
            nbd(p)=2
            l(p)=real(hydrau_coef_bounds(1),kind(l))! lower
            u(p)=real(hydrau_coef_bounds(2),kind(u))! upper
            
            nbd(routing_mesh%nb_nodes+p)=2
            l(routing_mesh%nb_nodes+p)=real(spreading_bounds(1),kind(l))! lower
            u(routing_mesh%nb_nodes+p)=real(spreading_bounds(2),kind(u))! upper
            
            p=p+1
        enddo
        
    end subroutine affect_bound
    
    subroutine linearise_control_vector(nmax,routing_mesh,routing_parameter,x)
        implicit none
        integer :: nmax ! maximum size of the control vector
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameter
        double precision,dimension(nmax), intent(inout) :: x
        
        x=real(0.,kind(x))
        x(1:routing_mesh%nb_nodes)=real(routing_parameter%hydraulics_coefficient,kind(x))
        x(routing_mesh%nb_nodes+1:2*routing_mesh%nb_nodes)=real(routing_parameter%spreading,kind(x))
        
    endsubroutine linearise_control_vector
    
    subroutine unlinearise_control_vector(nmax,routing_mesh,x,routing_parameter)
        implicit none
        integer :: nmax ! maximum size of the control vector
        type(type_routing_mesh), intent(in) :: routing_mesh
        double precision,dimension(nmax), intent(in) :: x
        type(type_routing_parameter), intent(inout) :: routing_parameter
        
        routing_parameter%hydraulics_coefficient=real(x(1:routing_mesh%nb_nodes),kind(routing_parameter%hydraulics_coefficient))
        routing_parameter%spreading=real(x(routing_mesh%nb_nodes+1:2*routing_mesh%nb_nodes),&
        &kind(routing_parameter%spreading))
        
    endsubroutine unlinearise_control_vector
    
    subroutine make_gradient_vector(nmax,routing_mesh,routing_parameterb,g)
        implicit none
        integer :: nmax ! maximum size of the control vector
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_parameter), intent(in) :: routing_parameterb
        double precision,dimension(nmax), intent(inout) :: g
        
        g=real(0.,kind(g))
        g(1:routing_mesh%nb_nodes)=real(routing_parameterb%hydraulics_coefficient,kind(g))
        g(routing_mesh%nb_nodes+1:2*routing_mesh%nb_nodes)=real(routing_parameterb%spreading,kind(g))
        
    endsubroutine make_gradient_vector
    
    
endsubroutine control

