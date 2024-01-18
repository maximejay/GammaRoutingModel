!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr


module mod_gamma_routing_setup
    implicit none
    
    
    type type_routing_setup
        integer :: npdt=1 !nb de pas de temps
        real :: dt = 3600. !time-step en seconds (s)
        real :: vmin = 0.1 !minimum velocity (m/s)
        real :: vmax = 10.0 !maximum velocity (m/s)
        real :: spreading_gamma_threshold=0.1 !objective value for the ratio Gamma(x+delta)/Gamma(x) where delta=spreading/dt and x the time-step position: spreading_gamma_threshold=0.1 is a good choice since that respect the differentiability condition and make the damping parameter the accurate for the modeller
        real :: elongation_factor=1.0 ! Time space deformation to reduce the computation time 1<=f<=2. Good values should not exeed 1.2
        real :: mode_discretization_step=0.1 !Discretization of the mode of the Gamma PDF
        real :: spreading_discretization_step=0.1 ! Discretization of the spreading of the Gamma PDF
        real,dimension(2) :: hydrau_coef_boundaries !depend on vmin,vmax
        real,dimension(2) :: spreading_boundaries !min_spreading=dt/dx, max_spreading=?
        real :: ponderation_regul=0.0 !factor of the regularisation
        real :: ponderation_cost=1.0 !factor of the regularisation
        integer :: pdt_start_optim=1 ! time step number to start the computation of the cost function
        integer :: auto_reg=0 ! approximate the lcurve, compute the best value for ponderation_regul
        integer :: hydraulics_coef_uniform=0
        integer :: spreading_uniform=0
        integer :: iter_max=100 ! maximum number of iteraration during the calibration
        character(3) :: velocity_computation="qmm" ! mm or m3/s : qmm or qm3
        character(4) :: criteria="rmse" ! mm or m3/s : qmm or qm3
        logical :: varying_spread=.false. !fixed or varying spreading coefficient
    end type type_routing_setup
    
    
    contains
    
    subroutine routing_setup_self_initialisation(routing_setup,npdt,dt,vmin,vmax,&
    &elongation_factor,mode_discretization_step,spreading_discretization_step,&
    &ponderation_regul,ponderation_cost,auto_reg,hydraulics_coef_uniform,spreading_uniform,&
    &iter_max,velocity_computation,criteria,varying_spread)
    
        ! Notes
        ! -----
        ! **routing_setup_self_initialisation** :
        !
        ! - Initialise the routing_setup derived type with user values
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (inout)
        ! ``npdt=1``                              number of time step (optional)
        ! ``dt=3600.``                            Time step in seconds (optional)
        ! ``vmin=0.1``                            Minimum flow speed in the networl (optional)
        ! ``vmax=10.``                            Maximum flow speed in the networl (optional)
        ! ``elongation_factor=1``                 Factor 1<=f<=2 to speed up the computation by deforming the time space (optional)
        ! ``mode_discretization_step=0.1``        Discretization step size for the mode of the Gamma PDF (optional)
        ! ``spreading_discretization_step=300.``  Discretization step size for the spreading (in seconds) of the Gamma PDF (optional)
        ! ``velocity_computation="qmm"``          Control how to compute the velocity. use the specific discharge "qmm" or the discharge "qm3" (optional) 
        ! ``criteria="nse"``                      Criteria for optimization nse or rmse (optional) 
        ! ``varying_spread=.false.``              Fixed or varying spreading coefficient
        ! =============================           ===================================
        
        implicit none
        
        type(type_routing_setup), intent(inout) :: routing_setup
        integer, optional :: npdt
        real, optional :: dt
        real, optional :: vmin
        real, optional :: vmax
        real, optional :: elongation_factor
        real, optional :: mode_discretization_step
        real, optional :: spreading_discretization_step
        real, optional :: ponderation_regul
        real, optional :: ponderation_cost
        integer, optional :: auto_reg
        integer, optional :: hydraulics_coef_uniform
        integer, optional :: spreading_uniform
        integer, optional :: iter_max
        character(3), optional :: velocity_computation
        character(4), optional :: criteria
        logical,optional :: varying_spread
        
        
        if (present(dt)) then
            routing_setup%npdt=npdt
        end if
        if (present(dt)) then
            routing_setup%dt=dt
        end if
        if (present(vmin)) then
            routing_setup%vmin=vmin
        end if
        if (present(vmax)) then
            routing_setup%vmax=vmax
        end if
        if (present(elongation_factor)) then 
            routing_setup%elongation_factor=elongation_factor
        end if
        if (present(mode_discretization_step)) then 
            routing_setup%mode_discretization_step=mode_discretization_step
        end if
        if (present(spreading_discretization_step)) then 
            routing_setup%spreading_discretization_step=spreading_discretization_step
            
            if (spreading_discretization_step>dt .and. routing_setup%varying_spread) then
                write(*,*) "spreading_discretization_step have to be lower then dt !"
                stop 0
            end if
            
        end if
        if (present(ponderation_regul)) then 
            routing_setup%ponderation_regul=ponderation_regul
        end if
        if (present(ponderation_cost)) then 
            routing_setup%ponderation_cost=ponderation_cost
        end if
        if (present(auto_reg)) then 
            routing_setup%auto_reg=auto_reg
        end if
        if (present(hydraulics_coef_uniform)) then 
            routing_setup%hydraulics_coef_uniform=hydraulics_coef_uniform
        end if
        if (present(spreading_uniform)) then 
            routing_setup%spreading_uniform=spreading_uniform
        end if
        if (present(iter_max)) then 
            routing_setup%iter_max=iter_max
        end if
        if (present(velocity_computation)) then 
            routing_setup%velocity_computation=velocity_computation
        end if
        if (present(criteria)) then 
            routing_setup%criteria=criteria
        end if
        if (present(varying_spread)) then 
            routing_setup%varying_spread=varying_spread
        end if
        
        !bounds auto calculation, can be changed after initialisation
        routing_setup%hydrau_coef_boundaries(1)=0.1
        routing_setup%hydrau_coef_boundaries(2)=10.
        routing_setup%spreading_boundaries(1)=0.5
!~         routing_setup%spreading_boundaries(2)=10800.
!~         routing_setup%spreading_boundaries(2)=routing_setup%dt+300.*sqrt(routing_setup%dt)
        routing_setup%spreading_boundaries(2)=5.
!~         routing_setup%spreading_boundaries(2)=routing_setup%dt/1000.+0.3*sqrt(routing_setup%dt)
        
    end subroutine routing_setup_self_initialisation
    
    
    subroutine routing_setup_clear(routing_setup)
        
        ! Notes
        ! -----
        ! **routing_setup_clear** :
        !
        ! - Clear the derived type routing_setup
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (inout)
        ! =============================           ===================================
        
        type(type_routing_setup), intent(inout) :: routing_setup
        
        type(type_routing_setup) :: routing_setup_new
        
        routing_setup=routing_setup_new
        
    endsubroutine routing_setup_clear
    
    
    subroutine routing_setup_copy(routing_setup,object_copy)
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_setup), intent(out) :: object_copy
        
        object_copy=routing_setup
    
    endsubroutine routing_setup_copy
    
end module mod_gamma_routing_setup
