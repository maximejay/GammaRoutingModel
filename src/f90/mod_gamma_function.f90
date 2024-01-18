!~ GammaRouting is a conceptual flow propagation model
!~ Copyright 2022, 2023 Hydris-hydrologie, Maxime Jay-Allemand

!~ This file is part of GammaRouting.

!~ GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

!~ GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

!~ You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

!~ Contact: maxime.jay.allemand@hydris-hydrologie.fr



module mod_gamma_function
    
    use mod_gamma_routing_setup
    use mod_gamma_routing_mesh
    use mod_gamma_routing_parameters
    use mod_gamma_routing_states
    
    implicit none
    
    
    contains
    
    
    subroutine compute_gamma_parameters(routing_setup,routing_mesh,routing_states)
        
        ! Notes
        ! -----
        ! **compute_gamma_parameters(routing_setup,routing_mesh,routing_states)** :
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
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        
        call compute_gamma_scale(routing_setup,routing_mesh,routing_states)
        call compute_gamma_window(routing_setup,routing_mesh,routing_states)
        call gamma_elongation_cells(routing_setup,routing_mesh,routing_states)
        
        if (allocated(routing_states%remainder)) then
            deallocate(routing_states%remainder)
        end if
        if (allocated(routing_states%states)) then
            deallocate(routing_states%states)
        end if
        allocate(routing_states%remainder(maxval(routing_states%window_length),routing_mesh%nb_nodes))
        allocate(routing_states%states(maxval(routing_states%window_length),routing_mesh%nb_nodes))
        
        !default value
        routing_states%remainder=0.
        routing_states%states=0.
        
        !make the tabulated gamma routing coefficient
        call tabulated_routing_coefficients_3D(routing_setup,routing_mesh,routing_states,"cdf")
        
    end subroutine compute_gamma_parameters
    
    
    subroutine compute_gamma_scale(routing_setup,routing_mesh,routing_states)
        
        ! Notes
        ! -----
        ! **compute_gamma_scale(routing_setup,routing_mesh,routing_states)** :
        !
        ! - Compute the scale coefficient of the Gamma PDF
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! =============================           ==================================
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
                
        !should be array of nbnodes
        real :: mu_min
        real :: mu_min_gamma
        real :: rate
        real :: delta_x
        integer :: i
        
        
        !étalement du rapport : Gamma(x+delta_x)/Gamma(x), deltat_x est en pas de temps, spread en  seconds/m et dt en seconds
        !valid for varying dx
        do i=1,size(routing_mesh%varying_dx)
        
            delta_x=routing_states%max_spreading*routing_mesh%varying_dx(i)/routing_setup%dt
            
            mu_min=routing_mesh%varying_dx(i)/(routing_setup%vmax*routing_setup%dt)
            
            mu_min_gamma=mu_min+1.0 !On caclul Gamma sur l'interval MU E [1:L+1]  et on calcul les coefficients sur Tau E [0:L]
            
            rate=(log(routing_setup%spreading_gamma_threshold)/(mu_min_gamma*log((mu_min_gamma+delta_x)/mu_min_gamma)-delta_x)) !Attention scale=1/Rate Rate=Beta
            
            if (1./rate < 0.1) then
                !this condition is reached for spreading <= dt
                write(*,*) "warnings : scale<0.1, scale=",1./rate
                write(*,*) "This could produce to calibration issues, such as a null value of the gradient."
            endif
        
            routing_states%scale_coef(i)=1./rate
        end do
        
        !valid for spreading cst and dx cst only
!~         delta_x=routing_states%max_spreading*minval(routing_mesh%dx)/routing_setup%dt
!~         !delta_x=routing_states%max_spreading*minval(routing_mesh%dx)**0.9/routing_setup%dt
        
!~         mu_min=minval(routing_mesh%dx)/(routing_setup%vmax*routing_setup%dt)
!~         mu_min_gamma=mu_min+1.0 !On caclul Gamma sur l'interval MU E [1:L+1]  et on calcul les coefficients sur Tau E [0:L]
!~         rate=(log(routing_setup%spreading_gamma_threshold)/(mu_min_gamma*log((mu_min_gamma+delta_x)/mu_min_gamma)-delta_x)) !Attention scale=1/Rate Rate=Beta
!~         !routing_states%scale=max(0.1,1./rate)
        
!~         if (1./rate < 0.1) then
!~             !this condition is reached for spreading <= dt
!~             write(*,*) "warnings : scale<0.1, scale=",1./rate
!~             write(*,*) "This could produce to calibration issues, such as a null value of the gradient."
!~         endif
        
!~         routing_states%max_scale=1./rate
        
    end subroutine compute_gamma_scale
    
    
    
    !computation of the scale coefficient
    subroutine generic_compute_gamma_scale(dx,dt,vmax,spread,epsilon,scale)
        
        ! Notes
        ! -----
        ! **generic_compute_gamma_scale(dx,dt,vmax,spread,epsilon,scale)** :
        !
        ! - Generic version of compute_gamma_scale : Compute the scale coefficient of the Gamma PDF
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``dx``                                  Maximum distance between 2 nodes, real (in)
        ! ``dt``                                  Time-step, real (in)
        ! ``vmax``                                Maximum velocity, real (in)
        ! ``spread``                              Maximum spreading coefficient (in)
        ! ``epsilon``                             Threshold so that Gamma pdf < epsilon for a the quantile mode+1 (in)
        ! ``scale``                               The computed scale coefficient (inout)
        ! =============================           ==================================
        
        implicit none
        
        real,intent(in) :: dx !discretization-space in meter (m)
        real,intent(in) :: dt !time-step en seconds (s)
        real,intent(in) :: vmax !maximum velocity (m/s)
        real,intent(in) :: spread !damping coefficient in seconds (s): spreading of the Gamma law
        real,intent(in) :: epsilon !objective value for the ratio Gamma(x+delta)/Gamma(x) where delta=spread/dt and x the time-step position: epsilon=0.1 is a good choice since that respect the differentiability condition and make the damping parameter the accurate for the modeller
        
        real,intent(out) :: scale
        
        real :: mu_min
        real :: mu_min_gamma
        real :: rate
        real :: delta_x
        
        !étalement du rapport : Gamma(x+delta_x)/Gamma(x), deltat_x est en pas de temps, spread en  seconds/m et dt en seconds
        delta_x=spread*dx/dt
        !delta_x=spread*dx**0.9/dt ! correction
        
        mu_min=dx/(vmax*dt)
        
        mu_min_gamma=mu_min+1.0 !On caclul Gamma sur l'interval MU E [1:L+1]  et on calcul les coefficients sur Tau E [0:L]
        rate=(log(epsilon)/(mu_min_gamma*log((mu_min_gamma+delta_x)/mu_min_gamma)-delta_x)) !Attention scale=1/Rate Rate=Beta
        
!~         if (1./rate < 0.1) then
!~             !this condition is reached for spreading <= dt
!~             write(*,*) "warnings : scale<0.1, scale=",1./rate
!~             write(*,*) "This could produce to calibration issues, such as a null value of the gradient."
!~         endif
        
!~         scale=max(0.1,1./rate)
        scale=1./rate
    
    end subroutine generic_compute_gamma_scale
    
    
    !computation of the window length to compute the gamma coefficient
    subroutine compute_gamma_window(routing_setup,routing_mesh,routing_states)
        
        ! Notes
        ! -----
        ! **compute_gamma_window(routing_setup,routing_mesh,routing_states)** :
        !
        ! - Compute the the Gamma window lenght, i.e the number of quantile fo calculate the Gamma Pdf/Cdf
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! =============================           ==================================
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        
        real :: spreading
        integer :: i,j
        integer :: index_scale_coef
        real :: scale
        
        !spreading=routing_states%max_spreading*minval(routing_mesh%dx)/routing_setup%dt !spread is the spreading of the gamma law in  seconds / m
        !spreading=routing_states%max_spreading*minval(routing_mesh%dx)**0.9/routing_setup%dt !spread is the spreading of the gamma law in seconds
        
        do j=1,routing_mesh%nb_nodes
        
            index_scale_coef=routing_mesh%index_varying_dx(j)
            
            spreading=routing_states%max_spreading*routing_mesh%dx(j)/routing_setup%dt !spread is the spreading of the gamma law in  seconds / m
            
            
            
            !we could compute the scale coef for each nodes, then the fun compute_gamma_scale will becom unuseful
            !------------------------------------------------------------------
            !call generic_compute_gamma_scale(routing_mesh%dx,routing_setup%dt,routing_setup%vmax,spreading,&
            !&routing_setup%spreading_gamma_threshold,scale)
            
            !routing_states%window_length(j)=ceiling( routing_mesh%dx(j)/(routing_setup%vmin*routing_setup%dt) &
            !& + 2/skweness_gamma( (routing_mesh%dx(j)/(routing_setup%vmin*routing_setup%dt))+1.0,&
            !&scale )* max(spreading,1.0) )  !Longueur de la fenetre d'itnertie du modèle en nb pdt
            !------------------------------------------------------------------
            
            
            routing_states%window_length(j)=ceiling( routing_mesh%dx(j)/(routing_setup%vmin*routing_setup%dt) &
            & + 2/skweness_gamma( (routing_mesh%dx(j)/(routing_setup%vmin*routing_setup%dt))+1.0,&
            &routing_states%scale_coef(index_scale_coef) )&
            & * max(spreading,1.0) )  !Longueur de la fenetre d'itnertie du modèle en nb pdt
            
        end do
        
!~         routing_states%window_shift=minval(routing_mesh%dx)/(routing_setup%vmax*routing_setup%dt)&
!~         &-floor(minval(routing_mesh%dx)/(routing_setup%vmax*routing_setup%dt))
        
        routing_states%window_shift=minval(routing_mesh%dx)/(routing_setup%vmax*routing_setup%dt)
        
        
        
        if (allocated(routing_states%quantile)) then
            deallocate(routing_states%quantile)
        endif
        allocate(routing_states%quantile(maxval(routing_states%window_length)))
        
        
        do i=1,maxval(routing_states%window_length)
            routing_states%quantile(i)=i
        end do
        
    end subroutine compute_gamma_window
    
    
    
    !computation of the window length to compute the gamma coefficient
    subroutine generic_compute_gamma_window(dx,dt,vmin,vmax,scale,spread,window_length,window_shift,quantile)
        
        ! Notes
        ! -----
        ! **generic_compute_gamma_window(dx,dt,vmax,spread,epsilon,scale)** :
        !
        ! - Generic version of compute_gamma_window : Compute the windows lenght, i.e the number of quantile where to compute the Gamma PDF/CDF
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``dx``                                  Maximum distance between 2 nodes, real (in)
        ! ``dt``                                  Time-step, real (in)
        ! ``vmin``                                Minimum velocity, real (in)
        ! ``vmax``                                Maximum velocity, real (in)
        ! ``scale``                               Scale coefficient of the Gamma pdf, real (in)
        ! ``spread``                              Maximum spreading coefficient, real (in)
        ! ``window_length``                       Window length equal to the number of the quantile, integer (out)
        ! ``window_shift``                        Window shift, shift the quantile so that the minimum quantile is equal to the minimum mode, real (out)
        ! ``quantile``                            The computed quantile, array(window_length), real (out)
        ! =============================           ==================================
        
        implicit none
        
        real,intent(in) :: dx !discretization-space in meter (m)
        real,intent(in) :: dt !time-step en seconds (s)
        real,intent(in) :: vmin !minimum velocity (m/s)
        real,intent(in) :: vmax !maximum velocity (m/s)
        real,intent(in) :: scale !scale parameter of the Gamma law: scale has to be computed according the spread coefficient
        real,intent(in) :: spread !damping coefficient in seconds (s/m): spreading of the Gamma law, used to compute the scale coefficient
        integer,intent(out) :: window_length !window lenght in time-step for the computation of the Gamma pdf
        !real,intent(out) :: shift !shift the 
        real,intent(out) :: window_shift !shift the window so that the 1st quantile correspond to the peak of the pdf when vmax is reached 
        real,dimension(:),allocatable,intent(out) :: quantile !quantile for the computation of the gamma function
        
        real :: spreading
        integer :: i
        
        spreading=spread/dt !spread is the spreading of the gamma law in seconds/m
        window_length=ceiling(dx/(vmin*dt) + 2/skweness_gamma((dx/(vmin*dt))+1,scale) * max(spreading,1.0) )  !Longueur de la fenetre d'itnertie du modèle en nb pdt
        
        !window_shift=dx/(vmax*dt)-floor(dx/(vmax*dt))
        window_shift=dx/(vmax*dt)
        
        allocate(quantile(window_length))
        do i=1,window_length
            quantile(i)=i
        end do
        
    end subroutine generic_compute_gamma_window
    
    
    !computation of the new_window_length, the number of memory cells and the adjusted_quantile according the elongation_coefficient
    subroutine gamma_elongation_cells(routing_setup,routing_mesh,routing_states)
        
        ! Notes
        ! -----
        ! **gamma_elongation_cells(routing_setup,routing_mesh,routing_states)** :
        !
        ! - recompute window length and update routing_states if the elongation coefficient > 1
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        Routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! =============================           ==================================
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        
        integer ::i,j
        integer :: new_windows_length
        real :: window_length
        real :: quantile
        
        do j=1,routing_mesh%nb_nodes
            
            new_windows_length=1
            window_length=0.0
            
            do i=1,routing_states%window_length(j)
            
                window_length=window_length+routing_setup%elongation_factor**(i-1)
                
                if (window_length>=routing_states%window_length(j)) then
                    new_windows_length=i
                    exit
                end if
                
            end do
            
            routing_states%window_length(j)=new_windows_length
            
        end do
        
        if (allocated(routing_states%quantile)) then
            deallocate(routing_states%quantile)
        endif
        allocate(routing_states%quantile(maxval(routing_states%window_length)))
        
        quantile=0.0
        do i=1,maxval(routing_states%window_length)
        
            quantile=quantile+routing_setup%elongation_factor**(i-1)
            routing_states%quantile(i)=quantile
            
        end do
        
    end subroutine gamma_elongation_cells
    
    
    !computation of the new_window_length, the number of memory cells and the adjusted_quantile according the elongation_coefficient
    subroutine generic_gamma_elongation_cells(initial_window_length,elongation_coefficient,new_windows_length,&
    &adjusted_quantile)
        
        ! Notes
        ! -----
        ! **generic_gamma_elongation_cells(initial_window_length,elongation_coefficient,new_windows_length,adjusted_quantile)** :
        !
        ! - Generic version of gamma_elongation_cells : recompute window length and update routing_states if the elongation coefficient > 1
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``initial_window_length``               initial_window_length, integer (in)
        ! ``elongation_coefficient``              elongation_coefficient, real (in)
        ! ``new_windows_length``                  new_windows_length, integer (in)
        ! ``adjusted_quantile``                   New quantile where to compute the Gamma pdf/cdf, array(new_windows_length), real (in)
        ! =============================           ==================================
        
        implicit none
        integer,intent(in) :: initial_window_length
        real,intent(in) :: elongation_coefficient
        integer,intent(out) :: new_windows_length
        real,dimension(:),allocatable,intent(out) :: adjusted_quantile
        
        integer ::i
        real :: window_length
        real :: quantile
        
        window_length=0
        do i=1,initial_window_length
            window_length=window_length+elongation_coefficient**(i-1)
            if (window_length>=initial_window_length) then
                new_windows_length=i
                exit
            end if
        end do
        
        allocate(adjusted_quantile(new_windows_length))
        
        quantile=0
        do i=1,new_windows_length
            quantile=quantile+elongation_coefficient**(i-1)
            adjusted_quantile(i)=quantile
        end do
    end subroutine generic_gamma_elongation_cells
    
    
    !computation of the skweness of the gamma function
    !The skweness is used to evaluate the window length
    !inputs:
    !mode: mode of the Gamma pdf
    !scale: scale coefficient of the gamma pdf
    !return:
    !skweness_gamma, a real
    function skweness_gamma(mode,scale)
        real :: skweness_gamma
        real :: mode
        real :: scale
        
        skweness_gamma=(2/sqrt(mode/scale + 1))
    end function skweness_gamma
    
    
    !computation of the gamma coefficient
    subroutine compute_gamma_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
        
        ! Notes
        ! -----
        ! **compute_gamma_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)** :
        !
        ! - compute the gamma coefficient for each quantile
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``scale``                               Scale coefficient of the Gamma pdf/cdf, real (in)
        ! ``mode``                                mode of the Gamma pdf/cdf, real (in)
        ! ``quatile``                             The quantile where to compute the gamma pdf/cdf, array(window_lenght), real (in)
        ! ``windows_shift``                       The shift value to shift the quantile so that the minimum mode equal to the minimum quantile
        ! ``density_function``                    use the 'pdf' or the 'cdf', string, prefer the cdf (numerical issues with the pdf)
        ! ``gamma_coefficient``                   The value of the Gamma pdf/cdf for each quantile, array(windows_length), real (in)
        ! =============================           ==================================
        
        implicit none
        real,intent(in) :: scale
        real,intent(in) :: mode
        real,dimension(:),intent(in) :: quantile
        real,intent(in) :: window_shift
        character(3),intent(in) :: density_function
        real,dimension(size(quantile)),intent(out) :: gamma_coefficient
        
        
        integer :: window_length
        integer :: i
        
        real*8 :: alpha,beta
        real :: center_left,center_right
        
        window_length=size(quantile)
        
        alpha=dble(1.0+(mode+1.0)/scale) !position parameter ot the pdf, here the pdf is centered on its mode. The mode is shifted +1
        beta=dble(1./scale) !rate coefficient of the gamma pdf, rate=1/scale
        
        gamma_coefficient=0.0
        do i=1,window_length
            if (density_function=="pdf") then
                gamma_coefficient(i)=GammaPDF(dble(quantile(i)+window_shift), alpha, beta) !here the quantile are shifted +shift to be sure that the  1st quantile correspond to the peak of the ditribution if vmax is reached
            end if
            
            if (density_function=="cdf") then
                if (i==1) then
                    center_left=1.0*(quantile(i))
                    center_right=0.5*(quantile(i+1)-quantile(i))
                end if
                if (i>1 .and. i<window_length) then
                    center_left=0.5*(quantile(i)-quantile(i-1))
                    center_right=0.5*(quantile(i+1)-quantile(i))
                endif
                if (i==window_length) then
                    center_left=0.5*(quantile(i)-quantile(i-1))
                    center_right=0.5*(1.5*(quantile(i)-quantile(i-1)))
                end if
                gamma_coefficient(i)=GammaCDF(dble(quantile(i)+window_shift+center_right), alpha, beta)-&
                &GammaCDF(dble(quantile(i)+window_shift-center_left), alpha, beta) !here the quantile are shifted +shift to be sure that the  1st quantile correspond to the peak of the ditribution if vmax is reached
            end if
        end do
        
    end subroutine compute_gamma_coefficient
    
    
    !computation of the routing coefficient
    subroutine compute_gamma_routing_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
        
        ! Notes
        ! -----
        ! **compute_gamma_routing_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)** :
        !
        ! - compute the gamma routing coefficient for each quantile
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``scale``                               Scale coefficient of the Gamma pdf/cdf, real (in)
        ! ``mode``                                mode of the Gamma pdf/cdf, real (in)
        ! ``quatile``                             The quantile where to compute the gamma pdf/cdf, array(window_lenght), real (in)
        ! ``windows_shift``                       The shift value to shift the quantile so that the minimum mode equal to the minimum quantile
        ! ``density_function``                    use the 'pdf' or the 'cdf', string, prefer the cdf (numerical issues with the pdf)
        ! ``gamma_coefficient``                   The value of the routing coefficient based on the Gamma pdf/cdf for each quantile, array(windows_length), real (in)
        ! =============================           ==================================
        
        
        implicit none
        real,intent(in) :: scale
        real,intent(in) :: mode
        real,dimension(:),intent(in) :: quantile
        real,intent(in) :: window_shift
        character(3),intent(in) :: density_function
        real,dimension(size(quantile)),intent(out) :: gamma_coefficient
        
        real,dimension(size(quantile)) :: gamma_values
        integer :: i
        integer :: window_length
        
        real*8 :: alpha,beta
        real :: center_left,center_right
        
        window_length=size(quantile)
        
        alpha=dble(1.0+(mode+1.0)/scale) !position parameter ot the pdf, here the pdf is centered on its mode. The mode is shifted +1
        beta=dble(1./scale) !rate coefficient of the gamma pdf, rate=1/scale
        
    
        gamma_values=0.0
        do i=1,window_length
            
            if (density_function=="pdf") then
                gamma_values(i)=GammaPDF(dble(quantile(i)+window_shift), alpha, beta) !here the quantile are shifted +shift to be sure that the  1st quantile correspond to the peak of the ditribution if vmax is reached
            end if
            if (density_function=="cdf") then
                if (i==1) then
                    center_left=1.0*(quantile(i))
                    center_right=0.5*(quantile(i+1)-quantile(i))
                end if
                if (i>1 .and. i<window_length) then
                    center_left=0.5*(quantile(i)-quantile(i-1))
                    center_right=0.5*(quantile(i+1)-quantile(i))
                endif
                if (i==window_length) then
                    center_left=0.5*(quantile(i)-quantile(i-1))
                    center_right=0.5*(1.5*(quantile(i)-quantile(i-1)))
                end if
                gamma_values(i)=GammaCDF(dble(quantile(i)+window_shift+center_right), alpha, beta)-&
                &GammaCDF(dble(quantile(i)+window_shift-center_left), alpha, beta) !here the quantile are shifted +shift to be sure that the  1st quantile correspond to the peak of the ditribution if vmax is reached
            end if
            
        end do
        
        gamma_coefficient=gamma_values/sum(gamma_values)
        
    end subroutine compute_gamma_routing_coefficient
    
    
    subroutine write_coefficients(coefficients,quantile,filename)
        real, dimension(:), intent(in) :: coefficients
        real, dimension(:), intent(in) :: quantile
        character(300),intent(in) :: filename
        
        integer :: i,n
        
        n=size(coefficients,1)
        open(10,file=trim(filename),status='replace')
            do i=1,n
                write(10,*) quantile(i),coefficients(i) 
            end do
        close(10)
        
    end subroutine write_coefficients
    
    
    subroutine tabulated_routing_coefficients_2D(routing_setup,routing_mesh,routing_states,density_function)
        
        ! Notes
        ! -----
        ! **tabulated_routing_coefficients_2D(routing_setup,routing_mesh,routing_states,density_function)** :
        !
        ! - Compute the tabulated routing coefficient for a range of mode
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! ``density_function``                    String 'pdf' or 'cdf, please use 'cdf' (in)
        ! =============================           ==================================
        
        implicit none
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        character(3),intent(in) :: density_function
        
        integer :: i,j
        real :: mode
        real, dimension(size(routing_states%quantile)) :: gamma_coefficient
        
        if (allocated(routing_states%tabulated_routing_coef)) then
            deallocate(routing_states%tabulated_routing_coef)
        end if
        
        allocate(routing_states%tabulated_routing_coef(size(routing_states%quantile),&
        &routing_states%nb_mode,1,size(routing_mesh%varying_dx)))
        routing_states%tabulated_routing_coef=0.0
        
        gamma_coefficient=0.
        mode=0
        
        do j=1,size(routing_mesh%varying_dx)
        
            do i=1,routing_states%nb_mode
                
                call compute_gamma_routing_coefficient(routing_states%scale_coef(j),mode,routing_states%quantile,&
                &routing_states%window_shift,density_function,gamma_coefficient)
                
                routing_states%tabulated_routing_coef(:,i,1,j)=gamma_coefficient
                mode=mode+routing_setup%mode_discretization_step
                
            end do
        
        end do
        
        
!~         if (allocated(routing_states%tabulated_routing_coef)) then
!~             deallocate(routing_states%tabulated_routing_coef)
!~         end if
!~         allocate(routing_states%tabulated_routing_coef(size(routing_states%quantile),routing_states%nb_mode,1))
!~         routing_states%tabulated_routing_coef=0.0
        
!~         gamma_coefficient=0.
!~         mode=0
!~         do i=1,routing_states%nb_mode
            
!~             call compute_gamma_routing_coefficient(routing_states%max_scale,mode,routing_states%quantile,&
!~             &routing_states%window_shift,density_function,gamma_coefficient)
            
!~             routing_states%tabulated_routing_coef(:,i,1)=gamma_coefficient
!~             mode=mode+routing_setup%mode_discretization_step
            
!~         end do
        
        call tabulated_delay_for_gamma(routing_states%nb_mode,&
        &routing_setup%mode_discretization_step,routing_states%tabulated_delay)
        
    end subroutine tabulated_routing_coefficients_2D
    
    
    subroutine generic_tabulated_routing_coefficients_2D(scale,quantile,window_shift,max_mode,frequency,&
    &density_function,tabulated_gamma_coefficient)
        real, intent(in) :: scale
        real,dimension(:), intent(in) :: quantile
        real, intent(in) :: window_shift
        real, intent(in) :: max_mode
        real, intent(in) :: frequency
        character(3),intent(in) :: density_function
        real, dimension(:,:,:), allocatable, intent(out) :: tabulated_gamma_coefficient
        
        integer :: nb_mode
        integer :: i
        real :: mode
        real, dimension(size(quantile)) :: gamma_coefficient
        
        nb_mode=ceiling((max_mode+1.0)/frequency)
        
        if (allocated(tabulated_gamma_coefficient)) then
            deallocate(tabulated_gamma_coefficient)
        end if
        allocate(tabulated_gamma_coefficient(size(quantile),nb_mode,1))
        tabulated_gamma_coefficient=0.0
        
        mode=0
        do i=1,nb_mode
            
            call compute_gamma_routing_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
            tabulated_gamma_coefficient(:,i,1)=gamma_coefficient
            mode=mode+frequency
            
        end do
        
    end subroutine generic_tabulated_routing_coefficients_2D
    
    
    subroutine tabulated_gamma_function_2D(routing_setup,routing_mesh,routing_states,&
    &density_function)
        
        ! Notes
        ! -----
        ! **tabulated_gamma_function_2D(routing_setup,routing_mesh,routing_states,density_function)** :
        !
        ! - Compute the tabulated Gamma function for a range of mode
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_mesh Derived Type (inout)
        ! ``density_function``                    String 'pdf' or 'cdf, please use 'cdf' (in)
        ! =============================           ==================================
        
         implicit none
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        character(3),intent(in) :: density_function
        
        integer :: i,j
        real :: mode
        real, dimension(size(routing_states%quantile)) :: gamma_coefficient
        
        if (allocated(routing_states%tabulated_routing_coef)) then
            deallocate(routing_states%tabulated_routing_coef)
        end if
        
        allocate(routing_states%tabulated_routing_coef(size(routing_states%quantile),&
        &routing_states%nb_mode,1,size(routing_mesh%varying_dx)))
        routing_states%tabulated_routing_coef=0.0
        
        gamma_coefficient=0.
        mode=0
        
        do j=1,size(routing_mesh%varying_dx)
        
            do i=1,routing_states%nb_mode
                
                call compute_gamma_coefficient(routing_states%scale_coef(j),mode,routing_states%quantile,&
                &routing_states%window_shift,density_function,gamma_coefficient)
                
                routing_states%tabulated_routing_coef(:,i,1,j)=gamma_coefficient
                mode=mode+routing_setup%mode_discretization_step
                
            end do
        
        end do
        
        
        call tabulated_delay_for_gamma(routing_states%nb_mode&
        &,routing_setup%mode_discretization_step,routing_states%tabulated_delay)
        
    end subroutine tabulated_gamma_function_2D
    
    
    subroutine generic_tabulated_gamma_function_2D(scale,quantile,window_shift,max_mode,frequency,&
    &density_function,tabulated_gamma_coefficient)
        
        real, intent(in) :: scale
        real,dimension(:), intent(in) :: quantile
        real, intent(in) :: window_shift
        real, intent(in) :: max_mode
        real, intent(in) :: frequency
        character(3),intent(in) :: density_function
        real, dimension(:,:,:), allocatable, intent(out) :: tabulated_gamma_coefficient
        
        integer :: nb_mode
        integer :: i
        real :: mode
        real, dimension(size(quantile)) :: gamma_coefficient
        
        nb_mode=ceiling((max_mode+1.0)/frequency)
        
        
        if (allocated(tabulated_gamma_coefficient)) then
            deallocate(tabulated_gamma_coefficient)
        end if
        allocate(tabulated_gamma_coefficient(size(quantile),nb_mode,1))
        tabulated_gamma_coefficient=0.0
        
        mode=0
        do i=1,nb_mode
            
            call compute_gamma_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
            tabulated_gamma_coefficient(:,i,1)=gamma_coefficient
            mode=mode+frequency
            
        end do
        
    end subroutine generic_tabulated_gamma_function_2D
    
    
    subroutine tabulated_delay_for_gamma(nb_mode,frequency,tabulated_delay)
        
        ! Notes
        ! -----
        ! **tabulated_delay_for_gamma(nb_mode,frequency,tabulated_delay)** :
        !
        ! - Return the tabulated delay corresponding to the tabulated gamma routing coefficient
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``nb_mode``                             Maximum number of discretized mode, real (in)
        ! ``frequency``                           Discretization step,real (in)
        ! ``tabulated_delay``                     The tabulated mode, array (in)
        ! =============================           ==================================
        
        implicit none
    
        integer, intent(in) :: nb_mode
        real, intent(in) :: frequency
        real, dimension(:), allocatable, intent(out) :: tabulated_delay
        
        integer :: i
        
        !nb_mode=ceiling((max_mode+1.0)/frequency)
        
        if (allocated(tabulated_delay)) then
            deallocate(tabulated_delay)
        end if
        allocate(tabulated_delay(nb_mode))
        tabulated_delay=0.0
        
        do i=1,nb_mode
            tabulated_delay(i)=real(i-1)*frequency
        end do
        
    end subroutine tabulated_delay_for_gamma
    
    
    subroutine write_tabulated_coefficients(tabulated_gamma_coefficient,quantile,filename)
        real, dimension(:,:,:), intent(in) :: tabulated_gamma_coefficient
        real, dimension(:), intent(in) :: quantile
        character(300),intent(in) :: filename
        
        integer :: i,j,n,m
        
        n=size(tabulated_gamma_coefficient,1)
        m=size(tabulated_gamma_coefficient,2)
        
        open(10,file=trim(filename),status='replace')
        
        do i=1,n
            write(10,*) quantile(i),(tabulated_gamma_coefficient(i,j,1),j=1,m)
        end do
        
        close(10)
        
    end subroutine write_tabulated_coefficients
    
    
    
    subroutine tabulated_routing_coefficients_3D(routing_setup,routing_mesh,routing_states,density_function)
        
        ! Notes
        ! -----
        ! **tabulated_routing_coefficients_3D(routing_setup,routing_mesh,routing_states,density_function)** :
        !
        ! - Compute the tabulated Gamma function for a range of mode and spreading
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``routing_setup``                       routing_setup Derived Type (in)
        ! ``routing_mesh``                        routing_mesh Derived Type (in)
        ! ``routing_states``                      Routing_states Derived Type (inout)
        ! ``density_function``                    String 'pdf' or 'cdf, please use 'cdf' (in)
        ! =============================           ==================================
        
        implicit none
        
        type(type_routing_setup), intent(in) :: routing_setup
        type(type_routing_mesh), intent(in) :: routing_mesh
        type(type_routing_states), intent(inout) :: routing_states
        character(3),intent(in) :: density_function
        
        integer :: i,j,k
        real :: mode,spread,scale
        real :: spreading_step
        real, dimension(size(routing_states%quantile)) :: gamma_coefficient
        
        
        if (allocated(routing_states%tabulated_routing_coef)) then
            deallocate(routing_states%tabulated_routing_coef)
        end if
        allocate(routing_states%tabulated_routing_coef(size(routing_states%quantile),&
        &routing_states%nb_mode,routing_states%nb_spreads,size(routing_mesh%varying_dx)))
        
        routing_states%tabulated_routing_coef=0.0
        
        scale=1.0
        spread=0.
        
        if (routing_setup%varying_spread .eqv. .false.) then
            spreading_step=routing_states%max_spreading
        end if
        if (routing_setup%varying_spread .eqv. .true.) then
            spreading_step=routing_setup%spreading_discretization_step
        end if
        
        write(*,*) "Tabulation of the Gamma routing coefficients ..."
        
        !loop on number of varying dx
        do k=1,size(routing_mesh%varying_dx)
            
            call wait_bar(k, size(routing_mesh%varying_dx))
            
            do j=1,routing_states%nb_spreads
            
                spread=real(j)*spreading_step
                
                call generic_compute_gamma_scale(routing_mesh%varying_dx(k),routing_setup%dt,routing_setup%vmax,&
                &spread,routing_setup%spreading_gamma_threshold,scale)
                            
                mode=0
                do i=1,routing_states%nb_mode
                    
                    !TODO: ici mode devrait tj être compris entre 1 et nb_mode
                    !et nb_mode=dx/vmin/dt - dx/vmax/dt
                    ! max_mode = dx/vmin/dt
                    ! min_mode = dx/vmax/dt
                    ! nb_mode = max_mode - min_mode
                    mode=real(i-1)*routing_setup%mode_discretization_step
                    
                    call compute_gamma_routing_coefficient(scale,mode,routing_states%quantile,&
                    routing_states%window_shift,density_function,gamma_coefficient)
                    
                    routing_states%tabulated_routing_coef(:,i,j,k)=gamma_coefficient
                    
                end do
                
            end do
            
        end do
        
        call tabulated_delay_for_gamma(routing_states%nb_mode,&
        &routing_setup%mode_discretization_step,routing_states%tabulated_delay)
        
        call tabulated_spreading_for_gamma(routing_states%nb_spreads,&
        &routing_setup%spreading_discretization_step,routing_states%tabulated_spreading)
        
    end subroutine tabulated_routing_coefficients_3D
    
    
    subroutine generic_tabulated_routing_coefficients_3D(dx,dt,vmax,epsilon,quantile,window_shift, max_mode, max_spreading, &
    &frequency_mode,frequency_spreading, density_function, tabulated_gamma_coefficient)
        real,intent(in) :: dx
        real,intent(in) :: dt
        real,intent(in) :: vmax
        real,intent(in) :: epsilon
        real,dimension(:), intent(in) :: quantile
        real, intent(in) :: window_shift
        real, intent(in) :: max_mode
        real, intent(in) :: max_spreading
        real, intent(in) :: frequency_mode
        real, intent(in) :: frequency_spreading
        character(3),intent(in) :: density_function
        real, dimension(:,:,:), allocatable, intent(out) :: tabulated_gamma_coefficient
        
        integer :: nb_mode,nb_spreads
        integer :: i,j
        real :: mode,spread,scale
        real, dimension(size(quantile)) :: gamma_coefficient
                    
        nb_mode=ceiling((max_mode+1.0)/frequency_mode)
        nb_spreads=ceiling(max_spreading/frequency_spreading)
        
        
        if (allocated(tabulated_gamma_coefficient)) then
            deallocate(tabulated_gamma_coefficient)
        end if
        allocate(tabulated_gamma_coefficient(size(quantile),nb_mode,nb_spreads))
        tabulated_gamma_coefficient=0.0
                    
        scale=1.0
        spread=0.
        
        do j=1,nb_spreads
        
            spread=real(j)*frequency_spreading
            call generic_compute_gamma_scale(dx=dx,dt=dt,vmax=vmax,spread=spread,epsilon=epsilon,scale=scale)
            
            mode=0
            do i=1,nb_mode
            
                mode=real(i-1)*frequency_mode
                call compute_gamma_routing_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
                tabulated_gamma_coefficient(:,i,j)=gamma_coefficient
                
            end do
            
        end do
    
    end subroutine generic_tabulated_routing_coefficients_3D
    
    
    
    subroutine generic_tabulated_gamma_function_3D(dx,dt,vmax,epsilon,quantile,window_shift, max_mode, max_spreading, &
    &frequency_mode,frequency_spreading, density_function, tabulated_gamma_coefficient)
        real,intent(in) :: dx
        real,intent(in) :: dt
        real,intent(in) :: vmax
        real,intent(in) :: epsilon
        real,dimension(:), intent(in) :: quantile
        real, intent(in) :: window_shift
        real, intent(in) :: max_mode
        real, intent(in) :: max_spreading
        real, intent(in) :: frequency_mode
        real, intent(in) :: frequency_spreading
        character(3),intent(in) :: density_function
        real, dimension(:,:,:), allocatable, intent(out) :: tabulated_gamma_coefficient
        
        integer :: nb_mode,nb_spreads
        integer :: i,j
        real :: mode,spread,scale
        real, dimension(size(quantile)) :: gamma_coefficient
                    
        nb_mode=ceiling((max_mode+1.0)/frequency_mode)
        nb_spreads=ceiling(max_spreading/frequency_spreading)
        
        
        if (allocated(tabulated_gamma_coefficient)) then
            deallocate(tabulated_gamma_coefficient)
        end if
        allocate(tabulated_gamma_coefficient(size(quantile),nb_mode,nb_spreads))
        tabulated_gamma_coefficient=0.0
                    
        scale=1.0
        spread=0.
        
        do j=1,nb_spreads
        
            spread=(j)*frequency_spreading
            call generic_compute_gamma_scale(dx=dx,dt=dt,vmax=vmax,spread=spread,epsilon=epsilon,scale=scale)
            
            mode=0
            do i=1,nb_mode
            
                mode=real(i-1)*frequency_mode
                call compute_gamma_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
                tabulated_gamma_coefficient(:,i,j)=gamma_coefficient
                
            end do
            
        end do
    
    end subroutine generic_tabulated_gamma_function_3D
    
    
    
    
    subroutine tabulated_spreading_for_gamma(nb_spreads,frequency_spreading,tabulated_spreading)
        
        ! Notes
        ! -----
        ! **tabulated_spreading_for_gamma(nb_spreads,frequency,tabulated_delay)** :
        !
        ! - Return the tabulated spreading corresponding to the tabulated gamma routing coefficient
        !        
        ! =============================           ===================================
        ! Parameters                              Description
        ! =============================           ===================================
        ! ``nb_spreads``                          numbers of spreading steps (in)
        ! ``frequency_spreading``                 Discretization step,real (in)
        ! ``tabulated_spreading``                 The tabulated spreading, array (out)
        ! =============================           ==================================
        
        implicit none
        integer, intent(in) :: nb_spreads
        real, intent(in) :: frequency_spreading
        real, dimension(:), allocatable, intent(out) :: tabulated_spreading
        
        integer :: i
        
        !nb_spreads=ceiling(max_spreading/frequency_spreading)+1
        
        if (allocated(tabulated_spreading)) then
            deallocate(tabulated_spreading)
        end if
        allocate(tabulated_spreading(nb_spreads))
        
        do i=1,nb_spreads
            tabulated_spreading(i)=real(i)*frequency_spreading
        end do
        
    end subroutine tabulated_spreading_for_gamma
    
    
    subroutine write_tabulated_coefficients_3D(tabulated_gamma_coefficient,quantile,filename)
        real, dimension(:,:,:), intent(in) :: tabulated_gamma_coefficient
        real, dimension(:), intent(in) :: quantile
        character(300),intent(in) :: filename
        
        integer :: i,j,k,l,n,m
        
        n=size(tabulated_gamma_coefficient,1)
        m=size(tabulated_gamma_coefficient,2)
        l=size(tabulated_gamma_coefficient,3)
        
        open(10,file=trim(filename),status='replace')
        
        do k=1,l
        
            do i=1,n
                write(10,*) quantile(i),(tabulated_gamma_coefficient(i,j,k),j=1,m)
            end do
            
            if (k<l) then
                write(10,*) ""
                write(10,*) ""
            end if
            
        end do
        
        close(10)
        
    end subroutine write_tabulated_coefficients_3D
    
    
    
    subroutine write_matrix(nline,ncols,matrix,filename)
        integer, intent(in) :: nline
        integer, intent(in) :: ncols
        real, dimension(nline,ncols), intent(in) :: matrix
        character(300),intent(in) :: filename
        
        integer :: i,j
        
        open(10,file=trim(filename),status='replace')
        do i=1,nline
            write(10,*) i,(matrix(i,j),j=1,ncols)
        end do
        close(10)
        
    end subroutine write_matrix
    
    
    subroutine wait_bar(iter, n_iter)
        
        implicit none
        
        integer, intent(in) :: iter, n_iter
        integer :: ch, per
        
        per = 100 * iter / n_iter
        
        if (per /= 100 * (iter - 1) / n_iter) then
             
             write(*,'(256a1)', advance='no') ( char(8), ch = 1, &
             & per/2 + 9 )
             write(*, '(2x,1i3,1a1,2x,1a1,256a1)', advance='no') per, &
             & '%', '|', ( '=', ch = 1, per / 2)
             
        end if
        
        if (iter == n_iter) write(*,'(a)') '|'
          
    end subroutine wait_bar


    !##############################################################################
    ! FUNCTION GammaCDF
    !
    !~ Following function came from :
    !~ https://github.com/fabiankindermann/ce-fortran/blob/main/installation/toolbox/toolbox.f90
    !
    ! Calculates cumulated Gamma distribution at point x.
    !##############################################################################
    function GammaCDF(x, alpha, beta)
     
        implicit none
     
     
        !##### INPUT/OUTPUT VARIABLES #############################################
     
        ! point where to calculate function
        real*8, intent(in) :: x
     
        ! shape parameter of the distribution
        real*8, optional :: alpha
     
        ! scale parameter of the distribution
        real*8, optional :: beta
     
        ! value of cumulated Gamma distribution at p
        real*8 :: gammaCDF
        
        
        !##### OTHER VARIABLES ####################################################
     
        real*8 :: alpha_c, beta_c
     
          
        !##### ROUTINE CODE #######################################################
     
        ! initialize parameters
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta
     
        ! check for validity of parameters
        if(alpha_c <= 0d0)then
            call error('GammaCDF','alpha has a non-positive value')
        endif
        if(beta_c <= 0d0)then
            call error('GammaCDF','beta has a non-positive value')
        endif
     
        if(x <= 0d0)then
            GammaCDF = 0d0
        else
            GammaCDF = incomplete_gamma(beta_c*x, alpha_c)
        endif
     
    end function gammaCDF
    
    
    !##############################################################################
    ! FUNCTION GammaPDF
    !
    ! Calculates Gamma density functions at point x.
    !##############################################################################
    function GammaPDF(x, alpha, beta)
     
        implicit none
     
     
        !##### INPUT/OUTPUT VARIABLES #############################################
     
        ! point where to calculate function
        real*8, intent(in) :: x
     
        ! shape parameter of the distribution
        real*8, optional :: alpha
     
        ! scale parameter of the distribution
        real*8, optional :: beta
     
        ! value of Gamma density at p
        real*8 :: gammaPDF
        
        
        !##### OTHER VARIABLES ####################################################
     
        real*8 :: alpha_c, beta_c
     
          
        !##### ROUTINE CODE #######################################################
     
        ! initialize parameters
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta
     
        ! check for validity of parameters
        if(alpha_c <= 0d0)then
            call error('GammaPDF','alpha has a non-positive value')
        endif
        if(beta_c <= 0d0)then
            call error('GammaPDF','beta has a non-positive value')
        endif
     
       
     
        if(x < 0d0 .or. alpha_c < 1d0 .and. x <= 0d0)then
            GammaPDF = 0d0
        else
            
            GammaPDF = beta_c**alpha_c/gamma_function(alpha_c)* &
                        x**(alpha_c-1d0)*exp(-beta_c*x)
                        
        endif
     
    end function GammaPDF
    
    
    !##############################################################################
    ! FUNCTION gamma_function
    !
    ! Calculates the gamma fuction.
    !##############################################################################
    function gamma_function(x)

        implicit none
     
     
        !##### INPUT/OUTPUT VARIABLES #############################################
     
        ! point where to calculate function
        real*8, intent(in) :: x
     
        ! value of log of the gamma function
        real*8 :: gamma_function
             
     
        !##### ROUTINE CODE #######################################################

        
        ! check validity of inputs solution
        if(x < 0d0)then
            call error('gamma_function','x is smaller than zero')
        endif
        
        gamma_function = my_log_gamma(x)
        gamma_function = exp(gamma_function)
        
    end function gamma_function
    


     !##############################################################################
    ! FUNCTION my_log_gamma
    !
    ! Calculates log of the gamma fuction.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm AS245 from 
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa245/asa245.html
    !
    !     REFERENCE: Macleod, A.J. (1989). Algorithm AS 245: A Robust and Reliable 
    !                Algorithm for the  Logarithm of the Gamma Function. Applied 
    !                Statistics, 38(2), 397-402.
    !##############################################################################
    function my_log_gamma(x_in)

        implicit none
     
     
        !##### INPUT/OUTPUT VARIABLES #############################################
     
        ! point where to calculate function
        real*8, intent(in) :: x_in
     
        ! value of log of the gamma function
        real*8 :: my_log_gamma
     
     
        !##### OTHER VARIABLES ####################################################
     
        real*8, parameter :: alr2pi = 9.18938533204673d-1
        real*8, parameter :: xlge = 5.10d6 
        real*8, parameter :: r1(9) = (/ -2.66685511495d0, &
                                        -2.44387534237d1, &
                                        -2.19698958928d1, & 
                                         1.11667541262d1, &
                                         3.13060547623d0, &
                                         6.07771387771d-1, &
                                         1.19400905721d1, &
                                         3.14690115749d1, &
                                         1.52346874070d1 /)
        real*8, parameter :: r2(9) = (/ -7.83359299449d1, &
                                        -1.42046296688d2, &
                                         1.37519416416d2, & 
                                         7.86994924154d1, &
                                         4.16438922228d0, &
                                         4.70668766060d1, &
                                         3.13399215894d2, & 
                                         2.63505074721d2, &
                                         4.33400022514d1 /)
        real*8, parameter :: r3(9) = (/ -2.12159572323d5, & 
                                         2.30661510616d5, &
                                         2.74647644705d4, &
                                        -4.02621119975d4, &
                                        -2.29660729780d3, &
                                        -1.16328495004d5, &
                                        -1.46025937511d5, &
                                        -2.42357409629d4, &
                                        -5.70691009324d2 /)
        real*8, parameter :: r4(5) = (/  2.79195317918525d-1, &
                                         4.917317610505968d-1, &
                                         6.92910599291889d-2, &
                                         3.350343815022304d0, &
                                         6.012459259764103d0 /)
        real*8 :: x, x1, x2, y
     
     
        !##### ROUTINE CODE #######################################################


        my_log_gamma = 0d0
        x = x_in
        
        ! check validity of inputs solution
        if(x < 0d0)then
            call error('my_log_gamma','x is smaller than zero')
        endif
        
        ! get solution for 0 < X < 0.5 and 0.5 <= x < 1.5
        if(x < 1.5d0)then
            if(x < 0.5d0)then
                my_log_gamma = -log(x)
                y = x + 1d0
                
                ! return if x is smaller than machine epsilon
                if(y <= 1d0)return
            else
                my_log_gamma = 0d0
                y = x
                x = (x-0.5d0) - 0.5d0
            endif
            my_log_gamma = my_log_gamma + x * ((((r1(5)*y + r1(4))*y + r1(3))*y &
                + r1(2))*y + r1(1)) / ((((y + r1(9))*y + r1(8))*y + r1(7))*y + r1(6))
        

        ! get solution for 1.5 <= x < 4.0
        elseif(x < 4.0d0)then
            y = (x - 1d0) - 1d0
            my_log_gamma = y * ((((r2(5)*x + r2(4))*x + r2(3))*x + r2(2))*x &
                + r2(1)) / ((((x + r2(9))*x + r2(8))*x + r2(7))*x + R2(6))

        ! get solution for 4.0 <= x < 12
        elseif(x < 12d0)then
            my_log_gamma = ((((r3(5)*x + r3(4))*x + r3(3))*x + r3(2))*x + r3(1)) / &
                ((((x + r3(9))*x + r3(8))*x + r3(7))*x + r3(6))
        else
            y = log(x)
            my_log_gamma = x * (y -1d0) - 0.5d0*y + alr2pi
            if(x <= XLGE)then
                x1 = 1d0/x
                x2 = x1 * x1
                my_log_gamma = my_log_gamma + x1 * ((r4(3)*x2 + r4(2))*x2 + r4(1)) / &
                    ((x2 + r4(5))*x2 + r4(4))
            endif
        endif
        
    end function my_log_gamma  
     
    
    
        !##############################################################################
    ! FUNCTION incomplete_gamma
    !
    ! Calculates the incomplete gamma integral.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm ASA239 from 
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa239/asa239.html
    !
    !     REFERENCE: Shea, B. (1988). Algorithm AS 239: Chi-squared and Incomplete Gamma Integral, Applied Statistics, 37(3), 466-473.
    !##############################################################################
    function incomplete_gamma(x, p)

        implicit none
        
        
        !##### INPUT/OUTPUT VARIABLES #############################################
     
        ! point where to calculate function
        real*8, intent(in) :: x
     
        ! parameter of the gamma function
        real*8, intent(in) :: p
     
        ! value of normal density at p
        real*8 :: incomplete_gamma
     
     
        !##### OTHER VARIABLES ####################################################
     
        real*8 :: a, an, arg, b, c, pn1, pn2, pn3, pn4, pn5, pn6, rn
        real*8, parameter :: elimit = - 88d0
        real*8, parameter :: oflo = 1d37
        real*8, parameter :: plimit = 1000d0
        real*8, parameter :: tol = 1d-14
        real*8, parameter :: xbig = 1d8
     
     
        !##### ROUTINE CODE #######################################################

        incomplete_gamma = 0d0
        
        ! check validity of inputs solution
        if(x < 0d0)then
            call error('incomeplete_gamma','x is smaller than zero')
        endif
        
        if(p <= 0d0)then
            call error('incomeplete_gamma','p is not positive')
        endif
        
        ! set incomeplete_gamma = 0 if x is zero
        if(x <= 0d0)then
            incomplete_gamma = 0d0
            return
        endif

        ! normal approximation for large values of p
        if(p > plimit)then
            pn1 = 3d0*sqrt(p)*((x/p)**(1d0/3d0) + 1d0/(9d0*p) - 1d0)
            incomplete_gamma = normalCDF(pn1)
            return
        endif
        
        ! set incomeplete_gamma = 1 if x is large
        if(x > xbig)then
            incomplete_gamma = 1d0
        endif
        
        ! Pearson series expansion
        if(x <= 1d0 .or. x < p)then
            arg = p*log(x) - x - my_log_gamma(p+1d0)
            c = 1d0
            incomplete_gamma = 1.0D+00
            a = p
        
            do
                a = a + 1d0
                c = c*x/a
                incomplete_gamma = incomplete_gamma + c
        
                if(c <= tol) then
                    exit
                endif
        
            enddo
        
            arg = arg + log (incomplete_gamma)
        
            if ( elimit <= arg ) then
              incomplete_gamma = exp ( arg )
            else
              incomplete_gamma = 0d0
            endif
            
        ! Use a continued fraction expansion.
        else
        
            arg = p*log (x) - x - my_log_gamma(p)
            a = 1d0 - p
            b = a + x + 1d0
            c = 0d0
            pn1 = 1d0
            pn2 = x
            pn3 = x + 1d0
            pn4 = x*b
            incomplete_gamma = pn3/pn4
        
            do
        
                a = a + 1d0
                b = b + 2d0
                c = c + 1d0
                an = a*c
                pn5 = b*pn3 - an*pn1
                pn6 = b*pn4 - an*pn2
        
                if ( abs(pn6) > 0d0 ) then
        
                    rn = pn5 / pn6
        
                    if ( abs(incomplete_gamma - rn) <= min(tol, tol*rn)) then
                        exit
                    endif
                    incomplete_gamma = rn
                endif
        
                pn1 = pn3
                pn2 = pn4
                pn3 = pn5
                pn4 = pn6
                
                !  Re-scale terms in continued fraction if terms are large.
                if (oflo <= abs(pn5) ) then
                    pn1 = pn1/oflo
                    pn2 = pn2/oflo
                    pn3 = pn3/oflo
                    pn4 = pn4/oflo
                endif
        
            enddo
        
            arg = arg + log (incomplete_gamma)
        
            if (arg >= elimit) then
              incomplete_gamma = 1d0 - exp (arg)
            else
              incomplete_gamma = 1d0
            endif
        endif
    
    end function incomplete_gamma
    
    
    !######################
    !########################################################
    ! FUNCTION normalCDF
    !
    ! Calculates cumulated normal distribution at point x.
    !
    ! PARTS OF THIS PROCEDURE WERE COPIED AND ADAPTED FROM:
    !     Fortran Code by John Burkardt available as Algorithm ASA066 from 
    !     https://people.sc.fsu.edu/~jburkardt/f_src/asa066/asa066.html
    !
    !     REFERENCE: Hill, D. (1973). Algorithm AS 66: The Normal Integral.
    !                Applied Statistics, 22(3), 424-427.
    !##############################################################################
    function normalCDF(x, mu, sigma)
     
        implicit none
     
     
        !##### INPUT/OUTPUT VARIABLES #############################################
     
        ! point where to calculate function
        real*8, intent(in) :: x
     
        ! expectation of distribution
        real*8, optional :: mu
     
        ! variance of distribution
        real*8, optional :: sigma
     
        ! value of the normal distribution at x
        real*8 :: normalCDF
     
     
        !##### OTHER VARIABLES ####################################################
     
        real*8 :: mu_c, sigma_c, xtrans, xabs, y
        real*8, parameter :: a0  = 0.5d0
        real*8, parameter :: a1  = 0.398942280444d0
        real*8, parameter :: a2  = 0.399903438504d0
        real*8, parameter :: a3  = 5.75885480458d0
        real*8, parameter :: a4  = 29.8213557808d0
        real*8, parameter :: a5  = 2.62433121679d0
        real*8, parameter :: a6  = 48.6959930692d0
        real*8, parameter :: a7  = 5.92885724438d0
        real*8, parameter :: b0  = 0.398942280385d0
        real*8, parameter :: b1  = 3.8052d-8
        real*8, parameter :: b2  = 1.00000615302d0
        real*8, parameter :: b3  = 3.98064794d-4
        real*8, parameter :: b4  = 1.98615381364d0
        real*8, parameter :: b5  = 0.151679116635d0
        real*8, parameter :: b6  = 5.29330324926d0
        real*8, parameter :: b7  = 4.8385912808d0
        real*8, parameter :: b8  = 15.1508972451d0
        real*8, parameter :: b9  = 0.742380924027d0
        real*8, parameter :: b10 = 30.789933034d0
        real*8, parameter :: b11 = 3.99019417011d0
     
     
        !##### ROUTINE CODE #######################################################
     
        ! initialize parameters
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))then
            if(sigma_c > 0d0)sigma_c = sqrt(sigma)
        endif
        
        if(sigma_c <= 0d0)then
            call error('normalCDF','sigma has zero or negative value')
        endif
     
        ! standardize evaluation point
        xtrans = (x - mu_c)/sigma_c
        
        ! calculate absolute value and quadratic
        xabs = abs(xtrans)
        y = a0*xtrans**2
        
        ! choose the right interval for calculation
        if(xabs <= 1.28d0)then
            normalCDF = a0-xabs*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
        elseif(xabs <= 12.7d0)then
            normalCDF = b0*exp(-y)/(xabs-b1+b2/(xabs+b3+b4/(xabs-b5+b6/(xabs+b7-b8/ &
                    (xabs+b9+b10/(xabs+b11))))))
        else
            normalCDF = 0d0
        endif
        
        ! transform if other side of the bell
        if(xtrans > 0d0)normalCDF = 1d0-normalCDF
     
    end function normalCDF
    
    
    
    !##############################################################################
    ! SUBROUTINE error
    !
    ! Throws error message and stops program.
    !##############################################################################
    subroutine error(routine, message)
    
    
        !##### INPUT/OUTPUT VARIABLES #############################################
        
        ! routine in which error occured
        character(len=*), intent(in) :: routine
        
        ! error message
        character(len=*), intent(in) :: message
        
        
        !##### ROUTINE CODE #######################################################
        
        ! write error message
        write(*,'(/a,a,a,a/)')'ERROR ',routine,': ',message
        
        ! stop program
        stop
    
    end subroutine error
    
    
    
    !Compute the tabulated IncompleteGamma function
    !With nonuniform memory cells
    !>@kmax - maximum delay step
    !>@kdl - number of delay nodes (depend on beta)
    !>@dlstep  - Computation step of the incompletegamma function
    !>@sp - Spreading Coefficient of the incompletegamma function
    !>beta - Memory elongation coeficient (beta>=1., no elongation for beta=1.)
    !>@dltable () - tabulated delays
    !>@gmtable () - incomplete gamma tabulated
!~     subroutine gamma_table(kmax,kdl,dlstep,sp,beta,dltable,coeftable)
!~         integer :: kmax
!~         integer :: kdl
!~         real, dimension(kdl) :: dltable
!~         real, dimension(kmax,kdl) :: coeftable
!~         real :: dlstep
!~         real :: sp
!~         real :: beta
        
!~         real :: x
!~         real :: x1
!~         real :: dl
!~         integer :: i,j
!~         real, dimension(kmax,kdl) :: gmtable
        
        
!~         do i=1,kdl
        
!~           if(i.eq.1) then
!~             dl=1.0
!~           else
!~             dl=dl+dlstep
!~           end if
          
!~           dltable(i)=dl
          
!~           do j=1,kmax
          
!~             if(j.eq.1) then 
!~               x=0.0
!~             else
!~               x=x+beta**real(j-2)
!~             end if
!~                 !write(*,*) x
!~             x1=x/sp
!~             gmtable(j,i)=gammp(dl,x1)
            
!~           end do
!~         end do
        
!~         !added by max 
!~         coeftable=0.
!~         do i=1,kdl
!~             do j=2,kmax
!~                 coeftable(j-1,i)=gmtable(j,i)-gmtable(j-1,i)
!~             end do
!~             coeftable(:,i)=coeftable(:,i)/sum(coeftable(:,i))
!~         end do
        
!~         return    
!~     end subroutine gamma_table


!
!  Some functions/subroutines to do with gamma functions. AN asterisk
!  means it is user callable, else its an internal subroutine.
!  These are coded from the Numerical Recipes algorithms.
!
! *  Errfun -  Find error function Erf(x)
! *  GammLn -  Find natural log of Gamma function gamma(x)
! *  GammP  -  Find incomplete Gamma function P(a,x)
! *  GammQ  -  Find incomplete Gamma function Q(a,x) = 1 - P(a,x)
!    GSer   -  Find P with series approximation
!    GCF    -  Find Q with continued fraction approximation
!
!
!  History:
!    nebk 03nov92 Original version.
!    nebk 27apr93 Rename erf to errfun
!    nebk 12nov93 Prevent some underflow messages on exponentials
!    nebk 06jan94 Call new ExpUN for underflow messages
!
!************************************************************************
!* GammP -- Find incomplete Gamma function P(a,x)
!& nebk
!: gamma function, utilities
!+
!~ real function gammp (a, x)
!~ !
!~       real a, x
!~ !
!~ ! Find the incomplete Gamma function P(a,x)
!~ ! Algorithm from Numerical Recipes p 161
!~ !
!~ !  Input
!~ !    a,x     r     Arguments
!~ !--
!~ !
!~ !-----------------------------------------------------------------------
!~       real gamser, gln, gammcf
!~ !-----------------------------------------------------------------------
!~       if (x.lt.0.0 .or. a.le.0.0) then
!~         write(*,*)'f','Argument out of range in GammP'
!~         stop 1
!~       end if    
!~ !
!~       if (x.lt.a+1.0) then
!~ !
!~ ! Use series representation
!~ !
!~         call gser (a, x, gamser, gln)
!~         gammp = gamser
!~       else
!~ !
!~ ! Use continued fraction representation
!~ !
!~         call gcf (a, x, gammcf, gln)
!~         gammp = 1.0 - gammcf
!~       end if
!~ !
!~ end function gammp
!~ !
!~ !************************************************************************
!~ !* GammQ -- Find incomplete Gamma function Q(a,x) = 1 - P(a,x)
!~ !& nebk
!~ !: gamma function, utilities
!~ !+
!~ real function gammq (a, x)
!~ !
!~       real a, x
!~ !
!~ ! Find the incomplete Gamma function Q(a,x) = 1 - P(a,x)
!~ ! Algorithm from Numerical Recipes p 162
!~ !
!~ !  Input
!~ !    a,x     r     Arguments
!~ !--
!~ !
!~ !-----------------------------------------------------------------------
!~       real gamser, gln, gammcf
!~ !-----------------------------------------------------------------------
!~       if (x.lt.0.0 .or. a.le.0.0) then
!~         write(*,*)'f','Argument out of range in GammQ'
!~         stop 1
!~       end if    
!~ !
!~       if (x.lt.a+1.0) then
!~ !
!~ ! Use series representation
!~ !
!~         call gser (a, x, gamser, gln)
!~         gammq = 1.0 - gamser
!~       else
!~ !
!~ ! Use continued fraction representation
!~ !
!~         call gcf (a, x, gammcf, gln)
!~         gammq = gammcf
!~       end if
!~ !
!~ end function gammq
!~ !

!~  !************************************************************************
!~ !* GammLn -- Find ln(Gamma(x))
!~ !& nebk
!~ !: gamma function, utilities
!~ !+
!~ function gammln(xx) result(gln)
!~ !     
!~       !real :: gammln
!~       real :: gln
!~       real :: xx
!~ !
!~ ! Find ln(Gamma(x)) where Gamma is the Gamma function
!~ ! Algorithm from Numerical Recipes p 157
!~ !
!~ !  Input
!~ !    xx     r     Argument
!~ !--
!~ !
!~ !-----------------------------------------------------------------------
!~       double precision,dimension(6) :: cof
!~       double precision :: stp, fpf, x, tmp, ser
!~       integer :: j
!~ !
!~       !save cof, stp
!~       cof(1)=76.18009173d0
!~       cof(2)=-86.50532033d0
!~       cof(3)=24.01409822d0
!~       cof(4)=-1.231739516d0
!~       cof(5)=0.120858003d-2
!~       cof(6)=-0.536382d-5
!~       !cof=[76.18009173d0, -86.50532033d0, 24.01409822d0,&
!~       !         &-1.231739516d0, 0.120858003d-2, -0.536382d-5]
!~       stp=2.50662827465d0
!~       fpf=5.5d0
!~       !data cof /76.18009173d0, -86.50532033d0,   24.01409822d0,&
!~       !         &-1.231739516d0,  0.120858003d-2, -0.536382d-5/
!~       !data stp /2.50662827465d0/
!~       !data fpf / 5.5d0/
!~ !-----------------------------------------------------------------------
!~       x = xx - 1.0
!~       tmp = x + fpf
!~       tmp = (x + 0.5)*log(tmp) - tmp
!~       ser = 1.0
!~       do j = 1, 6
!~         x = x + 1.0
!~         ser = ser + cof(j)/x
!~       end do
      
!~       gln = tmp + log(stp*ser)
!~       !gammln = real(tmp + log(stp*ser),kind(gammln))
!~ !
!~ end function gammln
!~ !


!~ !************************************************************************
!~ !
!~       subroutine gser (a, x, gamser, gln)
!~ !
!~       real a, x, gamser, gln
!~ !
!~ ! Find the incomplete Gamma function P(a,x) with its series
!~ ! represenation. Algorithm from Numerical Recipes p 162
!~ !
!~ !  Input
!~ !    a,x     r     Arguments
!~ !  Output
!~ !    gamser  r     P
!~ !    gln     r     ln(Gamma(a))
!~ !
!~ !-----------------------------------------------------------------------
!~       integer itmax
!~       real eps
!~       parameter (itmax = 100, eps = 3.0e-7)
!~ !
!~       real ap, sum, del
!~       integer i
!~ !
!~ !      real gammln, arg, expun
!~       real arg!, exp
!~ !-----------------------------------------------------------------------
!~       gln = gammln(a)
!~ !
!~       if (x.lt.0.0) then
!~         write(*,*)'f', 'Argument out of range in GSer'
!~         stop 1
!~       else if (x.eq.0.0) then
!~         gamser = 0.0
!~       else
!~         ap = a
!~         sum = 1.0 / a
!~         del = sum
!~         do i = 1, itmax
!~           ap = ap + 1.0
!~           del = del * x / ap
!~           sum = sum + del
!~           if (abs(del).lt.abs(sum)*eps) then 
!~             exit
!~           endif
!~         end do
        
!~         if (i>=itmax) then
!~             write(*,*)'f', 'No convergence in GSer'
!~             stop 1
!~         endif
        
!~         arg = -x + a*log(x) - gln
!~         gamser = sum*exp(arg)
        
!~       end if
!~ !
!~       end subroutine gser
!~ !
!~ !************************************************************************
!~ !
!~ subroutine gcf (a, x, gammcf, gln)
!~ !
!~       real a, x, gammcf, gln
!~ !
!~ ! Find the incomplete Gamma function Q(a,x) = 1 - P(a,x) with its
!~ ! continued fraction represenation. Algorithm from Numerical
!~ ! Recipes p 162
!~ !
!~ !  Input
!~ !    a,x     r     Arguments
!~ !  Output
!~ !    gammcf  r     Q(a,x)
!~ !    gln     r     ln(Gamma(a))
!~ !
!~ !-----------------------------------------------------------------------
!~       integer itmax
!~       real eps
!~       parameter (itmax = 100, eps = 3.0e-7)
!~ !
!~       real gold, a0, a1, b0, b1, fac, an, ana, anf, g
!~       integer i
!~ !
!~ !     real gammln, arg, expun
!~       real arg!, exp
!~ !-----------------------------------------------------------------------
!~       gln = gammln(a)
!~ !
!~       gold = 0.0
!~       a0 = 1.0
!~       a1 = x
!~       b0 = 0.0
!~       b1 = 1.0
!~       fac = 1.0
!~ !
!~       do i = 1, itmax
!~         an = real(i)
!~         ana = an - a
!~         a0 = (a1 + a0*ana) * fac
!~         b0 = (b1 + b0*ana) * fac
!~         anf = an * fac
!~         a1 = x*a0 + anf*a1
!~         b1 = x*b0 + anf*b1
!~         if (a1.ne.0.0) then
!~           fac = 1.0 / a1
!~           g = b1 * fac
!~           if (abs((g-gold)/g).lt.eps) then
!~             exit
!~           endif
!~           gold = g
!~         end if
!~       end do
      
!~       if (i>=itmax) then
!~         write(*,*)'f', 'No convergence in GCF', x
!~         stop 1
!~       endif
!~ !
!~       arg = -x + a*log(x) - gln
!~       gammcf = exp(arg) * g
!~ !
!~ end subroutine gcf
!~ !

!~ !
!~ !************************************************************************
!~ !* Errfun -- Find Erf(x)
!~ !& nebk
!~ !: error function, gamma function, utilities
!~ !+
!~       real function errfun (x)
!~ !
!~       real x
!~ !
!~ ! Find Erf(x) from Gamma function
!~ ! Algorithm from Numerical Recipes p 164
!~ !
!~ !  Input
!~ !    x     r     Argument
!~ !--
!~ !
!~ !-----------------------------------------------------------------------
!~ !-----------------------------------------------------------------------
!~       if (x.lt.0.0) then
!~         errfun = -gammp(0.5,x**2)
!~       else
!~         errfun =  gammp(0.5,x**2)
!~       end if
!~ !
!~       end function errfun



    
    
    
end module mod_gamma_function
