module module_generic_gamma_function

    use mod_gamma_function

    contains
    
!~     !computation of the skweness of the gamma function
!~     !The skweness is used to evaluate the window length
!~     !inputs:
!~     !mode: mode of the Gamma pdf
!~     !scale: scale coefficient of the gamma pdf
!~     !return:
!~     !skweness_gamma, a real
!~     function skweness_gamma(mode,scale_coef)
!~         real :: skweness_gamma
!~         real :: mode
!~         real :: scale_coef
        
!~         skweness_gamma=(2.0/sqrt(mode/scale_coef + 1.0))
!~     end function skweness_gamma
    
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
        
        scale=1./rate
    
    end subroutine generic_compute_gamma_scale
    
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
    
    subroutine generic_tabulated_gamma_function_3D(dx,dt,vmax,epsilon,quantile,window_shift, max_mode, max_sc, &
    &frequency_mode,frequency_spreading, density_function, tabulated_gamma_coefficient)
        real,intent(in) :: dx
        real,intent(in) :: dt
        real,intent(in) :: vmax
        real,intent(in) :: epsilon
        real,dimension(:), intent(in) :: quantile
        real, intent(in) :: window_shift
        real, intent(in) :: max_mode
        real, intent(in) :: max_sc
        real, intent(in) :: frequency_mode
        real, intent(in) :: frequency_spreading
        character(3),intent(in) :: density_function
        real, dimension(:,:,:), allocatable, intent(out) :: tabulated_gamma_coefficient
        
        integer :: nb_mode,nb_spreads
        integer :: i,j
        real :: mode,spread,scale
        real, dimension(size(quantile)) :: gamma_coefficient
                    
        nb_mode=ceiling((max_mode+1.0)/frequency_mode)
        nb_spreads=ceiling(max_sc/frequency_spreading)
        
        
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
                call compute_gamma_routing_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
                tabulated_gamma_coefficient(:,i,j)=gamma_coefficient
                
            end do
            
        end do
    
    end subroutine generic_tabulated_gamma_function_3D
    
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
            
            call compute_gamma_routing_coefficient(scale,mode,quantile,window_shift,density_function,gamma_coefficient)
            tabulated_gamma_coefficient(:,i,1)=gamma_coefficient
            mode=mode+frequency
            
        end do
        
    end subroutine generic_tabulated_gamma_function_2D
    
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
    
    subroutine generic_tabulated_routing_coefficients_3D(dx,dt,vmax,epsilon,quantile,window_shift, max_mode, max_sc, &
    &frequency_mode,frequency_spreading, density_function, tabulated_gamma_coefficient)
        real,intent(in) :: dx
        real,intent(in) :: dt
        real,intent(in) :: vmax
        real,intent(in) :: epsilon
        real,dimension(:), intent(in) :: quantile
        real, intent(in) :: window_shift
        real, intent(in) :: max_mode
        real, intent(in) :: max_sc
        real, intent(in) :: frequency_mode
        real, intent(in) :: frequency_spreading
        character(3),intent(in) :: density_function
        real, dimension(:,:,:), allocatable, intent(out) :: tabulated_gamma_coefficient
        
        integer :: nb_mode,nb_spreads
        integer :: i,j
        real :: mode,spread,scale
        real, dimension(size(quantile)) :: gamma_coefficient
                    
        nb_mode=ceiling((max_mode+1.0)/frequency_mode)
        nb_spreads=ceiling(max_sc/frequency_spreading)
        
        
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
    

end module module_generic_gamma_function

