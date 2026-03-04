module module_generic_gamma_routing

    contains

    pure subroutine generic_interpolated_routing_coefficients_bilinear(mode,spreading,index_dx,window_length,nb_mode,&
        &nb_spreads,nb_varying_dx,tabulated_delay,tabulated_spreading,tabulated_routing_coef,gamma_coefficient)
        implicit none
        real, intent(in) :: mode
        real, intent(in) :: spreading
        integer, intent(in) :: index_dx
        integer, intent(in) :: window_length
        integer, intent(in) :: nb_mode
        integer, intent(in) :: nb_spreads
        integer, intent(in) :: nb_varying_dx
        real, dimension(nb_mode),intent(in) :: tabulated_delay
        real, dimension(nb_spreads),intent(in) :: tabulated_spreading
        real, dimension(window_length,nb_mode,nb_spreads,nb_varying_dx),intent(in) :: tabulated_routing_coef
        real, dimension(window_length), intent(out) :: gamma_coefficient
        
        integer :: i
        integer :: ix1,iy1,ix2,iy2
        real :: x1,y1,x2,y2
        
        !search indices for delay x
        ix1=0
        ix2=1
        do i=1,nb_mode
            if (mode<tabulated_delay(i)) then
                ix1=i-1
                ix2=i
                exit
            endif
        end do
        
        !search indices for spreading y
        iy1=1
        iy2=2
        do i=2,nb_spreads
            if (spreading<tabulated_spreading(i)) then
                iy1=i-1
                iy2=i
                exit
            endif
        end do
        
        x1=tabulated_delay(ix1)
        x2=tabulated_delay(ix2)
        y1=tabulated_spreading(iy1)
        y2=tabulated_spreading(iy2)
        
        
        gamma_coefficient=&
        & ((mode-x2)/(x1-x2)) * ((spreading-y2)/(y1-y2)) * tabulated_routing_coef(:,ix1,iy1,index_dx)&
        &+((mode-x1)/(x2-x1)) * ((spreading-y2)/(y1-y2)) * tabulated_routing_coef(:,ix2,iy1,index_dx)&
        &+((mode-x2)/(x1-x2)) * ((spreading-y1)/(y2-y1)) * tabulated_routing_coef(:,ix1,iy2,index_dx)&
        &+((mode-x1)/(x2-x1)) * ((spreading-y1)/(y2-y1)) * tabulated_routing_coef(:,ix2,iy2,index_dx)
        
    end subroutine generic_interpolated_routing_coefficients_bilinear
    
    pure subroutine generic_interpolated_routing_coefficients_linear(delay,index_dx,window_length,nb_mode,&
    &nb_varying_dx,tabulated_delay,tabulated_gamma_coefficient,gamma_coefficient)
        implicit none
        real, intent(in) :: delay
        integer, intent(in) :: index_dx
        integer, intent(in) :: window_length
        integer, intent(in) :: nb_mode
        integer, intent(in) :: nb_varying_dx
        real, dimension(nb_mode),intent(in) :: tabulated_delay
        real, dimension(window_length,nb_mode,nb_varying_dx),intent(in) :: tabulated_gamma_coefficient
        real, dimension(window_length), intent(out) :: gamma_coefficient
        
        integer :: i,lim_sup_ind_mu,lim_inf_ind_mu
        real :: delay_inf,delay_sup
        
        lim_inf_ind_mu=0
        lim_sup_ind_mu=1
        !search indices
        do i=1,nb_mode
            if (delay<tabulated_delay(i)) then
                lim_inf_ind_mu=i-1
                lim_sup_ind_mu=i
                exit
            endif
        end do
        
        delay_inf=tabulated_delay(lim_inf_ind_mu)
        delay_sup=tabulated_delay(lim_sup_ind_mu)
        
        gamma_coefficient=tabulated_gamma_coefficient(:,lim_inf_ind_mu,index_dx)+&
        &(delay-delay_inf)*(tabulated_gamma_coefficient(:,lim_sup_ind_mu,index_dx)-&
        &tabulated_gamma_coefficient(:,lim_inf_ind_mu,index_dx))&
        & / (delay_sup-delay_inf)
            
    end subroutine generic_interpolated_routing_coefficients_linear
    
end module module_generic_gamma_routing
