module osc_mode_dark
    use IO
    implicit none
contains

!-------------------------------------------------------------------------------------------------------------------------------

    subroutine dark_frequency(MATRIX, in_guess, r_disc, gmode, fmode, pmode)
        
        real(DBPR), intent(in) :: r_disc
        real(DBPR), allocatable, intent(in) :: MATRIX(:,:)
        real(DBPR), intent(out) :: gmode, fmode, pmode
        integer :: j, flag, exit_flag
        real(DBPR) :: i, step, f_guess, p_guess, g_guess, temp
        real(DBPR) :: nodes, f_node, p_node, g_node, in_guess(:), NaN = 0.0_DBPR
        real(DBPR), allocatable :: h_in(:), r(:), m1(:), m2(:), p1(:), p2(:), e1(:), e2(:), &
                                dedp1(:), dedp2(:), lamb(:), phi(:), dphi(:), dE_dP_total(:), guesses(:)
        
        call read_TOV_data(MATRIX, h_in, r, m1, m2, p1, p2, e1, e2, dedp1, dedp2, lamb, phi, dphi)
        
        flag = 0

        allocate(dE_dP_total(size(r)))
        do j= 1, size(r)
            if ((p1(j)==0 .and. p2(j).ne.0) .or. (isnan(dedp1(j)))) then
                dE_dP_total(j)= dedp2(j) 
            else if ((p1(j).ne.0 .and. p2(j)==0) .or. (isnan(dedp2(j)))) then
                dE_dP_total(j)= dedp1(j)
            else
                dE_dP_total(j)= dedp1(j)*(p1(j)+e1(j))/(p1(j)+e1(j)+p2(j)+e2(j))+ &
                                dedp2(j)*(p2(j)+e2(j))/(p1(j)+e1(j)+p2(j)+e2(j))
            end if
        end do
        
        ! This do loop finds the first (exit if flag_sign== n), n points where the function cuts the x axis
        ! These points will be the initial guesses for f, p and g modes. 
        ! We will then do newton raphson for these initial guesses. 

        exit_flag = 0

        222 if (in_guess(2)<0) then
            
            i    = 1
            step = 0.1
            call find_guess(MATRIX, i, step, 3, r_disc, guesses, nodes)
            
            f_guess = guesses(1)
            p_guess = guesses(2)
            g_guess = guesses(3)
            
            call find_freq(h_in, f_guess, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, f_node, r_disc, fmode)
            call find_freq(h_in, p_guess, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, p_node, r_disc, pmode)
            
        else
            
            f_guess = in_guess(1)
            p_guess = in_guess(2)
            g_guess = in_guess(3)
            
            call find_freq(h_in, f_guess, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, f_node, r_disc, fmode)
            call find_freq(h_in, p_guess, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, p_node, r_disc, pmode)

        end if  

        if (((p_node .ne. 1.0_DBPR).or.(f_node .ne. 0.0_DBPR)).and.exit_flag<1) then
            in_guess = (/-10.0_DBPR,-10.0_DBPR,-10.0_DBPR/)
            exit_flag = exit_flag+1
            goto 222
        end if

        if (r_disc>0.) then
            if (g_guess<0) then
                i    = 1
                step = 0.1
                call find_guess(MATRIX, i, step, 1, r_disc, guesses, nodes)
                g_guess = guesses(1)
            end if
            call find_freq(h_in, g_guess, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, g_node, r_disc, gmode)
            gmode = gmode*47.71314928e-3_DBPR
        else
            gmode = -10.0_DBPR
        end if
        
        fmode = fmode*47.71314928e-3_DBPR
        pmode = pmode*47.71314928e-3_DBPR

        if (gmode>fmode) then
            ! gmode has the highest frequency and is actually the pmode. fmode has the lowest frequency. So fix the labels. 
            temp  = gmode
            gmode = fmode
            fmode = pmode
            pmode = temp
        end if

        if ((f_node.ne.0.0_DBPR).or.(p_node.ne.1.0_DBPR).or.(gmode>0.and.g_node.ne.1)) then
        
            if (flag == 1) then
                if (f_node.ne.0.0_DBPR) fmode = NaN/NaN
                if (p_node.ne.1.0_DBPR) pmode = NaN/NaN
                if (g_node.ne.1.0_DBPR) gmode = NaN/NaN
            else
                in_guess = (/-10.0_DBPR,-10.0_DBPR,-10.0_DBPR/)
                flag = flag+1
                goto 222
            end if

        end if

        ! print*, 'f-node= ', f_node, 'p1-node= ',p_node, 'g-node= ', g_node
        
    end subroutine dark_frequency

!-------------------------------------------------------------------------------------------------------------------------------

    subroutine find_guess(MATRIX, i, step, stop, r_disc, guesses, nodes)
        
        integer, intent(in) :: stop
        real(DBPR), intent(in) :: r_disc
        real(DBPR), allocatable, intent(in) :: MATRIX(:,:)
        real(DBPR), intent(out) :: nodes
        real(DBPR), allocatable, intent(out) ::  guesses(:)
        integer :: j, flag_sign
        real(DBPR) :: i, step, func_val, sign_prev, sign_
        real(DBPR), allocatable :: h_in(:), r(:), m1(:), m2(:), p1(:), p2(:), e1(:), e2(:), &
                                dedp1(:), dedp2(:), lamb(:), phi(:), dphi(:), dE_dP_total(:)
        
        allocate(guesses(0))

        call read_TOV_data(MATRIX, h_in, r, m1, m2, p1, p2, e1, e2, dedp1, dedp2, lamb, phi, dphi)

        allocate(dE_dP_total(size(r)))
        do j= 1, size(r)
            if ((p1(j)==0 .and. p2(j).ne.0) .or. (isnan(dedp1(j)))) then
                dE_dP_total(j)= dedp2(j) 
            else if ((p1(j).ne.0 .and. p2(j)==0) .or. (isnan(dedp2(j)))) then
                dE_dP_total(j)= dedp1(j)
            else
                dE_dP_total(j)= dedp1(j)*(p1(j)+e1(j))/(p1(j)+e1(j)+p2(j)+e2(j))+ &
                                dedp2(j)*(p2(j)+e2(j))/(p1(j)+e1(j)+p2(j)+e2(j))
            end if
        end do
        
        flag_sign = 0
        call mode_rk4(h_in, i, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, r_disc, func_val, nodes)
        sign_ = sign(real(1, DBPR), func_val)
        
        666 do

            sign_prev = sign_
            call mode_rk4(h_in, i, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, r_disc, func_val, nodes)
            sign_     = sign(real(1, DBPR), func_val)
            
            if (sign_/= sign_prev) then 
                guesses = [guesses, i-step]
                flag_sign = flag_sign +1
                goto 666
            end if
            
            i=i+step

            if (flag_sign== stop) then 
                exit
            end if

        end do

    end subroutine find_guess

!-------------------------------------------------------------------------------------------------------------------------------

    subroutine find_freq(h_in, guess, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, nodes, r_disc, root)

        real(DBPR), intent(in) :: h_in(:), r(:), p1(:), p2(:), e1(:), e2(:), &
                                dE_dP_total(:), lamb(:), phi(:), dphi(:), r_disc
        real(DBPR), intent(out) ::  root, nodes
        real(DBPR) :: tol, error, h, original_guess, guess, new_guess, f, fh, fdash

        111 original_guess = guess

        h     = 1e-10_DBPR     ! This is the order of error on derivative for Newton-Raphson method
        tol   = 1e-10_DBPR     ! Tolerance for closeness of values in Newton-Raphson
        error = 10             ! Temporary initial value for error
        
        f = 0
        fh = 0
        fdash = 0
        new_guess = 0

        do while (error>tol)
        
            call mode_rk4(h_in, guess, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, r_disc, f, nodes)
            call mode_rk4(h_in, guess-h, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, r_disc, fh, nodes)

            fdash     = (f- fh)/h
            new_guess = guess- (f/fdash)
            error     = abs(new_guess-guess)
            guess     = new_guess
            
        end do
        
        root = guess

        if (root<0.0_DBPR .or. isnan(root)) then
            guess = original_guess+0.01_DBPR
            goto 111
        end if
    
    end subroutine find_freq

!-------------------------------------------------------------------------------------------------------------------------------

    subroutine mode_rk4(h_in, omega, r, p1, p2, e1, e2, dE_dP_total, lamb, phi, dphi, r_disc, fomega, nodes) 

        real(DBPR), intent(in) :: h_in(:), omega, r(:), p1(:), p2(:), e1(:), e2(:), &
                                dE_dP_total(:), lamb(:), phi(:), dphi(:), r_disc
        real(DBPR), intent(out) :: fomega, nodes
        integer :: i
        real(DBPR) :: l, W, V, Wm, Vm, Wp, Vp, c1, c2, x(2), k1(2), k2(2), k3(2), k4(2), sign_W_p, sign_V_p, sign_W, sign_V

        l = 2.0_DBPR
        W = r(1)**(l+1)
        V = -(r(1)**l)/l
        x = [W, V]
        
        nodes= 0
        sign_W_p= sign(real(1, DBPR), W)
        sign_V_p= sign(real(1, DBPR), V)
        
        do i= 1, size(r)-1
            
            Wm = x(1)
            Vm = x(2)
            
            k1 = h_in(i)* mode_func(r(i), x, omega, dE_dP_total(i), lamb(i), phi(i), dphi(i))
            
            k2 = h_in(i)* mode_func(r(i)+0.5_DBPR*h_in(i), x+0.5_DBPR*k1, omega, &
                (dE_dP_total(i)+dE_dP_total(i+1))*0.5_DBPR, (lamb(i)+lamb(i+1))*0.5_DBPR, (phi(i)+phi(i+1))*0.5_DBPR, &
                (dphi(i)+dphi(i+1))*0.5_DBPR)
            
            k3 = h_in(i)* mode_func(r(i)+0.5_DBPR*h_in(i), x+0.5_DBPR*k2, omega, &
                (dE_dP_total(i)+dE_dP_total(i+1))*0.5_DBPR, (lamb(i)+lamb(i+1))*0.5_DBPR, (phi(i)+phi(i+1))*0.5_DBPR, &
                (dphi(i)+dphi(i+1))*0.5_DBPR)
            
            k4 = h_in(i)* mode_func(r(i)+h_in(i), x+k3, omega, dE_dP_total(i+1), lamb(i+1), phi(i+1), dphi(i+1))
            
            x  = x+ (1.0_DBPR/6.0_DBPR)*(k1+ 2*k2+ 2*k3+ k4)
            
            if (abs(r(i)-r_disc)<1E-10) then

                Wp   = x(1)
                Vp   = x(2)
                c1   = exp(lamb(i)-2*phi(i))*(omega**2)*(r_disc**2)                          
                c2   = ((e1(i)+e2(i)+ p1(i)+p2(i))/(e1(i+1)+e2(i+1)+ p1(i+1)+p2(i+1)))*((c1*Vm) + (dphi(i)*Wm))
                x(1) = Wm
                x(2) = (c1**(-1.))*(c2 - dphi(i)*Wm)
                
            end if

            sign_W= sign(real(1, DBPR), x(1))
            sign_V= sign(real(1, DBPR), x(2))
            if ((sign_W.ne.sign_W_p) .and. (sign_V.ne.sign_V_p)) then
                nodes= nodes+1
                sign_W_p = sign_W
                sign_V_p = sign_V
            end if
            
        end do

        W = x(1)
        V = x(2)
        
        fomega = (omega**2)*exp(lamb(i))*exp(-2*phi(i))*V+ (1/(r(i)**2))*dphi(i)*W
        
    end subroutine mode_rk4

!-------------------------------------------------------------------------------------------------------------------------------

    function mode_func(r, x, omega, dE_dP_total, lamb, phi, dphi) result(dW_dV)

        real(DBPR), intent(in) :: r, omega, dE_dP_total, lamb, phi, dphi
        real(DBPR), intent(in) :: x(:) 
        real(DBPR) :: dW_dv(2), l, W, V, dW, dV

        W = x(1)
        V = x(2)
        l = 2.
        
        dW = dE_dP_total* ((omega**2)* (r**2)* exp(lamb)* exp(-2*phi)* V+ (dphi* W))- (l*(l+1)* exp(lamb)* V)  
        dV = 2* dphi* V- (exp(lamb)* W)/ (r**2)
        
        dW_dV = [dW, dV]
        
    end function mode_func

!-------------------------------------------------------------------------------------------------------------------------------

end module osc_mode_dark