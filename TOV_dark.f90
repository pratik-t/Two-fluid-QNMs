module TOV_dark
    use IO
    use odepack_mod
    use osc_mode_dark
    implicit none
    
contains

! TOV ODEs--------------------------------------------------------------------------------------------------------------------

function TOV_dark_func(e1, e2, r, x) result(dm_dp_dphi)

    real(DBPR), intent(in) :: e1, e2, r
    real(DBPR), intent(in) :: x(:)
    real(DBPR) :: dm_dp_dphi(size(x))
    real(DBPR) :: M, p1, p2, phi, dm1, dm2, dp1, dp2, dphi
    real(DBPR), parameter :: PI= 3.14159265358979323846_DBPR
    
    if (size(x)==3) then

        M   = x(1)
        p1  = x(2)
        phi = x(3)

        dm1  = 4*PI*e1*r**2
        dphi = (M + 4*PI*p1*r**3)/(r*(r - 2*M))
        dp1  = -(e1 + p1)*dphi

        dm_dp_dphi= [dm1, dp1, dphi]

    else

        M   = x(1) + x(2)
        p1  = x(3)
        p2  = x(4)
        phi = x(5)

        dm1  = 4*PI*e1*r**2
        dm2  = 4*PI*e2*r**2
        dphi = (M + 4*PI*(p1+p2)*r**3)/(r*(r - 2*M))
        dp1  = -(e1 + p1)*dphi
        dp2  = -(e2 + p2)*dphi
        
        dm_dp_dphi= [dm1, dm2, dp1, dp2, dphi]

    end if
    
end function TOV_dark_func

! Tidal Deformability ODEs----------------------------------------------------------------------------------------------------

function tidal_func(r, y, e, p, temp, lamb, dphi) result(dy) 

    real(DBPR), intent(in) :: y, r, e, p, temp, lamb, dphi
    real(DBPR) :: Q, dy
    real(DBPR), parameter :: PI= 3.14159265358979323846_DBPR

    Q  = 4*PI*exp(2*lamb)*(5*e+ 9*p+ temp)- (6*exp(2*lamb))/(r**2)- 4*(dphi)**2
    dy = (-1./r)*(y**2+ y*exp(2*lamb)*(1.+ 4*PI*(r**2)*(p- e))+ Q*r**2)
    
end function tidal_func

!Solving TOV------------------------------------------------------------------------------------------------------------------

subroutine SolveTOV_dark(h_in, p1_data, p2_data, e1_data, e2_data, ec1, ec2, r_lim, &
                    e1_disc, r1_fin, r2_fin, m1_fin, m2_fin, yR, in_guess, gnu, fnu, pnu) 
        
    real(DBPR), intent(in) :: h_in(3), p1_data(:), p2_data(:), e1_data(:), e2_data(:), ec1, ec2, r_lim, e1_disc
    real(DBPR), intent(in), optional :: in_guess(:)
    real(DBPR), intent(out) :: m1_fin, m2_fin, r1_fin, r2_fin
    real(DBPR), intent(out), optional :: yR, gnu, fnu, pnu
    integer :: i, flag, sz, dc
    real(DBPR), parameter :: PI= 3.14159265358979323846_DBPR
    real(DBPR) :: h, h_twof, h_norm, h_dark, r0, e1, e2, p1, p2, p1_limit, p2_limit, &
                mc, m1, m2, pc1, pc2, phi, phic, phi_fin, phi_far, err_twof, &
                err_norm, err_dark, r_disc, m1_disc, m2_disc, p1_disc, p2_disc, e2_disc, phi_disc
    real(DBPR), allocatable :: x0_prev(:), x0(:), x1(:), h_array(:), r_array(:), m1_array(:), m2_array(:), &
                            p1_array(:), p2_array(:), e1_array(:), e2_array(:), de_dp_1(:), de_dp_2(:), &
                            phi_array(:), dphi_array(:), lamb_array(:), MATRIX(:,:)

    
    call interpolate(e1_data, p1_data, ec1, pc1)
    call interpolate(e2_data, p2_data, ec2, pc2)
    
    if (ec1==0) pc1 = 0
    if (ec2==0) pc2 = 0

    mc   = 0.0_DBPR
    phic = 0.0_DBPR
    
    p1_limit = p1_data(1)         ! The limit replacing p<0
    p2_limit = p2_data(1)
    m1_fin   = 0.0_DBPR
    m2_fin   = 0.0_DBPR
    r1_fin   = 0.0_DBPR
    r2_fin   = 0.0_DBPR

    flag = 0
    
    r0   = 1E-15_DBPR                ! Starting point of radius. This increases by h in each step.
    x0   = [mc, mc, pc1, pc2, phic]  ! [Initial mass, Initial pressure, Initial phi]
    p1   = x0(3)                     ! p= initial pressure. This will change in the loop.
    p2   = x0(4)
    e1   = ec1                       ! e= initial energy density. This will change in the loop.
    e2   = ec2                       
    
    err_twof = 1E-8_DBPR
    err_norm = 1E-8_DBPR
    err_dark = 1E-8_DBPR

    h_twof   = h_in(1)
    h_norm   = h_in(2)
    h_dark   = h_in(3)

    allocate(m1_array(0), m2_array(0))

    r_array   = [r0]
    p1_array  = [pc1]
    p2_array  = [pc2]
    e1_array  = [ec1]
    e2_array  = [ec2]
    phi_array = [phic]

    if(ec1==0.0_DBPR) then
        x0   = [mc, pc2, phic]
        flag = 1
    else if(ec2==0.0_DBPR) then
        x0   = [mc, pc1, phic]
        flag = 1
    end if 
    
!Main RK4 Loop----------------------------------------------------------------------------------------------------------------

    do while(p1>=0.0_DBPR .or. p2>=0.0_DBPR)
        
        if (p1>p1_limit .and. p2>p2_limit) then
            
111         x0_prev = x0
            h       = h_twof
            m1_fin  = x0(1)
            m2_fin  = x0(2)
            r1_fin  = r0
            r2_fin  = r0
            phi_fin = x0(5)
            
            call rk4_(h, e1, e2, r0, x0, x1)
            
            if ((x1(3)<p1_limit.and.abs(x0(3)-p1_limit)>err_twof).or.(x1(4)<p2_limit.and.abs(x0(4)-p2_limit)>err_twof)) then
                h_twof = h_twof/10
                goto 111
            else 
                x0  = x1
                r0  = r0+h
                m1  = x0(1)
                m2  = x0(2)
                p1  = x0(3)
                p2  = x0(4)
                phi = x0(5)
                call interpolate(p1_data, e1_data, p1, e1)
                call interpolate(p2_data, e2_data, p2, e2)
            end if
    
        else if (p1>p1_limit .and. p2<p2_limit) then
           
222         if (flag==0) then
                flag = 1
                if (allocated(x1)) then
                    deallocate(x0, x1)
                    allocate(x0(3), x1(3))
                end if
                call interpolate(p1_data, e1_data, x0_prev(3), e1)
                m2_fin = x0_prev(2)
                x0 = [x0_prev(1)+x0_prev(2), x0_prev(3), x0_prev(5)]
            end if
            
            h       = h_norm
            m1_fin  = x0(1)-m2_fin
            r1_fin  = r0
            phi_fin = x0(3)
            e2      = 0.0_DBPR
            
            call rk4_(h, e1, e2, r0, x0, x1)
            
            if (x1(2)<p1_limit.and.abs(x0(2)-p1_limit)>err_norm) then
                h_norm = h_norm/10
                goto 222
            else
                x0  = x1
                r0  = r0+h
                m1  = x0(1)
                m2  = m2_fin
                p1  = x0(2)
                p2  = 0.0_DBPR
                phi = x0(3)    
                call interpolate(p1_data, e1_data, p1, e1)
            end if
            
        else if (p1<p1_limit .and. p2>p2_limit) then

333         if (flag==0) then                
                flag = 1
                if (allocated(x1)) then
                    deallocate(x0, x1)
                    allocate(x0(3), x1(3))
                end if
                call interpolate(p2_data, e2_data, x0_prev(4), e2)
                m1_fin = x0_prev(1)
                x0 = [x0_prev(1)+x0_prev(2), x0_prev(4), x0_prev(5)]
            end if
            
            h       = h_dark
            m2_fin  = x0(1)-m1_fin
            r2_fin  = r0
            phi_fin = x0(3)
            e1      = 0.0_DBPR
            
            call rk4_(h, e2, e1, r0, x0, x1)

            if (x1(2)<p2_limit.and.abs(x0(2)-p2_limit)>err_dark) then
                h_dark = h_dark/10
                goto 333
            else
                x0 = x1
                r0  = r0+h
                m1  = m1_fin
                m2  = x0(1)
                p1  = 0.0_DBPR
                p2  = x0(2)
                phi = x0(3)
                call interpolate(p2_data, e2_data, p2, e2)
            end if

        else if (p1<p1_limit .and. p2<p2_limit) then
            exit
        end if  
        
        if (r1_fin*1000.0_DBPR-r_lim > 1.0_DBPR) then
            r1_fin = r_lim+1000.0_DBPR
            r2_fin = r_lim+1000.0_DBPR
            exit
        end if

        ! if(present(yR)) print*, r0, e1, e2

        if ((p1>p1_limit.and.p2>p2_limit).or.(p1==0.0_DBPR.and.p2>p2_limit).or.(p2==0.0_DBPR.and.p1>p1_limit)) then
            r_array   = [r_array, r0]
            m1_array  = [m1_array, m1_fin]
            m2_array  = [m2_array, m2_fin]
            p1_array  = [p1_array, p1]
            p2_array  = [p2_array, p2]
            e1_array  = [e1_array, e1]
            e2_array  = [e2_array, e2]
            phi_array = [phi_array, phi]
        end if
    
    end do
    
    m1_array  = [m1_array, m1_fin]
    m2_array  = [m2_array, m2_fin]
    
!φ array---------------------------------------------------------------------------------------------------------------------
        
    phi_far   = (1./2.)*log(1.- (2.*(m1_fin+m2_fin)/max(r1_fin, r2_fin)))
    phi_array = phi_array+ (phi_far-phi_fin)

!Discontinuity Inclusion-----------------------------------------------------------------------------------------------------

    if (e1_disc<ec1 .and. e1_disc>0.0_DBPR) then

        sz = size(r_array)

        call interpolate(e1_array(sz:1:-1), r_array(sz:1:-1), e1_disc, r_disc)
        call interpolate(r_array, m1_array,   r_disc, m1_disc)
        call interpolate(r_array, m2_array,   r_disc, m2_disc)
        call interpolate(r_array, p1_array,   r_disc, p1_disc)
        call interpolate(r_array, p2_array,   r_disc, p2_disc)
        call interpolate(r_array, e2_array,   r_disc, e2_disc)
        call interpolate(r_array, phi_array,  r_disc, phi_disc)

        do i= 1, size(r_array)
            if (e1_array(i)-e1_disc<0) exit
        end do
        dc= i

        r_array    = [r_array(1:dc-1),    r_disc,       r_array(dc:sz)]
        m1_array   = [m1_array(1:dc-1),   m1_disc,      m1_array(dc:sz)]
        m2_array   = [m2_array(1:dc-1),   m2_disc,      m2_array(dc:sz)]
        p1_array   = [p1_array(1:dc-1),   p1_disc,      p1_array(dc:sz)]
        p2_array   = [p2_array(1:dc-1),   p2_disc,      p2_array(dc:sz)]
        e1_array   = [e1_array(1:dc-1),   e1_disc,      e1_array(dc:sz)]
        e2_array   = [e2_array(1:dc-1),   e2_disc,      e2_array(dc:sz)]
        phi_array  = [phi_array(1:dc-1),  phi_disc,     phi_array(dc:sz)]

    else
        r_disc = -10.0_DBPR

    end if

    allocate(h_array(size(r_array)), de_dp_1(size(r_array)), de_dp_2(size(r_array)), &
    dphi_array(size(r_array)), lamb_array(size(r_array)))

    do i= 1, size(r_array)-1
        h_array(i) = r_array(i+1)-r_array(i)
        de_dp_1(i) = (e1_array(i+1)-e1_array(i))/(p1_array(i+1)-(p1_array(i)))
        de_dp_2(i) = (e2_array(i+1)-e2_array(i))/(p2_array(i+1)-(p2_array(i)))
    end do

    h_array(size(r_array)) = h_array(size(r_array)-1)
    de_dp_1(size(r_array)) = de_dp_1(size(r_array)-1)
    de_dp_2(size(r_array)) = de_dp_2(size(r_array)-1)

    dphi_array = ((m1_array+m2_array) + 4*PI*(p1_array+p2_array)*(r_array**3))/&
        (r_array*(r_array - 2*(m1_array+m2_array)))

    lamb_array = log( (1.- 2.*(m1_array+m2_array)/r_array)**(-1./2.) )
    
    if (present(yR)) then
        do i= 1, size(m1_array)
            write(11,*) r_array(i), p1_array(i), p2_array(i), m1_array(i), m2_array(i), phi_array(i)
        end do
    end if
    
!Tidal Deformability---------------------------------------------------------------------------------------------------------
    
    ! if (present(yR)) then

    !     y0 = 2.0_DBPR
    !     y_array = [y0]
        
    !     do i= 2, size(r_array)
            
    !         if (abs(p2_array(i-1))<1E-10 .and. abs(e2_array(i-1))<1E-10) then
    !             call rk4_tide(h_array(i-1), y0, r_array(i-1), e1_array(i-1), p1_array(i-1), de_dp_1(i-1), &
    !                         lamb_array(i-1), dphi_array(i-1))
    !             y_array = [y_array, y0]
                
    !         else if (abs(p1_array(i-1))<1E-10 .and. abs(e1_array(i-1))<1E-10) then
    !             call rk4_tide(h_array(i-1), y0, r_array(i-1), e2_array(i-1), p2_array(i-1), de_dp_2(i-1), &
    !                         lamb_array(i-1), dphi_array(i-1))
    !             y_array = [y_array, y0] 
                
    !         else 
    !             call rk4_tide(h_array(i-1), y0, r_array(i-1), e1_array(i-1), p1_array(i-1), de_dp_1(i-1), &
    !                         lamb_array(i-1), dphi_array(i-1), e2_array(i-1), p2_array(i-1), de_dp_2(i-1))
    !             y_array = [y_array, y0]
                
    !         end if

    !     end do
    !     yR = y_array(size(r_array))
        
    ! end if
    
!Non Radial Modes-------------------------------------------------------------------------------------------------------------

    if (present(gnu)) then
        
        allocate(MATRIX(13,size(r_array)))
        MATRIX(1,:)  = h_array
        MATRIX(2,:)  = r_array
        MATRIX(3,:)  = m1_array
        MATRIX(4,:)  = m2_array
        MATRIX(5,:)  = p1_array
        MATRIX(6,:)  = p2_array
        MATRIX(7,:)  = e1_array
        MATRIX(8,:)  = e2_array
        MATRIX(9,:)  = de_dp_1
        MATRIX(10,:) = de_dp_2
        MATRIX(11,:) = lamb_array
        MATRIX(12,:) = phi_array
        MATRIX(13,:) = dphi_array
        
        call dark_frequency(MATRIX, in_guess, r_disc, gnu, fnu, pnu)
        
    end if

!Returning--------------------------------------------------------------------------------------------------------------------

    m1_fin = m1_fin*1000.0_DBPR/1.4766_DBPR    ! Output mass in units of solar mass
    m2_fin = m2_fin*1000.0_DBPR/1.4766_DBPR
    r1_fin = r1_fin*1000.0_DBPR                ! Output radius in kilometers
    r2_fin = r2_fin*1000.0_DBPR
    
end subroutine SolveTOV_dark

!RK4-------------------------------------------------------------------------------------------------------------------------

subroutine rk4_(h_, e1_, e2_, r_, x0_, x1_)

    real(DBPR), intent(in) :: h_, e1_, e2_, r_
    real(DBPR), intent(in) :: x0_(:)
    real(DBPR), allocatable, intent(out) :: x1_(:)
    real(DBPR) :: k1_(size(x0_)), k2_(size(x0_)), k3_(size(x0_)), k4_(size(x0_))

    k1_ = h_* TOV_dark_func (e1_, e2_, r_, x0_)
    k2_ = h_* TOV_dark_func (e1_, e2_, r_+0.5_DBPR*h_, x0_+0.5_DBPR*k1_)
    k3_ = h_* TOV_dark_func (e1_, e2_, r_+0.5_DBPR*h_, x0_+0.5_DBPR*k2_)
    k4_ = h_* TOV_dark_func (e1_, e2_, r_+h_, x0_+k3_)
    
    x1_ = x0_+ (1.0_DBPR/6.0_DBPR)*(k1_+2.0_DBPR*k2_+2.0_DBPR*k3_+k4_)

end subroutine rk4_

!RK4 tide--------------------------------------------------------------------------------------------------------------------

subroutine rk4_tide(h_, y_, r_, e_, p_, temp_, lamb_, dphi_)

    real(DBPR), intent(in) :: h_, r_, e_, p_, temp_, lamb_, dphi_
    real(DBPR), intent(inout) :: y_
    real(DBPR) :: k1_, k2_, k3_, k4_

    k1_ = h_* tidal_func (r_, y_, e_, p_, temp_, lamb_, dphi_)
    k2_ = h_* tidal_func (r_+0.5_DBPR*h_, y_+0.5_DBPR*k1_, e_, p_, temp_, lamb_, dphi_)
    k3_ = h_* tidal_func (r_+0.5_DBPR*h_, y_+0.5_DBPR*k2_, e_, p_, temp_, lamb_, dphi_)
    k4_ = h_* tidal_func (r_+h_, y_+k3_, e_, p_, temp_, lamb_, dphi_)

    y_  = y_+ (1.0_DBPR/6.0_DBPR)*(k1_+2.0_DBPR*k2_+2.0_DBPR*k3_+k4_)
    
end subroutine rk4_tide

!----------------------------------------------------------------------------------------------------------------------------

subroutine newton_raphson(h_in, p1_data, p2_data, e1_data, e2_data, e1, e2, r_lim, e1_disc, dark_frac, root_e2, &
                        R1, R2, M1, M2, yR, in_guess, gnu, fnu, pnu)

    real(DBPR), intent(in) :: h_in(:), p1_data(:), p2_data(:), e1_data(:), e2_data(:), e1, e2, r_lim, e1_disc, dark_frac
    real(DBPR), intent(in), optional :: in_guess(:)
    real(DBPR), intent(out) ::  root_e2, R1, R2, M1, M2
    real(DBPR), intent(out), optional :: yR, gnu, fnu, pnu
    integer :: iter, flag
    real(DBPR) :: h, tol, error, prev_guess, guess, new_guess, f_pre, f, fh, fdash

    h     = 1e-10_DBPR        ! This is the order of error on derivative for Newton-Raphson method

    123 tol   = 1e-7_DBPR     ! Tolerance for closeness of values in Newton-Raphson
    error = 10                ! Temporary initial value for error
    guess = e2
    iter  = 0
    flag  = 0

    do while (error>tol)

        call SolveTOV_dark(h_in, p1_data, p2_data, e1_data, e2_data, e1, guess-h, r_lim, -10.0_DBPR, R1, R2, M1, M2)
        fh = abs(M2/(M1+M2)- dark_frac)

        call SolveTOV_dark(h_in, p1_data, p2_data, e1_data, e2_data, e1, guess, r_lim, -10.0_DBPR, R1, R2, M1, M2)
        f = abs(M2/(M1+M2)- dark_frac)
  
        prev_guess = guess
        fdash     = (f- fh)/h
        
        new_guess = guess- (f/fdash)
        error     = f
        root_e2   = guess
        guess     = new_guess
        iter      = iter+1
        
        ! print*, h, guess/1.323423_DBPR, fdash, f, M2/(M1+M2), iter

        if (abs(fdash)==0 .and. h<1E-7_DBPR) then
            h = h*10.0_DBPR
            goto 123
        end if

        if (fdash==0.0_DBPR .or. iter>=50 .or. guess<0) then
            flag = 1
            exit
        end if

    end do

    if (flag==1) then
        
        iter  = 0
        guess = e2
        h     = 1000000
        
        777 if (guess-h <0 ) then
            h= h/10 
            goto 777
        end if

        do while (error>tol)
            
            new_guess = guess-h
            call SolveTOV_dark(h_in, p1_data, p2_data, e1_data, e2_data, e1, new_guess, r_lim, -10.0_DBPR, R1, R2, M1, M2)
            f_pre = f
            f = (M2/(M1+M2))-dark_frac

            ! print*, h, new_guess/1.323423_DBPR, R1, R2, M2/(M1+M2), iter
            
            if (f<0 .or. new_guess-h<0) then
                h= h/10
            else if (f>1E-12) then
                guess= new_guess
            end if
            
            error= abs(f-f_pre)
            iter= iter+1
        end do
        
        root_e2 = new_guess

    end if

    ! print*, '>>>', root_e2/1.323423_DBPR
    
    if (present(yR)) call SolveTOV_dark(h_in, p1_data, p2_data, e1_data, e2_data, e1, root_e2, r_lim, e1_disc, &
                                        R1, R2, M1, M2, yR= yR, in_guess= in_guess, gnu= gnu, fnu= fnu, pnu= pnu)
    
end subroutine newton_raphson

!----------------------------------------------------------------------------------------------------------------------------

subroutine SolveTOV_dark_lsoda(h, p1_data, p2_data, e1_data, e2_data, ec1, ec2, r_lim, r1_fin, r2_fin, m1_fin, m2_fin, &
                            yR, in_guess, gnu, fnu, pnu)
    
    type(lsoda_class) :: lsoda_two, lsoda_one
    real(DBPR), intent(in) :: h, ec1, ec2, r_lim
    real(DBPR), intent(in) :: p1_data(:), p2_data(:), e1_data(:), e2_data(:)
    real(DBPR), intent(out) :: m1_fin, m2_fin, r1_fin, r2_fin
    real(DBPR), intent(in), optional :: in_guess(:)
    real(DBPR), intent(out), optional :: yR, gnu, fnu, pnu
    integer :: i, r_flag, p1_flag, p2_flag, istate, itask, neq2, neq1
    real(DBPR) :: r, rout, x2(5), x1(3)
    real(DBPR) :: e1, e2, p1, p2, p1_limit, p2_limit, mc, pc1, pc2, phic, phi_fin, phi_far, y0, r_disc
    real(DBPR), allocatable :: h_array(:), r_array(:), m1_array(:), m2_array(:), &
                            p1_array(:), p2_array(:), e1_array(:), e2_array(:), de_dp_1(:), de_dp_2(:), &
                            phi_array(:), dphi_array(:), lamb_array(:), MATRIX(:,:), y_array(:), &
                            p_tot(:), e_tot(:), temp_array(:)
    real(DBPR) :: rtol, atol
    real(DBPR), parameter :: PI= 3.14159265358979323846_DBPR
    
    !------------------------------------------------------------------------------------------------------------------------
    
    neq2 = size(x2)
    neq1 = size(x1)
    call lsoda_two%initialize(TOV_two_lsoda, neq2, g=surface_two, ng=1, istate=istate,mxstep= 500000000)
    call lsoda_one%initialize(TOV_one_lsoda, neq1, g=surface_one, ng=1, istate=istate,mxstep= 500000000)

    rtol   = 5E-7_DBPR
    atol   = 1E-8_DBPR   
    itask  = 1
    istate = 1
    
    !------------------------------------------------------------------------------------------------------------------------

    call interpolate(e1_data, p1_data, ec1, pc1)
    call interpolate(e2_data, p2_data, ec2, pc2)
    
    if (ec1==0) pc1 = 0
    if (ec2==0) pc2 = 0

    mc   = 0.0_DBPR
    phic = 0.0_DBPR
    
    p1_limit = p1_data(1)         ! The limit replacing p<0
    p2_limit = p1_data(1)
    m1_fin   = 0.0_DBPR
    m2_fin   = 0.0_DBPR
    r1_fin   = 0.0_DBPR
    r2_fin   = 0.0_DBPR
    
    r    = 1E-15_DBPR                ! Starting point of radius. This increases by h in each step.
    x2   = [mc, mc, pc1, pc2, phic]  ! [Initial mass, Initial pressure, Initial phi]
    p1   = x2(3)                     ! p= initial pressure. This will change in the loop.
    p2   = x2(4)
    e1   = ec1                       ! e= initial energy density. This will change in the loop.
    e2   = ec2                       

    allocate(m1_array(0), m2_array(0))

    r_array   = [r]
    m1_array  = [mc]
    m2_array  = [mc]
    p1_array  = [pc1]
    p2_array  = [pc2]
    e1_array  = [ec1]
    e2_array  = [ec2]
    phi_array = [phic]

    if(ec1==0.0_DBPR) then
        x1   = [mc, pc2, phic]
        p1_flag = 1
    else if(ec2==0.0_DBPR) then
        x1   = [mc, pc1, phic]
        p2_flag = 1
    end if

    r_flag = 0
    p1_flag = 0
    p2_flag = 0
    rout = r+h    

    r_disc = -100.0_DBPR
    
    if (p1_flag==0 .and. p2_flag==0) then
            
        do
            call lsoda_two%integrate(x2, r, rout, rtol, [atol], itask, istate)
            
            r1_fin  = r
            r2_fin  = r
            m1_fin  = x2(1)
            m2_fin  = x2(2)
            p1      = x2(3)
            p2      = x2(4)
            phi_fin = x2(5)       
            
            ! print*, r, p1, p2, p1_limit
            
            if (p1<p1_limit .and. p1>0) p1_flag = 1
            if (p2<p2_limit .and. p2>0) p2_flag = 1
            
            if (p1<0) then
                p1      = p1_limit
                p1_flag = 1
            end if

            if (p2<0) then
                p2      = p2_limit
                p2_flag = 1
            end if
            
            call interpolate(p1_data, e1_data, p1, e1)
            call interpolate(p2_data, e2_data, p2, e2)
            
            if (present(yR)) then
                r_array   = [r_array, r]
                m1_array  = [m1_array, m1_fin]
                m2_array  = [m2_array, m2_fin]
                p1_array  = [p1_array, p1]
                p2_array  = [p2_array, p2]
                e1_array  = [e1_array, e1]
                e2_array  = [e2_array, e2]
                phi_array = [phi_array, phi_fin]
            end if
        
            if (p1_flag>0 .or. p2_flag>0) exit

            if ( istate<0 ) then
                write (6,12) istate
                12 format (///' Error halt in TOV... ISTATE = ', i3)
                stop 1
            else
                rout = rout + h
            end if
            
        end do

        if ((istate < 0 .or. istate /= 3 .or. lsoda_two%jroot(1) /= -1) .and. (r_flag==0.and.p1_flag==0.and.p2_flag==0)) then
            print*, 'ISTATE = ', istate
            error stop 'Root finding Failed'
        endif

        if (p2_flag==1) then
            x1 = [m1_fin+m2_fin, p1, phi_fin] !Dark Matter finished
            e2 = 0
            p2 = 0
        else if (p1_flag==1) then
            x1 = [m2_fin+m1_fin, p2, phi_fin] !Normal Matter finished
            e1 = 0
            p1 = 0
        end if
        
        if (r==rout) then
            rout=r+h
        end if

    end if
   
    itask  = 1
    istate = 1

    do 
        call lsoda_one%integrate(x1, r, rout, rtol, [atol], itask, istate)
        
        if (p2_flag==1) then
            
            r1_fin  = r
            m1_fin  = x1(1)-m2_fin
            p1      = x1(2)
            if (abs(p1-p1_limit)<atol .and. p1>0) p1_flag = 1
            if (p1<0) then
                p1      = p1_limit
                p1_flag = 1
            end if
            call interpolate(p1_data, e1_data, p1, e1)

        else if (p1_flag==1) then
            
            r2_fin  = r
            m2_fin  = x1(1)-m1_fin
            p2      = x1(2)
            if (abs(p2-p2_limit)<atol .and. p2>0) p2_flag = 1
            if (p2<0) then
                p2      = p2_limit
                p2_flag = 1
            end if
            call interpolate(p2_data, e2_data, p2, e2)

        end if
        
        ! print*, r, p1, p2
        
        phi_fin = x1(3)       
        
        if (r1_fin*1000 > r_lim) then
            ! Exit if radius is greater than r_lim since these stars won't be considered anyways.
            ! Flag for this condition is r_fin= r_lim+ 1000
            r1_fin = r_lim+1000.0_DBPR
            r2_fin = r_lim+1000.0_DBPR
            r_flag  = 1
            exit
        end if
        
        if (present(yR)) then
            r_array   = [r_array, r]
            m1_array  = [m1_array, m1_fin]
            m2_array  = [m2_array, m2_fin]
            p1_array  = [p1_array, p1]
            p2_array  = [p2_array, p2]
            e1_array  = [e1_array, e1]
            e2_array  = [e2_array, e2]
            phi_array = [phi_array, phi_fin]
        end if

        if (p1_flag==1 .and. p2_flag==1) exit
        
        if ( istate<0 ) then
            write (6,13) istate
            13 format (///' Error halt in TOV... ISTATE = ', i3)
            stop 1
        else
            rout = rout + h
        end if

    end do
        
    !Tidal Deformability-----------------------------------------------------------------------------------------------------
    
    if (present(yR)) then
        
        !φ array-------------------------------------------------------------------------------------------------------------
        
        phi_far   = (1./2.)*log(1.- (2.*(m1_fin+m2_fin)/max(r1_fin, r2_fin)))
        phi_array = phi_array+ (phi_far-phi_fin)

        allocate(h_array(size(r_array)), de_dp_1(size(r_array)), de_dp_2(size(r_array)), &
        dphi_array(size(r_array)), lamb_array(size(r_array)))

        do i= 1, size(r_array)-1
            h_array(i) = r_array(i+1)-r_array(i)
            de_dp_1(i) = (e1_array(i+1)-e1_array(i))/(p1_array(i+1)-(p1_array(i)))
            de_dp_2(i) = (e2_array(i+1)-e2_array(i))/(p2_array(i+1)-(p2_array(i)))
        end do

        h_array(size(r_array)) = h_array(size(r_array)-1)
        de_dp_1(size(r_array)) = de_dp_1(size(r_array)-1)
        de_dp_2(size(r_array)) = de_dp_2(size(r_array)-1)

        dphi_array = ((m1_array+m2_array) + 4*PI*(p1_array+p2_array)*(r_array**3))/&
            (r_array*(r_array - 2*(m1_array+m2_array)))

        lamb_array = log( (1.- 2.*(m1_array+m2_array)/r_array)**(-1./2.) )
        
        allocate(temp_array(size(r_array)), e_tot(size(r_array)), p_tot(size(r_array)))
        e_tot = e1_array+e2_array
        p_tot = p1_array+p2_array
        do i=1, size(r_array)
            if (isnan(de_dp_1(i))) then
                temp_array(i) = (e2_array(i)+p2_array(i))*de_dp_2(i)
            else if (isnan(de_dp_2(i))) then
                temp_array(i) = (e1_array(i)+p1_array(i))*de_dp_1(i)
            else
                temp_array(i) = (e1_array(i)+p1_array(i))*de_dp_1(i)+(e2_array(i)+p2_array(i))*de_dp_2(i)
            end if
        end do
        
        y0 = 2.0_DBPR
        y_array = [y0]

        do i=2,size(r_array)
            call rk4_tide(h_array(i-1), y0, r_array(i-1), e_tot(i-1), p_tot(i-1), temp_array(i-1), &
                        lamb_array(i-1), dphi_array(i-1))
            y_array = [y_array, y0]            
        end do
        yR = y_array(size(r_array))
        
    end if
    
    ! if (present(yR)) then
        do i=1, size(r_array)
            write(99, *) e1_array(i), p1_array(i), e2_array(i), p2_array(i), r_array(i), m1_array(i), m2_array(i)
        end do
    ! end if

    !Non Radial Modes---------------------------------------------------------------------------------------------------------

    if (present(gnu)) then
        
        allocate(MATRIX(13,size(r_array)))
        MATRIX(1,:)  = h_array
        MATRIX(2,:)  = r_array
        MATRIX(3,:)  = m1_array
        MATRIX(4,:)  = m2_array
        MATRIX(5,:)  = p1_array
        MATRIX(6,:)  = p2_array
        MATRIX(7,:)  = e1_array
        MATRIX(8,:)  = e2_array
        MATRIX(9,:)  = de_dp_1
        MATRIX(10,:) = de_dp_2
        MATRIX(11,:) = lamb_array
        MATRIX(12,:) = phi_array
        MATRIX(13,:) = dphi_array
        
        call dark_frequency(MATRIX, in_guess, r_disc, gnu, fnu, pnu)
        
    end if

    !Returning---------------------------------------------------------------------------------------------------------------

    m1_fin = m1_fin*1000.0_DBPR/1.4766_DBPR    ! Output mass in units of solar mass
    m2_fin = m2_fin*1000.0_DBPR/1.4766_DBPR
    r1_fin = r1_fin*1000.0_DBPR                ! Output radius in kilometers
    r2_fin = r2_fin*1000.0_DBPR
    
!----------------------------------------------------------------------------------------------------------------------------

contains

    subroutine TOV_two_lsoda(self_, Neq_, r_, x2_, Ydash_, ierr_)
        implicit none
        class(lsoda_class), intent(inout) :: self_
        integer, intent(in) :: Neq_
        real(DBPR), intent(in) :: r_
        real(DBPR), intent(in) :: x2_(Neq_)
        real(DBPR), intent(out) :: Ydash_(Neq_)
        integer, intent(out) :: ierr_
        real(DBPR) :: e1_, p1_, e2_, p2_
        ierr_ = 0
        self_%iwork = self_%iwork
        
        if(x2_(3)> p1_data(size(p1_data))) then
            p1_ = p1_data(size(p1_data))
        else
            p1_ = x2_(3)
        end if

        if(x2_(4)> p2_data(size(p2_data))) then
            p2_ = p2_data(size(p2_data))
        else
            p2_ = x2_(4)
        end if

        call interpolate(p1_data, e1_data, p1_, e1_)
        call interpolate(p2_data, e2_data, p2_, e2_)
        
        Ydash_ = TOV_dark_func(e1_, e2_, r_, x2_)
        
    end subroutine TOV_two_lsoda

    !------------------------------------------------------------------------------------------------------------------------

    subroutine TOV_one_lsoda(self_, Neq_, r_, x1_, Ydash_, ierr_)
        implicit none
        class(lsoda_class), intent(inout) :: self_
        integer, intent(in) :: Neq_
        real(DBPR), intent(in) :: r_
        real(DBPR), intent(in) :: x1_(Neq_)
        real(DBPR), intent(out) :: Ydash_(Neq_)
        integer, intent(out) :: ierr_
        real(DBPR) :: e_, p_
        ierr_ = 0
        self_%iwork = self_%iwork
        
        if (p2_flag==1) then
            if(x1_(2)> p1_data(size(p1_data))) then
                p_ = p1_data(size(p1_data))
            else
                p_ = x1_(2)
            end if
            call interpolate(p1_data, e1_data, p_, e_)
        else if (p1_flag==1) then
            if(x1_(2)> p2_data(size(p2_data))) then
                p_ = p2_data(size(p2_data))
            else
                p_ = x1_(2)
            end if
            call interpolate(p2_data, e2_data, p_, e_)
        end if

        Ydash_ = TOV_dark_func(e_, e_, r_, x1_)
        
        
    end subroutine TOV_one_lsoda
    
    !------------------------------------------------------------------------------------------------------------------------

    subroutine surface_two(self_, neq_, r_, x2_, ng, gout, ierr_) !The root function, ng is number of roots
        implicit none
        class(lsoda_class), intent(inout) :: self_
        integer, intent(in) :: neq_, ng
        real(DBPR), intent(in) :: r_, x2_(neq_)
        integer, intent(out) :: ierr_
        real(DBPR), intent(out) :: gout(ng)
        real(DBPR) :: temp_
        ierr_ = 0
        temp_ = r_
        self_%iwork = self_%iwork

        gout(1) = min(x2_(3)-p1_limit,  x2_(4)-p2_limit)
        
    end subroutine surface_two

    !------------------------------------------------------------------------------------------------------------------------

    subroutine surface_one(self_, neq_, r_, x1_, ng, gout, ierr_) !The root function, ng is number of roots
        implicit none
        class(lsoda_class), intent(inout) :: self_
        integer, intent(in) :: neq_, ng
        real(DBPR), intent(in) :: r_, x1_(neq_)
        integer, intent(out) :: ierr_
        real(DBPR), intent(out) :: gout(ng)
        real(DBPR) :: temp_
        ierr_ = 0
        temp_ = r_
        self_%iwork = self_%iwork

        if (p2_flag==1) then
            gout(1) = x1_(2)-p1_limit
        else if (p1_flag==1) then
            gout(1) = x1_(2)-p2_limit
        end if
        
    end subroutine surface_one

    !------------------------------------------------------------------------------------------------------------------------

end subroutine SolveTOV_dark_lsoda

!----------------------------------------------------------------------------------------------------------------------------

subroutine newton_raphson_lsoda(p1_data, p2_data, e1_data, e2_data, e1, e2, r_lim, dark_frac, root_e2, &
                        R1, R2, M1, M2, yR, in_guess, gnu, fnu, pnu)

    real(DBPR), intent(in) :: p1_data(:), p2_data(:), e1_data(:), e2_data(:), e1, e2, r_lim, dark_frac
    real(DBPR), intent(in), optional :: in_guess(:)
    real(DBPR), intent(out) ::  root_e2, R1, R2, M1, M2
    real(DBPR), intent(out), optional :: yR, gnu, fnu, pnu
    integer :: iter
    real(DBPR) :: h, tol, error, guess, new_guess, f_pre, f
    
    tol   = 1e-4_DBPR     ! Tolerance for closeness of values in Newton-Raphson
    error = 10                ! Temporary initial value for error
    guess = e2
    iter  = 0
        
    iter  = 0
    guess = e2
    h     = 1000000
    
    777 if (guess-h <0 ) then
        h= h/10 
        goto 777
    end if
    
    do while (error>tol)
        
        new_guess = guess-h
        call SolveTOV_dark_lsoda(5E-4_DBPR, p1_data, p2_data, e1_data, e2_data, e1, new_guess, r_lim, R1, R2, M1, M2)
        
        f_pre = f
        f = (M2/(M1+M2))-dark_frac
        
        ! print*, h, new_guess/1.323423_DBPR, R1, R2, M2/(M1+M2), iter
        
        if (f<0 .or. new_guess-h<0) then
            h= h/10
        else if (f>1E-12) then
            guess= new_guess
        end if
        
        error= abs(f-f_pre)
        iter= iter+1
    end do
    
    root_e2 = new_guess

    ! print*, '>>>', root_e2/1.323423_DBPR
    
    if (present(yR)) call SolveTOV_dark_lsoda(1E-5_DBPR, p1_data, p2_data, e1_data, e2_data, e1, root_e2, r_lim, &
                                        R1, R2, M1, M2, yR= yR, in_guess= in_guess, gnu= gnu, fnu= fnu, pnu= pnu)
    
end subroutine newton_raphson_lsoda


end module TOV_dark