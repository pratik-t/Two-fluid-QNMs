module IO
    implicit none
    integer, parameter:: DBPR=kind(0.d0) !double precision
contains

!-----------------------------------------------------------------------------------------------------------------------------

subroutine read_eos(ID, filename, a_data, b_data, c_data)

    character(len=*), intent(in) :: filename
    integer, intent(in) :: ID
    character(len=1000):: line
    integer :: n_cols, pos, iline
    real(DBPR), allocatable, intent(out) :: a_data(:), b_data(:)
    real(DBPR), allocatable, intent(out), optional :: c_data(:)
    real(DBPR) :: a_point, b_point, c_point
        
    allocate(a_data(0), b_data(0))
    if (present(c_data)) allocate(c_data(0))

    open(unit= ID+7, file= filename, status= 'old', action= 'read')
    read(ID+7,'(a)') line
    pos= 1
    n_cols= 0

    ! Finds number of columns
    do
        iline = verify(line(pos:), ' ')  ! Find next non-blank.
        if (iline == 0) then 
            exit                         ! No word found.
        end if
        n_cols = n_cols + 1              ! Found something.
        pos = pos + iline - 1            ! Move to start of the word.
        iline = scan(line(pos:), ' ')    ! Find next blank.
        if (iline == 0) then
            exit                         ! No blank found.
        end if
        pos = pos + iline - 1            ! Move to the blank.
    end do 
    close(ID+7)    

    if (n_cols<2) then
        print*, 'The input files contain only one column. Make sure the columns are space separated and not tab separated'
        stop
    end if

    open(unit= ID+7, file= filename, status= 'old', action= 'read')
    
    do
       
        if (n_cols==3) then
            read(ID+7,*, end=100) a_point, b_point, c_point

        else if (n_cols==2) then
            read(ID+7,*, end=100) a_point, b_point
            c_point = -10.0_DBPR
        end if

        a_data = [a_data, a_point]
        b_data = [b_data, b_point]
        if (present(c_data)) c_data = [c_data, c_point]

    end do

    100 close(ID+7) 
        
end subroutine read_eos

!-------------------------------------------------------------------------------------------------------------------------------

subroutine make_dark_eos(n_points, step, m, cs, cv, emax, p, e)
    
    integer, intent(in) :: n_points
    real(DBPR), intent(in) :: step, m, cs, cv, emax
    real(DBPR), allocatable, intent(out) :: p(:), e(:)
    integer :: i
    real(DBPR) :: kf, meff, meff_guess, z, n, p_point, e_point, p_point_prev
    real(DBPR), parameter :: PI= 3.14159265358979323846_DBPR, hc3= (197.33_DBPR)**3

    allocate(e(0), p(0))
    
    p_point = -99999

    kf= 0
    meff_guess = m
    
    do i= 1, n_points

        kf = kf+ i/step
        
        if (cs==0) then 
            meff = m
        else
            call find_meff(m, meff_guess, cs, kf, meff)
            meff_guess = meff
        end if
       
        ! print*, m, meff
       
        z = kf/meff
        n = (kf**3)/(3*(PI**2)*hc3)

        if (cs==0) then
            e_point = ((meff**4/(8*pi**2))*(((2*z**3+ z)*sqrt(1+ z**2))- asinh(z))/hc3)+ (0.5_DBPR*(cv**2)*(n**2)*hc3)
            p_point = ((meff**4/(24*pi**2))*(((2*z**3- 3*z)*sqrt(1+ z**2))+ 3*asinh(z))/hc3)+ (0.5_DBPR*(cv**2)*(n**2)*hc3)
        else

            p_point_prev = p_point

            e_point = ((meff**4/(8*pi**2))*(((2*z**3+ z)*sqrt(1+ z**2))- asinh(z))/hc3)+ (0.5_DBPR*(cv**2)*(n**2)*hc3)+ &
                        (((m-meff)**2)/(2*(cs**2)*hc3)) 
            p_point = ((meff**4/(24*pi**2))*(((2*z**3- 3*z)*sqrt(1+ z**2))+ 3*asinh(z))/hc3)+ (0.5_DBPR*(cv**2)*(n**2)*hc3)- &
                        (((m-meff)**2)/(2*(cs**2)*hc3)) 

            ! print*, p_point, ((meff**4/(24*pi**2))*(((2*z**3- 3*z)*sqrt(1+ z**2))+ 3*asinh(z))/hc3), (0.5_DBPR*(cv**2)*(n**2)*hc3),&
            !         (((m-meff)**2)/(2*(cs**2)*hc3)), m, meff
            ! print*, e_point, p_point, n

            if ((p_point<0 .and. p_point_prev>0).and. n>1E-8_DBPR) then
                p= [0]
                e= [0]
                exit
            end if
            

        end if
       
        if (p_point>=1E-15_DBPR) then
            p= [p, p_point]
            e= [e, e_point]
            ! print*, '>', e_point, p_point, n
        end if

        if (e_point>= emax) then
            exit
        end if

    end do
    
end subroutine make_dark_eos

!-------------------------------------------------------------------------------------------------------------------------------

subroutine find_meff(m, meff_guess, cs, kf, meff)

    real(DBPR), intent(in) :: m, meff_guess, cs, kf
    real(DBPR), intent(out) :: meff
    integer :: iter, flag
    real(DBPR) :: h, tol, error, prev_guess, guess, new_guess, f_pre, f, fh, fdash

    h     = 1e-10_DBPR        ! This is the order of error on derivative for Newton-Raphson method

    123 tol   = 1e-7_DBPR     ! Tolerance for closeness of values in Newton-Raphson
    error = 10                ! Temporary initial value for error
    guess = meff_guess
    iter  = 0
    flag  = 0

    do while (error>tol)

        fh = meff_func(guess-h, m, cs, kf)
        f  = meff_func(guess, m, cs, kf)
        
        prev_guess = guess

        fdash     = (f- fh)/h
        new_guess = guess- (f/fdash)
        error     = f
        meff      = guess
        guess     = new_guess
        iter      = iter+1
        
        ! print*, m, guess, f, iter

        if (abs(fdash)==0 .and. h<1E-4_DBPR) then
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
        guess = meff_guess
        h     = 1000000
        
        777 if (guess-h <0 ) then
            h= h/10 
            goto 777
        end if

        do while (error>tol)
            
            new_guess = guess-h
            f_pre = f
            f = meff_func(new_guess, m, cs, kf)

            ! print*,'>', m, new_guess, f, iter
            
            if (f<0 .or. new_guess-h<0) then
                h= h/10
            else if (f>1E-12) then
                guess= new_guess
            end if
            
            error= abs(f-f_pre)
            iter= iter+1
        end do
        
        meff = new_guess

    end if

contains

    function meff_func(meff_, m_, cs_, kf_) result(trance) 

        real(DBPR), intent(in) :: meff_, m_, cs_, kf_
        real(DBPR) :: z_, trance
        real(DBPR), parameter :: PI= 3.14159265358979323846_DBPR
        
        z_ = kf_/meff_
        trance = meff_- m_+ (cs_**2)*((meff_**3)/(2*PI**2))*(z_*sqrt(1+ z_**2)- asinh(z_))
    
    end function meff_func
    
end subroutine find_meff

!-------------------------------------------------------------------------------------------------------------------------------
    
subroutine make_filename(fn, ID1, ID2, ID3, dir_in, dir_out, out_text, input_fn, output_fn)
    character(len=*), intent(in) :: fn, ID1, ID2, ID3, dir_in, dir_out, out_text
    character(len=*), intent(out) ::  input_fn, output_fn
    logical :: exists
    character(len=128) :: in_path, out_path
    
    INQUIRE( FILE = 'DATA\'//dir_in(1:len_trim(dir_in))//'/.', EXIST = exists )
    if (exists .eqv. .true.) then
        in_path = dir_in(1:len_trim(dir_in))
    else 
        print*, "Directory `"// dir_in // "' does not exist. Kindly ensure that all equation of state directories" // &
                " are present in the 'DATA' directory"
        stop
    end if
    input_fn = 'DATA\'//in_path(1:len_trim(in_path))//'\'//fn(1:len_trim(fn))

    ! Making directory DATA if it doesn't exist, and specifying the out_path to it where the data will be stored.

    INQUIRE( FILE = 'DATA\'//dir_out(1:len_trim(dir_out))//'\.', EXIST = exists )
    if (exists .eqv. .false.) then
        call system("mkdir "//'DATA\'//dir_out(1:len_trim(dir_out)))
    end if
    out_path = dir_out(1:len_trim(dir_out))//'\'

    output_fn = 'DATA\'//out_path(1:len_trim(out_path))//fn(1:len_trim(fn)-4)//&
            '_'//ID1(1:len_trim(ID1))//'_'//ID2(1:len_trim(ID2))//'_'//ID3(1:len_trim(ID3))//out_text(1:len_trim(out_text))//'.dat'

end subroutine make_filename
    
!-------------------------------------------------------------------------------------------------------------------------------

subroutine read_TOV_data(MATRIX, h, r, m1, m2, p1, p2, e1, e2, de_dp_1, de_dp_2, lamb, phi, dphi)
    
    real(DBPR), allocatable, dimension(:,:), intent(in) :: MATRIX
    real(DBPR), allocatable, dimension(:), intent(out) :: h, r, m1, m2, p1, p2, e1, e2, &
                                                        de_dp_1, de_dp_2, lamb, phi, dphi

    h       = MATRIX(1,:)
    r       = MATRIX(2,:) 
    m1      = MATRIX(3,:) 
    m2      = MATRIX(4,:) 
    p1      = MATRIX(5,:) 
    p2      = MATRIX(6,:) 
    e1      = MATRIX(7,:) 
    e2      = MATRIX(8,:)
    de_dp_1 = MATRIX(9,:)
    de_dp_2 = MATRIX(10,:)
    lamb    = MATRIX(11,:) 
    phi     = MATRIX(12,:) 
    dphi    = MATRIX(13,:)             
    
end subroutine read_TOV_data

!-------------------------------------------------------------------------------------------------------------------------------

subroutine interpolate(xdata, ydata, xval, yval) 

    real(DBPR), intent(in) :: xval
    real(DBPR), intent(in) :: xdata(:), ydata(:)
    real(DBPR), intent(out) :: yval
    integer :: data_size, start, mid, finish, range, low, high
    real(DBPR) :: slope

    data_size = size(xdata)
    
    start  = 1
    finish = data_size
    range  = finish-start
    mid    = (start+finish)/2
    
    ! Do a binary search of data to find the bounds for given xval
    
    if (xval>xdata(finish) .and. xval-xdata(finish)>1E-10) then
        print*, 'Value out of range', xval, '>', xdata(finish)
        
    else 
        do while (xdata(mid) /= xval .and. range>0)
            if (xval> xdata(mid)) then
                start= mid+1
            else 
                finish= mid-1
            end if
            range= finish-start
            mid= (start+finish)/2
        end do
        
        if (xdata(mid-1)<xval .and. xval<xdata(mid)) then
            low= mid-1
        else
            low= mid
        end if

        high= low+1
        
        ! print*, 'Value', xval, 'is between: ', xdata(low), ', ', xdata(high)
        
        if (xdata(low)==xval) then
            yval= ydata(low)
            
        else
            ! slope = (log10(ydata(high))-log10(ydata(low)))/(log10(xdata(high))-log10(xdata(low)))
            ! yval  = log10(ydata(low))+ slope*(log10(xval)- log10(xdata(low)))
            ! yval  = 10**yval
            slope = (ydata(high)-ydata(low))/(xdata(high)-xdata(low))
            yval  = ydata(low)+ slope*(xval- xdata(low))    
        end if
    
    end if
        
end subroutine interpolate

!-------------------------------------------------------------------------------------------------------------------------------

subroutine continuity_check(x_data, z_data, disc_start, disc_end)

    real(DBPR), allocatable, intent(inout) :: x_data(:), z_data(:)
    integer, intent(out) :: disc_start, disc_end
    integer :: i
        
    do i= 1, size(x_data)-1
        if (z_data(i)==0) then
            disc_start  = i
            disc_end    = i+1
            x_data(i+1) = x_data(i)+1E-10_DBPR
            exit
        else
            disc_start = -10
            disc_end   = -10    
        end if       
    end do

end subroutine continuity_check

!-------------------------------------------------------------------------------------------------------------------------------


end module IO
