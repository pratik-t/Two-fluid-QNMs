program Nonradial_modes
    use IO
    use TOV_dark
    use omp_lib
    implicit none
    
!-------------------------------------------------------------------------------------------------------------------------------!
    
    !---------------------------------------------------------------------------------------!
    !                                                                                       !
    ! EOS file should have the first column as pressure and second one as energy density    !
    ! third column should be number density                                                 !
    !                                                                                       !
    ! 1 mev/fm^3 = 1.7827 × 10^12 g/cm^3                                                    !
    ! 1 mev/fm^3 = 1.6022 × 10^33 dyne/cm^3                                                 !
    ! 1 g/cm^3= 7.4237 × 10^(-19) /km^2                                                     !
    ! 1 dyne/cm^3= 8.2601 × 10^(-40) /km^2                                                  !
    !                                                                                       !
    ! From these we have the following relations: (/Mm^2= 10^6 /km^2)                       !
    ! 1 mev/fm^3= 1.323423 /Mm^2 (energy density in Megameters= 10^3 km)                    !
    ! 1 mev/fm^3= 1.323423 /Mm^2 (pressure in Megameters)                                   !
    ! Mass of sun in km= 1.4766 km                                                          !
    !                                                                                       !
    !---------------------------------------------------------------------------------------!

    !------------------------------------------!
    !                                          !
    ! ω [1/Mm]= 299.7905585 ω [1/sec]          !
    ! ν= ω/2π                                  !
    ! ν [kHz]= 299.7905585/2000π ω [1/Mm]      !
    ! ν [kHz]= 47.71314928e-3 ω [1/Mm]         !
    !                                          !
    !------------------------------------------!

    !---------------------------------------------------------------!
    !                                                               !
    ! Compactness   = Mass*1.4766/Radius                            !
    ! Surf_redshift = ((1.- (2*Mass*1.4766/Radius))**(-1./2.))-1.   !
    ! tidal_def     = (2./3.)*k2*(Radius/(Mass*1.4766))**5          !
    !                                                               !
    !---------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------------
    
    !$ integer :: nThreads
    character(len=1) :: creturn
    character(len= 128) :: ID_frac, ID_cs, ID_cv, ID_mass, filename, input_filename, output_filename
    integer :: t_start, t_stop, clock_max, clock_rate, i, j, start1, end1, start2, end2, counter, disc_start, disc_end, dmf, cs, cv
    integer, allocatable :: cs_vals(:), cv_vals(:)
    real(DBPR), allocatable :: DM_fracs(:), p1_data(:), p2_data(:), e1_data(:), e2_data(:), &
                                temp_data(:), MATRIX_out(:,:)
    real(DBPR) :: h_in(3), h_twof, h_norm, h_dark, DM_mass, dark_frac, frac, root, r_lim, &
                M1, M2, R1, R2, yR, in_guess(3), gnu, fnu, pnu, e1_disc, p1_disc, NaN=0.0_DBPR

!-------------------------------------------------------------------------------------------------------------------------------
    !$ nThreads = omp_get_max_threads()
    !$ call omp_set_num_threads(nThreads)
    
    creturn = achar(13)  ! Carriage return
    
    write(*, '(/a)', advance='no') 'Enter EOS filename: '  
    read(*, *) filename

    write(*, '(/a)', advance='no') 'Enter mass of dark matter particle [GeV]: '
    read(*, *) ID_mass
    read(ID_mass, *) DM_mass

! Finding if directory EOS exists, making output directory DATA_MR, and making input and output filenames (with path) with the output filenames ending with mass of DM particle. 

    call make_filename(filename,'','','','EOS', 'DDATA_'//ID_mass(1:len_trim(ID_mass))//'_GeV', &
                        '_MR', input_filename, output_filename)
    
!-------------------------------------------------------------------------------------------------------------------------------
    
    call read_eos(1, input_filename, p1_data, e1_data, temp_data)

    ! Converting from Mev/fm^3 to /Mm^2
    p1_data   = p1_data*1.323423_DBPR
    e1_data   = e1_data*1.323423_DBPR

!------------------------------------------------------------------------------------------------------------------------------ 
! Rest of the inputs 
    
    if (temp_data(1)>0) then
        call continuity_check(p1_data, temp_data, disc_start, disc_end)
        if (disc_start>0 .and. disc_end>0) then
            e1_disc = e1_data(disc_end)
            p1_disc = p1_data(disc_end)
        else
            e1_disc = -10.0_DBPR
            p1_disc = -10.0_DBPR
        end if
    end if
    
    if (e1_disc>0 .and. p1_disc>0) then

        write(*, '(/a, f10.5, a, f10.5, a, f10.5, a)') 'EOS has density discontinuity at: Energy Density= ', &
        e1_data(disc_start), ' <-> ', e1_data(disc_end),' [MeV/fm^3] and Pressure= ', p1_disc, ' [MeV/fm^3]'

        e1_disc = e1_disc*1.323423
        p1_disc = p1_disc*1.323423

    end if

    ! The desired limiting radius in [Km] beyond which MR won't be calculated
    r_lim= 60

    write(*, '(/a)', advance='no') 'Press "Enter" to start: '
    read(*,*)
    write(*, *)

    call system_clock (t_start, clock_rate, clock_max)
    
!-----------------------------------------------------------------------------------------------------------------------------

    do i= 1, size(e1_data)
        if (e1_data(i)/1.323423>120.) then 
            start1= i
            exit
        end if
    end do
    
    end1   = start1+1!size(e1_data)

!-----------------------------------------------------------------------------------------------------------------------------
    
    h_twof    = 5E-6_DBPR !5E-6
    h_norm    = 5E-6_DBPR !1E-6
    h_dark    = 5E-6_DBPR !5E-6
    h_in      = [h_twof, h_norm, h_dark]
    in_guess  = (/-10, -10, -10/)
    
    !DM Fraction cannot be 0 or 1. 

    ! DM_fracs = [0.01_DBPR, 0.02_DBPR, 0.03_DBPR, 0.04_DBPR, 0.05_DBPR, 0.06_DBPR, 0.07_DBPR, 0.08_DBPR, 0.09_DBPR, 0.10_DBPR, &
    ! 0.11_DBPR, 0.12_DBPR, 0.13_DBPR, 0.14_DBPR, 0.15_DBPR, 0.16_DBPR, 0.17_DBPR, 0.18_DBPR, 0.19_DBPR, 0.20_DBPR]
    ! cs_vals  = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50]
    ! cv_vals  = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50]

    DM_fracs = [0.00_DBPR]
    cs_vals  = [0]
    cv_vals  = [0]


    !$OMP PARALLEL DO DEFAULT(PRIVATE) COLLAPSE(3) SHARED(filename, DM_mass, ID_mass, DM_fracs, cs_vals, cv_vals, &
    !$OMP& h_in, start1, end1, p1_data, e1_data, e1_disc, r_lim) SCHEDULE(DYNAMIC)

    do dmf= 1, size(DM_fracs)
    do cs = 1, size(cs_vals)
    do cv = 1, size(cv_vals)

        counter = 10000*dmf+100*cs+cv
        
        print*, counter
        dark_frac = DM_fracs(dmf)
        write(ID_frac, '(f4.2)') dark_frac
        write(ID_cs, '(i3.3)') cs_vals(cs)
        write(ID_cv, '(i3.3)') cv_vals(cv)

        call make_filename(filename, ID_frac, ID_cs, ID_cv, 'EOS', 'DDATA_'//ID_mass(1:len_trim(ID_mass))//'_GeV', &
                    '_data', input_filename, output_filename)
            
        call make_dark_eos(500000, 100.0_DBPR, DM_mass*1000, real(cs_vals(cs)/1000.0, kind=8), &
                            real(cv_vals(cv)/1000.0, kind=8), 5E10_DBPR, p2_data, e2_data)
        print*, '>>>>>', e2_data(size(e2_data))
        if (size(e2_data)==1) then
            if (allocated(MATRIX_out)) deallocate(MATRIX_out)
            allocate(MATRIX_out(1,10))
            MATRIX_out(1, :)= [-1, -1, -1, -1, -1, -1, -1]
            goto 789
        end if

        ! Converting from Mev/fm^3 to /Mm^2
        p2_data   = p2_data*1.323423_DBPR
        e2_data   = e2_data*1.323423_DBPR

        start2 = 10
        end2   = size(e2_data)

        if (allocated(MATRIX_out)) deallocate(MATRIX_out)
        if (dark_frac==1) then
            allocate(MATRIX_out(end2-start2,10))
        else
            allocate(MATRIX_out(end1-start1,10))
        end if

        do i= start1+1, end1 
        
            do j= start2+1, end2
                
                ! call SolveTOV_dark([1E-4_DBPR, 1E-4_DBPR, 1E-3_DBPR], p1_data, p2_data, e1_data, e2_data, &
                !                     e1_data(i), e2_data(j), r_lim, -10.0_DBPR, R1, R2, M1, M2)

                call SolveTOV_dark_lsoda(5E-4_DBPR, p1_data, p2_data, e1_data, e2_data, e1_data(i), e2_data(j), r_lim, &
                                        R1, R2, M1, M2)
                
                frac = M2/(M1+M2)
                print*, '>', dark_frac, e1_data(i)/1.323423_DBPR, e2_data(j)/1.323423_DBPR, frac

                if (frac-dark_frac>=0) then

                    call newton_raphson_lsoda(p1_data, p2_data, e1_data, e2_data, e1_data(i), e2_data(j), &
                                        r_lim, dark_frac, root, R1, R2, M1, M2, yR, in_guess, gnu, fnu, pnu)
                    
                    if (frac<dark_frac+0.01_DBPR) start2 = j
                    
                    if (R1 > r_lim) then
                        MATRIX_out(i-start1, :)= [-1, -1, -1, -1, -1, -1, -1]
                    else

                        if (isnan(fnu) .or. isnan(pnu) .or. isnan(gnu)) then
                            in_guess = [-10, -10, -10]
                        else
                            in_guess = [fnu, pnu, gnu]/47.71345159236942258889E-3_DBPR
                        end if
                        
                        if (gnu<0.0_DBPR) gnu = NaN/NaN

                        MATRIX_out(i-start1, :)= [e1_data(i)/1.323423_DBPR, root/1.323423_DBPR, &
                                                R1, R2, M1, M2, yR, gnu, fnu, pnu]
                        
                        print*, output_filename(1:len_trim(output_filename)), e1_data(i)/1.323423_DBPR, &
                        root/1.323423_DBPR, M2/(M1+M2), M1+M2, fnu, pnu

                    end if
                    
                    exit  

                end if
                
            end do
        
        end do
        
        789 open(unit= counter, file= output_filename, status= 'replace', action= 'write')
        
        write(counter, '(10a)') 'Central ε Normal [MeV/fm^3] :: ', 'Central ε dark [MeV/fm^3] :: ', &
        'R Normal [KM] :: ', 'R Dark [Km] :: ', 'M Normal [M_sun] :: ', 'M Dark [M_sun] :: ', 'yR :: ', &
        'g mode [KHz] :: ', 'f mode [KHz] :: ', 'p1 mode [KHz]'
        
        ! write(counter, '(2a)') 'ε dark [MeV/fm^3] :: ', 'p dark [MeV/fm^3]'
        ! do i= 1, size(e2_data)
        !     write(counter,'(2(E20.10))') e2_data(i), p2_data(i)
        ! end do

        do i= 1, size(MATRIX_out, 1)
            if (MATRIX_out(i,1).ne.-1) write(counter, '(10(E20.10))') Matrix_out(i, :)
        end do

    end do
    end do
    end do
    
    !$OMP END PARALLEL DO
    
    !-------------------------------------------------------------------------------------------------------------------------
    
    call system_clock (t_stop, clock_rate, clock_max)

    write(*,*) 
    write(*,*) 'Elapsed time is: ', real(t_stop-t_start)/real(clock_rate)

end program Nonradial_modes