!------------------------------------------------------------------------
!  To calculate the Raman tensor using finite element analysis
!  following the procedure as suggested by Umari and Pasquarello,
!  Diamond & Related Materials 14, 1255 (2005).
!  required files: (1) dXidR.dat containing the double derivative of
!                      of the Force vector with respect to
!                      to Electric field Ei and Ej
!                      Here we consider only the zz component of the
!                      dXidR tensor  (The code can be used for each
!                      component of the tensor separately )
!                  Format: 1 1
!                          atomicspecies  dXidR(1) dXidR(1) dXidR(1)
!                          1 2
!                          atomicspecies  dXidR(1) dXidR(1) dXidR(1) 
!                          .. and so on for all components
!                  (2) the file containing the eigenvectors of the
!                      Dynamical matrix ( use the fileig key in the
!                      dynmat.x input)
!                  (3) input file containing:
!                       line # 1:  Number of Atoms, no of atomic species
!                       line # 2:  Mass(i), i = 1,N_species
!------------------------------------------------------------------------

      program main
      implicit none
      REAL, allocatable :: dXidR ( :, :, :, :)
      REAL, allocatable :: amass (:)
      REAL*8, allocatable :: thzfreq(:), cmfreq(:)
      INTEGER, allocatable :: ityp (:)
      COMPLEX*8, allocatable :: eig( :, :, :)

      REAL*8, allocatable :: RamTens(:,:,:), DiffRamCross(:)
      
      REAL :: omega, omega_au
      INTEGER :: nat, ntype, nmodes, modnum, comp, rsp
      REAL*8 :: tfreq, cfreq 
      INTEGER :: i, j, k, at_ind, alpha, beta, mod_ind
      COMPLEX*8 :: e1, e2, e3
      REAL :: q (3), q1, q2, q3

      REAL*8 :: lam, factor, Iraman, laserfreq, const
      REAL*8 :: A_n, B_n2, depol_ratio

      INTEGER ::  ibin, Nbin
      REAL ::  f, fval, famp, dbin, degauss, freq_min, freq_max, dosval
      REAL, allocatable :: dos(:)
      
      INTEGER :: cmp1, cmp2
    
      CHARACTER (len=100) :: tmp, line   

!---------------------------------------------!
      INTERFACE 

      FUNCTION lorentz (e,eval,dbin,degauss)
      REAL :: lorentz
      REAL, INTENT(IN) :: e,eval,dbin,degauss
      END FUNCTION

      END INTERFACE
!---------------------------------------------!

      lam = 512.15   ! 512.15 nm - wave length of incident light in nm (Sodium D line)
      laserfreq = 2.99792*10**5/lam   ! in THz                                     

      comp = 3       ! All 3 directions
      degauss = 10    ! Lorentzian broadening in cm-1
      dbin = 0.2     ! binwidth of Raman spectra in cm-1

!---------------------------------------------!      
!     OPEN INPUT FILE      
      open (21, file='input.in', status='unknown')     
!     no of atomic species
      read (21,*) nat, ntype
!     atomic mass of each species
      allocate (amass (ntype))      
      read (21,*) (amass(i), i=1,ntype)
      close(21)
!---------------------------------------------!
!     OPEN  dXidR.dat file 
      allocate (ityp (nat))
      allocate (dXidR ( comp, comp, nat, 3))
      open (22, file='dXidR.dat', status='unknown')

      
      do i = 1, comp
       do j = 1, comp
        read (22,*) cmp1, cmp2
        do at_ind = 1, nat
         read (22,*) ityp(at_ind), (dXidR(i,j,at_ind,k), k=1,3)
        enddo
       enddo
      enddo

      close(22)
  
!     OPEN fileout.eig 
      open (23, file='vector.eig', action='READ', status='OLD')
      read (23, '(a50)') tmp

      read (23, '(a50)') tmp

      read (23, '(a100)') line
      call get_qvec(line, q(1), q(2), q(3))
!      write (6,'(1x,''q = '',3f12.4)') q(1), q(2), q(3)

      read (23, '(a50)') tmp

      nmodes=3*nat
      allocate (thzfreq(nmodes),cmfreq(nmodes))
      allocate (eig(nmodes,nat,3))
      do i = 1, nmodes
      read (23, '(a100)') line
      call get_mode(line, modnum, tfreq, cfreq)
        thzfreq(i)=tfreq
        cmfreq(i)=cfreq
!        write (6,9010) i, thzfreq(i), cmfreq(i)
        do at_ind = 1, nat
          read (23,'(a100)') line
          call get_eigvec(line, e1, e2, e3)
          eig(i,at_ind,1) = e1
          eig(i,at_ind,2) = e2
          eig(i,at_ind,3) = e3
        enddo
      enddo
      read (23, '(a100)') tmp
      close (23)
      !---------------------------------------------!

      allocate (RamTens(nmodes,comp,comp))
      allocate (DiffRamCross(nmodes))
      RamTens(:,:,:)=0.00
 
      do mod_ind = 1, nmodes
      do alpha = 1, comp
      do beta = 1, comp
       do at_ind = 1, nat
        j = ityp(at_ind)
        do k = 1, 3
         RamTens(mod_ind, alpha, beta) = RamTens(mod_ind, alpha, beta) &
       - ( (dXidR(alpha,beta,at_ind,k)*REAL(eig(mod_ind,at_ind,k))) & 
           /sqrt(amass(j)) &
          )
        enddo
       enddo    
!       RamTens(mod_ind,alpha,beta) = RamTens(mod_ind,alpha,beta) & 
!                                     / omega 
      enddo
      enddo
      enddo

      !---------------------------------------------!

      open (24, file='RAMAN_TENSOR.OUT',action='WRITE')
      do mod_ind = 1, nmodes
       write (24,*)'---------------------------------------'
       write (24,'(3x,I3,3x,f15.6,1x)') mod_ind, cmfreq(mod_ind)
       write (24,*)'---------------------------------------'
       do alpha = 1, comp
        write (24,'(3x,3f10.6,3x)') & 
          (RamTens(mod_ind,alpha,beta), beta=1,comp)
       enddo
       write (24,*)'---------------------------------------'
      enddo
      close (24)
       
      open (25, file='DIFF_RAMAN_CROSS-SECTION.OUT',action='WRITE')
      write (25,*) ' #   Freq (cm-1)  Freq(THz)   & 
                     Raman Activity   Depol. ratio'
      do mod_ind = 1, nmodes

!	const =  2*PI^2 * h / 45.0*c  !hw unit:kg*cm  pi=3.14159265359
        const = 0.9695017991974687 * 3.711397914 / 1.66053   !hw  unit: 10**-33 cm*kg

        factor = const*((laserfreq - thzfreq(mod_ind))/100)**4 &
         /(thzfreq(mod_ind)*(1 - exp(-1.0*thzfreq(mod_ind)/6.209376)) ) !hw (0.1/1.66) is used to convert unit of Raman activity.
                                                                        !from A^4/amu to cm^4/kg

        Iraman = 0.0
        A_n = 0.0
        B_n2 = 0.0

        A_n = ( RamTens(mod_ind,1,1) + RamTens(mod_ind,2,2) & 
                      + RamTens(mod_ind,3,3) )/3.0
        B_n2 =(  (RamTens(mod_ind,1,1) - RamTens(mod_ind,2,2))**2 &
                +(RamTens(mod_ind,2,2) - RamTens(mod_ind,3,3))**2 &
                +(RamTens(mod_ind,3,3) - RamTens(mod_ind,1,1))**2 &  
          + 6.0*( RamTens(mod_ind,1,2)**2 &
                + RamTens(mod_ind,2,3)**2 &
                + RamTens(mod_ind,3,1)**2 ) )/2.0
        Iraman = 45.0*A_n**2 + 7.0*B_n2

        depol_ratio = 3.0*B_n2/(45.0*A_n**2 + 4.0*B_n2) 

        DiffRamCross(mod_ind) = factor * Iraman
        write (25,'(1x,2f15.6,3x,f12.6,3x,f12.6)') &
        cmfreq(mod_ind), thzfreq(mod_ind), Iraman, depol_ratio

      enddo
      close (25)

!---------------------------------------------!
!     BROADENING with Lorentzian     
!---------------------------------------------!
      freq_min = 0.0 
      freq_max = cmfreq(nmodes) + 200
      
      Nbin = (freq_max - freq_min )/ dbin
      allocate (dos(Nbin))
      dos (:) = 0.0

!     Sort the frequency in the bins
!     Exclude first the first 6 bands 
      do mod_ind = 7, nmodes
         do ibin = 1, Nbin
            f = freq_min + (dbin*ibin)
            fval= cmfreq(mod_ind)
            famp = DiffRamCross(mod_ind)      
            dosval = famp*lorentz(f,fval,dbin,degauss)
            dos(ibin) = dos(ibin) + dosval
         enddo 
      enddo
!     done

     open(26,file='Spectra_Raman_lorentz.OUT',action='WRITE')
      do ibin =1,Nbin
         f = freq_min + (ibin*dbin) 
         write(26,2000) f, dos(ibin)
      enddo
      close(26)
 
      2000 format(f12.3, 5X,f12.6)
!---------------------------------------------!
      end program main
!---------------------------------------------!
    
        subroutine get_qvec(l, q1, q2, q3)
            implicit none
            character(len=*), intent(in) :: l
            real, intent(out) :: q1, q2, q3
            integer :: idx

            ! Search for "number#"
            idx = index(l, 'q =') + len('q =')

            ! Get the integer after that word
            read(l(idx:idx+12), '(3f12.4)') q1, q2, q3

        end subroutine get_qvec

        subroutine get_mode(l, mnum, tfreq, cfreq)
            implicit none
            character(len=*), intent(in) :: l
            real*8, intent(out) :: tfreq, cfreq
            integer, intent(out) :: mnum
            integer :: idx

            ! Search for "number#"
            idx = index(l, 'freq (') + len('freq (')

            ! Get the integer after that word
            read(l(idx:idx+5), '(I5)') mnum

            idx = index(l, ') =') + len(') =')
            ! Get the integer after that word
            read(l(idx:idx+15), '(f15.6)') tfreq

            idx = index(l, '[THz] =') + len('[THz] =')
            ! Get the integer after that word
            read(l(idx:idx+15), '(f15.6)') cfreq

        end subroutine get_mode

        subroutine get_eigvec(l, e1, e2, e3)
            implicit none
            character(len=*), intent(in) :: l
            complex*8, intent(out) :: e1, e2, e3
            integer :: idx

            ! Search for "number#"
            idx = index(l, '(') + len('(')

            ! Get the integer after that word
            read(l(idx:idx+72), '(3 (f10.6,1x,f10.6,3x))') e1, e2, e3

        end subroutine get_eigvec
       !-----------------------------------------------------------------
       !  Define the broadening function here (Lorentzian Broadening)
       !-----------------------------------------------------------------
        FUNCTION lorentz (e,eval,dbin,degauss)
            IMPLICIT NONE
            REAL :: lorentz
            REAL :: sigma, A, B
            REAL, INTENT(IN) :: e, eval,  dbin, degauss
      
            ! Declare local constant Pi
            REAL, PARAMETER :: Pi = 3.1415927
      
            sigma = degauss
            A = 1.0/ (Pi)
            B = (e - eval)**2 + sigma**2
            lorentz = A * (sigma / B)
!             lorentz = 1.0/ B
        END FUNCTION lorentz
       !-----------------------------------------------------------------

