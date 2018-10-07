!> \file green_harmonics.f90
!!    It contains the \ref green_harmonics module.
!*************************************************************
!  green_harmonics
!
!>    This module calculates the U matrix and the tensor moment
!!    matrices.
!
!*************************************************************
module seedsig_harmonics
   implicit none

   private
   !> array of wigner3j(j1,0,j2,0,j3,0) used in the u4 subroutine. Notice that j2
   !!             must be .ge. j1+j3 as it is used as the sumation parameter in u4. 
   real(KIND=8)           :: wigner3j000(0:3,0:6,0:3)
   !> array of (n!)^(-1/2)
   real(KIND=8), parameter :: invsqrtfac(0:15) = (/1d0,1d0,0.7071067811865475d0,0.4082482904638630d0,0.2041241452319315d0,&
         0.912870929175276d-1,0.37267799624996494d-1,0.14085904245475276d-1,0.4980119205559973d-2,0.1660039735186657d-2,&
         0.524950656957260d-3,0.15827857841616381d-3,0.45691089927761735d-4,0.12672428274337012d-4,0.338684892d-5,&
         0.874480631d-6/)
   !> array of (n/2)!
   integer(KIND=SELECTED_INT_KIND(17)),parameter :: fac2(0:28) = & 
         (/1_8,0_8,1_8,0_8,2_8,0_8,6_8,0_8,24_8,0_8,120_8,0_8,720_8,0_8,5040_8,0_8,40320_8,0_8,362880_8,0_8,3628800_8,0_8, &
         39916800_8,0_8,479001600_8,0_8,6227020800_8,0_8,87178291200_8/)

   public :: harmonic_u4,harmonic_setup,harmonic_gamma,uj2slater

contains


!*************************************************************
!  facapprox
!
!>    calculates the approximate value of n!/(m1!m2!m3!)
!!
!!    n! = exp( p(n) + log(2*pi)/2 - n + (n+1/2)*log(n))
!!    where p(n) is a continued fraction with the parameters
!!    acoeff (see code) as coefficients.
!!
!*************************************************************
   pure function facapprox2(twon,twom1,twom2,twom3)
   implicit none
   real(KIND=8)               :: facapprox2     ! function 
   integer,intent(in)         :: twon           !< 2n
   integer,optional,intent(in):: twom1          !< 2m1
   integer,optional,intent(in):: twom2          !< 2m2
   integer,optional,intent(in):: twom3          !< 2m3

   !> coefficients for the expansion (see Wikipedia)
   real(KIND=8),parameter     :: acoeff(0:8) = (/ 1d0/12d0,1d0/30d0,53d0/210d0,195d0/371d0,22999d0/22737d0,29944523d0/19773142d0,&
                                               109535241009d0/48264275462d0,29404527905795295658d0/9769214287853155785d0,&
                                               455377030420113432210116914702d0/113084128923675014537885725485d0 /)
   real(KIND=8),parameter  :: ln2pi2=.918938533204672741d0     !< log(2pi)/2
   real(KIND=8)            :: pn,pm1,pm2,pm3 ! continued fractions
   real(KIND=8)            :: n,m1,m2,m3     ! real valued m values
   real(KIND=8)            :: dexpcm1,dexpcm2,dexpcm3 ! temp variables
   real(KIND=8)            :: c              ! temp 
   integer                 :: i              ! loop index
   logical                 :: m1flag         ! flag for calculating the expansion
   logical                 :: m2flag         ! flag for calculating the expansion
   logical                 :: m3flag         ! flag for calculating the expansion
   
   ! Deal with the optional arguments.
   m1flag = .false.
   m2flag = .false.
   m3flag = .false.
   dexpcm1 =1d0
   dexpcm2 =1d0
   dexpcm3 =1d0
   if (present(twom1)) then
      if (twom1 .ge. 14) then
         m1flag = .true.
         m1 = twom1*5d-1
         pm1 = 1d0
      else
         dexpcm1 = fac2(twom1)
      endif
   endif
   if (present(twom2)) then
      if (twom2 .ge. 14) then
         m2flag = .true.
         m2 = twom2*5d-1
         pm2 = 1d0
      else
         dexpcm2 = fac2(twom2)
      endif
   endif
   if (present(twom3)) then
      if (twom3 .ge. 14) then
         m3flag = .true.
         m3 = twom3*5d-1
         pm3 = 1d0
      else
         dexpcm3 = fac2(twom3)
      endif
   endif

   n = twon*0.5d0
   pn=1d0
   ! Continued fraction
   do i=8,1,-1
      if (twon .ge. 14) pn =  n+acoeff(i)/pn
      if (m1flag) pm1 =  m1+acoeff(i)/pm1
      if (m2flag) pm2 =  m2+acoeff(i)/pm2
      if (m3flag) pm3 =  m3+acoeff(i)/pm3
   enddo
   
   c = 0d0
   if (twon .ge. 14) c = acoeff(0)/pn+ln2pi2-n+(n+0.5d0)*dlog(n)
   if (m1flag) c = c - ( acoeff(0)/pm1+ln2pi2-m1+(m1+0.5d0)*dlog(m1) )
   if (m2flag) c = c - ( acoeff(0)/pm2+ln2pi2-m2+(m2+0.5d0)*dlog(m2) )
   if (m3flag) c = c - ( acoeff(0)/pm3+ln2pi2-m3+(m3+0.5d0)*dlog(m3) )

   facapprox2=dexp(c)

   if (twon .lt. 14) facapprox2 = facapprox2*fac2(twon)
   if (.not. m1flag) facapprox2 = facapprox2/dexpcm1
   if (.not. m2flag) facapprox2 = facapprox2/dexpcm2
   if (.not. m3flag) facapprox2 = facapprox2/dexpcm3

   end function


!*************************************************************
!  
!  harmonic_setup
!
!>    calculates and stores the Wigner 3j symbols for m=0.
!!
!*************************************************************
   subroutine harmonic_setup()
      implicit none
      integer  :: j1,j2,j3    ! loop counters

      do j3=lbound(wigner3j000,3),ubound(wigner3j000,3)
         do j2=lbound(wigner3j000,2),ubound(wigner3j000,2)
            do j1=lbound(wigner3j000,1),ubound(wigner3j000,1)
               wigner3j000(j1,j2,j3) = wigner3j2(2*j1,0,2*j2,0,2*j3,0)
            enddo
         enddo
      enddo

      return
   end subroutine



!*************************************************************
!  
!  wigner3j2
!
!>    gives wigner3j(twoj1/2,twom1/2,twoj2/2,twom2/2,twoj3/2,twom3/2)
!!       to be able to use integer arguments.
!!
!!       The summation formulas were copied 
!!       from NIST Digital Library of Mathematical Functions (http://dlmf.nist.gov/),
!!       but they are just Clebsch-Gordan coefficients formulas with the extra transformation:
!!       wigner3j(j1,m1,j2,m2,j3,m3) = (-1)^(j1-j2-m3)/sqrt(2*j3+1)*<j1,m1,j2,m2|j3,-m3>
!!
!!       @author PT   
!!
!!       Nb: The fixed size parameter arrays makes this routine limited to j1+j2+j3 .le. 12.
!!       This can easily be extended by increasing the invsqrtfac array. One may also 
!!       modify the routine to use a factorial function/array that takes both a nominator
!!       and a denominator, and distribute the Delta part or the nomenator of the symmetric part
!!       to balance out the other two terms i.e. x^(1/2)*y^(1/2)*1/z = (x/z)*(y/x)^(1/2)
!!
!*************************************************************
   pure function wigner3j2(twoj1,twom1,twoj2,twom2,twoj3,twom3)
      implicit none
      real(KIND=8)         :: wigner3j2   ! wigner 3j coeff
      integer,intent(in)   :: twoj1       !< 2*j1
      integer,intent(in)   :: twom1       !< 2*m1
      integer,intent(in)   :: twoj2       !< 2*j2
      integer,intent(in)   :: twom2       !< 2*m2
      integer,intent(in)   :: twoj3       !< 2*j3
      integer,intent(in)   :: twom3       !< 2*m3
      integer              :: sumstart,sumstop,twos      ! variables for the summation range
      real(KIND=8)         :: sumc,phase                 ! summation temp variables

         ! check to avoid unreasonable values
         wigner3j2 = huge(wigner3j2) ! set a super-unreasonable value to return if the numerical checks are violated
         ! Basic check: twoj1 + twoj2 + twoj3 .le. 28, twoj1,twoj2,twoj3 non-negative, and m in {-j,-j+1,...,j-1,j}.
         if (twoj1 + twoj2 + twoj3 .gt. 28 .or. twoj1 .lt. 0 .or. twoj2 .lt. 0 .or.  twoj3 .lt. 0 .or.&
             mod(twoj1+twom1,2) .ne. 0 .or. mod(twoj2+twom2,2) .ne. 0 .or. mod(twoj3+twom3,2) .ne. 0 .or.&
             abs(twom1) .gt. twoj1 .or. abs(twom2) .gt. twoj2 .or. abs(twom3) .gt. twoj3) return

         wigner3j2 = 0d0
         ! Check that the sum of the Jz quantum number is zero, and that the triangle conditions is satisfied.
         if (twom1 + twom2 + twom3 .ne. 0 .or. abs(twoj1-twoj2) .gt. twoj3 .or. twoj1+twoj2 .lt. twoj3) return

         ! Phase
         if (mod((twoj1-twoj2-twom3)/2,2) .eq. 0) then
            wigner3j2 = 1d0
         else
            wigner3j2 = -1d0
         endif

         ! Delta(j1j2j3)
         wigner3j2 = wigner3j2 * &
          sqrt(dble(fac2(twoj1+twoj2-twoj3)) * &
               dble(fac2(twoj1-twoj2+twoj3)) * &
               dble(fac2(-twoj1+twoj2+twoj3))) * &
                    invsqrtfac((twoj1+twoj2+twoj3)/2+1)

         ! Symmetric   
         wigner3j2 = wigner3j2 * &
          sqrt(dble(fac2(twoj1+twom1))*dble(fac2(twoj1-twom1)) * &
               dble(fac2(twoj2+twom2))*dble(fac2(twoj2-twom2)) * &
               dble(fac2(twoj3+twom3))*dble(fac2(twoj3-twom3)))

         ! Summation
         sumc = 0d0
         sumstart = -min(0,twoj3-twoj2+twom1,twoj3-twoj1-twom2) ! the factorials must get non-negative arguments.
         sumstop  = min(twoj1+twoj2-twoj3,twoj1-twom1,twoj2+twom2)
         do twos = sumstart,sumstop,2
            if (mod(twos/2,2) .eq. 0) then
               phase = 1d0
            else
               phase = -1d0
            endif
            sumc = sumc + phase/(dble(fac2(twos)) * &
                                 dble(fac2(twoj1+twoj2-twoj3-twos)) * & 
                                 dble(fac2(twoj1-twom1-twos)) * &
                                 dble(fac2(twoj2+twom2-twos)) * &
                                 dble(fac2(twoj3-twoj2+twom1+twos)) *& 
                                 dble(fac2(twoj3-twoj1-twom2+twos)))
         enddo
         wigner3j2 = wigner3j2*sumc

      return
   end function

!*************************************************************
!  
!  wigner6j2
!
!>    gives wigner6j(twoj1/2,twom1/2,twoj2/2,twom2/2,twoj3/2,twom3/2)
!!       to be able to use integer arguments.
!!
!!       The summation formulas were copied 
!!       from NIST Digital Library of Mathematical Functions (http://dlmf.nist.gov/),
!!       Should be rewritten using eq. 34.4.2 when speed is of importance (now 
!!       eq. 34.4.1 is used).
!!
!!       @author OG
!!
!!       Nb: The fixed size parameter arrays makes this routine limited to j1+j2+j3 .le. 12.
!!       This can easily be extended by increasing the invsqrtfac array. One may also 
!!       modify the routine to use a factorial function/array that takes both a nominator
!!       and a denominator, and distribute the Delta part or the nomenator of the symmetric part
!!       to balance out the other two terms i.e. x^(1/2)*y^(1/2)*1/z = (x/z)*(y/x)^(1/2)
!!
!*************************************************************
   pure function wigner6j2(twoj1,twol1,twoj2,twol2,twoj3,twol3)
      implicit none
      real(KIND=8)         :: wigner6j2   ! wigner 3j coeff
      integer,intent(in)   :: twoj1       !< 2*j1
      integer,intent(in)   :: twol1       !< 2*l1
      integer,intent(in)   :: twoj2       !< 2*j2
      integer,intent(in)   :: twol2       !< 2*l2
      integer,intent(in)   :: twoj3       !< 2*j3
      integer,intent(in)   :: twol3       !< 2*l3
      integer              :: twom1,twom2,twom3          ! 2*m1, 2*m2, 2*,3
      integer              :: twom1p,twom2p,twom3p       ! 2*m1',2*m2',2*m3'
      integer              :: sumstart,sumstop,twos      ! variables for the summation range
      real(KIND=8)         :: sumc,phase                 ! summation temp variables

      wigner6j2 = 0d0

      do twom1 = -twoj1,twoj1,2
         do twom2 = -twoj2,twoj2,2
            do twom3 = -twoj3,twoj3,2
               do twom1p = -twol1,twol1,2
                  do twom2p = -twol2,twol2,2
                     do twom3p = -twol3,twol3,2
                        wigner6j2 = wigner6j2 + (-1)**((twol1+twom1p+twol2+twom2p+twol3+twom3p)/2) &
                                  * wigner3j2(twoj1,twom1,twoj2,twom2,twoj3,twom3) & 
                                  * wigner3j2(twoj1,twom1,twol2,twom2p,twol3,-twom3p) & 
                                  * wigner3j2(twol1,-twom1p,twoj2,twom2,twol3,twom3p) & 
                                  * wigner3j2(twol1,twom1p,twol2,-twom2p,twoj3,twom3) 
                     end do
                  end do
               end do
            end do
         end do
      end do

      return
   end function



!*************************************************************
!
!  clebschgordan2
!
!>    calculates Clebsh-Gordan coefficients of half the arguments,
!!    i.e. <twoj1/2,twom1/2,twoj2/2,twom2/2 | twoj3/2,twom3/2>
!
!*************************************************************
   pure function clebschgordan2(twoj1,twom1,twoj2,twom2,twoj3,twom3)
      implicit none
      real(KIND=8)         :: clebschgordan2          ! output Clebsh-Gordan
      integer,intent(in)   :: twoj1       !< 2*j1
      integer,intent(in)   :: twom1       !< 2*m1
      integer,intent(in)   :: twoj2       !< 2*j2
      integer,intent(in)   :: twom2       !< 2*m2
      integer,intent(in)   :: twoj3       !< 2*j3
      integer,intent(in)   :: twom3       !< 2*m3
   
      if (mod(twoj1-twoj2+twom3,2) .eq. 0) then ! phase = -1
         clebschgordan2 = sqrt(dble(twoj3+1))*wigner3j2(twoj1,twom1,twoj2,twom2,twoj3,-twom3)
      else ! phase = 1
         clebschgordan2 = -sqrt(dble(twoj3+1))*wigner3j2(twoj1,twom1,twoj2,twom2,twoj3,-twom3)
      endif
      return

   end function

!*************************************************************
!
!  harmonic_u4
!
!>    calculates the 4-index U-matrix for a single set of orbitals (j1 = j2 = j3 = j4 = l).
!!
!!    Definition:
!!    u(j1,m1,j2,m2,j3,m3,j4,m4) = sum_{k,q} (-1)**(m1+m2+q)*sqrt((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))*
!!           wigner3j(j1,-m1,k,-q,j3,m3)*wigner3j(j1,0,k,0,j3,0)*wigner3j(j2,0,k,0,j4,0)*wigner3j(j2,-m2,k,q,j4,m4)*
!!           slater(k,j1,j2,j3,j4)*delta(spin1,spin3)*delta(spin2,spin4)
!!
!!    This expression is obtained from making a multipole expansion of 1/|r-r'|, and identifying the angular
!!    integrals with the slater integrals. The slater integrals is finally written in terms of the
!!    Wigner 3j symbols.
!!    See for example: http://www.cond-mat.de/events/correl16/manuscripts/eder.pdf
!!    Or: http://www.cond-mat.de/events/correl11/manuscripts/pavarini.pdf
!
!*************************************************************
   pure subroutine harmonic_u4(slater,u)
   implicit none
   real(KIND=8),intent(in)       :: slater(:)         !< Slater parameters
   complex(KIND=8),intent(out)   :: u(:,:,:,:)        !< 4-index U-matrix
   integer                       :: l                 !< l-quantum number
   integer                       :: nspin             !< number of spins
   integer                       :: m1,m2,m3,m4,k,spin1,spin2     ! loop counters
   integer                       :: msize,twol,twok,twom1,twom2,twom3,twom4,twoq,im1,im2,im3,im4,spinoffset,loffset
   real(KIND=8)                  :: c,a

      u = (0d0,0d0)
      if (size(u,1) .ne. size(u,2) .or. size(u,2) .ne. size(u,3) .or. size(u,3) .ne. size(u,4)) return

      msize = size(u,1)
      if (mod(msize,2) .eq. 0) then
         l = (size(u,1)/2-1)/2
         nspin = 2
      else
         l = (size(u,1)-1)/2
         nspin = 1
      endif
      spinoffset = 2*l+1
      loffset = l+1
      twol = 2*l

      ! Definition:
      ! u(j1,m1,j2,m2,j3,m3,j4,m4) = sum_kq (-1)**(m1+m2+q)*sqrt((2*j1+1)*(2*j2+1)*(2*j3+1)*(2*j4+1))*
      !           wigner3j(j1,-m1,k,-q,j3,m3)*wigner3j(j1,0,k,0,j3,0)*wigner3j(j2,0,k,0,j4,0)*wigner3j(j2,-m2,k,q,j4,m4)*
      !           slater(k,j1,j2,j3,j4)*delta(spin1,spin3)*delta(spin2,spin4)

      do k=0,2*l,2 ! as k = j2 =< j1+j3 = l + l
         ! m independent part: slater(k,l,l,l,l)*wigner3j(l,0,k,0,l,0)*wigner3j(l,0,k,0,l,0)*sqrt((2*l+1)*(2*l+1)*(2*l+1)*(2*l+1))
         twok = 2*k
         c = slater(k/2+1)*(wigner3j000(l,k,l)**2)*(2*l+1)**2
         do m4 = -l,l
            twom4 = 2*m4
            do m3 = -l,l
               twom3 = 2*m3
               do m2 = max(-l,-l+m3+m4,-k+m4),min(l,l+m3+m4,k+m4) ! -l < m1 < l, -k < q < k
                  twom2 = 2*m2
                  m1 = m3+m4-m2 ! m1+m2 = m3+m4
                  twom1 = 2*m1
                  twoq = (twom2 - twom4) ! wigner3j only non-zero if q-m2+m4 = m3-m1-q = 0
                  if (mod(m1+2*m2-m4,2) .eq. 0) then
                     a =  c*wigner3j2(twol,-twom1,twok,-twoq,twol,twom3)*wigner3j2(twol,-twom2,twok,twoq,twol,twom4)
                  else
                     a = -c*wigner3j2(twol,-twom1,twok,-twoq,twol,twom3)*wigner3j2(twol,-twom2,twok,twoq,twol,twom4)
                  endif
                  do spin1=1,nspin
                     do spin2=1,nspin 
                        im1 = spinoffset*(spin1-1)+loffset+m1
                        im2 = spinoffset*(spin2-1)+loffset+m2
                        im3 = spinoffset*(spin1-1)+loffset+m3
                        im4 = spinoffset*(spin2-1)+loffset+m4
                        u(im1,im2,im3,im4) = u(im1,im2,im3,im4) + a
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
   end subroutine harmonic_u4




!*************************************************************
!
!  facfac2
!
!>    returns the double factorial (n!!)
!!
!*************************************************************
   pure function facfac2(twon)
   implicit none
   real(KIND=8)       :: facfac2       ! n!!
   integer,intent(in) :: twon          !< 2* n

      facfac2 = 0
      if (mod(twon,2) .ne. 0) return

      if (mod(twon/2,2) .eq. 0) then
         ! Even
         facfac2 = 2**(twon/4)*fac2(twon/2) ! (2n)!! = n!*2^n
      else
         ! Odd
         facfac2 = facapprox2(twon,twon/2-1)/2**(twon/2-1)  ! (2n+1)!! = (2n+1)!/(n!*2^n)
      endif         
      return
   end function facfac2



!*************************************************************
!
!  harmonic_nbar
!
!>    returns the tensor moment normalization constant coeff.
!!
!!    See definition in eqn. 27 in Phys. Rev. B 80, 035121 (2009)
!*************************************************************
   pure function harmonic_nbar(twok,twop,twor)
   implicit none
   complex(KIND=8)    :: harmonic_nbar             ! constant
   integer,intent(in) :: twok       !< 2* index k
   integer,intent(in) :: twop       !< 2* index p
   integer,intent(in) :: twor       !< 2* index r
   integer            :: twog,twogk,twogp,twogr    ! temporary indices
   real(KIND=8)       :: gfacfac                   ! coefficient (see below in the equations)

      twog = twok+twop+twor
      twogk = twop+twor-twok
      twogp = twok+twor-twop
      twogr = twok+twop-twor
      if ( twok+twop .lt. twor .or. abs(twok-twop) .gt. twor .or. mod(twok,2) .ne. 0 .or. mod(twor+twop,2) .ne. 0) then
         harmonic_nbar = 0
      else 
         !Make gfacfac = g!! / ((g-2k)!! (g-2p)!! (g-2r)!!)
         if (twop .eq. 0) then
            ! gfacfac = 1/((g-2k)!!(g-2r)!!) = 1/((g-2k)!!)^2 
            gfacfac = 1d0/facfac2(twogk)**2
         else
            ! gfacfac = g/((g-2k)!!(g-2r)!!)
            gfacfac = twog/(2*facfac2(twogk)*facfac2(twogr))
         endif
         
         gfacfac = gfacfac*dsqrt(dble(fac2(twogk)*fac2(twogp)*fac2(twogr)))*invsqrtfac(twog/2+1)
         ! Multiply with i^(k+p+r)
         select case(mod(twog/2,4))
         case (1)
            harmonic_nbar = cmplx(0d0,gfacfac,KIND=8)
         case (2)
            harmonic_nbar = -gfacfac
         case (3)
            harmonic_nbar = -cmplx(0d0,gfacfac,KIND=8)
         case default
            harmonic_nbar = gfacfac
         end select
      endif

      return
   end function harmonic_nbar


!*************************************************************
!
!  harmonic_n
!
!>    returns the tensor moment normalization factor.
!!
!!    See definition in eqn. 21 in Phys. Rev. B 80, 035121 (2009)
!*************************************************************
   pure function harmonic_n(twol,twok)
   real(KIND=8)         :: harmonic_n     ! normalization factor
   integer,intent(in)   :: twok           !< 2* k index
   integer,intent(in)   :: twol           !< 2* l index

      harmonic_n = fac2(2*twol)/dsqrt(dble(fac2(2*twol-twok)*fac2(2*twol+twok+2)))
      return
   end function harmonic_n

!*************************************************************
!
!  harmonic_gamma
!
!>    generates the gamma matrices.
!!
!!    See equation 5.37 in Multipoles in Correlated Electorn
!!    Materials, PhD Thesis, F. Cricchio
!!
!!    @todo Finish commenting the input (at least!!!). Maybe
!!       Oscar is the best person for looking at this.
!*************************************************************
   subroutine harmonic_gamma(M,k,p,r,rho,w2)
      implicit none
      integer,intent(in)            :: M  !< size of the matrix
      integer,intent(in)            :: k  !< k - even numbers density, odd current
      integer,intent(in)            :: p  !< p - charge or spin
      integer,intent(in)            :: r  !< r - coupling of charge and spin, odd breaks time reversal symmetry
      integer,allocatable           :: qn(:,:)    !< conserved quantum numbers
      complex(KIND=8),intent(in)    :: rho(:,:)   !< density matrix
      real(KIND=8),intent(out)      :: w2  !< w - multipoles
      complex(KIND=8),allocatable   :: gamma_matrix(:,:,:) !< gamma matrices (M,M,t)
      complex(KIND=8),allocatable   :: val(:)
      real(KIND=8)                  :: c     !< temp variable
      integer                       :: t  !< t - coupling index
      integer                       :: i,a,b,twomla,twomlb,twomsa,twomsb,twol,twox,twoy,phase,er
     

      allocate(gamma_matrix(M,M,-r:r),val(-r:r),qn(2,M),stat=er)
      gamma_matrix = (0d0,0d0)
      val = (0d0,0d0)

      ! quantum numbers
      if (M .eq. 10) then
         qn(1,:) = (/ -2,-1, 0, 1, 2,-2,-1, 0, 1, 2 /)
         qn(2,:) = (/ -1,-1,-1,-1,-1, 1, 1, 1, 1, 1 /)
      elseif (M .eq. 14) then
         qn(1,:) = (/ -3,-2,-1, 0, 1, 2, 3,-3,-2,-1, 0, 1, 2, 3 /)
         qn(2,:) = (/ -1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1  /)
      endif

      ! Get the l quantum number from the size of gamma_matrix
      twol = maxval(qn(1,:))*2

      ! Try a quick exit
      if (k+p .lt. r .or. abs(k-p) .gt. r .or. p .gt. 1 .or. k .gt. twol) return

      ! Calculate the matrix elements of gamma
      do t=-r,r
         do a=1,M
            do b=1,M
               twomla = qn(1,a)*2 ! 2*ml
               twomlb = qn(1,b)*2 ! 2*ml
               twomsa = qn(2,a)   ! 2*ms
               twomsb = qn(2,b)   ! 2*ms
               c = 0
               ! Sum over x and y (no phase factor since -x-y = t)
               ! sum_{x,y} wigner3j(l,-ml1,k,x,l,ml2)*wigner3j(s,-ms1,p,y,s,ms2)*wigner3j(k,x,r,t,p,y)
               do twox=-2*k,2*k,2
                  twoy = 2*t-twox ! t-x-y = 0
                  if (abs(twoy) .le. 2*p .and. twomsb+twoy .eq. twomsa .and. -twomla+twox .eq. -twomlb) then
                     c = c + wigner3j2(twol,-twomla,2*k,twox,twol,twomlb)*&
                          wigner3j2(1,-twomsa,2*p,twoy,1,twomsb)*&
                          wigner3j2(2*k,-twox,2*r,2*t,2*p,-twoy)
                  endif
               enddo
               ! phase factor: (-1)^(-ml1+l-ms1+s+k+p+t)
               phase = (-twomla+twol-twomsa+1+2*k+2*p+2*t)/2
               if (mod(phase,2) .eq. 0) then
                  gamma_matrix(a,b,t) = c
               else
                  gamma_matrix(a,b,t) = -c
               endif
            enddo
         enddo

         ! Multiply with the normalization
         gamma_matrix(:,:,t) = (gamma_matrix(:,:,t)/harmonic_nbar(2*k,2*p,2*r))/harmonic_n(twol,2*k)/harmonic_n(1,2*p)

      enddo

!      write(*,*) size(val)
      do t=-r,r
         val(t)=sum((/(dot_product(gamma_matrix(:,i,t),rho(:,i)),i=1,M)/))
         write(*,'(1x,99f14.8)') val(t)
      enddo

      w2 = 0
      do t=-r,r
         if (mod(t,2) .eq. 0) then
            w2 = w2 + val(-t)*val(t)
         else
            w2 = w2 - val(-t)*val(t)
         endif
      enddo
      write(*,'(1x,A,99f14.8)') "w2",w2

      deallocate(gamma_matrix,val,qn)
      return

   end subroutine harmonic_gamma

!*******************************************************************
!
!>    Reads U and J, gives back an array of Slater parameters
!!
!*******************************************************************

   subroutine uj2slater(M,uval,jval,slater)
      integer,intent(in)                    :: M           !< d or f orbital
      real(KIND=8),intent(in)               :: uval, jval  !< size of the matrix
      real(KIND=8),intent(out)              :: slater(:)   !< slater parameters

      if(M .eq. 10) then

         ! d-orbitals => F0, F2, F4
         ! F0 = U
         ! F4/F2 = 0.625 Empirical
         ! J= (F2 + F4)/14
         ! -> F2 = 14*J/(1+0.625)
         slater(1) = uval
         slater(2) = jval*14.d0/(1d0+0.625d0)
         slater(3) = 0.625d0*slater(2)

      elseif(M .eq. 14) then

         ! f-orbitals => F0, F2, F4, F6
         ! F0 = U
         ! F4/F2 = 0.668 Empirical
         ! F6/F2 = 0.494 Empirical
         ! J = (286F2 + 195F4 + 250F6)/6435
         ! F2 = 6435*J/(286+195*0.668+250*0.494)
         slater(1) = uval
         slater(2) = jval*6435/(286+195*0.668d0+250*0.494d0)
         slater(3) = 0.668d0*slater(2)
         slater(4) = 0.494d0*slater(2)

      endif

   end subroutine

end module seedsig_harmonics

!
