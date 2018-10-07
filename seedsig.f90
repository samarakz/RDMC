    program seedsig
    !------------------------------------------------------------------------------
    ! This program calculates the static self-energy via a given density.         |
    !------------------------------------------------------------------------------
    use seedsig_harmonics
    use seedsig_selfenergy
    implicit none

    character(len=24),parameter   :: inpname = 'seedsig.inp', outname='seedsig.log', initso='so.log'
    character(len=24)             :: label               ! cluster label
    complex(kind=8),allocatable   :: u(:,:,:,:)          ! U matrix will be read from input
    complex(kind=8),allocatable   :: rho(:,:)            ! the density matrix
    complex(kind=8),allocatable   :: sig(:,:)            ! the static part of the self energy for one atom
    complex(kind=8),allocatable   :: dc(:,:)             ! the double counting DC
    complex(kind=8),allocatable   :: sig_dc(:,:)         ! the self energy minus DC
    complex(kind=8),allocatable   :: sig_dc_afm(:,:)     ! the self energy minus DC
    complex(kind=8),allocatable   :: D(:,:)              ! occupation matrix, mostly diagonal
    complex(kind=8),allocatable   :: u_rand(:,:)         ! the random hermitian matrix
    complex(KIND=8),allocatable   :: Sx(:,:)             ! Sx
    complex(KIND=8),allocatable   :: Sy(:,:)             ! Sy
    complex(KIND=8),allocatable   :: Sz(:,:)             ! Sz
    complex(KIND=8),allocatable   :: Lp(:,:)             ! L+
    complex(KIND=8),allocatable   :: Lx(:,:)             ! Lx
    complex(KIND=8),allocatable   :: Ly(:,:)             ! Ly
    complex(KIND=8),allocatable   :: Lz(:,:)             ! Lz
    complex(KIND=8),allocatable   :: SO(:,:)             ! Spin-orbit matrix
    complex(KIND=8),allocatable   :: Jz(:,:)             ! Jz
    complex(KIND=8),allocatable   :: iden(:,:)           ! Identity matrix
    complex(KIND=8),allocatable   :: Trans(:,:)          ! Transformation matrix
    complex(kind=8),allocatable   :: work(:)             ! temporary work space
    complex(kind=8),allocatable   :: rot(:,:)            ! rotation matrix
    real(kind=8),allocatable      :: rwork(:), W(:)      ! temporary matrices
    real(kind=8)                  :: uval, jval          ! parametrs U and J
    real(KIND=8)                  :: slater(4)           ! slater parameters
    real(KIND=8)                  :: mag(3)              ! magnetization Tr(rho.S)
    real(KIND=8)                  :: lmag(3)             ! magnetization Tr(rho.L)
    real(kind=8)                  :: x,y, tmp_1,tmp_2    ! temporary values
    real(kind=8)                  :: Uavg, Javg          ! U and J parameters
    real(kind=8)                  :: nel                 ! Occupation
    real(kind=8)                  :: alpha               ! mixing
    real(kind=8)                  :: w2                  ! multipoles to power 2
    integer                       :: k                   ! even density, odd current 
    integer                       :: p                   ! 0 charge, 1 spin 
    integer                       :: r                   ! coupling of charge and spin, odd breaks time reversal symmetry 
    integer                       :: id(5,1)             ! the id of the cluster
    integer                       :: id_fm2(5,1)         ! the id of the afm cluster
    integer                       :: id_afm1(5,1)         ! the id of the afm cluster
    integer                       :: id_afm2(5,1)         ! the id of the afm cluster
    integer                       :: seed, er            ! random number and error message
    integer                       :: n, a, i, j, l       ! loop counters
    integer                       :: u_size              ! number of nonzero elements in u matrix
    logical                       :: socflag             ! spin-orbit flag
    logical                       :: magflag             ! spin-polariz DC flag
    logical                       :: AFMflag             ! spin-orbit flag
    logical                       :: diag_occ            ! give only diagonal occupation matrix
    integer                       :: ntot                ! number of clusters
    integer                       :: M, lwork, info
    real, parameter               :: ev2ry = 1d0/13.605693009d0 ! to convert eV to Ry
    real, parameter               :: pi = 3.14159265

    call system_clock(seed)      ! to avoid generating the same random numbers.
    call srand(seed)
 
    open(unit=3, file=outname,status='REPLACE',action='WRITE',iostat=er)
    write(3,'(1X, "The input file is read from file:  ", A24)') inpname

    open(unit=1,file=inpname,status='OLD',action='READ',iostat=er)

    if(er>0) then
       write(*,"(1X,' ERROR: File ',A,' does not exist!', I5)") inpname, er
       return
    endif

    read(1,*,iostat=er) ntot, socflag, magflag, AFMflag, diag_occ, alpha
    read(1,*,iostat=er) uval, jval
    read(1,*,iostat=er) k,p,r

    ! U and J are scaled from eV to Ry
    uval = ev2ry*uval
    jval = ev2ry*jval

    ! Initialize the harmonic module
    call harmonic_setup()

    main: do n=1, ntot
            ! Read the input file to get the id and the number of orbitals
             id = 0
             id_fm2 = 0
             id_afm1 = 0
             id_afm2 = 0
             read(1,*,iostat=er) id(:,1)
             write(3,"(1X,'The id of the cluster:',5i4)") id(:,1)
             if (AFMflag) then
                read(1,*,iostat=er) id_fm2(:,1)
                write(3,"(1X,'The id of the cluster:',5i4)") id_fm2(:,1)
                read(1,*,iostat=er) id_afm1(:,1)
                write(3,"(1X,'The id of the cluster:',5i4)") id_afm1(:,1)
                read(1,*,iostat=er) id_afm2(:,1)
                write(3,"(1X,'The id of the cluster:',5i4)") id_afm2(:,1)
             endif
             ! Make the label 
             call lda_id2label(id(:,1),label)
             write(3,*) "label  ",label

             ! Get the size from the l quantum number: id = (t, l, e, site, basis)
             M = 2*(2*id(2,1) + 1)

            ! Allocate and initialize the arrays
             allocate(u(M,M,M,M),u_rand(M,M),rho(M,M),D(M,M),iden(M,M),Sx(M,M),Sy(M,M),Sz(M,M),SO(M,M), &
                      Lp(M,M),Lx(M,M),Ly(M,M),Lz(M,M),Jz(M,M),Trans(M,M),rot(M,M),sig(M,M),dc(M,M),     &
                      sig_dc(M,M),sig_dc_afm(M,M),W(M),stat=er)

             ! u      = (0d0,0d0)
             D      = (0d0,0d0)
             sig    = (0d0,0d0)
             dc     = (0d0,0d0)
             sig_dc = (0d0,0d0)
             sig_dc_afm = (0d0,0d0)
             u_rand = (0d0,0d0)
             Sx     = (0d0,0d0)
             Sy     = (0d0,0d0)
             Sz     = (0d0,0d0)
             Lp     = (0d0,0d0)
             Lx     = (0d0,0d0)
             Ly     = (0d0,0d0)
             Lz     = (0d0,0d0)
             Jz     = (0d0,0d0)
             SO     = (0d0,0d0)
             Trans  = (0d0,0d0)
             rot    = (0d0,0d0)
             iden   = (0d0,0d0)

             ! constructing Sx, Sy and Sz operators.
             forall(i=1:M/2) Sx(i,i+M/2) = dcmplx(0.5d0,0d0)
             forall(i=1:M/2) Sx(i+M/2,i) = dcmplx(0.5d0,0d0)
             forall(i=1:M/2) Sy(i+M/2,i) = dcmplx(0d0,-0.5d0)
             forall(i=1:M/2) Sy(i,i+M/2) = dcmplx(0d0, 0.5d0)
             forall(i=1:M/2) Sz(i,i)     = dcmplx(-0.5d0,0d0)
             forall(i=1:M/2) Sz(i+M/2,i+M/2) = dcmplx(0.5d0,0d0)
             forall(i=1:M) iden(i,i) = dcmplx(1d0,0d0)

             ! constructing L+
             do i=2,M/2
                Lp(i,i-1) = sqrt(dble(id(2,1)*(id(2,1)+1) - (-id(2,1)+i-2)*(-id(2,1)+i-2+1)))
             enddo

             forall(i=2:M/2) Lp(i+M/2,i+M/2-1) = Lp(i,i-1)

             Lx = (Lp+transpose(dconjg(Lp)))/2d0
             Ly = (Lp-transpose(dconjg(Lp)))/dcmplx(0d0,2d0)

             forall(i=1:M/2) Lz(i,i) = dcmplx(-id(2,1)+i-1,0d0)
             forall(i=1:M/2) Lz(i+M/2,i+M/2) = Lz(i,i)

             Jz = Sz + Lz

             do i=1,M
                SO = matmul(Lx,Sx) + matmul(Ly,Sy) + matmul(Lz,Sz)
             enddo

       !      write(3,*) ""
       !      write(3,*) "SO: "
       !      write(3,*) "Real:"
       !      do i=1, M
       !         write(3,'(1x,99f12.8)') dble(SO(i,:))
       !      enddo

       !      write(3,*) ""
       !      write(3,*) "Imaginary:"
       !      do i=1, M
       !         write(3,'(1x,99f12.8)') dimag(SO(i,:))
       !      enddo


            ! Read the suggested orbital occupation 
             read(1,*) (W(i),i=1,size(D,1))
             if (er /= 0) write(*,*) "ERROR: reading occupation."
             do i=1, size(D,1)
                D(i,i) = dcmplx(W(i),0d0)
             enddo
             deallocate(W)

             write(3,*) " Diagonal part of D:"
             write(3,'(1x,99f12.7)') (dble(D(i,i)),i=1,M)
             write(3,*) " "

             ! U and J are mapped to slater parameters
             call uj2slater(M,uval,jval,slater)
             write(3,'(1x,99f12.7)') slater(:)
             write(*,*) "Slater parameters:", slater(:)

             ! Construct the U matrix file
             call harmonic_u4(slater,u)
             

            ! generating a M x M random Hermitian matrix.
             do i=1, M
                do j=i, M
                   x = rand()
                   y = rand()
                   u_rand(i,j) = dcmplx(x,y)
                   if (i .eq. j) u_rand(i,j) = dcmplx(x,0d0)
                   u_rand(j,i) = dconjg(u_rand(i,j))
                enddo
             enddo

             ! in case of scalar-relativistic.
             do i=1, M/2
                do j=M/2+1,M
                   u_rand(i,j) = dcmplx(0d0,0d0)
                   u_rand(j,i) = dcmplx(0d0,0d0)
                enddo
             enddo
  

             write(3,*) "Real part of the random matrix:"
             do i=1, M
                write(3,'(1x,99f12.7)') dble(u_rand(i,:))
             enddo

             write(3,*) ""
             write(3,*) "Imaginary part of the random matrix:"
             do i=1, M
                write(3,'(1x,99f12.7)') dimag(u_rand(i,:))
             enddo

             Trans = u_rand
     !       Trans = (1-alpha)*u_rand + alpha*D

             ! to get the eigenvectors of u_rand matrix.
             lwork = 1024         ! 64*M-32
             allocate(work(lwork), rwork(3*M-2), W(M),stat=er)
             call zheev('V','L',M/2,Trans(1:M/2,1:M/2),M/2,W(1:M/2),work,lwork,rwork,info)
             call zheev('V','L',M/2,Trans(M/2+1:M,1+M/2:M),M/2,W(1+M/2:M),work,lwork,rwork,info)


             if (socflag) then
                ! construct the rotation matrix
                x = pi*rand()
                y = 2*pi*rand()

                do i=1,M/2
                   rot(i,i) = cos(pi*x)
                   rot(i+M/2,i+M/2) = cos(pi*x)
                   rot(i,i+M/2) = -1*exp(dcmplx(0d0,-1d0)*pi*y)*sin(pi*x)
                   rot(i+M/2,i) =    exp(dcmplx(0d0, 1d0)*pi*y)*sin(pi*x)
                enddo

                Trans(1:M,1:M) = matmul(rot,Trans)
!               Trans(1:M,1:M) = matmul(Trans,transpose(dconjg(Trans)))  ! should be unitary
             endif

             write(3,*) ""
             write(3,*) "Transformation matrix after diagonalization and rotation: "
             write(3,*) "Real:"
             do i=1, M
                write(3,'(1x,99f12.7)') dble(Trans(i,:))
             enddo

             write(3,*) ""
             write(3,*) "Imaginary:"
             do i=1, M
                write(3,'(1x,99f12.7)') dimag(Trans(i,:))
             enddo


             ! if diagonal occupation is favored
             if (diag_occ) then
                Trans = 0.0d0
                do i=1,M
                   Trans(i,i) = dcmplx(1.0d0,0d0)
                enddo
             endif


             ! constructing the initial density
             rho(1:M,1:M) = matmul(Trans,matmul(D,transpose(dconjg(Trans))))
             rho = (1-alpha)*rho + alpha*D

        !     rho(7,9) = -0.1d0
        !     rho(9,7) = -0.1d0

             write(3,*) ""
             write(3,*) "Real part of the density: "
             do i=1, M
                write(3,'(1x,99f12.7)') dble(rho(i,:))
             enddo
             write(3,*) "Imaginary part of the density: "
             do i=1, M
                write(3,'(1x,99f12.7)') dimag(rho(i,:))
             enddo

             ! generates the gamma matrices 
             call harmonic_gamma(M,k,p,r,rho,w2)
             write(3,*) ""
             write(3,'(1x," ID:  ",A10,"   w(k=",i1,",p=",i1,",r=",i1,")^2:",15x,es14.6)') &
                   trim(label),k,p,r,w2
             
             write(*,'(1x," ID:  ",A10,"   w(k=",i1,",p=",i1,",r=",i1,")^2:",15x,es14.6)') &
                   trim(label),k,p,r,w2

             ! constructing the static part of the self energy
             do i=1, M
                do j=1, M
                   sig(1:M,1:M) = sig(1:M,1:M) + (u(:,i,:,j)-u(:,i,j,:))*rho(j,i)
                enddo
             enddo

             write(3,*) ""
             write(3,*) "Real part of the self energy:"
             do i=1, M
                write(3,'(1x,99f12.7)') dble(sig(i,:))
             enddo

             write(3,*) ""
             write(3,*) "Imaginary part of the self energy:"
             do i=1, M
                write(3,'(1x,99f12.7)') dimag(sig(i,:))
             enddo

             ! Get the average U and J values
             Uavg = 0
             Javg = 0
             do i=1,M
                do j=1,M
                   ! <U> and <J>
                   Uavg = Uavg + u(i,j,i,j)
                   Javg = Javg + u(i,j,j,i)
                enddo
             enddo
             Javg = (Javg-Uavg/M)/(M*(M/2-1))
             Uavg = Uavg/M**2
             write(3,*) "Average U and J:",Uavg,Javg

             ! Get the spherical averages
             mag = 0d0
             lmag = 0d0
             nel = 0d0
             do i=1,M
                nel     = nel     + dot_product(iden(1:M,i),rho(1:M,i))     ! Tr(rho)
                mag(1)  = mag(1)  + dot_product(rho(1:M,i),Sx(1:M,i))       ! Tr(rho*.Sx)
                mag(2)  = mag(2)  + dot_product(rho(1:M,i),Sy(1:M,i))       ! Tr(rho*.Sy)
                mag(3)  = mag(3)  + dot_product(rho(1:M,i),Sz(1:M,i))       ! Tr(rho*.Sz)
                lmag(1) = lmag(1) + dot_product(rho(1:M,i),Lx(1:M,i))       ! Tr(rho*.Lx)
                lmag(2) = lmag(2) + dot_product(rho(1:M,i),Ly(1:M,i))       ! Tr(rho*.Ly)
                lmag(3) = lmag(3) + dot_product(rho(1:M,i),Lz(1:M,i))       ! Tr(rho*.Lz)
             enddo


             open(unit=4,file=initso,position='append',status='old',action='readWRITE',iostat=er)

             mag = 2*mag
             write(3,*) ""
             write(3,*) "Number of electrons: ", nel
             write(3,*) "Spin moments:"
             write(3,'(1x,4f12.7)') (mag(i), i=1,3), sqrt(dot_product(mag,mag))
             write(*,'(1x,4f12.7)') (mag(i), i=1,3), sqrt(dot_product(mag,mag))

             write(3,*) "Orbital moments:"
             write(3,'(1x,4f12.8)') (lmag(i), i=1,3), sqrt(dot_product(lmag,lmag))

             tmp_1 = dot_product(mag,lmag)/sqrt(dot_product(mag,mag))/sqrt(dot_product(lmag,lmag))
             tmp_2 = (180d0/3.14d0)*acos(sign(min(abs(tmp_1),1d0),tmp_1))

!            write S, L, S.L/|S||L|,angle and the w2 to the so.log file
!             write(4,*) "  Spin   Orbital   S.L/|S||L|   Angle  W2"
             write(4,'(1x,4f10.3,f14.5)') sqrt(dot_product(mag,mag)), sqrt(dot_product(lmag,lmag)), tmp_1, tmp_2, w2

             close(4)

             ! lda instead of lsda for average double counting
             if(.not. magflag) then
                mag = 0d0
             endif

             ! DC potential
             dc(1:M,1:M) = 0.5d0*(Uavg*(2*nel-1) - Javg*(nel-1))*iden(1:M,1:M) + &
                            (-Javg)*mag(1)*Sx(1:M,1:M) +&
                            (-Javg)*mag(2)*Sy(1:M,1:M) +&
                            (-Javg)*mag(3)*Sz(1:M,1:M)

             write(3,*) ""
             write(3,*) "Real part of DC:"
             do i=1, M
                write(3,'(1x,99f12.7)') dble(dc(i,:))
             enddo

             write(3,*) ""
             write(3,*) "Imaginary part of DC:"
             do i=1, M
                write(3,'(1x,99f12.7)') dimag(dc(i,:))
             enddo

            
             ! Sigma - DC
             sig_dc(1:M,1:M) = sig(1:M,1:M) - dc(1:M,1:M)

             write(3,*) ""
             write(3,*) "Real part of the self energy after DC:"
             do i=1, M
                write(3,'(1x,99f12.7)') dble(sig_dc(i,:))
             enddo

             write(3,*) ""
             write(3,*) "Imaginary part of the self energy after DC:"
             do i=1, M
                write(3,'(1x,99f12.7)') dimag(sig_dc(i,:))
             enddo

             call selfenergy_write(sig_dc,id)
             call selfenergy_write(sig_dc,id_fm2)

             ! Create the self energy for the other cluster in a AFM system.
             if(AFMflag) then
                 ! Use dc as a temporary array
                 do i=1, M
                    do j=1,M
                       dc(i,j)=(-1)**(i-j)*sig_dc(M+1-j,M+1-i)
                    enddo
                 enddo
                 sig_dc_afm = dc

                 write(3,*) ""
                 write(3,*) "Real part of the self energy of the AFM atom:"
                 do i=1, M
                    write(3,'(1x,99f12.7)') dble(sig_dc_afm(i,:))
                 enddo

                 write(3,*) ""
                 write(3,*) "Imaginary part of the self energy of AFM atom:"
                 do i=1, M
                    write(3,'(1x,99f12.7)') dimag(sig_dc_afm(i,:))
                 enddo

                 call selfenergy_write(sig_dc_afm,id_afm1)
                 call selfenergy_write(sig_dc_afm,id_afm2)

             endif

             deallocate(u,u_rand,rho,D,iden,Sx,Sy,Sz,SO,Lp,Lz,Jz,Trans,rot,sig,dc,sig_dc,sig_dc_afm,work,rwork,W)

          enddo  main

         ! Finally, close the input file
          close(1)

    end program

