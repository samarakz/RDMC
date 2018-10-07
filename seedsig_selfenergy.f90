module seedsig_selfenergy
implicit none

contains
!*******************************************************************************************
!
!  selfenergy_write
!
!>    it writes the static self-energy to file.
!
!     @TODO: Add the dynamical part of the self-energy to the input.
!
!*******************************************************************************************

subroutine selfenergy_write(sig_static,id)
implicit none
complex(KIND=8),intent(in)       :: sig_static(:,:)   !< The static self-energy
integer,intent(in)               :: id(:,:)           !< the id array
character(len=40)                :: label             ! cluster label to be written
integer,parameter                :: sig_unit=742      ! Binary file "sig" unit
integer                          :: msize             ! size of the matrices
integer                          :: line, n           ! loop counters
integer                          :: nmsl              ! total number of frequencies in the linear mesh
integer                          :: emeshsize         ! size of the real energy mesh
integer                          :: nfermi            ! position of the fermi level within emesh
logical                          :: sig_nonzeroreal   ! write sig_real to file
integer,parameter                :: sig_version = 3   ! Current version/format of the sig file
integer                          :: nwisdom           ! size of wisdom array
integer                          :: lda_id_size0      ! the lda_id_size
integer                          :: er                ! error message
real(KIND=8)                     :: nel               ! number of electrons for this sig
real(KIND=8)                     :: nel_old           ! number of electrons for the previous sig
real(KIND=8)                     :: sigdiff           ! convergency parameter of the self-energy
real(KIND=8)                     :: green_mu          ! chemical potential
real(KIND=8)                     :: green_mu_old      ! old chemical potential
real(KIND=8)                     :: TTT               ! Temperature in Ry
real(KIND=8)                     :: emeshdiff         ! spacing between the points of the real energy mesh
real(KIND=8)                     :: etot_dccorr       !
logical                          :: hsig_updated      ! hsig_updated flag for this iteration
logical                          :: sig_green_mu      ! sigma was generated with the current my
logical                          :: sig_gave_lda      ! sigma generated the current lda hamiltonian
logical                          :: cfflag            !

   ! Get the label
   call lda_id2label(id(:,1),label)

   ! Open the self-energy file as binary file
   open(sig_unit,file="sig-"//TRIM(ADJUSTL(label)),form='unformatted',status='replace',iostat=er)

   ! Prepare the header of the sig binary file
   n = 1
   nmsl = 1                  
   nwisdom = 1                 
   TTT = 0d0                  
   emeshsize = 1              
   nfermi = 1                 
   emeshdiff = 0d0           
   green_mu = 0d0             
   green_mu_old = 0d0        
   nel = 0d0                 
   nel_old = 0d0             
   sigdiff = 1d0             
   hsig_updated = .false.     ! Tells the subroutine that sigma or H_LDA was updated
   sig_green_mu = .false.     ! Tells the subroutine that sigma was calculated using this green_mu
   sig_gave_lda=.false.       ! sigma was used to generate the current LDA calculation


      write (sig_unit) sig_version      
      write (sig_unit) n,nmsl,nwisdom,TTT,emeshsize,nfermi,emeshdiff,green_mu,green_mu_old,hsig_updated,nel_old,nel,sigdiff,&
                   sig_green_mu,sig_gave_lda

   ! Write the cluster size
   n = size(id,2)  ! Should be 1
   lda_id_size0 = size(id,1)
   msize = size(sig_static,1)
   cfflag = .false.
   etot_dccorr = 0d0

   ! Check if sig_real should be written to file
   sig_nonzeroreal = .false.
   ! write the cluster header
   write(sig_unit) lda_id_size0,n,cfflag,sig_nonzeroreal,etot_dccorr
   write(sig_unit) id       ! The id array gives a unique label of the cluster

   ! write the self-energy arrays
   write(sig_unit) sig_static ! wisdom
   write(sig_unit) sig_static ! sig_static
   write(sig_unit) sig_static ! sig
   if (sig_nonzeroreal) write(sig_unit) sig_static ! sig_real

   ! Close file
   close(sig_unit)

   return

end subroutine selfenergy_write

!*******************************************************************
! id2label
!
!>    generates the cluster label from an id array
!!
!*******************************************************************
   pure subroutine lda_id2label(id,label)
   implicit none
   integer,intent(in)            :: id(:)   !< identification array
   integer                       :: lda_id_size
   character(len=*),intent(out)  :: label             !< cluster label
   integer                       :: j                 !< loop counter
   character(len=3)              :: label_element     !< label for the cluster

      lda_id_size = size(id)
      label = ""
      do j=1,lda_id_size
         if (abs(id(j)) .ge. 100) then
            write(label_element,'(i3)') abs(id(j))
         else
            write(label_element,'(i2.2)') abs(id(j))
         endif
         label = TRIM(label) // TRIM(ADJUSTL(label_element))
      enddo
   end subroutine

end module
