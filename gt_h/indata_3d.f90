MODULE indata

  USE static

  IMPLICIT NONE

  SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Input parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER :: nband_v,ni,nf,ncell,nk1,nkplus,mm1,mplus,mm2,nm2,nmod,NKGt,NMODES
  INTEGER :: nx,ny,nz,Ndx,Ndy,Ndz,nkx,nkyz,nky,nkz,ngt,Ncy,Ncz
  INTEGER :: NEK
  INTEGER :: Nomp
  
  real(dp) :: dx, dy, dz, t0x, t0y, t0z
  REAL(DP) :: ac,ac1,ac2,ac3,Ecut,omega,delta_gap

  REAL(DP), ALLOCATABLE :: Ev(:)
  REAL(DP), ALLOCATABLE :: Kyz(:)
  REAL(DP), ALLOCATABLE :: KGt_kyz(:,:,:)
  complex(dp), allocatable ::  HL(:,:,:), TL(:,:,:), U_LCBB(:,:,:), U_PSI(:,:,:,:,:)

  
  LOGICAL :: refine, gap_corr
  LOGICAL :: domag, lspinorb, okvan
  
  CHARACTER(LEN=100) :: outdir, indir2
  CHARACTER(LEN=100) :: input_file_DFT
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  NAMELIST /indata_dimensionality/             &
       &  ncell,                               &
       &  ncy,                                 &
       &  ncz,                                 &
       &  Nomp,                                &
       &  Nky,                                 &
       &  Nkz
  NAMELIST /indata_basis/                      &
       &  nband_v,                             &
       &  ni,                                  &
       &  nf,                                  &
       &  refine,                              &
       &  nkx,                                 &
       &  nkplus,                              &
       &  Mplus,                               &
       &  Ecut,                                &
       &  delta_gap
  NAMELIST /indata_basis2/                     &
       &  nm2,                                 &
       &  mm2,                                 &   
       &  indir2
  NAMELIST /indata_cell/                       &
       &  ac1,                                 &
       &  ac2,                                 &
       &  ac3
 NAMELIST /indata_inout/                       &
       &  outdir,                              &
       &  input_file_DFT

CONTAINS


  SUBROUTINE indata_grid()
    IMPLICIT NONE

    integer :: i, j

    ac=ac1

    omega=ac1*ac2*ac3
    write(*,*)'unit cell volume',omega

    WRITE(*,*)'NKY=',NKY
    WRITE(*,*)'NKZ=',NKZ

    Nkyz=nky*nkz
    write(*,*)'Nkyz',Nkyz

  END SUBROUTINE indata_grid


  SUBROUTINE indata_readinput()

    IMPLICIT NONE    
    delta_gap = 0.0_dp
    READ(*,NML=indata_dimensionality)
    READ(*,NML=indata_basis)

    gap_corr=.false.
    if(delta_gap /= 0.0_dp)  gap_corr=.true.
    write(*,*)'gap_corr',gap_corr,delta_gap
    nk1=nkx
    nek=0
    
    if(nek > 0)then
       allocate(Ev(nEk))
       read(*,*)Ev
       write(*,*)'ev',ev
    end if
    
    if(ncell==2)READ(*,NML=indata_basis2)

    READ(*,NML=indata_cell)
    write(*,*)  ac1
    write(*,*)  ac2
    write(*,*)  ac3

    READ(*,NML=indata_inout)
    write(*,*)'reading input file',input_file_DFT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  END SUBROUTINE indata_readinput


    
SUBROUTINE invertR(A,n,np)
  
  IMPLICIT NONE
  
  INTEGER :: info,n,np
  INTEGER, DIMENSION(np) :: ipiv
  REAL(DP), dimension(np,np) :: A
  REAL(DP), dimension(np*np) :: work
  
  CALL dgetrf(n,n,A,n,ipiv,info)
  CALL dgetri(n,A,n,ipiv,work,np*np,info)
  
END SUBROUTINE invertR



SUBROUTINE SUB_DEF_Z0(Mi,Mf,ny,A,subband)
implicit none
integer :: ny,mi,mf
real(dp) :: subband(1:(Mf-Mi+1))
!real(dp) :: WE(1:(ny))
complex(dp) :: A(1:NY,1:NY),Uii(1:NY,1:Mf-Mi+1)
integer :: INFO!,LWORK
integer, allocatable :: iwork(:), supp(:)
complex(dp), allocatable :: work(:)
real(dp), allocatable :: rwork(:)
REAL(DP), EXTERNAL :: DLAMCH


allocate(WORK(20*ny))
allocate(RWORK(24*ny))
allocate(IWORK(10*ny))
allocate(Supp(2*ny))

call ZHEEVR('N','I','U',ny,A,ny,0.0,0.0,mi,mf,2*DLAMCH('S'),Mf-Mi+1,subband,Uii,ny,SUPP,WORK,20*ny,RWORK,24*ny,IWORK,10*ny,INFO)

deallocate(work)
deallocate(rwork)
deallocate(supp)
deallocate(iwork)
!write(*,*)'info=',info
if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: SUB_DEF HAS FAILED. INFO=',INFO
   stop
endif

END SUBROUTINE SUB_DEF_Z0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function stringa(ii) result(rr)
    integer, intent(in) :: ii
    character(:), allocatable :: rr
    character(range(ii)+2) :: tt
    write(tt,'(i0)') ii
    rr = trim(tt)
  end function stringa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

END MODULE indata

