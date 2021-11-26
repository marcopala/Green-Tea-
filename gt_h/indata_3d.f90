MODULE indata

  USE static

  IMPLICIT NONE

  SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Energy integral parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(DP), ALLOCATABLE :: en_global(:) !Vettore globale delle energie in 0-EOp
  REAL(DP), ALLOCATABLE :: w(:)  !Vettore globale dei pesi
  REAL(DP) :: emin, emax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER :: NUMEL_3D !number of elements 3D geometry
  INTEGER :: NUMN_3D  !number of nodes 3D geometry
  INTEGER :: NTOT_X   !total number of nodes in x direction
  INTEGER :: NTOT_Y   !total number of nodes in y direction
  INTEGER :: NTOT_Z   !total number of nodes in z direction
  INTEGER :: NUMBOUND_3D  !number of diricleth boundary nodes
  INTEGER :: NUMELECT_3D  !number of electrodes
  INTEGER :: LWORK_3D     !number of unknows in Poisson solution
  INTEGER :: nnz_3D       !number of non-zero entrues of the Poisson laplacian 
 
  REAL(DP), ALLOCATABLE :: coord_3D(:,:)     !nodes coordinates (number in progressive order)
  REAL(DP), ALLOCATABLE :: coord_3D_ord(:,:) !nodes coordinates ("unknows first" order)
  INTEGER, ALLOCATABLE :: lista_3D(:,:)              !element nodes (number in progressive order)
  INTEGER, ALLOCATABLE :: lista_3D_ord(:,:)          !element nodes ("unknows first" order)
  INTEGER, ALLOCATABLE :: whichkind_3D(:)
  INTEGER, ALLOCATABLE :: whichkind_3D_ord(:)
  INTEGER, ALLOCATABLE :: type_3D(:,:,:)

  REAL(DP), ALLOCATABLE :: dop_vec(:)
  REAL(DP), ALLOCATABLE :: Ev(:)
  INTEGER, ALLOCATABLE :: trap_vec(:,:)
  INTEGER, ALLOCATABLE :: coul(:)

  REAL(DP) :: potelectr,potelectr0

  REAL(DP), ALLOCATABLE ::   epsilon_3D(:)  ! dielectric constant for each element
  CHARACTER(LEN=12), ALLOCATABLE ::  nomelectr(:) ! nome degli elettrodi
  INTEGER, ALLOCATABLE :: nodelectr(:) !numero di nodi in ogni elettrodo
  
  INTEGER, ALLOCATABLE :: map_3D(:) !Transformation map between the two orders
  
  INTEGER, ALLOCATABLE :: connect_3D(:,:) !Connectivity for the unknows

  INTEGER, PARAMETER :: NMAX = 2048
  REAL(DP), PARAMETER :: LMAX = 102.4d-9 !in metres

  REAL(DP) :: c0,c1,c2,c3,c4,c5

  real(4) :: t1,t2,t3,t4,tt1,tt2

  INTEGER :: Nop, NEK

  INTEGER :: NUMVG,NUMVD,NUMBZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Input parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER :: source_len
  INTEGER :: drain_len 
  INTEGER :: gate_len

  REAL(DP) :: source_dop_val
  REAL(DP) :: channel_dop_val
  REAL(DP) :: drain_dop_val
   
  INTEGER :: tox_h_res
  INTEGER :: tox_w_res
  INTEGER :: tox_h_ch
  INTEGER :: tox_w_ch

  INTEGER :: tsi_h_res
  INTEGER :: tsi_w_res
  INTEGER :: tsi_h_ch
  INTEGER :: tsi_w_ch

  INTEGER :: spacer

  LOGICAL :: top_gate
  LOGICAL :: bot_gate
  LOGICAL :: sx_gate
  LOGICAL :: dx_gate

  CHARACTER(3) :: Tdir
  CHARACTER(1) :: chtype  
  REAL(DP) :: T_xx       
  REAL(DP) :: T_yy       
  REAL(DP) :: T_zz       
  REAL(DP) :: T_yz       
  REAL(DP) :: T_xz       
  REAL(DP) :: T_xy       
  INTEGER :: Ndeltax, Ndeltay, Ndeltaz
  REAL(DP) :: deltax
  REAL(DP) :: deltay
  REAL(DP) :: deltaz

  INTEGER :: nsolv,nsolc
  CHARACTER(10), allocatable :: material(:)
  REAL(DP), allocatable :: mole(:)
  INTEGER, allocatable :: NX_REG(:)
  INTEGER :: NRG

  REAL(DP) :: ERROR_INNER
  INTEGER :: MAX_ITER_INNER
  REAL(DP) :: ERROR_OUTER
  INTEGER :: MAX_ITER_OUTER
  REAL(DP) :: alphapot       
  REAL(DP) :: TEMP
  INTEGER :: NKT

!!! EPM 
  INTEGER :: nband,nband_v,ni,nf,ncell,nk1,nkplus,mm1,mplus,mm2,nm2,nkx,nkx_d,Nk,nkyz,nky,nkz,nmod,Ncut,nx,ny,nz,ns,ns_d,NKGt,NMODES
  INTEGER :: max_g,max_gx,max_gy,max_gz,Nd,Ndx,Ndy,Ndz,np,ng,ngv,ngx,ngt,npol,Ncy,Ncz,Ncx_D
  INTEGER, ALLOCATABLE :: indicione(:), indicino(:), indgv(:,:,:), ijg(:,:)!, Nijg(:,:) !,ig(:),igg(:)

  integer, allocatable :: iky(:),ikz(:),Nkgt_kyz(:),Nk_kyz(:), Nez(:), Ney(:), Nex(:)

  real(dp) :: dx, dy, dz, t0x, t0y, t0z, Kmax
  REAL(DP) :: b1(3),b2(3),b3(3),T(3),R(3),A1(3),A2(3),A3(3)

  REAL(DP), ALLOCATABLE :: kx(:), kx_d(:), ky(:), kz(:), Kyz(:),  KGt(:,:), hkl(:,:), hklv(:,:), G(:,:), G2(:), Gv(:,:)
  REAL(DP), ALLOCATABLE :: deg_ky(:), deg_kz(:), deg_kyz(:), kv_max(:), kc_min(:), Utc(:,:), Utv(:,:), KGt_kyz(:,:,:)
  complex(dp), allocatable ::  HL(:,:,:), TL(:,:,:), U_LCBB(:,:,:), U_PSI(:,:,:,:,:)

  REAL(DP), ALLOCATABLE :: form_factor(:,:,:)

!!! MATERIAL
  REAL(DP) :: E_GAP,ref,ref_ev,ref_ec,ac,ac1,ac2,ac3,g_spin,Ecut,omega,zeta,R0_a,R0_c,A0_a,B0_a,A2_a,R2_a,A0_c,B0_c,A2_c,R2_c,mu,eta,delta_gap
  LOGICAL :: refine, gap_corr

  CHARACTER :: domag, lspinorb, okvan

!!! ROUGHNESS
  REAL(DP) :: rms
  REAL(DP) :: c_len
  INTEGER :: seed1
  INTEGER :: seed2
  INTEGER :: seed3
  INTEGER :: seed4

  REAL(DP) :: Dt
  REAL(DP) :: Ebound
  INTEGER :: trap_d
  INTEGER :: trap_seed

  REAL(DP) :: L_imp
  REAL(DP) :: Dit
  INTEGER :: charge_type
  INTEGER :: coul_seed
  REAL(DP) :: Hk
  LOGICAL :: VXC

  REAL(DP) :: Eop
  INTEGER :: Nop_g
  INTEGER :: Nop_f2
  INTEGER :: NOP_f3 
  REAL(DP) :: Dop_f2
  REAL(DP) :: Dop_f3
  REAL(DP) :: Dop_g
  INTEGER :: Nsub, Nomp
  REAL(DP) :: Dac_e
  REAL(DP) :: Dac_h
  REAL(DP) :: SCBA_tolerance
  INTEGER :: SCBA_max_iter
  REAL(DP), ALLOCATABLE :: Dop_g_mat(:)
  REAL(DP), ALLOCATABLE :: Dac_e_mat(:)
  REAL(DP), ALLOCATABLE :: Dac_h_mat(:)

  REAL(DP) :: VGMIN
  REAL(DP) :: VGMAX
  REAL(DP) :: DELTAVG
  REAL(DP) :: VDMIN
  REAL(DP) :: VDMAX
  REAL(DP) :: DELTAVD
  REAL(DP) :: BZMIN
  REAL(DP) :: BZMAX
  REAL(DP) :: DELTABZ
  REAL(DP) :: mus
  REAL(DP) :: mud
  REAL(DP) :: workgate
  REAL(DP) :: DIEL_SC
  REAL(DP) :: DIEL_OX
  CHARACTER(LEN=100) :: outdir, indir2
  CHARACTER(LEN=100) :: input_file_DFT
           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Auxiliary variables
  INTEGER  :: NUMBOUNDOLD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!All integer quantities are referred as NUMBER OF ELEMENTS!!!!!!!!!!!!!!!!!!  
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


    ncx_d=(source_len+drain_len+gate_len+2*spacer) 
    write(*,*)'NCX_D=',NCX_D
    
!    ncy = tsi_w_ch
    write(*,*)'NCY_D=',NCY

!    ncz= tsi_h_ch 
    write(*,*)'NCZ_D=',NCZ
 
    
    write(*,*)'Ndeltax=',Ndeltax
    write(*,*)'Ndeltay=',Ndeltay
    write(*,*)'Ndeltaz=',Ndeltaz

    !nky=1+Ncy/2!2*ncy
    !nkz=1+Ncz/2!2*ncz

    WRITE(*,*)'NKY=',NKY
    WRITE(*,*)'NKZ=',NKZ

    allocate(ky(nky))
    allocate(kz(nkz))
    


    allocate(deg_ky(nky),deg_kz(nkz))
    deg_ky=2.0_dp    
    deg_kz=2.0_dp

    do i=1,nky
       ky(i) = dble(i-1)/dble(Ncy) !1.0_dp*dble(i-ncy-1)/dble(ncy) ! 2pi/a0 units
       if(i==1 .or. i==nky) deg_ky(i)=1.0_dp
       write(*,*)'Ky',i,ky(i),deg_ky(i)
    end do
    do j=1,nkz
       kz(j) = dble(j-1)/dble(Ncz) !1.0_dp*dble(j-ncz-1)/dble(ncz) ! 2pi/a0 units
       if(j==1 .or. j==nkz) deg_kz(j)=1.0_dp
       write(*,*)'Kz',j,kz(j),deg_kz(j)
    end do

    Nkyz=nky*nkz
    write(*,*)'Nkyz',Nkyz
    allocate(deg_kyz(Nkyz))
    deg_kyz=deg_ky*deg_kz

    forall(j=1:nkz, i=1:nky) deg_kyz(i+(j-1)*nky)=deg_ky(i)*deg_kz(j)

    do j=1,nkz
       do i=1,nky 
          write(*,*)'deg_kyz',i+(j-1)*nky,deg_kyz(i+(j-1)*nky)
       end do
    end do
!stop

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


  
 SUBROUTINE map_potential(potA,potB,which,map,nx,ny,nz,lwork)
 
 REAL(DP), INTENT(OUT) :: potA(1:nx,1:ny,1:nz)

 REAL(DP), INTENT(IN) :: potB(0:lwork-1)
 INTEGER, INTENT(IN) :: which(0:nx*ny*nz-1)
 INTEGER, INTENT(IN) :: map(0:nx*ny*nz-1)
 INTEGER, INTENT(IN) :: nx
 INTEGER, INTENT(IN) :: ny
 INTEGER, INTENT(IN) :: nz
 INTEGER, INTENT(IN) :: lwork

 INTEGER :: ii,x_index,y_index,z_index

 potA=0.0_dp
 
 DO ii=0, nx*ny*nz-1

 x_index=ii/(ny*nz) +1
 y_index=mod(mod(ii,ny*nz),ny) +1
 z_index=mod(ii,ny*nz)/ny +1
 
 IF(which(ii).ge.1)THEN
 !Gated region added to the potential defined on the workspace only
 potA(x_index,y_index,z_index)=potelectr!(1)
 ELSE
 potA(x_index,y_index,z_index)=potB(map(ii))
 END IF
 
 END DO

 END SUBROUTINE map_potential


  
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

