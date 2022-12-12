! Copyright or © or Copr. Marco Pala (February 24, 2022)

! e-mail: marco.pala@cnrs.fr

! This software is a computer program whose purpose is 
! to perform self-consistent simulations of nanosystems with a full ab initio approach
! by using the density functional theory and the non-equilibrium Green's function method.

! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 

! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited liability. 

! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 

! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.

MODULE indata

  USE static

  IMPLICIT NONE

  SAVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Energy integral parameters!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(DP), ALLOCATABLE :: en_global(:) !Energy vector in the interval [0 Eop]
  REAL(DP), ALLOCATABLE :: w(:)  !weights
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
  INTEGER :: nnz_3D       !number of non-zero entries of the laplacian
  INTEGER :: NUMBOUNDOLD 
 
  REAL(DP), ALLOCATABLE :: coord_3D(:,:)      !nodes coordinates (number in progressive order)
  REAL(DP), ALLOCATABLE :: coord_3D_ord(:,:)  !nodes coordinates ("unknows first" order)
  INTEGER,  ALLOCATABLE :: list_3D(:,:)       !element nodes (number in progressive order)
  INTEGER,  ALLOCATABLE :: list_3D_ord(:,:)   !element nodes ("unknows first" order)
  INTEGER,  ALLOCATABLE :: whichkind_3D(:)
  INTEGER,  ALLOCATABLE :: whichkind_3D_ord(:)
  INTEGER,  ALLOCATABLE :: type_3D(:,:,:)

  REAL(DP), ALLOCATABLE :: dop_vec(:)
  INTEGER,  ALLOCATABLE :: coul(:)
  REAL(DP)              :: potelectr, potelectr0
  REAL(DP), ALLOCATABLE ::   epsilon_3D(:)  ! dielectric constant for each element
  INTEGER, ALLOCATABLE  :: nodelectr(:) 
  INTEGER, ALLOCATABLE  :: map_3D(:) !Transformation map between the two orders
  INTEGER, ALLOCATABLE  :: connect_3D(:,:) !Connectivity for the unknows

  INTEGER  :: NUMVG,NUMVD,NUMBZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Input parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER  :: source_len
  INTEGER  :: drain_len 
  INTEGER  :: gate_len

  REAL(DP) :: source_dop_val
  REAL(DP) :: channel_dop_val
  REAL(DP) :: drain_dop_val
   
  INTEGER  :: tox_lft, to2_lft
  INTEGER  :: tox_rgt, to2_rgt
  INTEGER  :: tox_top, to2_top
  INTEGER  :: tox_bot, to2_bot

  INTEGER  :: tsc_h
  INTEGER  :: tsc_w

  INTEGER  :: spacer

  LOGICAL  :: top_gate
  LOGICAL  :: bot_gate
  LOGICAL  :: lft_gate
  LOGICAL  :: rgt_gate
  LOGICAL  :: schottky_source
  LOGICAL  :: onlyT
  LOGICAL  :: in_pot
  LOGICAL  :: magnetic

  CHARACTER(2) :: updw
  CHARACTER(1) :: chtype  
  INTEGER  :: Ndeltax, Ndeltay, Ndeltaz
  REAL(DP) :: deltax
  REAL(DP) :: deltay
  REAL(DP) :: deltaz

  INTEGER  :: nsolv,nsolc
  
  REAL(DP) :: ERROR_INNER
  INTEGER  :: MAX_ITER_INNER
  REAL(DP) :: ERROR_OUTER
  INTEGER  :: MAX_ITER_OUTER
  REAL(DP) :: alphapot       
  REAL(DP) :: TEMP
  INTEGER  :: NKT

  INTEGER  :: nband,nband_c,nband_v,nk1,MM,ni,nf,nkx,nkx_d,Nk,nkyz,nky,nkz,nmod,nmodes_ph
  INTEGER  :: nx, ny, nz, NKGt, NMODES
  INTEGER  :: Nrx,Nry,Nrz,ngt,npol,Ncy,Ncz,Ncx_D
  integer,  allocatable :: Nez(:), Ney(:), Nex(:), nband_val(:), off_k_nvb(:)
  real(dp)              :: dx, dy, dz
  REAL(DP), ALLOCATABLE :: kx(:), ky(:), kz(:), Kyz(:), KGt(:,:), k_vec(:,:), KGt_kyz(:,:,:), kq_vec(:,:) 
  REAL(DP), ALLOCATABLE :: deg_ky(:), deg_kz(:), deg_kyz(:), kv_max(:,:), kc_min(:,:), off_set(:)

  LOGICAL,  ALLOCATABLE :: k_selec(:)
  LOGICAL               :: phonons, dfpt, vec_field_new, vec_field_old
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type I_BLOCKS
     integer, allocatable :: I(:)
  end type I_BLOCKS
  type(I_blocks),allocatable :: ind_kx(:,:), ind_bnd(:,:)
  
  type H_blocks
     complex(dp), allocatable :: H(:,:)
  end type H_blocks
  type(H_blocks),allocatable :: HL(:,:), TL(:,:), ULCBB(:,:)
  type(H_blocks),allocatable :: Si_m05(:,:), Si(:,:)
  
  type M_blocks
     COMPLEX(DP), ALLOCATABLE :: M(:,:,:)
  end type M_blocks
  type(M_blocks),allocatable :: el_ph_mtrx(:,:,:,:)

  type K_tensor
     complex(dp), allocatable :: K(:,:,:,:)
  end type K_tensor
  type(k_tensor),allocatable :: U_PSI(:,:), AJ(:,:), BJ(:,:), CJ(:,:)
  

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  REAL(DP)  :: E_GAP,ref,ac,ac1,ac2,ac3,g_spin,Ecut,mu,eta

  CHARACTER :: domag, lspinorb, okvan
  
  INTEGER               :: num_mat, num_reg, num_het, nqs, nqmodes
  INTEGER,  ALLOCATABLE :: imat(:), ihet(:), nc_reg(:), typ_mat(:), ind_q(:,:)
  INTEGER,  ALLOCATABLE :: nm_mat(:), mat_1(:), mat_0(:),ihh(:,:)
  REAL(DP), ALLOCATABLE :: ref_ev(:),ref_ec(:)
  REAL(DP), ALLOCATABLE :: omega_q(:,:),x_q(:,:)

  REAL(DP) :: Eop
  INTEGER  :: Nop_g
  REAL(DP) :: Dop_g
  INTEGER  :: Nsub,Nomp
  REAL(DP) :: Dac, Dac_e
  REAL(DP) :: Dac_h
  REAL(DP) :: SCBA_alpha
  REAL(DP) :: SCBA_tolerance
  INTEGER  :: SCBA_max_iter

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
  REAL(DP) :: DIEL_O2
  
  CHARACTER(LEN=60) :: outdir
  CHARACTER(LEN=60) :: inputdir
  CHARACTER(LEN=100) :: input_file_DFPT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!All integer quantities are referred as NUMBER OF ELEMENTS!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  NAMELIST /indata_mesh/                       &
       &  Ndeltax,                             &
       &  Ndeltay,                             &
       &  Ndeltaz  
  NAMELIST /indata_lenghts/                    &
       &  source_len,                          &
       &  drain_len,                           &
       &  gate_len,                            &
       &  spacer
  NAMELIST /indata_oxide/                      &
       &  tox_top,                             &
       &  tox_bot,                             &
       &  tox_lft,                             &
       &  tox_rgt,                             &
       &  to2_top,                             &
       &  to2_bot,                             &
       &  to2_lft,                             &
       &  to2_rgt
  NAMELIST /indata_channel/                    &
       &  tsc_h,                               &
       &  tsc_w  
  NAMELIST /indata_regions/                    &
       &  num_mat,                             &
       &  num_reg
  NAMELIST /indata_heterojunctions/            &
       &  num_het                              
  NAMELIST /indata_doping/                     &
       &  source_dop_val,                      &
       &  drain_dop_val,                       &
       &  channel_dop_val
  NAMELIST /indata_gate/                       &
       &  top_gate,                            &
       &  bot_gate,                            &
       &  lft_gate,                            &
       &  rgt_gate                            
  NAMELIST /indata_device/                     &
       &  chtype,                              &
       &  DIEL_SC,                             &      
       &  DIEL_OX,                             &     
       &  DIEL_O2	   
  NAMELIST /indata_dimensionality/             &
       &  Ncy,                                 &
       &  Ncz,                                 &
       &  nkx,                                 &
       &  Nky,                                 &
       &  Nkz
  NAMELIST /indata_basis/                      &
       &  nsolv,                               &
       &  nsolc,                               &
       &  magnetic,                            &
       &  updw,                                &  
       &  npol,                                &
       &  g_spin,                              &
       &  Ecut
  NAMELIST /indata_cell/                       &
       &  ac1,                                 &
       &  ac2,                                 &
       &  ac3
  NAMELIST /indata_convergence/                &
       &  ERROR_INNER,                         &
       &  MAX_ITER_INNER,                      &
       &  ERROR_OUTER,                         &
       &  MAX_ITER_OUTER,                      &
       &  alphapot,                            &
       &  Nomp
  NAMELIST /indata_energy/                     &
       &  dfpt,                                &
       &  nqs,                                 &
       &  nqmodes,                             &
       &  Eop,                                 &
       &  Nop_g,                               &
       &  Dop_g,                               &
       &  Dac,                                 &
       &  Nsub,                                &
       &  SCBA_alpha,                          &
       &  SCBA_tolerance,                      &
       &  SCBA_max_iter,                       &
       &  TEMP,                                &
       &  NKT,                                 &
       &  eta                            
  NAMELIST /indata_stimulus/                   &
       &  in_Pot,                              &
       &  vec_field_new,                       &
       &  vec_field_old,                       &
       &  onlyT,                               &
       &  VGMIN,                               &
       &  VGMAX,                               &
       &  DELTAVG,                             & 
       &  VDMIN,                               &
       &  VDMAX,                               &
       &  DELTAVD,                             &
!       &  BZMIN,                               &
!       &  BZMAX,                               &
!       &  DELTABZ,                             &
       &  workgate                         
 NAMELIST /indata_inout/                       &
       &  outdir, inputdir,                    &
       &  input_file_DFPT


CONTAINS
  
  SUBROUTINE indata_grid()
    
    IMPLICIT NONE
    
    ac=ac1
    
    write(*,*)'unit cell volume',ac1*ac2*ac3
    
    ncx_d=(source_len+drain_len+gate_len+2*spacer)
    
    write(*,*)'NCX_D=',NCX_D
    write(*,*)'NCY_D=',NCY
    write(*,*)'NCZ_D=',NCZ
    
    write(*,*)'Ndeltax=',Ndeltax
    write(*,*)'Ndeltay=',Ndeltay
    write(*,*)'Ndeltaz=',Ndeltaz
    
    WRITE(*,*)'NKX=',NKx
    WRITE(*,*)'NKY=',NKY
    WRITE(*,*)'NKZ=',NKZ
    
    allocate(ky(1:nky))
    allocate(kz(1:nkz))
    
    allocate(deg_ky(nky),deg_kz(nkz))
    deg_ky=2.0_dp    
    deg_kz=2.0_dp  
    
  END SUBROUTINE indata_grid


  SUBROUTINE indata_readinput()

    IMPLICIT NONE    
    INTEGER       :: i,j,l,nc,iy,iz
    CHARACTER(20) :: comment
    
!!! DEFAULT VALUES !!!
    vec_field_new =.FALSE.
    vec_field_old =.FALSE.
    onlyT         =.FALSE.
    in_pot        =.FALSE.
    dfpt          =.FALSE.
    tox_top=0
    tox_bot=0
    tox_lft=0
    tox_rgt=0
    to2_top=0
    to2_bot=0
    to2_lft=0
    to2_rgt=0
    tsc_w=1
    tsc_h=1 
    spacer=0
    num_mat=1
    num_reg=1 
    num_het=0 
    top_gate = .false.
    bot_gate = .false.
    lft_gate = .false.
    rgt_gate = .false.     
    DIEL_SC=1.0d0
    DIEL_OX=1.0d0
    DIEL_O2=1.0d0
    ncy=1
    ncz=1
    nkx=1
    nky=1
    nkz=1
    source_dop_val  = 0.0d20
    drain_dop_val   = 0.0d20
    channel_dop_val = 0.0d20
    ERROR_INNER = 1.0d-10
    MAX_ITER_INNER=30
    ERROR_OUTER = 1.0d-3
    MAX_ITER_OUTER=30
    alphapot=1.0_dp
    Nomp=1
    Dop_g=0.0_dp
    Nop_g=1
    Dac=0.0_dp
    Nsub=3
    SCBA_alpha=1.0_dp
    SCBA_max_iter=50
    SCBA_tolerance=1.0d-3
    eta=1.0d-6
    TEMP=300.0_dp
    NUMBZ=1
    BZMIN=0.0_dp
    BZMAX=0.0_dp
!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    READ(*,NML=indata_mesh)
    READ(*,NML=indata_lenghts)
    READ(*,NML=indata_oxide)
    READ(*,NML=indata_channel)

    if(tox_lft+tox_rgt+tsc_w .ne. Ndeltay)then
       write(*,*)'structure error along the lateral direction (y-axis)'
       write(*,*)tox_lft,tox_rgt,tsc_w 
       stop
    end if
    if(tox_bot+tox_top+tsc_h.ne. Ndeltaz)then
       write(*,*)'structure error along the vertical direction (z-axis)'
       write(*,*)tox_bot,tox_top,tsc_h 
       stop
    end if

    READ(*,NML=indata_regions)             !        num_mat,    num_reg
    allocate(nc_reg(num_reg),typ_mat(num_reg))
    READ(*,*)comment
    write(*,*)comment
    do i=1,num_reg
       read(*,*)j,nc_reg(i),typ_mat(i)
    end do
    !check coherence of the inputs
    do i=1,num_reg
       if(typ_mat(i) > num_mat)stop
    end do
    nc=0
    do i=1,num_reg
       nc=nc+nc_reg(i)
    end do
    if( nc /= source_len+drain_len+gate_len+2*spacer )then
       write(*,*)'Error! nc_reg wrong...'
       stop
    end if
    allocate(imat(1:nc),ihet(1:nc+1))
    imat=0
    ihet=0
    l=0

    do i=1,num_reg
       do j=1,nc_reg(i)
          l=l+1
          imat(l)=typ_mat(i)
          !write(*,*)l,'imat',imat(l)
       end do
    end do
    READ(*,NML=indata_heterojunctions)
    if(num_het>0)then
       READ(*,*)comment
       write(*,*)comment
    allocate(mat_0(num_het),mat_1(num_het))
       allocate(ihh(num_mat,num_mat))
       ihh=0
       do i=1,num_het
          read(*,*)j,mat_0(i),mat_1(i)
          ihh(mat_1(i),mat_0(i))=j
       end do

    end if
    ihet(1)=imat(1)
    do i=2,nc
       if(imat(i) == imat(i-1))then
          ihet(i)=imat(i)
       else
          ihet(i)= num_mat+ihh(imat(i),imat(i-1))
       end if
    end do
    ihet(nc+1)=imat(nc)
    
    allocate(ref_ec(num_mat),ref_ev(num_mat))
    !!!stop

    READ(*,NML=indata_doping)
    READ(*,NML=indata_gate)
    READ(*,NML=indata_device)
    READ(*,NML=indata_dimensionality)
    READ(*,*)comment
    write(*,*)comment
    NKyz=nky*nkz
    write(*,*)'Nkyz',Nkyz
    allocate(k_vec(3,Nkyz))
    allocate(off_k_nvb(Nkyz))
    allocate(k_selec(NKyz))
    off_k_nvb=0
    do iz=1,nkz
       do iy=1,nky
          l = iy + (iz-1)*nky
          read(*,*) k_vec(2,l), k_vec(3,l), k_selec(l)  !!!!!, off_k_nvb(l)
       end do
    end do

    allocate(kq_vec(3,nkx*Nkyz))
    do iz=1,nkz
       do iy=1,nky
          l = iy + (iz-1)*nky
          if(nkx == 8)then
             kq_vec(1,1+(l-1)*nkx)=-0.375_dp
             kq_vec(2,1+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,1+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,2+(l-1)*nkx)=-0.25_dp
             kq_vec(2,2+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,2+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,3+(l-1)*nkx)=-0.125_dp
             kq_vec(2,3+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,3+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,4+(l-1)*nkx)= 0.0_dp
             kq_vec(2,4+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,4+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,5+(l-1)*nkx)= 0.125_dp
             kq_vec(2,5+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,5+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,6+(l-1)*nkx)= 0.25_dp
             kq_vec(2,6+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,6+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,7+(l-1)*nkx)= 0.375_dp
             kq_vec(2,7+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,7+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,8+(l-1)*nkx)= 0.5_dp
             kq_vec(2,8+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,8+(l-1)*nkx)= k_vec(3,l)
          else          if(nkx == 4)then
             kq_vec(1,1+(l-1)*nkx)=-0.25_dp
             kq_vec(2,1+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,1+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,2+(l-1)*nkx)= 0.0_dp
             kq_vec(2,2+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,2+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,3+(l-1)*nkx)= 0.25_dp
             kq_vec(2,3+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,3+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,4+(l-1)*nkx)= 0.5_dp
             kq_vec(2,4+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,4+(l-1)*nkx)= k_vec(3,l)
          else          if(nkx == 2)then
             kq_vec(1,1+(l-1)*nkx)= 0.0_dp
             kq_vec(2,1+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,1+(l-1)*nkx)= k_vec(3,l)
             kq_vec(1,2+(l-1)*nkx)= 0.5_dp
             kq_vec(2,2+(l-1)*nkx)= k_vec(2,l)
             kq_vec(3,2+(l-1)*nkx)= k_vec(3,l)
          else
             write(*,*)'pb with kq_vec definition: unexpected nkx',nkx
             stop
          end if
       end do
    end do
    
    magnetic='F'
    updw='ni'
    READ(*,NML=indata_basis)
    if(magnetic)then
       select case (updw)
       case( 'ni')
          write(*,*)'error: magnetization direction must be either up or dw'
          stop
       case('up')
           write(*,*)'magnetic polarization is up'
       case('dw')
           write(*,*)'magnetic polarization is dw'
       end select
    end if
    READ(*,*)comment
    write(*,*)comment
    allocate(nm_mat(num_mat),nband_val(num_mat),off_set(num_mat))
    do i=1,num_mat
       read(*,*)j,nm_mat(i),nband_val(i),off_set(i)
    end do

    READ(*,NML=indata_cell)

    write(*,*)  'ac1 =',ac1
    write(*,*)  'ac2 =',ac2
    write(*,*)  'ac3 =',ac3
    
    allocate(AJ(Nkyz,num_mat)) 
    allocate(BJ(Nkyz,num_mat)) 
    allocate(CJ(Nkyz,num_mat))  
    allocate(HL(Nkyz,num_mat))    
    allocate(TL(NKyz,num_mat+num_het))
    allocate(ULCBB(NKyz,num_mat))
    allocate(U_PSI(NKyz,num_mat))
    allocate(Si(Nkyz,num_mat))   
    allocate(Si_m05(Nkyz,num_mat))  

    
    deltax=ac1/dble(Ndeltax)
    deltay=ac2/dble(Ndeltay)
    deltaz=ac3/dble(Ndeltaz)
    write(*,*)'deltax=',deltax
    write(*,*)'deltay=',deltay
    write(*,*)'deltaz=',deltaz

    READ(*,NML=indata_convergence)
    READ(*,NML=indata_energy)
    READ(*,NML=indata_stimulus)
    READ(*,NML=indata_inout)

    write(*,*) 'Channel height  (no ox)=',tsc_h*deltaz
    write(*,*) 'Channel width (no ox)=',tsc_w*deltay 
    write(*,*) 'Channel length =',gate_len*deltax
     
    write(*,*) 'Drain height  (no ox)=',tsc_h*deltaz
    write(*,*) 'Drain width (no ox)=',tsc_w*deltay 
    write(*,*) 'Drain lenght  =',drain_len*deltax    
 
    write(*,*) 'Source height  (no ox)=',tsc_h*deltaz
    write(*,*) 'Source width (no ox)=',tsc_w*deltay 
    write(*,*) 'Source length =',source_len*deltax    

    write(*,*) 'Spacer length =',spacer*deltax

    write(*,*) 'Total length =', (gate_len+source_len+drain_len+2*spacer)*Ndeltax*deltax 
    write(*,*) 'Total height =', ac3+(to2_top+to2_bot)*deltaz
    write(*,*) 'Total width =', ac2+(to2_lft+to2_rgt)*deltay

    write(*,*) 'Damping factor= ', alphapot
    write(*,*) 'NKT', NKT

    
  
    phonons=.FALSE.
    if(dac > 0.0_dp .or. dop_g > 0.0_dp) phonons=.TRUE.
    if(dfpt)  phonons=.TRUE.

    if(phonons) then
       write(*,*)'PHONONS ON'
    else
       write(*,*)'PHONONS OFF, balistic simulation'
    end if

    if(phonons) allocate(el_ph_mtrx(NKyz,Nkx,NKyz,num_mat))
    
    if(phonons)then
       allocate(ind_kx(NKyz,num_mat))
       allocate(ind_bnd(NKyz,num_mat))
       if (dfpt) allocate(ind_q(nkx,nkyz))
    end if

    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF((mod(tsc_h-tsc_h,2).ne.0).or.(mod(tsc_w-tsc_w,2).ne.0))STOP "Error in structure declaration"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  END SUBROUTINE indata_readinput

SUBROUTINE indata_structure_init()

IMPLICIT NONE

INTEGER :: icc,ii,jj,offset,vertical_offset,lateral_offset,plane_index,inplane_y,inplane_z,type

icc=0

!!! THIS CONVERTS the number of unit cells along x in number of elements!!!!
source_len=source_len*Ndeltax
drain_len=drain_len*Ndeltax
gate_len=gate_len*Ndeltax
spacer=spacer*Ndeltax

!!! Number of NODES in the macrogrid
NTOT_X=(source_len+drain_len+gate_len+2*spacer)+1
NTOT_Y=(tsc_w+tox_lft+tox_rgt+to2_lft+to2_rgt)+1
NTOT_Z=(tsc_h+tox_top+tox_bot+to2_top+to2_bot)+1

write(*,*)'NTOT_X=',NTOT_X
write(*,*)'NTOT_Y=',NTOT_Y
write(*,*)'NTOT_Z=',NTOT_Z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NUMN_3D=NTOT_X*NTOT_Y*NTOT_Z  
NUMEL_3D=(NTOT_X-1)*(NTOT_Y-1)*(NTOT_Z-1) 

write(*,*)'Total number of nodes=',NUMN_3D
write(*,*)'Total number of elements=',NUMEL_3D

lateral_offset=1 !(((tsc_w+2*tox_w)-(tsc_w+2*tox_w_ch))/2 + 1)
vertical_offset=1 !(((tsc_h+2*tox_h)-(tsc_h+2*tox_h_ch))/2 + 1)

NUMBOUND_3D=0
NUMBOUNDOLD=0

!!! change june 2021
IF(schottky_source)THEN
   NUMBOUND_3D=NUMBOUND_3D+source_len*NTOT_X*NTOT_Y*NTOT_Z   
END IF
!!! 
IF(lft_gate)THEN
NUMBOUND_3D=NUMBOUND_3D+(gate_len)*NTOT_Z*lateral_offset
icc=icc+1
END IF
NUMBOUNDOLD=NUMBOUND_3D
IF(rgt_gate)THEN
NUMBOUND_3D=NUMBOUND_3D+(gate_len)*NTOT_Z*lateral_offset
icc=icc+1
END IF
NUMBOUNDOLD=NUMBOUND_3D
IF(top_gate)THEN
NUMBOUND_3D=NUMBOUND_3D+(gate_len)*NTOT_Y*vertical_offset
icc=icc+1
END IF
NUMBOUNDOLD=NUMBOUND_3D
IF(bot_gate)THEN
NUMBOUND_3D=NUMBOUND_3D+(gate_len)*NTOT_Y*vertical_offset
icc=icc+1
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!
!Removing double counting
!!!!!!!!!!!!!!!!!!!!!!!!!
IF(lft_gate.and.top_gate)THEN
NUMBOUND_3D=NUMBOUND_3D-vertical_offset*lateral_offset*(gate_len)
END IF
IF(lft_gate.and.bot_gate)THEN
NUMBOUND_3D=NUMBOUND_3D-vertical_offset*lateral_offset*(gate_len)
END IF
IF(rgt_gate.and.top_gate)THEN
NUMBOUND_3D=NUMBOUND_3D-vertical_offset*lateral_offset*(gate_len)
END IF
IF(rgt_gate.and.bot_gate)THEN
NUMBOUND_3D=NUMBOUND_3D-vertical_offset*lateral_offset*(gate_len)
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!

NUMELECT_3D=icc
write(*,*)'Number of gates',icc
write(*,*)'Total number of gate nodes',NUMBOUND_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Allocation of arrays!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ALLOCATE(coord_3D(1:3,0:NUMN_3D-1))
    coord_3D(:,:)=0.0_dp
    ALLOCATE(whichkind_3D(0:NUMN_3D-1))
    whichkind_3D(:)=0
    ALLOCATE(list_3D(1:8,0:NUMEL_3D-1))
    list_3D(:,:)=0
    ALLOCATE(epsilon_3D(0:NUMEL_3D-1))
    epsilon_3D(:)=0.0_dp 
    IF(NUMELECT_3D.gt.0)THEN
    potelectr=0.0_dp
    END IF
    ALLOCATE(type_3D(1:NTOT_X-1,1:NTOT_Y-1,1:NTOT_Z-1))
    type_3D=0

LWORK_3D = NUMN_3D - NUMBOUND_3D

write(*,*)'Number of unknows in Poisson solution=', LWORK_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!Definition of the real coordinates and!!!!!!!!!
!!!!!!!!!!!!of the grid point qualifiers!!!!!!!!!!!!! 
!! -1 Unknonw for Poisson & Schroedinger 2D
!!  0 Unknonw for Poisson
!!  1 Not an unknown
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

whichkind_3D=-1 !Default for unknows nodes
jj=0
DO ii=0, NUMN_3D-1

plane_index=ii/(NTOT_Y*NTOT_Z)
inplane_y=mod(ii-plane_index*(NTOT_Y*NTOT_Z),NTOT_Y)
inplane_z=(ii-plane_index*(NTOT_Y*NTOT_Z))/NTOT_Y

coord_3D(1,ii)=plane_index*deltax
coord_3D(2,ii)=inplane_y*deltay
coord_3D(3,ii)=inplane_z*deltaz

IF((inplane_y.eq.0).or.(inplane_z.eq.0).or.(inplane_y.eq.NTOT_Y-1).or.(inplane_z.eq.NTOT_Z-1))whichkind_3D(ii)=0 !Boundary node. Not un unknown for Schroedinger 2D

!!! change june 2021
IF(schottky_source)THEN
   IF((plane_index.le.(source_len)))THEN
      whichkind_3D(ii)=11 !set nodes
   END IF
END IF

IF((plane_index.gt.(source_len+spacer)).and.(plane_index.le.(source_len+spacer+gate_len)))THEN
!Gated area
IF(rgt_gate)THEN
IF(inplane_y.ge.((lateral_offset-1)+tox_lft+tox_rgt+tsc_w+to2_lft+to2_rgt))THEN
whichkind_3D(ii)=1 !Gate node
END IF
END IF
IF(lft_gate)THEN
IF(inplane_y.le.(lateral_offset-1))THEN
whichkind_3D(ii)=1 !Gate node
END IF
END IF
IF(bot_gate)THEN
IF(inplane_z.le.(vertical_offset-1))THEN
whichkind_3D(ii)=1 !Gate node
END IF
END IF
IF(top_gate)THEN
IF(inplane_z.ge.((vertical_offset-1)+tox_top+tox_bot+tsc_h+to2_top+to2_bot))THEN
whichkind_3D(ii)=1 !Gate node
END IF
END IF
END IF

!!!IF(plane_index.eq.0)whichkind_3D(ii)=1 !! THE DIRICHLECT NODES ARE LOCATED AT THE SOURCE (FIRST SLICE)

END DO
write(*,*)'number of zero whichkind',jj
!stop
lateral_offset  =0!((tsc_w+tox_lft+tox_rgt)-(tsc_w+2*tox_w_ch))/2
vertical_offset =0!((tsc_h+2*tox_h)-(tsc_h+2*tox_h_ch))/2

DO ii=0, NUMEL_3D -1

plane_index=ii/(NUMEL_3D/(NTOT_X-1))
offset=plane_index*(NUMEL_3D/(NTOT_X-1))
inplane_y=mod(ii-offset,(NTOT_Y-1))
inplane_z=(ii-offset)/(NTOT_Y-1)

list_3D(1,ii)=plane_index*NTOT_Y*NTOT_Z+mod(ii,(NTOT_Y-1)*(NTOT_Z-1))+inplane_z
list_3D(2,ii)=list_3D(1,ii)+1
list_3D(3,ii)=list_3D(1,ii)+1+NTOT_Y
list_3D(4,ii)=list_3D(1,ii)+NTOT_Y
list_3D(5,ii)=(plane_index+1)*NTOT_Y*NTOT_Z+mod(ii,(NTOT_Y-1)*(NTOT_Z-1))+inplane_z
list_3D(6,ii)=list_3D(5,ii)+1
list_3D(7,ii)=list_3D(5,ii)+1+NTOT_Y
list_3D(8,ii)=list_3D(5,ii)+NTOT_Y


epsilon_3D(ii)=DIEL_0*DIEL_O2     !Default as external oxide
type=2

IF((plane_index.lt.source_len).or.(plane_index.ge.(source_len+2*spacer+gate_len)))THEN
!RESERVOIR REGION!
IF((inplane_z.ge.(to2_bot)).and.(inplane_z.lt.(to2_top+tox_bot+tsc_h+tox_top)).and. &
   (inplane_y.ge.(to2_lft)) .and.(inplane_y.lt.(to2_lft +tox_lft +tsc_w+tox_rgt)))THEN
IF((inplane_z.ge.(tox_bot+to2_bot)).and.(inplane_z.lt.(tox_bot+to2_bot+tsc_h)).and. &
   (inplane_y.ge.(tox_lft +to2_lft)) .and.(inplane_y.lt.(tox_lft +to2_lft +tsc_w)))THEN
epsilon_3D(ii)=DIEL_0*DIEL_SC  !Semiconductor
type=1
ELSE
epsilon_3D(ii)=DIEL_0*DIEL_OX  !Internal Oxide
type=2
END IF
END IF
ELSE
IF((plane_index.gt.(source_len+spacer)).or.(plane_index.le.(source_len+spacer+gate_len)))THEN
!GATED REGION!
IF((inplane_z.ge.(to2_bot+vertical_offset)).and.(inplane_z.lt.(to2_bot+tox_bot+tsc_h+vertical_offset+tox_top)).and. &
   (inplane_y.ge.(to2_lft +lateral_offset)) .and.(inplane_y.lt.(to2_lft +tox_lft +tsc_w+lateral_offset+tox_rgt)))THEN
IF((inplane_z.ge.tox_bot +to2_bot+vertical_offset).and.(inplane_z.lt.tox_bot+to2_bot+tsc_h+vertical_offset).and. &
   (inplane_y.ge.tox_lft  +to2_lft+lateral_offset)  .and.(inplane_y.lt.tox_lft +to2_lft +tsc_w+lateral_offset))THEN
epsilon_3D(ii)=DIEL_0*DIEL_SC !Semiconductor
type=1
ELSE
epsilon_3D(ii)=DIEL_0*DIEL_OX !Internal Oxide
type=2
END IF
END IF
IF((inplane_z.lt.vertical_offset).and.top_gate)epsilon_3D(ii)=DIEL_0*DIEL_METAL !Gate
IF((inplane_z.ge.(tox_bot+tox_top+tsc_h+to2_bot+to2_top+vertical_offset)).and.bot_gate)epsilon_3D(ii)=DIEL_0*DIEL_METAL !Gate
IF((inplane_y.lt.lateral_offset).and.rgt_gate)epsilon_3D(ii)=DIEL_0*DIEL_METAL !Gate
IF((inplane_y.ge.(tox_lft+tox_rgt+tsc_w+to2_lft+to2_rgt+lateral_offset)).and.lft_gate)epsilon_3D(ii)=DIEL_0*DIEL_METAL !Gate
ELSE
!SPACER REGION!
IF((inplane_z.ge.(to2_bot+vertical_offset)).and.(inplane_z.lt.(to2_bot+tsc_h+vertical_offset+tox_bot+tox_top)).and. &
   (inplane_y.ge.(to2_lft +lateral_offset)) .and.(inplane_y.lt.(to2_lft +tsc_w+lateral_offset+tox_lft+tox_rgt)))THEN
IF((inplane_z.ge.(tox_bot+to2_bot+vertical_offset)).and.(inplane_z.lt.(tox_bot+to2_bot+tsc_h+vertical_offset)).and. &
   (inplane_y.ge.(tox_lft +to2_lft +lateral_offset)) .and.(inplane_y.lt.(tox_lft +to2_lft +tsc_w+lateral_offset)))THEN
epsilon_3D(ii)=DIEL_0*DIEL_SC !Semiconductor
type=1
ELSE
epsilon_3D(ii)=DIEL_0*DIEL_OX !Internal  Oxide
type=2
END IF
END IF
END IF
END IF

type_3D(plane_index+1,inplane_y+1,inplane_z+1)=type

END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   Reordering of the nodes:    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 1. Unknowns Schroedinger+Poisson  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 2. Unknowns Poisson           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 3. Gate nodes                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ALLOCATE(coord_3D_ord(1:3,0:NUMN_3D-1))
    coord_3D_ord(:,:)=0.0_dp
    ALLOCATE(whichkind_3D_ord(0:NUMN_3D-1))
    whichkind_3D_ord(:)=0
    ALLOCATE(list_3D_ord(1:8,0:NUMEL_3D-1))
    list_3D_ord(:,:)=0
    ALLOCATE(map_3D(0:NUMN_3D-1))
    map_3D(:)=0
    
jj=0
    !//// Unknows Schro+Poisson////
    DO ii=0,NUMN_3D-1
       IF(whichkind_3D(ii).eq.-1) THEN
          map_3D(ii)=jj
          jj=jj+1
       END IF
    END DO
    !//// Unknowns Poisson //////
    DO ii=0,NUMN_3D-1   
       IF(whichkind_3D(ii).eq.0) THEN
          map_3D(ii)=jj
          jj=jj+1
       END IF
    END DO
    !///// gate ////////////
    DO ii=0,NUMN_3D-1
       IF(whichkind_3D(ii).eq.1) THEN
          map_3D(ii)=jj
          jj=jj+1
       END IF
    END DO
    
    !///// schottky ////////////
    DO ii=0,NUMN_3D-1
       IF(whichkind_3D(ii).eq.11) THEN
          map_3D(ii)=jj
          jj=jj+1
       END IF
    END DO


    IF (jj.ne.NUMN_3D) STOP 'indata_struttura_init: errore in reordering'

    DO ii=0,NUMN_3D-1
       whichkind_3D_ord(map_3D(ii))=whichkind_3D(ii)
       coord_3D_ord(1,map_3D(ii))=coord_3D(1,ii)
       coord_3D_ord(2,map_3D(ii))=coord_3D(2,ii)
       coord_3D_ord(3,map_3D(ii))=coord_3D(3,ii)
    END DO
    DO ii=0,NUMEL_3D-1
       list_3D_ord(1,ii)=map_3D(list_3D(1,ii))
       list_3D_ord(2,ii)=map_3D(list_3D(2,ii))
       list_3D_ord(3,ii)=map_3D(list_3D(3,ii))
       list_3D_ord(4,ii)=map_3D(list_3D(4,ii))
       list_3D_ord(5,ii)=map_3D(list_3D(5,ii))
       list_3D_ord(6,ii)=map_3D(list_3D(6,ii))
       list_3D_ord(7,ii)=map_3D(list_3D(7,ii))
       list_3D_ord(8,ii)=map_3D(list_3D(8,ii))
    END DO


    write(*,*) 'Calling INDATA_CONNECTIVITY...'
    CALL INDATA_CONNECTIVITY()
    write(*,*) '...done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!Non zero entries of the Poisson laplacian!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nnz_3D=0
    DO ii=0,LWORK_3D-1
       nnz_3D=nnz_3D+(connect_3D(ii,7)+1) 
    END DO

    write(*,*) "nnz_3D=",nnz_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!Doping and coulumb auxiliary vec!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ALLOCATE(dop_vec(0:LWORK_3D-1))
   ALLOCATE(coul(0:LWORK_3D-1))
   dop_vec=0.0_dp
   coul=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!Allocation of global energy and integral weights arrays!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    ALLOCATE(en_global(1:Nsub))
    ALLOCATE(w(1:Nsub))

    call gaussianlegendre(0.0_dp, 1.0_dp, en_global, w, Nsub)
    
    en_global=en_global*Eop
    w=w*Eop
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END SUBROUTINE indata_structure_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE indata_CONNECTIVITY()
  
  IMPLICIT NONE
  INTEGER  :: ii,jj,nel,icc,found
  REAL(DP) :: xi,yi,xj,yj,zi,zj
  INTEGER  :: checked(0:NUMN_3D-1)    
  
  ALLOCATE(connect_3D(0:LWORK_3D-1,1:7)) 
  
  connect_3D(:,:)=0
  
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !\\\\ internal nodes connectivity \\
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  icc=0
  found=0
  checked(:)=0
  
  DO nel=0,NUMEL_3D-1
     
     DO ii=1,8
        
        IF ((whichkind_3D_ord(list_3D_ord(ii,nel)).le.0)) THEN 
           
           xi=coord_3D_ord(1,list_3D_ord(ii,nel))
           yi=coord_3D_ord(2,list_3D_ord(ii,nel))
           zi=coord_3D_ord(3,list_3D_ord(ii,nel))
           
           DO jj=1,8
              
              found=0
              IF (whichkind_3D_ord(list_3D_ord(jj,nel)).le.0) THEN 
                 
                 IF(list_3D_ord(jj,nel).ne.list_3D_ord(ii,nel))THEN
                    
                    xj=coord_3D_ord(1,list_3D_ord(jj,nel))
                    yj=coord_3D_ord(2,list_3D_ord(jj,nel))
                    zj=coord_3D_ord(3,list_3D_ord(jj,nel))   
                    
                    IF (((xi.ne.xj).and.((yi.eq.yj).and.(zi.eq.zj))).or.&
                         ((xi.eq.xj).and.((yi.ne.yj).and.(zi.eq.zj))).or.&
                         ((xi.eq.xj).and.((yi.eq.yj).and.(zi.ne.zj))))THEN
                       DO icc=1,connect_3D(list_3D_ord(ii,nel),7)
                          IF(connect_3D(list_3D_ord(ii,nel),icc).eq.list_3D_ord(jj,nel)) found = 1
                       END DO
                       IF(found .ne. 1) THEN
                          connect_3D(list_3D_ord(ii,nel),7)=connect_3D(list_3D_ord(ii,nel),7)+1  
                          connect_3D(list_3D_ord(ii,nel),connect_3D(list_3D_ord(ii,nel),7))=list_3D_ord(jj,nel)
                       END IF
                    END IF
                    
                 END IF
                 
              END IF
           END DO
           
        END IF
     END DO
  END DO
  
END SUBROUTINE indata_CONNECTIVITY



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE indata_build_doping(doping,coul,s_val,d_val,c_val,n,lwork)

  REAL(DP), INTENT(OUT) :: doping(0:lwork-1)
  INTEGER,  INTENT(OUT) :: coul(0:lwork-1)

  REAL(DP), INTENT(IN)  :: s_val
  REAL(DP), INTENT(IN)  :: d_val
  REAL(DP), INTENT(IN)  :: c_val
  INTEGER,  INTENT(IN)  :: lwork
  INTEGER,  INTENT(IN)  :: n

  INTEGER :: ii, plane_index, inplane_y, inplane_z
  REAL(DP) :: dop
   
  doping=0.0_dp
  coul=0

  write(*,*)'tox_lft',tox_lft,'tsc_w',tsc_w
  write(*,*)'tox_bot',tox_bot,'tsc_h',tsc_h,'tox_top',tox_top

  open(unit=99,file='doping_profile.dat',status='unknown')
  DO ii=0, n-1

     plane_index=ii/(NTOT_Y*NTOT_Z)
     inplane_y=mod(ii-plane_index*(NTOT_Y*NTOT_Z),NTOT_Y)
     inplane_z=(ii-plane_index*(NTOT_Y*NTOT_Z))/NTOT_Y

     dop=0.0_dp

     IF((inplane_y .ge.tox_lft +to2_lft) .and.(inplane_y.le.(tox_lft +to2_lft +tsc_w)).and.  &
     &  (inplane_z .ge.tox_bot+to2_bot).and.(inplane_z.le.(tox_bot+to2_bot+tsc_h)) )THEN
        IF( plane_index .le.  source_len)                     dop = s_val !Source  doped node
        IF( plane_index .ge. (source_len+gate_len+2*spacer))  dop = d_val !Drain   doped node
        IF((plane_index .gt.  source_len) .and. &
           (plane_index .lt. (source_len+gate_len+2*spacer))) dop = c_val !Channel doped node
        
        if(whichkind_3D(ii).le.0) doping(map_3D(ii))=ELCH*dop

        write(99,*)(plane_index-1)*deltax*1.0d7,dop
        
     END IF
     
  END DO
  close(99)
  
  END SUBROUTINE indata_build_doping

 SUBROUTINE map_potential(potA,potB,which,map,nx,ny,nz,lwork)
 
 REAL(DP), INTENT(OUT) :: potA(1:nx,1:ny,1:nz)

 REAL(DP), INTENT(IN)  :: potB(0:lwork-1)
 INTEGER,  INTENT(IN)  :: which(0:nx*ny*nz-1)
 INTEGER,  INTENT(IN)  :: map(0:nx*ny*nz-1)
 INTEGER,  INTENT(IN)  :: nx
 INTEGER,  INTENT(IN)  :: ny
 INTEGER,  INTENT(IN)  :: nz
 INTEGER,  INTENT(IN)  :: lwork
 INTEGER               :: ii,x_index,y_index,z_index

 potA=0.0_dp
 
 DO ii=0, nx*ny*nz-1

 x_index=ii/(ny*nz) +1
 y_index=mod(mod(ii,ny*nz),ny) +1
 z_index=mod(ii,ny*nz)/ny +1
 
 IF(which(ii).ge.1)THEN
 potA(x_index,y_index,z_index)=potelectr
 ELSE
 potA(x_index,y_index,z_index)=potB(map(ii))
 END IF
 
!!!! IMPOSE NEUMANN conditions for the potential seen by GREEN
  pota(1:NX,1,1:NZ)  = pota(1:NX,2,1:NZ)
  pota(1:NX,NY,1:NZ) = pota(1:NX,NY-1,1:NZ)
  pota(1:NX,1:NY,1)  = pota(1:NX,1:NY,2)
  pota(1:NX,1:NY,NZ) = pota(1:NX,1:NY,NZ-1)

 END DO

 END SUBROUTINE map_potential

 SUBROUTINE inverse_map_charge(rhoA,rhoB,which,map,nx,ny,nz,lwork)
 
 REAL(DP), INTENT(OUT) :: rhoA(1:nx,1:ny,1:nz)

 REAL(DP), INTENT(IN)  :: rhoB(0:lwork-1)
 INTEGER,  INTENT(IN)  :: which(0:nx*ny*nz-1)
 INTEGER,  INTENT(IN)  :: map(0:nx*ny*nz-1)
 INTEGER,  INTENT(IN)  :: nx
 INTEGER,  INTENT(IN)  :: ny
 INTEGER,  INTENT(IN)  :: nz
 INTEGER,  INTENT(IN)  :: lwork
 INTEGER               :: ii,x_index,y_index,z_index

 rhoA=0.0_dp
 
 DO ii=0, nx*ny*nz-1

 x_index=ii/(ny*nz) +1
 y_index=mod(mod(ii,ny*nz),ny) +1
 z_index=mod(ii,ny*nz)/ny +1
 
 IF(which(ii).le.0)THEN
 rhoA(x_index,y_index,z_index)=rhoB(map(ii))
 END IF
 
 END DO

 END SUBROUTINE inverse_map_charge

 SUBROUTINE map_charge(rhoA,rhoB,which,map,nx,ny,nz,lwork)
 
 REAL(DP), INTENT(OUT)  :: rhoA(0:lwork-1)
 REAL(DP), INTENT(IN)   :: rhoB(1:nx,1:ny,1:nz)
 INTEGER,  INTENT(IN)   :: which(0:nx*ny*nz-1)
 INTEGER,  INTENT(IN)   :: map(0:nx*ny*nz-1) 
 INTEGER,  INTENT(IN)   :: nx
 INTEGER,  INTENT(IN)   :: ny
 INTEGER,  INTENT(IN)   :: nz
 INTEGER,  INTENT(IN)   :: lwork
 INTEGER                :: ii,x_index,y_index,z_index

 rhoA=0.0_dp

 DO ii=0, nx*ny*nz-1

 x_index=ii/(ny*nz) +1
 y_index=mod(mod(ii,ny*nz),ny) +1
 z_index=mod(ii,ny*nz)/ny +1
 
 IF(which(ii).le.0)THEN
    rhoA(map(ii))=rhoB(x_index,y_index,z_index)
 END IF
 
 END DO

 END SUBROUTINE map_charge


 SUBROUTINE shift_potential(potA,potB,shift,econd,which,map,list_ord,which_ord,&
      epsilon3D,nx,ny,nz,numel,lwork)

 REAL(DP), INTENT(OUT)  :: potA(0:lwork-1)
 REAL(DP), INTENT(IN)   :: potB(0:lwork-1)
 REAL(DP), INTENT(IN)   :: shift
 REAL(DP), INTENT(IN)   :: econd
 INTEGER,  INTENT(IN)   :: which(0:nx*ny*nz-1)
 INTEGER,  INTENT(IN)   :: map(0:nx*ny*nz-1)
 INTEGER,  INTENT(IN)   :: list_ord(1:8,0:numel-1)
 INTEGER,  INTENT(IN)   :: which_ord(0:nx*ny*nz-1)
 REAL(DP), INTENT(IN)   :: epsilon3D(0:numel-1)
 INTEGER,  INTENT(IN)   :: nx
 INTEGER,  INTENT(IN)   :: ny
 INTEGER,  INTENT(IN)   :: nz
 INTEGER,  INTENT(IN)   :: lwork
 INTEGER,  INTENT(IN)   :: numel
 INTEGER                :: ll, nn, ii
 INTEGER                :: checked(0:lwork-1)

 checked=0

 DO ii=0, nx*ny*nz-1
    IF(which(ii).le.0)THEN
       potA(map(ii))=potB(map(ii))+econd 
    endif
 enddo

 DO ll=0, numel-1
    
    IF((epsilon3D(ll).eq.DIEL_0*DIEL_OX).or.(epsilon3D(ll).eq.DIEL_0*DIEL_O2))THEN
       DO nn=1,8
          
          IF(which_ord(list_ord(nn,ll)).le.0)THEN
             
             IF(checked(list_ord(nn,ll)).ne.1)THEN
                potA(list_ord(nn,ll))=potA(list_ord(nn,ll))+shift
                checked(list_ord(nn,ll))=1
             END IF
             
          END IF
          
       END DO
    END IF
 END DO
 
END SUBROUTINE shift_potential

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gaussianlegendre(x1, x2, x, w, n) ! Gauss-Legendre integrals

  implicit none
  
  INTEGER,  intent(in)  :: n
  REAL(DP), intent(in)  :: x1, x2
  REAL(DP), intent(out) :: x(n), w(n)
  INTEGER               :: i, j, m
  REAL(DP)              :: p1, p2, p3, pp, xl, xm, z, z1
  REAL(DP), parameter   :: eps=3.0d-14
      
  m = (n+1)/2
  xm = 0.5_dp*(x2+x1)
  xl = 0.5_dp*(x2-x1)
  do i=1,m
    z = cos(3.141592654_dp*(i-0.25_dp)/(n+0.5_dp))
    z1 = 0.0
    do while(abs(z-z1) .gt. eps)
      p1 = 1.0_dp
      p2 = 0.0_dp
      do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
      end do
      pp = n*(z*p1-p2)/(z*z-1.0_dp)
      z1 = z
      z = z1 - p1/pp
    end do
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = (2.0_dp*xl)/((1.0_dp-z*z)*pp*pp)
    w(n+1-i) = w(i)
  end do

end subroutine gaussianlegendre

  
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
  integer                  :: ny,mi,mf
  real(dp)                 :: subband(1:(Mf-Mi+1))
  complex(dp)              :: A(1:NY,1:NY),Uii(1:NY,1:Mf-Mi+1)
  integer                  :: INFO
  integer, allocatable     :: iwork(:), supp(:)
  complex(dp), allocatable :: work(:)
  real(dp), allocatable    :: rwork(:)
  REAL(DP), EXTERNAL       :: DLAMCH


  allocate(WORK(20*ny))
  allocate(RWORK(24*ny))
  allocate(IWORK(10*ny))
  allocate(Supp(2*ny))
  
  call ZHEEVR('N','I','U',ny,A,ny,0.0,0.0,mi,mf,2*DLAMCH('S'),&
       Mf-Mi+1,subband,Uii,ny,SUPP,WORK,20*ny,RWORK,24*ny,IWORK,10*ny,INFO)
  
  deallocate(work)
  deallocate(rwork)
  deallocate(supp)
  deallocate(iwork)
  
  if (INFO.ne.0)then
     write(*,*)'SEVERE WARNING: SUB_DEF HAS FAILED. INFO=',INFO
     stop
  endif

END SUBROUTINE SUB_DEF_Z0

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function stringa(ii) result(rr)
  integer,       intent(in) :: ii
  character(:), allocatable :: rr
  character(range(ii)+2)    :: tt
  
  write(tt,'(i0)') ii
  rr = trim(tt)

end function stringa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer function ind_Kyz(v)
   implicit none

   INTEGER  :: iz,iy,l
   real(dp) :: v(2),vv(2)

   ind_kyz=0

   vv=v

   if(abs(V(1))>0.5_dp+1.0d-3) vv(1)=abs(abs(v(1))-1.0_dp)
   if(abs(V(2))>0.5_dp+1.0d-3) vv(1)=abs(abs(v(2))-1.0_dp)
   
   do iz=1,nkz
      do iy=1,nky
         l = iy + (iz-1)*nky
         if(  abs(abs(vv(1))-abs(k_vec(2,l))) < 1.0d-3 .and. &
              abs(abs(vv(2))-abs(k_vec(3,l))) < 1.0d-3 )then
            ind_kyz=l
            exit
         end if
      end do
   end do
   if( ind_kyz == 0 )then
      write(*,*)'pb w ind_kyz', ind_kyz,v
      stop
   end if
   
 end function ind_Kyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END MODULE indata

