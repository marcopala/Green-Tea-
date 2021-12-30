MODULE Pseudopot_so_gen

  USE static
  USE indata

  IMPLICIT NONE

  SAVE

!!! Variables  

CONTAINS


subroutine read_QE_output
implicit none

real(dp) :: a_1(3),a_2(3),a_3(3)
real(dp) :: b_1(3),b_2(3),b_3(3)
integer(k15), allocatable :: miller(:,:), ind_miller(:,:,:),N_projs(:), N_PW(:), ind_kG(:,:), ityp(:), miller_2D(:,:), ind_miller_2D(:,:),ind_cube(:)
real(dp), allocatable :: x_at(:), y_at(:), z_at(:), Gvec(:,:), R_at(:,:)
character(len=6), allocatable :: specie(:)
real(dp), allocatable :: k_vec(:,:),coeff(:),kx1(:),kx2(:)
character(len=80) :: comment
integer(k15) :: i,j,n,nnpw,npw,N_bands
integer :: ikx,l,k,m,nn,nnn,mmm,ivec(3),nrx,nrx0,nry,nrz,n1,n2,n3,m1,igt,jgt,ix,iy,iz,N_spin, nspin
integer(k15) :: N_at, N_typ_at, N_Gvec, N_k_pts, N_beta, max_mill_1,max_mill_2,max_mill_3, min_mill_1,min_mill_2,min_mill_3
real(dp) :: aux1,aux2,vec(3),t0,a0
real(dp) :: Ecutoff,ref
complex(dp) :: tmp
real(4) :: t1,t2
real(dp), allocatable :: E(:), KGt(:,:), Gx(:)
complex(dp), allocatable :: beta_fun(:,:),KS_fun(:,:),Psi_mod(:,:),Psi_mod_1(:,:),A(:,:),B(:,:),C(:,:),D(:,:),Q(:,:),Uh(:,:),U(:,:),betafunc(:,:)
complex(dp), allocatable :: Deeq(:,:,:),Deeq_so(:,:,:,:),HV(:,:),HLL(:,:),TLL(:,:),HLLL(:,:),TLLL(:,:),HLLLL(:,:),TLLLL(:,:)

complex(dp), allocatable :: Vloc(:,:),HCC(:,:)
complex(dp), allocatable :: dal(:),dbe(:),work(:),rwork(:),kval(:),id(:,:)
integer :: ncf,nm,nplus,iyz,ii,jj
real(dp), allocatable :: hkl(:,:)
integer :: is,ic,ip,jp,nt,na,ikb,info
integer(k15), allocatable :: i_h(:,:)
complex(dp) :: dummy(1,1)
character(256) :: nome

logical :: en_select
!LOGICAL :: skip_diag = .false.


if(.not. refine ) write(*,*)'WRITING HAMILTONIANS'
if( refine )  write(*,*)'READING HAMILTONIANS'


t0=hbareV*hbareV/2.0_dp/m0 !eV cm**2
a0=ac ! cm

call omp_set_num_threads(Nomp) !!! this set the environement variable OMP_NUM_THREADS = Nomp
write(*,*)
write(*,'(a,I3,a)')' This simulation will use',Nomp,' threads'
write(*,*)
open(unit=10,file=TRIM(input_file_DFT),status='unknown')
read(10,'(a)') comment
read(10,*) a_1(:), a_2(:), a_3(:)

!write(*,*)a_1(1)*bohr
!write(*,*)a_2(2)*bohr
!write(*,*)a_3(3)*bohr

read(10,'(a)') comment
write(*,*) comment
read(10,*) b_1(:), b_2(:), b_3(:)
write(*,*)
!write(*,*)dot_product(a_1,b_1)/(2.0_dp*pi),dot_product(a_1,b_2),dot_product(a_1,b_3)
!write(*,*)dot_product(a_2,b_1),dot_product(a_2,b_2)/(2.0_dp*pi),dot_product(a_2,b_3)
!write(*,*)dot_product(a_3,b_1),dot_product(a_3,b_2),dot_product(a_3,b_3)/(2.0_dp*pi)
!write(*,*)

b_1=b_1/(2.0_dp*pi/a0)/bohr !b1 in units of 2pi/a0
b_2=b_2/(2.0_dp*pi/a0)/bohr !b2 in units of 2pi/a0
b_3=b_3/(2.0_dp*pi/a0)/bohr !b3 in units of 2pi/a0
write(*,*) 'b1', b_1
write(*,*) 'b2', b_2
write(*,*) 'b3', b_3


read(10,'(a)') comment
read(10,*) N_typ_at
read(10,'(a)') comment
read(10,*) N_at
read(10,'(a)') comment
allocate(specie(N_at))
allocate(ityp(N_at))
allocate(x_at(N_at),y_at(N_at),z_at(N_at))
allocate(R_at(3,N_at))
do i=1,N_at
   read(10,*) ityp(i)
   read(10,'(a,3e25.15)') specie(i), R_at(:,i)
end do

read(10,'(a)') comment
read(10,*) N_Gvec
read(10,'(a)') comment
write(*,*) comment

allocate(miller(3,N_Gvec))
READ(10,'(3i8)') (miller(:,i), i=1,N_Gvec)



allocate(ind_miller(minval(miller(1,:)):maxval(miller(1,:)),&
     minval(miller(2,:)):maxval(miller(2,:)),&
     minval(miller(3,:)):maxval(miller(3,:))))
ind_miller=0
do i=1,N_Gvec
   ind_miller(miller(1,i),miller(2,i),miller(3,i))=i
end do


allocate(Gvec(3,N_Gvec))
do i=1,N_Gvec
   Gvec(1:3,i)=miller(1,i)*b_1(1:3)+miller(2,i)*b_2(1:3)+miller(3,i)*b_3(1:3)
end do
Gvec=dble(nint(Gvec)) !!! Gvec in units of 2pi/a0

write(*,*)'N_Gvec =',N_Gvec
read(10,'(a)') comment
write(*,*) comment
read(10,*) N_spin
write(*,*)'NSPIN =',N_spin

IF ( n_spin <= 2 ) THEN
   npol = 1
   domag=.FALSE.
   lspinorb=.FALSE.
ELSE
   READ(10,'(a)') comment
   READ(10,*) domag, lspinorb
   npol = 2
ENDIF
write(*,*) 'domag =',domag
write(*,*) 'lspinorb =',lspinorb
write(*,*) 'NPOL =',npol

allocate(Vloc(N_Gvec,N_spin))
write(*,*) 'shape vloc',shape(vloc)
   read(10,'(a)') comment
   write(*,*) comment
do is=1,N_spin
   read(10,'(a)') comment
   write(*,*) comment
   do i=1,N_Gvec     
      read(10,'(2e25.15)') Vloc(i,is)!aux1, aux2
      !Vloc(i,is)=dcmplx(aux1,aux2)
   end do
end do

read(10,'(a)') comment
read(10,*) okvan
if (okvan) then
write(*,*)'Error, US PP not implemented. The simulation will be stopped!'
stop
end if


read(10,'(a)') comment
allocate(N_projs(N_typ_at))
read(10,*) N_projs(1:N_typ_at)
write(*,*)'N_projs =',N_projs

nn=maxval(N_projs(:))

IF ( n_spin /= 4 ) THEN
   allocate(hkl(nn,nn))
   allocate(Deeq(nn,nn,N_typ_at))
ELSE
   allocate(hkl(nn,nn))
   allocate(Deeq_so(nn,nn,n_spin,N_typ_at))
END IF


read(10,'(a)') comment
do j=1,N_typ_at
   read(10,'(a)') comment
   write(*,*)
   read(10,*) i,N_projs(j)
   IF ( n_spin /= 4 ) THEN
      read(10,*)Deeq(1:N_projs(j),1:N_projs(j),j)
   ELSE
      do is=1,N_spin
         read(10,*)comment
         read(10,*)Deeq_so(1:N_projs(j),1:N_projs(j),is,j)
      end do
   END IF
end do

deallocate(hkl)
read(10,'(a)') comment
read(10,*) N_beta

read(10,'(a)') comment
write(*,*) comment
read(10,*) N_k_pts
write(*,*) N_k_pts
write(*,*) 'N_k_pts =',N_k_pts

allocate(N_PW(N_k_pts))

read(10,'(a)') comment
write(*,*) comment
read(10,*) N_PW(1:N_k_pts)

nkx=n_k_pts/Nkyz
if(mod(n_k_pts,Nkyz) /=0)then
   write(*,*)'Fatal error in the number of k points!', n_k_pts, Nkyz
   stop
end if

if(N_SPIN == 2) then
   if(mod(nkx,2) /=0)then
      write(*,*)'error in the number of k-points'
      stop
   else
      nkx = nkx/N_SPIN      
   end if
end if

write(*,*) 'NK points =',N_k_pts
write(*,*) 'NKx points =',Nkx
write(*,*) 'NKyz points =',Nkyz

allocate(k_vec(3,N_k_pts))
allocate(ind_kG(N_Gvec,N_k_pts))
ind_kG=0

!stop
Ecutoff=ryd*Ecut !!!!/3.0_dp

nrx=0
do n2=-1000,1000!minval(miller(1,:)),maxval(miller(1,:))
   vec(1:3)=n2*b_1(1:3)
   if(t0*vec(1)**2*(2.0_dp*pi/a0)**2 < 4.0_dp*Ecutoff )nrx=nrx+1
end do
nrx=(nrx-1)
write(*,*)'NRX =',nrx
!nrx=maxval(miller(1,:))-minval(miller(1,:))+1
nry=0
do n2=minval(miller(2,:)),maxval(miller(2,:))
   vec(2:3)=n2*b_2(2:3)
   if(t0*vec(2)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff )nry=nry+1
end do
!nry=nry-1
write(*,*)'NRY =',nry
nrz=0
do n3=minval(miller(3,:)),maxval(miller(3,:))
   vec(2:3)=n3*b_3(2:3)
   if(t0*vec(3)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff )nrz=nrz+1
end do
!nrz=nrz-1
write(*,*)'NRZ =',nrz

ndx=nrx
ndy=nry
ndz=nrz

ncf=min((nrx/2-2),8)
allocate(coeff(0:ncf))

call coefficienti (2*ncf,coeff)
write(*,*)'Discretization order =',ncf
do i=0,ncf
   write(*,*)i,coeff(i)
end do
!stop

Ngt=0
do n2=-nry/2,nry/2!-maxval(miller(2,:)),maxval(miller(2,:))
do n3=-nrz/2,nrz/2!maxval(miller(3,:)),maxval(miller(3,:))
vec(2:3)=n2*b_2(2:3)+n3*b_3(2:3)
if(t0*dot_product(vec(2:3),vec(2:3))*(2.0_dp*pi/a0)**2 <= Ecutoff) then ! eV
!if(t0*vec(2)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff .and. t0*vec(3)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff)then
   Ngt=Ngt+1
end if

end do
end do

max_mill_1=maxval(miller(1,:))
max_mill_2=maxval(miller(2,:))
max_mill_3=maxval(miller(3,:))
min_mill_1=minval(miller(1,:))
min_mill_2=minval(miller(2,:))
min_mill_3=minval(miller(3,:))

write(*,*)'Ngt =',Ngt

allocate(miller_2D(3,Ngt))
miller_2D=0

j=0
do n2=-nry/2,nry/2!-maxval(miller(2,:)),maxval(miller(2,:))
do n3=-nrz/2,nrz/2!-maxval(miller(3,:)),maxval(miller(3,:))
vec(2:3)=n2*b_2(2:3)+n3*b_3(2:3) !!!! THIS IMPLIES THAT b_1 is orthonal to the transverse plane
if(t0*dot_product(vec(2:3),vec(2:3))*(2.0_dp*pi/a0)**2 <= Ecutoff) then ! eV
!if(t0*vec(2)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff .and. t0*vec(3)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff)then 
   j=j+1
   miller_2D(2,j)=n2
   miller_2D(3,j)=n3
end if

end do
end do
allocate(ind_miller_2D(-nry/2:nry/2,-nrz/2:nrz/2))
ind_miller_2D=0

do i=1,NGt
   ind_miller_2D(miller_2D(2,i),miller_2D(3,i))=i
end do

!ncell=1!2!4!1!3
nrx0=nrx!/ncell
if(mod(nrx,nrx0)/=0)then
   write(*,*)'nrx0problem',nrx,nrx0,ecut
   stop
end if
write(*,*)'ncell =',ncell
write(*,*)'max miller_2D 2',maxval(miller_2D(2,:))
write(*,*)'max miller_2D 3',maxval(miller_2D(3,:))

allocate(Gx(nrx))
Gx=0.0_dp
do i=-nrx/2,nrx/2-1!maxval(miller(1,:)),maxval(miller(1,:))!-1
   m1=i
   if ( m1 < 0 ) then
     m1 = m1 + nrx + 1
   else
      m1 = m1 + 1
   end if
   Gx(m1)=dble(i)*b_1(1)
end do

Dx=a_1(1)*bohr/dble(nrx)
Dy=a_2(2)*bohr/dble(nry)
Dz=a_3(3)*bohr/dble(nrz)

write(*,*)'a0/Dx =',a0/Dx,Dx
write(*,*)'a0/Dy =',a0/Dy,Dy
write(*,*)'a0/Dz =',a0/Dz,Dz

    nx=nrx-1
    ny=nry-1
    nz=nrz-1

    allocate(kc_min(Nkyz))
    allocate(kv_max(Nkyz))

    IF ( N_SPIN == 2 ) then
       write(*,*)'MAGNETIZATION ON'
       nspin=2
    else
       nspin=1
    END IF
    write(*,*)
    nnpw=maxval(N_PW(:))
    write(*,*)'nnpw =',nnpw
    write(*,*)
    
do is=1,nspin ! LOOP OVER THE SPIN INDEX IN CASE NSPIN==2 (collinear magnetic system)   
do iyz=1,Nkyz  ! Loops over the transverse k vectors
    
do ikx=1,Nkx
   write(*,*)'ikx =',ikx,'is =',is
   read(10,'(a)') comment
   if(NSPIN==2) then
      read(10,*) i
      write(*,*) 'Collinear magnetic system, current_spin = ', i
      if (i < 1 .or. i > 2 ) then
         write(*,*)'Error: mismatch in spin index'
         stop
      end if
   end if
   read(10,*) k_vec(:,ikx+(iyz-1)*nkx)
   k_vec(:,ikx+(iyz-1)*nkx)=k_vec(:,ikx+(iyz-1)*nkx)/bohr/(2.0_dp*pi/a0)
   write(*,*) k_vec(:,ikx+(iyz-1)*nkx)

   read(10,'(a)') comment
   read(10,*) N_PW(ikx+(iyz-1)*nkx)
   npw=N_PW(ikx+(iyz-1)*nkx)
   
   read(10,'(a)') comment
   
   do i=1,npw
      read(10,*)ind_kG(i,ikx+(iyz-1)*nkx)
   end do
   if(.not.(allocated(beta_fun)))then
      allocate(beta_fun(nnpw,N_beta*nkx))
   end if
   do l=1,N_beta
      read(10,'(a)') comment
      READ(10,'(2e25.15)') beta_fun(1:npw,l+(ikx-1)*N_beta)
   end do

   read(10,'(a)') comment
   read(10,*) N_bands
   write(*,'(a, i4, a)')' Using',N_bands,' bands from the DFT simulation'


 ! KS EIGEN-FUNCTIONS
   if(.not.(allocated(KS_fun)))then
      allocate(KS_fun(nnpw*npol,N_bands*nkx))
      KS_FUN=0.0_dp
   end if
   do j=1,N_bands
      read(10,'(a)') comment
      read(10,'(a)') comment
      IF ( n_spin /= 4 ) THEN
         read(10,'(2e25.15)') KS_fun(1:npw,j+(ikx-1)*N_bands)
      ELSE
         read(10,'(2e25.15)') KS_fun(1:npw,j+(ikx-1)*N_bands)
         read(10,'(2e25.15)') KS_fun(1+nnpw:npw+nnpw,j+(ikx-1)*N_bands)
      END IF
   end do

end do


allocate(Psi_mod(nkx*nrx*Ngt*npol,N_bands*nkx))
Psi_mod=0.0_dp
allocate(A(nkx*nrx*Ngt*npol,N_bands*nkx))
A=0.0_dp
do ikx=1,nkx
   npw=N_PW(ikx+(iyz-1)*nkx)
   
   allocate(B(Nrx*Ngt,N_bands))
   B=0.0_dp
   do jp=1,npol
      do j=1,npw
         if(abs(miller(1,ind_kG(j,ikx+(iyz-1)*nkx)))<nrx/2 .and.  &
              abs(miller(2,ind_kG(j,ikx+(iyz-1)*nkx)))<nry/2 .and. &
              abs(miller(3,ind_kG(j,ikx+(iyz-1)*nkx)))<nrz/2)then
            m1=miller(1,ind_kG(j,ikx+(iyz-1)*nkx))
            if ( m1 < 0 ) then
               m1 = m1 + nrx + 1
            else
               m1 = m1 + 1
            end if
            jgt=ind_miller_2D(miller(2,ind_kG(j,ikx+(iyz-1)*nkx)),miller(3,ind_kG(j,ikx+(iyz-1)*nkx)))
            if(jgt.ne.0)then
               B(m1+(jgt-1)*nrx,1:N_bands)=KS_fun(j+(jp-1)*nnpw,1+(ikx-1)*N_bands:ikx*N_bands)
            end if
         end if
      end do
      do jgt=1,Ngt
         do m1=1,Nrx
            A(m1+(ikx-1)*nrx+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt,1+(ikx-1)*N_bands:ikx*N_bands)=B(m1+(jgt-1)*nrx,1:N_bands)
         end do
      end do
   end do
   deallocate(B)
end do


allocate(Uh(Nkx*nrx,Nkx*nrx))
do ikx=1,nkx
   
   do m1=1,nrx
      do n=1,Nkx*nrx
         Uh(m1+(ikx-1)*nrx,n)=&
              exp(dcmplx(0.0_dp,1.0_dp)*(k_vec(1,ikx+(iyz-1)*nkx)+Gx(m1))*(2.0_dp*pi)*dble(n-1)/dble(nrx))/sqrt(dble(Nkx*nrx))
      end do
   end do
end do
allocate(B(nkx*nrx,nkx*N_bands))
do jp=1,npol
   do jgt=1,Ngt
      call zgemm ('c','n',nkx*nrx,nkx*N_bands,nkx*nrx,alpha,Uh,nkx*nrx,&
           A(1+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt:jgt*nkx*nrx+(jp-1)*nkx*nrx*ngt,1:nkx*N_Bands),nkx*nrx,beta,B(1:nkx*nrx,1:nkx*n_bands),nkx*nrx)
      do n=1,nkx*nrx
         Psi_mod(n+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt,1:Nkx*N_bands)=B(n,1:Nkx*N_bands)
      end do
   end do
end do
deallocate(A,B)

M=N_bands
NM=nkx*N_bands
NMODES=NM
write(*,*)'Total number of bands from DFT simulation =',NM

if(.not. refine)then
write(*,*)
write(*,'(a,I3,a)')' WARNING: by using Nomp = ',Nomp,' the required memory is approximatively '
write(*,'(F8.1,a)')     dble(nrx*Ngt)*dble(N_beta*nkx)*16.0d-9+& !betafunc
     dble(Nomp*nkx*nrx*npol)*dble(nkx*nrx*Ngt*npol)*16.0d-9+&  ! HCC
     dble(Nomp*2*Nrx*npol)*dble(nrx*Ngt*npol)*16.0d-9+&  ! HLL,TLL
     dble(Nomp*nrx)*dble(nrx*ngt)*16.0d-9+&  ! Q
     dble(Nomp*2*nkx*nrx)*dble(nkx*nrx)*16.0d-9+&  ! B, U
     dble(2*nm)*dble(nrx*Ngt*npol)*16.0d-9, ' Gb'
write(*,*)
write(*,'(a)')' PLEASE CHECK THAT YOUR SYSTEM CAN MANAGE THIS SIMULATION!'
write(*,'(a)')' Also note that, by decreasing Nomp, the allocated memory can be reduced up to '
write(*,'(F8.1,a)')    dble(nrx*Ngt)*dble(N_beta*nkx)*16.0d-9+& !betafunc
     dble(nkx*nrx*npol)*dble(nkx*nrx*Ngt*npol)*16.0d-9+&  ! HCC
     dble(2*Nrx*npol)*dble(nrx*Ngt*npol)*16.0d-9+&  ! HLL,TLL
     dble(nrx)*dble(nrx*ngt)*16.0d-9+&  ! Q
     dble(2*nkx*nrx)*dble(nkx*nrx)*16.0d-9+&  ! B, U
     dble(2*nm)*dble(nrx*Ngt*npol)*16.0d-9, ' Gb'
end if

allocate(C(Nrx*Ngt*npol,NM))
C=0.0_dp
do ikx=1,nkx
   do jp=1,npol
      do jgt=1,Ngt
         C(1+(jgt-1)*nrx+(jp-1)*nrx*ngt:jgt*nrx+(jp-1)*nrx*ngt,1+(ikx-1)*M:ikx*M)=&
              Psi_mod(1+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt:nrx+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt,1+(ikx-1)*N_bands:ikx*N_bands)
      end do
   end do
end do

!!!!!
if(.not.refine)then
if(allocated(KGt))deallocate(KGt)
allocate(KGt(4,1*Ngt))
KGt=0.0_dp
do j=1,Ngt
   kgt(2,j)=miller_2D(2,j)*b_2(2)+miller_2D(3,j)*b_3(2) + k_vec(2,1+(iyz-1)*nkx)
   kgt(3,j)=miller_2D(2,j)*b_2(3)+miller_2D(3,j)*b_3(3) + k_vec(3,1+(iyz-1)*nkx)
   kgt(4,j)=(kgt(2,j)**2+kgt(3,j)**2)
end do
allocate(U(NGt*npol,(nry)*(nrz)))
do iy=1,NRY
   do iz=1,NRZ
      j=iy+(iz-1)*(NRY)
      do jp=1,npol
      U(1+(jp-1)*NGt:jp*NGt,j)=exp( cmplx(0.0_dp,-1.0_dp)*KGt(2,1:NGt)*2.0_dp*pi/a0*dble(iy)*Dy+&
           cmplx(0.0_dp,-1.0_dp)*KGt(3,1:NGt)*2.0_dp*pi/a0*dble(iz)*Dz )/sqrt(dble((nry)*(nrz)))
   end do
   end do
end do
allocate(A(NM,(nry)*(nrz)))
allocate(B(NM,ngt*npol))

do ix=1,nrx
   open(unit=2000+ix,file='Nyz_nkyz_'//TRIM(STRINGA(iyz))//'_x_'//TRIM(STRINGA(ix))//'.dat',status='unknown')
   do ip=1,npol
      do jgt=1,Ngt
         do i=1,NM
            B(i,jgt+(ip-1)*ngt)=conjg(C(ix+(jgt-1)*NRX+(ip-1)*nrx*ngt,i))
         end do
      end do
   end do
   
!!!!   B(1:NM,1:ngt)=transpose(conjg(C(1+(ix-1)*NGt:ix*NGt,1:NM)))
   call ZGEMM('n','n',NM,(nry)*(nrz),NGt*npol,alpha,B,NM,U,NGt*npol,beta,A,NM)
   do iy=1,nry
      do iz=1,nrz
         tmp=0.0d0
         do i=1,nband_v!1,NM
            tmp=tmp+(A(i,iy+(iz-1)*(nry))*dconjg( A(i,iy+(iz-1)*(nry)) ))
         end do
         write(2000+ix,*)(iy-1)*dy*1e8,(iz-1)*dz*1e8,dble(tmp)
               end do
      write(2000+ix,*)
   end do
   close(2000+ix)
end do
deallocate(A,B,U)

end if
!!

deallocate(Psi_mod)
allocate(PSI_MOD(Nrx*Ngt*npol,NM))
Psi_mod=0.0_dp

if(.not. refine)then
write(*,*)
t1=SECNDS(0.0)
write(*,*)'Starting the orthonormalization...'
ic=1
call MGS(nrx0*Ngt*npol,NM,C(1:nrx0*ngt*npol,1:NM))
do jp=1,npol
   do jgt=1,Ngt
      do ix=1,nrx0
         PSI_MOD(ix+(ic-1)*Nrx0+(jgt-1)*nrx+(jp-1)*nrx*ngt,1+(ic-1)*NM:ic*NM)=C(ix+(jgt-1)*Nrx0+(jp-1)*nrx0*ngt,1:NM)      
      end do  
   end do
end do
do j=1,NM
   do i=1,Nrx*Ngt*npol
      if(C(i,j)/=C(i,j))then
         write(*,*)'prob with the MGS',j
         write(*,*)'stopping the simulaton'
         stop
         exit
      end if
   end do
end do
t2=SECNDS(t1)
write(*,'(a,F7.3,a)')' Orthonormalization done in ',t2,' s'
write(*,*)
end if

deallocate(C)
deallocate(ks_fun)


if(ncell==1)then
   ic=1
if(nspin==1) then
   open(unit=13,file=TRIM(outdir)//'PPsi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
else
if(is==1)   open(unit=13,file=TRIM(outdir)//'PPsi_Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
if(is==2)   open(unit=13,file=TRIM(outdir)//'PPsi_Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
end if
if(.not. refine)then
   write(*,*)
   write(*,*)'Writing PPSI files...'
   do j=1,NM
      do i=1,Nrx*Ngt*npol
         write(13,*)PSI_MOD(i,j)
      end do
   end do
   write(*,*)'done'
else if(refine)then
   write(*,*)
   write(*,*)'reading PPSI files...'
   do j=1,NM
      do i=1,Nrx*Ngt*npol
         read(13,*)PSI_MOD(i,j)
      end do
   end do
   write(*,*)'done'
end if
close(13)
end if
write(*,*)

if(ncell==2)then
   ic=1 
   if(nspin==1) then
      nome=trim(indir2)//'PPsi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat'
   else
      if(is==1)nome=trim(indir2)//'PPsi_Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat'
      if(is==2)nome=trim(indir2)//'PPsi_Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat'
   end if
   write(*,*)nome
   allocate(PSI_MOD_1(Nrx*Ngt*npol,NM2))
   open(unit=13,file=nome,status='unknown')
      do j=1,NM2
      do i=1,Nrx*Ngt*npol
         read(13,*)PSI_MOD_1(i,j)
      end do
   end do
   close(13)
end if


if(.not. refine)then
   
allocate(ind_cube(nrx*nkx*ngt))
ind_cube=0
do ikx=1,nkx
   npw=N_PW(ikx+(iyz-1)*nkx)
   do j=1,npw
      if(abs(miller(1,ind_kG(j,ikx+(iyz-1)*nkx)))<nrx/2 .and. &
           abs(miller(2,ind_kG(j,ikx+(iyz-1)*nkx)))<nry/2 .and. &
           abs(miller(3,ind_kG(j,ikx+(iyz-1)*nkx)))<nrz/2)then
         m1=miller(1,ind_kG(j,ikx+(iyz-1)*nkx))
         if ( m1 < 0 ) then
            m1 = m1 + nrx + 1
         else
            m1 = m1 + 1
         end if
         jgt=ind_miller_2D(miller(2,ind_kG(j,ikx+(iyz-1)*nkx)),miller(3,ind_kG(j,ikx+(iyz-1)*nkx)))
         if(jgt/=0) ind_cube(m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=j
      end if
   end do
end do

allocate(betafunc(nrx*Ngt,N_beta*nkx))
betafunc=0.0_dp
do igt=1,Ngt
   do ikx=1,nkx
      do n1=1,nrx
         i=ind_cube(n1+(ikx-1)*nrx+(igt-1)*nrx*nkx)
         if(i/=0)then
            betafunc(n1+(igt-1)*nrx,1+(ikx-1)*N_beta:ikx*N_beta)=beta_fun(i,1+(ikx-1)*N_beta:ikx*N_beta)
         end if
      end do
   end do
end do
deallocate(beta_fun)

allocate(A(nrx*Ngt*npol,NM))
allocate(C(nrx*Ngt*npol,NM))
A=0.0_dp
C=0.0_dp

write(*,*)
write(*,*)'Transforming H in the reduced basis (this can be long) ...'
write(*,*)
t1=SECNDS(0.0)
!$omp parallel default(none) private(igt,jgt,ikx,ikb,nt,na,m1,n1,n,nn,l,ivec,&
!$omp i,j,ip,jp,tmp,vec,B,U,HLL,TLL,HCC,Q,D)&
!$omp shared(nomp,is,ind_kg,nkx,nrx,iyz,npol,nm,ngt,ncf,N_projs,ityp,N_typ_at,N_spin,&
!$omp N_at,N_beta,ind_cube,ind_miller,miller,miller_2D,min_mill_1,max_mill_1,b_2,b_3,&
!$omp a0,t0,Dx,Dy,Dz,coeff,k_vec,Uh,betafunc,Deeq,Deeq_so,PSI_MOD,A,C,vloc)

allocate(HCC(nkx*nrx*npol,nkx*nrx*Ngt*npol))
allocate(B(nkx*nrx,nkx*nrx),U(nkx*nrx,nkx*nrx))
allocate(Q(nrx,nrx*ngt))
allocate(HLL(Nrx*npol,nrx*Ngt*npol))
allocate(TLL(Nrx*npol,nrx*Ngt*npol))

!$omp do
do igt=1,Ngt
HCC=0.0_dp

do ikx=1,nkx

IF (N_SPIN /= 4)THEN
   ikb=0
   DO nt=1,N_typ_at
      allocate(D(nrx,N_projs(nt)))
      DO na=1,n_at
         IF ( nt == ityp(na) ) THEN
            call zgemm('n','n',nrx,N_projs(nt),N_projs(nt),alpha,&
                 betafunc(1+(igt-1)*nrx:igt*nrx,ikb+1+(ikx-1)*N_beta:ikb+N_projs(nt)+(ikx-1)*N_beta),&
                 nrx,Deeq(1:N_projs(nt),1:N_projs(nt),nt),N_projs(nt),beta,D,nrx)
            call zgemm('n','c',nrx,nrx*ngt,N_projs(nt),alpha,D(1:nrx,1:N_projs(nt)),nrx,&
                 betafunc(1:nrx*ngt,ikb+1+(ikx-1)*N_beta:ikb+N_projs(nt)+(ikx-1)*N_beta),&
                 nrx*ngt,beta,Q(1:nrx,1:nrx*ngt),nrx)
            !Hkk(1:nrx,1:nrx*ngt,ikx) = Hkk(1:nrx,1:nrx*ngt,ikx)+ryd*(Q(1:nrx,1:nrx*ngt))
            do jgt=1,ngt
               HCC(1+(ikx-1)*nrx:ikx*nrx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=&
               HCC(1+(ikx-1)*nrx:ikx*nrx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Q(1:nrx,1+(jgt-1)*nrx:jgt*nrx))
            end do
            ikb=ikb+N_projs(nt)
         end IF
      end do
      deallocate(D)
   end do
ELSE
do l=1,4
      ikb=0
   DO nt=1,N_typ_at
      allocate(D(nrx,N_projs(nt)))
      DO na=1,n_at
         IF ( nt == ityp(na) ) THEN
            call zgemm('n','n',nrx,N_projs(nt),N_projs(nt),alpha,&
                 betafunc(1+(igt-1)*nrx:igt*nrx,ikb+1+(ikx-1)*N_beta:ikb+N_projs(nt)+(ikx-1)*N_beta),&
                 nrx,Deeq_so(1:N_projs(nt),1:N_projs(nt),l,nt),N_projs(nt),beta,D,nrx)
            call zgemm('n','c',nrx,nrx*ngt,N_projs(nt),alpha,D(1:nrx,1:N_projs(nt)),nrx,&
                 betafunc(1:nrx*ngt,ikb+1+(ikx-1)*N_beta:ikb+N_projs(nt)+(ikx-1)*N_beta),&
                 nrx*ngt,beta,Q(1:nrx,1:nrx*ngt),nrx)
            ! Hkk(1:nrx,1:nrx*ngt,ikx+(l-1)*nkx) = Hkk(1:nrx,1:nrx*ngt,ikx+(l-1)*nkx)+ryd*(Q(1:nrx,1:nrx*ngt))
            do jgt=1,ngt
               if(l==1)   HCC(1+(ikx-1)*nrx:ikx*nrx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=&
                    HCC(1+(ikx-1)*nrx:ikx*nrx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Q(1:nrx,1+(jgt-1)*nrx:jgt*nrx))
               if(l==2)   HCC(1+(ikx-1)*nrx:ikx*nrx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
                    HCC(1+(ikx-1)*nrx:ikx*nrx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Q(1:nrx,1+(jgt-1)*nrx:jgt*nrx))
               if(l==3)   HCC(1+(ikx-1)*nrx+nrx*nkx:ikx*nrx+nrx*nkx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=&
                    HCC(1+(ikx-1)*nrx+nrx*nkx:ikx*nrx+nrx*nkx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Q(1:nrx,1+(jgt-1)*nrx:jgt*nrx))
               if(l==4)   HCC(1+(ikx-1)*nrx+nrx*nkx:ikx*nrx+nrx*nkx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
                    HCC(1+(ikx-1)*nrx+nrx*nkx:ikx*nrx+nrx*nkx,1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt:nrx+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Q(1:nrx,1+(jgt-1)*nrx:jgt*nrx))
            end do
            ikb=ikb+N_projs(nt)
         end IF
      end do
      deallocate(D)
   end do
end do
END IF

do m1=1,nrx
   do jgt=1,Ngt
      j=ind_cube(m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)
      if(j/=0)then
         do n1=1,nrx
            i=ind_cube(n1+(ikx-1)*nrx+(igt-1)*nrx*nkx)
            if(i/=0)then
               ivec=miller(:,ind_kG(i,ikx+(iyz-1)*nkx))-miller(:,ind_kG(j,ikx+(iyz-1)*nkx))
               n=ind_miller(ivec(1),ivec(2),ivec(3))

IF (N_SPIN /=4) THEN
  if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*Vloc(n,is)
ELSE
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Vloc(n,1)+Vloc(n,4))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
        Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Vloc(n,2)-im*Vloc(n,3))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=&
        Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Vloc(n,2)+im*Vloc(n,3))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
        Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Vloc(n,1)-Vloc(n,4))
END IF
if(ivec(1) > 0 )then
   if(ivec(1)-nrx>=min_mill_1 .and. ivec(1)-nrx<=max_mill_1 )then
      n=ind_miller(ivec(1)-nrx,ivec(2),ivec(3))
IF (N_SPIN /=4) THEN
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*Vloc(n,is)
ELSE
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Vloc(n,1)+Vloc(n,4))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
        Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Vloc(n,2)-im*Vloc(n,3))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=&
        Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Vloc(n,2)+im*Vloc(n,3))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
        Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Vloc(n,1)-Vloc(n,4))
END IF
end if
end if
if(ivec(1) < 0 )then
   if(ivec(1)+nrx>=min_mill_1 .and. ivec(1)+nrx<=max_mill_1 ) then
      n=ind_miller(ivec(1)+nrx,ivec(2),ivec(3))
IF (N_SPIN /=4) THEN
   if( n /= 0 )  Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*Vloc(n,is)
ELSE
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Vloc(n,1)+Vloc(n,4))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
        Hcc(n1+(ikx-1)*nrx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Vloc(n,2)-im*Vloc(n,3))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)=&
        Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx)+ryd*(Vloc(n,2)+im*Vloc(n,3))
   if( n /= 0 ) Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)=&
        Hcc(n1+(ikx-1)*nrx+nrx*nkx,m1+(ikx-1)*nrx+(jgt-1)*nrx*nkx+nrx*nkx*ngt)+ryd*(Vloc(n,1)-Vloc(n,4))
END IF
end if
end if
end if
end do
end if
end do
end do


end do

HLL=0.0_dp
TLL=0.0_dp

   
do ip=1,npol
   do jp=1,npol
      do jgt=1,Ngt
         U(1:nkx*nrx,1:nkx*nrx)=HCC(1+(ip-1)*nkx*nrx:nkx*nrx+(ip-1)*nkx*nrx,1+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt:jgt*nkx*nrx+(jp-1)*nkx*nrx*ngt)
         call zgemm ('n','n',nkx*nrx,nkx*nrx,nkx*nrx,alpha,U,nkx*nrx,Uh,nkx*nrx,beta,B,nkx*nrx)
         call zgemm ('c','n',nkx*nrx,nkx*nrx,nkx*nrx,alpha,Uh,nkx*nrx,B,nkx*nrx,beta,U,nkx*nrx)
         HCC(1+(ip-1)*nkx*nrx:nkx*nrx+(ip-1)*nkx*nrx,1+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt:jgt*nkx*nrx+(jp-1)*nkx*nrx*ngt)=U(1:nkx*nrx,1:nkx*nrx)
      end do
   end do
end do


if(igt <= (Ngt/nomp) .and. mod(igt,10)==0) write(*,*)floor(dble(igt-1)/dble(Ngt/nomp)*100), '% done'

do ip=1,npol
do jp=1,npol
   do jgt=1,Ngt
      do n=1,Nrx
         do nn=1,Nrx
            HLL(n+(ip-1)*nrx,nn+(jgt-1)*nrx+(jp-1)*nrx*ngt)=HCC(n+(ip-1)*nkx*nrx,nn    +(jgt-1)*Nkx*nrx+(jp-1)*nkx*nrx*ngt)
            TLL(n+(ip-1)*nrx,nn+(jgt-1)*nrx+(jp-1)*nrx*ngt)=HCC(n+(ip-1)*nkx*nrx,nn+nrx+(jgt-1)*Nkx*nrx+(jp-1)*nkx*nrx*ngt)
         end do
      end do
   end do
end do

if(nkx==2)then
   do jp=1,npol
      do jgt=1,Ngt
         do n=1,Nrx
            do nn=n+1,Nrx
               TLL(n+(ip-1)*nrx,nn+(jgt-1)*nrx+(jp-1)*nrx*ngt)=0.0_dp
            end do
            do nn=n,n
               TLL(n+(ip-1)*nrx,nn+(jgt-1)*nrx+(jp-1)*nrx*ngt)=&
               TLL(n+(ip-1)*nrx,nn+(jgt-1)*nrx+(jp-1)*nrx*ngt)/2.0_dp
            end do
         end do
      end do
   end do
end if

end do


   vec(2:3)= ( miller_2D(2,igt)*b_2(2:3) + miller_2D(3,igt)*b_3(2:3) + k_vec(2:3,1+(iyz-1)*nkx) )*(2.0_dp*pi/a0) !!! WARNING! We suppose that b_1 is along x
   tmp=2.0_dp*t0/(Dx**2)*coeff(0)
   do j=0,ncf
      tmp=tmp+2.0_dp*t0/(Dy**2)*coeff(j)*cos(vec(2)*dble(j)*Dy)&
           +  2.0_dp*t0/(Dz**2)*coeff(j)*cos(vec(3)*dble(j)*Dz)
   end do

do ip=1,npol
   do jp=1,npol
      if(ip==jp)then
         do n=1,Nrx
            HLL(n+(ip-1)*nrx,n+(igt-1)*nrx+(jp-1)*nrx*ngt)=HLL(n+(ip-1)*nrx,n+(igt-1)*nrx+(jp-1)*nrx*ngt)+tmp
         end do
      
         do j=1,ncf
            do n=1,Nrx-j
               HLL(n+(ip-1)*nrx,n+j+(igt-1)*Nrx+(jp-1)*nrx*ngt)=&
                    HLL(n+(ip-1)*nrx,n+j+(igt-1)*Nrx+(jp-1)*nrx*ngt)+coeff(j)*t0/(Dx**2)
               HLL(n+j+(ip-1)*nrx,n+(igt-1)*Nrx+(jp-1)*nrx*ngt)=&
                    HLL(n+j+(ip-1)*nrx,n+(igt-1)*Nrx+(jp-1)*nrx*ngt)+coeff(j)*t0/(Dx**2)
            end do
         end do
         
         do j=1,ncf
            do n=1,j
               TLL(nrx-n+1+(ip-1)*nrx,j-n+1+(igt-1)*Nrx+(jp-1)*nrx*ngt)=&
                    TLL(nrx-n+1+(ip-1)*nrx,j-n+1+(igt-1)*Nrx+(jp-1)*nrx*ngt)+coeff(j)*t0/(Dx**2)
            end do
         end do
      end if
   end do

end do


do ip=1,npol

   call ZGEMM('n','n',nrx,NM,nrx*Ngt*npol,alpha,HLL(1+(ip-1)*nrx:nrx+(ip-1)*nrx,:),nrx,&
        PSI_MOD,nrx*Ngt*npol,beta,A(1+(igt-1)*nrx+(ip-1)*nrx*ngt:igt*nrx+(ip-1)*nrx*ngt,:),nrx)
   call ZGEMM('n','n',nrx,NM,nrx*Ngt*npol,alpha,TLL(1+(ip-1)*nrx:nrx+(ip-1)*nrx,:),nrx,&
        PSI_MOD,nrx*Ngt*npol,beta,C(1+(igt-1)*nrx+(ip-1)*nrx*ngt:igt*nrx+(ip-1)*nrx*ngt,:),nrx)

end do

end do
!$omp end do

deallocate(B,U,Q)
deallocate(HLL,TLL)
deallocate(HCC)
!$omp end parallel

deallocate(ind_cube)
deallocate(betafunc)

t2=SECNDS(t1)
write(*,'(a,F7.3,a)')'100% done in ',t2,' s'

allocate(HLLL(NM,NM),TLLL(NM,NM))
call ZGEMM('c','n',NM,NM,nrx*Ngt*npol,alpha,PSI_MOD,nrx*Ngt*npol,A,nrx*Ngt*npol,beta,HLLL,NM)
call ZGEMM('c','n',NM,NM,nrx*Ngt*npol,alpha,PSI_MOD,nrx*Ngt*npol,C,nrx*Ngt*npol,beta,TLLL,NM)
if(ncell==2)then
   allocate(TLL(NM2,NM))
   call ZGEMM('c','n',NM2,NM,nrx*Ngt*npol,alpha,PSI_MOD_1,nrx*Ngt*npol,C,nrx*Ngt*npol,beta,TLL,NM2)
end if
deallocate(A,C)

write(*,*)'Computing the Hamiltonian blocks...'


if(iyz==1)then
   allocate(A(NM,NM))
   allocate(E(NM))
   allocate(B(NM,NM))

   A=HLLL+TLLL+transpose(dconjg(TLLL))!+TLLLL+transpose(dconjg(TLLLL))
   call SUB_DEF_Z(1,NM,NM,A,E,B)
   ref=E(nband_v)
   write(*,*)'TOP VB AT GAMMA',ref

   deallocate(A)
   deallocate(B)
   deallocate(E)
end if

forall (i = 1:NM) HLLL(i,i)=HLLL(i,i)-ref !!! THIS sets the zero energy point at the top of EV(Gamma)


if(gap_corr)then

   if(nm <= nband_v)then
      write(*,*)'errore in gap_corr'
      stop
   end if
   
   n=1000
   allocate(A(NM,NM),B(NM,NM),C(NM,NM))
   allocate(E(NM))
   allocate(TLL(NM,NM),HLL(NM,NM))
   HLL=0.0_dp
   TLL=0.0_dp

   do ikx=1,n+1
      write(*,*)'ikx',ikx,dble(ikx-1)/dble(n)*2.0_dp*pi
      A=HLLL+TLLL*exp(dcmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
           transpose(dconjg(TLLL))*exp(dcmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
      call SUB_DEF_Z(1,NM,NM,A,E,B)
      
      A=HLLL+TLLL*exp(dcmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
           transpose(dconjg(TLLL))*exp(dcmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
      call ZGEMM('c','n',NM,NM,NM,alpha,B,NM,A,NM,beta,C,NM)
      call ZGEMM('n','n',NM,NM,NM,alpha,C,NM,B,NM,beta,A,NM)

      forall (i = nband_v+1:NM) A(i,i)=A(i,i)+delta_gap
!      A=0.0_dp
!      do k=1,nband_v
!         A(k,k)=E(k)
!      end do
!      do k=nband_v+1,NM
!         A(k,k)=E(k)+delta_gap
      !      end do
      
      call ZGEMM('n','n',NM,NM,NM,alpha,B,NM,A,NM,beta,C,NM)
      call ZGEMM('n','c',NM,NM,NM,alpha,C,NM,B,NM,beta,A,NM)
      
      HLL=HLL+A/dble(n+1)
      TLL=TLL+A/dble(n+1)*exp(dcmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
      
   end do
   
   HLLL=HLL
   TLLL=TLL

   deallocate(TLL,HLL)
   deallocate(A,B,C)
   deallocate(E)  
   
end if


end if
if(refine)then
   allocate(HLLL(NM,NM),TLLL(NM,NM))
end if

ic=1

open(unit=13,file=TRIM(outdir)//'HH00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
if(.not. refine)then
   do i=1,nm
      do j=1,nm
         write(13,*)HLLL(i,j)
      end do
   end do
else if(refine)then
   do i=1,nm
      do j=1,nm
         read(13,*)HLLL(i,j)
      end do
   end do
end if
close(13)


!allocate(A(nm,nm))
!A=transpose(dconjg(TLLL))
!ic=1
open(unit=13,file=TRIM(outdir)//'HH01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
if(.not. refine)then
do i=1,nm
   do j=1,nm
      write(13,*)TLLL(i,j)
   end do
end do
else if(refine)then
do i=1,nm
   do j=1,nm
      read(13,*)TLLL(i,j)
   end do
end do
end if
close(13)
!deallocate(A)


deallocate(Uh)


n=40
M=min(12*nband_v/10,NM)
!write(*,*)'M',M
allocate(hkl(n+1,M))
allocate(E(M))
allocate(A(NM,NM))
do ikx=1,n+1
   A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
        transpose(dconjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
   call SUB_DEF_Z0(1,M,NM,A,E)
   do k=1,M
      hkl(ikx,k)=e(k)
   end do
end do
deallocate(A)
deallocate(E)
do k=1,M!1,min(10,M)
   ii=k!nband_v-min(10,M)/2+k
   !write(*,*)k,ii  
   if(nspin == 1)then
      open(unit=300+k,file=TRIM(outdir)//'Evdisp_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   else
      if(is==1)open(unit=300+k,file=TRIM(outdir)//'Evdisp_up_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
      if(is==2)open(unit=300+k,file=TRIM(outdir)//'Evdisp_dw_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   end if
   do ikx=1,n+1
      write(300+k,*)dble(ikx-1-n/2)/dble(n),hkl(ikx,ii)
   end do
   close(300+k)
end do
deallocate(hkl)

!!!!! STARTING the REFINEMENT


allocate(kx1(nk1))
if(nk1==6)kx1=(/-3.0_dp/6.0_dp, -2.0_dp/6.0_dp, -1.0_dp/6.0_dp, 0.0_dp, 1.0_dp/6.0_dp,  2.0_dp/6.0_dp, 3.0_dp/6.0_dp/)
if(nk1==5)kx1=(/-2.0_dp/5.0_dp, -1.0_dp/5.0_dp, 0.0_dp, 1.0_dp/5.0_dp, 2.0_dp/5.0_dp/)
if(nk1==4)kx1=(/-0.25_dp, 0.0_dp, 0.25_dp, 0.5_dp/)
if(nk1==3)kx1=(/-1.0_dp/3.0_dp, 0.0_dp, 1.0_dp/3.0_dp/)
if(nk1==2)kx1=(/0.0_dp, 0.5_dp/)
if(nk1==1)kx1=(/0.0_dp/)
nnn=nkplus
allocate(kx2(nnn))
if(nnn==2)kx2=(/ 0.0, 0.5/)
if(nnn==3)kx2=(/ -0.3333, 0.0, 0.3333/)
if(nnn==4)kx2=(/-0.25, 0.0, 0.25, 0.5/)
if(nnn==8)kx2=(/-0.375, -0.25, -0.125, 0.0, 0.125, 0.25, 0.375, 0.5/)
nplus=mplus
nn=nplus*nnn
!!!!! inizio
mm1=nf-ni+1
mmm=nk1*mm1+nn
!write(*,*)'reduced basis size=',nk1*mm1,nn,mmm
allocate(B(NM,mmm))
allocate(A(NM,NM))
allocate(E(nf-ni+1))
allocate(C(NM,nf-ni+1))
do ikx=1,nk1
   A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*kx1(ikx)*2.0_dp*pi)+transpose(dconjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*kx1(ikx)*2.0_dp*pi)
   call SUB_DEF_Z(ni,nf,NM,A,E,C)
   B(:,1+(ikx-1)*mm1:ikx*mm1)=C(:,1:mm1)
end do
deallocate(C,E)

if(nplus>0)then
   allocate(E(nplus))
   allocate(C(NM,nplus))
   do ikx=1,nnn
      write(*,*)'adding kx',kx2(ikx)
      A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*kx2(ikx)*2.0_dp*pi)+&
           transpose(dconjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*kx2(ikx)*2.0_dp*pi)
      call SUB_DEF_Z(nf+1,nf+nplus,NM,A,E,C)
      B(:,nk1*mm1+1+(ikx-1)*nplus:nk1*mm1+ikx*nplus)=C(:,1:nplus)
      write(*,*)nf+1,nf+nplus,nk1*mm1+1+(ikx-1)*nplus,nk1*mm1+ikx*nplus
   end do
   deallocate(C,E)
end if
deallocate(A)
write(*,*)'Reduced basis order =',mmm
call MGS(NM,mmm,B) !!! Modified Gram-Schmidt to orthonormalize the Psi functions at different kx

!!!! MY CHANGE
NMODES=mmm
!write(*,*)'NUMBER OF BLOCH STATES=',NMODES


   ic=1
   allocate(C(Nrx*Ngt*npol,NModes))
   call ZGEMM('n','n',Nrx*Ngt*npol,NModes,NM,alpha,psi_mod,Nrx*Ngt*npol,B,NM,beta,C,Nrx*Ngt*npol)

   if(ncell==1)then
      if( nspin == 1 )then
         open(unit=13,file=TRIM(outdir)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
      else
         if(is==1)open(unit=13,file=TRIM(outdir)//'Psi_Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
         if(is==2)open(unit=13,file=TRIM(outdir)//'Psi_Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
      end if
      do j=1,NModes
         do i=1,Nrx*Ngt*npol
            write(13,*)C(i,j)
         end do
      end do
      close(13)
      ic=1
      if( nspin == 1 )then
         open(unit=13,file='Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
      else
         if(is==1) open(unit=13,file='Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
         if(is==2) open(unit=13,file='Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(ic))//'.dat',status='unknown')    
      end if
      do i=1,NM
         do j=1,mmm
            write(13,*)B(i,j)
         end do
      end do
      close(13)
   end if


   deallocate(C)

   if(ncell==2)then
      ic=1
      if( nspin == 1 )then
         nome=trim(indir2)//'Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(ic))//'.dat'
      else
         if(is==1) nome=trim(indir2)//'Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(ic))//'.dat'
         if(is==2) nome=trim(indir2)//'Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(ic))//'.dat'
      end if
      write(*,*)nome
      allocate(C(NM2,MM2))
      open(unit=13,file=nome,status='unknown')
      do i=1,NM2
         do j=1,MM2
            read(13,*)C(i,j)
         end do
      end do
      close(13)
      allocate(HV(MM2,mmm))
      allocate(A(NM2,mmm))
      call ZGEMM('n','n',NM2,mmm,NM,alpha,TLL,NM2,B,NM,beta,A,NM2)
      call ZGEMM('c','n',MM2,mmm,NM2,alpha,C,NM2,A,NM2,beta,HV,MM2)
      deallocate(C,A)
      deallocate(TLL)
      
   end if
   

allocate(HLL(mmm,mmm),TLL(mmm,mmm))
allocate(A(NM,mmm))
call ZGEMM('n','n',NM,mmm,NM,alpha,HLLL,NM,B,NM,beta,A,NM)
call ZGEMM('c','n',mmm,mmm,NM,alpha,B,NM,A,NM,beta,HLL,mmm)
call ZGEMM('n','n',NM,mmm,NM,alpha,TLLL,NM,B,NM,beta,A,NM)
call ZGEMM('c','n',mmm,mmm,NM,alpha,B,NM,A,NM,beta,TLL,mmm)
deallocate(A)
deallocate(kx1)
if(allocated(kx2))deallocate(kx2)


n=40
allocate(hkl(n+1,mmm))
allocate(E(mmm))
allocate(A(mmm,mmm))
do ikx=1,n+1
   A=HLL+TLL*exp(cmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
        transpose(dconjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
   call SUB_DEF_Z0(1,mmm,mmm,A,E)
   do k=1,mmm
      hkl(ikx,k)=e(k)
   end do
end do
deallocate(A)
deallocate(E)
do k=1,nf-ni+nplus+10!min(10,M)
   ii=k!nband_v-min(10,M)/2+k
   !write(*,*)k,i
   if( nspin == 1 )then
      open(unit=300,file=TRIM(outdir)//'Edisp_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   else
      if(is==1)  open(unit=300,file=TRIM(outdir)//'Edisp_up_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
      if(is==2)  open(unit=300,file=TRIM(outdir)//'Edisp_dw_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   end if
   do ikx=1,n+1
      write(300,*)dble(ikx-1-n/2)/dble(n),hkl(ikx,ii)
   end do
   close(300)
end do
deallocate(hkl)


!stop



!!!!
en_select=.false.
if(en_select)then
   deallocate(HLL,TLL)
!   ne=8
!   allocate(Ev(ne))
!   Ev(1)=-1.1 
!   Ev(2)=-0.1 
!   Ev(3)= 1.3 
!   Ev(4)= 1
!   Ev(5)= 1
!   Ev(6)= 2 
!   Ev(7)= 2.5
!   Ev(8)= 5
   
   
   write(*,*)nek,'Ev',Ev
   
   allocate(A(2*NM,2*NM),C(2*NM,2*NM), U(2*NM,2*NM), id(NM,NM))
   allocate(DAL(2*NM))
   allocate(DBE(2*NM))
   allocate(work(4*NM),rwork(16*NM))
   allocate(kval(10000),Q(NM,10000))
   id=0.0_dp
   forall (i=1:NM) id(i,i)=1.0_dp
   j=0
   do i=1,nEK
      A=0.0_dp
      C=0.0_dp
      A(1:NM,1:NM)=HLLL(1:NM,1:NM)-Ev(i)*id(1:NM,1:NM)
      A(1:NM,1+NM:2*NM)=TLLL(1:NM,1:NM)
      A(1+NM:2*NM,1:NM)=id(1:NM,1:NM)
      C(1:NM,1:NM)=-transpose(conjg(TLLL(1:NM,1:NM)))
      C(1+NM:2*NM,1+NM:2*NM)=id(1:NM,1:NM)
      write(*,*)2*nm,size(A),size(C)
      call ZGGEV('n', 'v', 2*NM, A, 2*NM, C, 2*NM, DAL, DBE, dummy, 1, U, 2*NM, WORK, 4*NM, rwork, INFO)
            
      write(*,*)i,'INFO = ',info
      do n=1,2*NM
         if(abs(DBE(n))>0)then
            if(abs(abs(DAL(n)/DBE(n))-1.0_dp)<1.0d-10)then
               !write(*,*)ev(i),DAL(n)/DBE(n),abs(aimag(log(DAL(n)/DBE(n))))/(2.0*pi)
              ! if(abs(abs(aimag(log(DAL(n)/DBE(n))))/(2.0*pi)-Kv(i))<1.0d-2 )then
               !   write(*,*)ev(i),log(DAL(n)/DBE(n))/(2.0*pi),kv(i)
                  j=j+1
                  kval(j)=DAL(n)/DBE(n)
                  Q(1:NM,j)=U(1:NM,n)
               !end if
            end if
         end if
      end do
   end do
   write(*,*)j
   do n=1,j
      write(*,*)n,aimag(log(kval(n))),aimag(log(kval(n)))/(2.0*pi)
   end do
   deallocate(work,rwork,id,A,C,U,dal,dbe)

      
   mmm=j
   write(*,*)'en-selected-basis size =',mmm
   allocate(C(NM,mmm))
   C(1:NM,1:mmm)=Q(1:NM,1:mmm)
   deallocate(Q)
   call MGS(NM,mmm,C) !!! Modified Gram-Schmidt to orthonormalize the Psi functions at different kx
  
allocate(HLLLL(mmm,mmm),TLLLL(mmm,mmm))
allocate(A(NM,mmm))
call ZGEMM('n','n',NM,mmm,NM,alpha,HLLL,NM,C,NM,beta,A,NM)
call ZGEMM('c','n',mmm,mmm,NM,alpha,C,NM,A,NM,beta,HLLLL,mmm)
call ZGEMM('n','n',NM,mmm,NM,alpha,TLLL,NM,C,NM,beta,A,NM)
call ZGEMM('c','n',mmm,mmm,NM,alpha,C,NM,A,NM,beta,TLLLL,mmm)
deallocate(A,C)


n=40
allocate(E(mmm))
allocate(A(mmm,mmm),hkl(n+1,mmm))
do ikx=1,n+1
   A=HLLLL+TLLLL*exp(cmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+transpose(dconjg(TLLLL))*exp(cmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
   call SUB_DEF_Z0(1,mmm,mmm,A,E)
   hkl(ikx,1:mmm)=e(1:mmm)
end do
do ii=1,mmm
   if(nspin == 1)then
      open(unit=300,file=TRIM(outdir)//'disp_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   else
      if(is==1)open(unit=300,file=TRIM(outdir)//'disp_up_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
      if(is==2)open(unit=300,file=TRIM(outdir)//'disp_dw_'//TRIM(STRINGA(ii))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   end if
   do ikx=1,n+1
      write(300,*)dble(ikx-1-n/2)/dble(n),hkl(ikx,ii)
   end do
   close(300)
end do
deallocate(A)
deallocate(E,hkl)
stop
!!!!
else 
NM=mmm
NMODES=NM
write(*,*)'NUMBER OF BLOCH STATES=',NMODES
end if
!!!! END

if(ncell==1)then

   do ic=1,ncell
      if( nspin == 1 )then
         open(unit=13,file=TRIM(outdir)//'H00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
      else
         if(is==1) open(unit=13,file=TRIM(outdir)//'H00_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
         if(is==2) open(unit=13,file=TRIM(outdir)//'H00_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
      end if
      do i=1,nm
         do j=1,nm
            write(13,*)HLL(i,j)
         end do
      end do
      close(13)
   end do

   allocate(A(nm,nm))
   A=transpose(dconjg(TLL))
   ic=1
   if( nspin == 1 )then
      open(unit=13,file=TRIM(outdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
   else
      if(is==1) open(unit=13,file=TRIM(outdir)//'H01_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
      if(is==2) open(unit=13,file=TRIM(outdir)//'H01_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
   end if
   do i=1,nm
      do j=1,nm
         write(13,*)A(i,j)
      end do
   end do
   close(13)
   deallocate(A)
end if

if(ncell==2)then
   allocate(A(nm,mm2))
   A=transpose(dconjg(HV))
   ic=1
   if( nspin == 1 )then
      open(unit=13,file=TRIM(outdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nhet_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
   else
      if (is==1) open(unit=13,file=TRIM(outdir)//'H01_up_nkyz_'//TRIM(STRINGA(iyz))//'_nhet_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
      if (is==2) open(unit=13,file=TRIM(outdir)//'H01_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nhet_'//TRIM(STRINGA(ic))//'.dat',status='unknown')
   end if
   do i=1,nm
      do j=1,mm2
         write(13,*)A(i,j)
      end do
   end do
   close(13)
   deallocate(A)
   deallocate(HV)
end if


!stop

if(allocated(KGt))deallocate(KGt)
allocate(KGt(4,1*Ngt))
KGt=0.0_dp
do j=1,Ngt
   kgt(2,j)=miller_2D(2,j)*b_2(2)+miller_2D(3,j)*b_3(2) + k_vec(2,1+(iyz-1)*nkx)
   kgt(3,j)=miller_2D(2,j)*b_2(3)+miller_2D(3,j)*b_3(3) + k_vec(3,1+(iyz-1)*nkx)
   kgt(4,j)=(kgt(2,j)**2+kgt(3,j)**2)
end do

if(.not.allocated(kgt_kyz))then
   allocate(kgt_kyz(4,1:1*Ngt,Nkyz))
   kgt_kyz=0
end if

KGt_kyz(1:4,1:Ngt,iyz)=KGt(1:4,1:1*Ngt)


if(.not.allocated(HL))then
   allocate(HL(NM,NM,Nkyz))
   HL=0.0_dp
end if
if(.not.allocated(TL))then
   allocate(TL(NM,NM,Nkyz))
   TL=0.0_dp
end if
if(.not.allocated(U_LCBB))then
   allocate(U_LCBB(Nrx*NGt*npol,NM,Nkyz))
   U_LCBB=0.0_dp
end if

!!!! MY CHANGE
HL(:,:,iyz)=HLL(:,:)
TL(:,:,iyz)=TLL(:,:)
allocate(A(Nrx*NGt*npol,NM))
call zgemm('n','n',Nrx*NGt*npol,NM,nkx*M,alpha,PSI_MOD,Nrx*NGt*npol,B,nkx*M,beta,A,Nrx*NGt*npol)
U_LCBB(:,:,iyz)=A(:,:)

deallocate(A,B)
deallocate(PSI_MOD)
if(allocated(PSI_MOD_1))deallocate(PSI_MOD_1)
deallocate(HLLL,TLLL)


deallocate(HLL,TLL)

write(*,*)
write(*,*)'END Kyz',iyz
write(*,*)
write(*,*)
end do !end do iyz
end do ! fine end do is
read(10,'(a)') comment

close(10)
write(*,*)' END of DFT INPUT FILE'

write(*,*)'reduced basis size=',nmodes
IF(N_SPIN/=4)THEN
   deallocate(Deeq)
ELSE
   deallocate(Deeq_so)
END IF



!stop

end subroutine read_QE_output



SUBROUTINE SUB_DEF_Z(Mi,Mf,ny,A,subband,U)
implicit none
integer :: ny,mi,mf
real(dp) :: subband(1:(Mf-Mi+1))
!real(dp) :: WE(1:(ny))
complex(dp) :: A(1:NY,1:NY),U(1:NY,1:(Mf-Mi+1))
integer :: INFO!,LWORK
integer, allocatable :: iwork(:), supp(:)
complex(dp), allocatable :: work(:)
real(dp), allocatable :: rwork(:)
REAL(DP), EXTERNAL :: DLAMCH


allocate(WORK(20*ny))
allocate(RWORK(24*ny))
allocate(IWORK(10*ny))
allocate(Supp(2*ny))

call ZHEEVR('V','I','U',ny,A,ny,0.0,0.0,mi,mf,2*DLAMCH('S'),(Mf-Mi+1),subband,U,ny,SUPP,WORK,20*ny,RWORK,24*ny,IWORK,10*ny,INFO)
!!call ZHEEVR('N','I','L',ny,A,ny,0.0_dp,0.0_dp,1,ny,2*DLAMCH('S'),(ny),subband,U,ny,SUPP,WORK,20*ny,RWORK,24*ny,IWORK,10*ny,INFO)

deallocate(work)
deallocate(rwork)
deallocate(supp)
deallocate(iwork)
!write(*,*)'info=',info
if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: SUB_DEF HAS FAILED. INFO=',INFO
   stop
endif     

END SUBROUTINE SUB_DEF_Z


subroutine coefficienti (p,cv)
implicit none
integer :: n,d,p,i,imin,imax,ioffset,noffset
real(dp), allocatable :: A(:,:), b(:), c(:)
real(dp) :: dfatt, cv(0:p/2)

d=2
if(mod(p,2)/=0)then
write(*,*)'p must be even! Stopping the calculation in coefficienti','p=',p
stop
end if

allocate(A(d+p-1,d+p-1))
imin=-p/2
imax=p/2
do i=imin,imax
   ioffset = 1-imin
   noffset=1
   do n=0,p
      A(n+noffset,i+ioffset)=i**n
   end do
end do
allocate(b(p+1),c(p+1))
b=0.0_dp 
b(3)=1.0_dp 
dfatt = 2.0_dp

call invertR(A,p+1,p+1)

c=-dfatt*matmul(A,b)
cv(0:p/2)=c(p/2+1:p+1)
cv(0)=cv(0)/2.0_dp

deallocate(A)
deallocate(b,c)

end subroutine coefficienti

 
 subroutine MGS(np,nm,Q)
 !!!! MODIFIED GRAM SCHMIDT ALGORITHM
 implicit none
 
 integer, intent(IN) :: np, nm
 complex(dp), intent(INOUT) ::  Q(np,nm)
 integer :: j,k
 complex(dp), allocatable :: X(:,:),R(:,:)!,v(:)
 
 allocate(X(np,nm))  

 X=Q
 
 allocate(R(nm,nm)) 
 R=0.0_dp
 Q=0.0_dp
 R(1,1)=sqrt(dble(dot_product(X(:,1),X(:,1))))
 Q(:,1)=X(:,1)/R(1,1)
  
 do k=2,nm
    do j=1,k-1
       R(j,k)=dot_product(Q(:,j),X(:,k))
       X(:,k) =  X(:,k) - R(j,k)*Q(:,j)
    end do
    R(k,k) = sqrt(dble(dot_product(X(:,k),X(:,k))))
    Q(:,k) = X(:,k)/R(k,k) 
 end do

deallocate(R,X)

end subroutine MGS
 

END MODULE Pseudopot_so_gen
