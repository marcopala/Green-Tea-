! Copyright or © or Copr. Marco Pala (February 24, 2022)

! e-mail: marco.pala@c2n.upsaclay.fr ;  marco.pala@cnrs.fr

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

MODULE Pseudopot_so_gen

  USE static
  USE indata

  IMPLICIT NONE

  SAVE

!!! Variables  

CONTAINS


subroutine read_QE_output
  implicit none
  
character(len=80) :: comment
character(len=6), allocatable :: specie(:)

integer :: N_at, N_typ_at, ikx,l,k,m,nn,nnn,mmm,nrx,nrx0,nry,nrz,n1,n2,n3,m1,igt,jgt,ix,iy,iz
integer :: NPOL, N_spin, n_bands, nspin, N_beta
integer :: ivec(3),ncf,nm,nplus,iyz,ii,jj
integer :: is,ic,ip,jp,nt,na,ikb,info
integer(k15) :: N_Gvec, N_k_pts,i,j,n,nnpw,npw
integer(k15) :: max_mill_1, max_mill_2,max_mill_3, min_mill_1, min_mill_2, min_mill_3
integer(k15), allocatable :: miller(:,:), ind_miller(:,:,:), miller_2D(:,:), ind_miller_2D(:,:)
integer(k15), allocatable :: N_projs(:), N_PW(:), ityp(:), ind_cube(:), ind_kG(:,:)
integer, allocatable :: ind_kx(:), ind_bnd(:), inds(:)

real(dp) :: a_1(3),a_2(3),a_3(3)
real(dp) :: b_1(3),b_2(3),b_3(3)
real(dp) :: vec(3),t0,a0
real(dp) :: Ecutoff,ref,tmp1,tmp2
real(4) :: t1,t2

real(dp), allocatable :: x_at(:), y_at(:), z_at(:), Gvec(:,:), R_at(:,:)
real(dp), allocatable :: k_vec(:,:),coeff(:),kx1(:),kx2(:), in_kx(:)
real(dp), allocatable :: E(:), KGt(:,:), Gx(:), hkl(:,:)

complex(dp) :: tmp, dummy(1,1)
complex(dp), allocatable :: beta_fun(:,:),KS_fun(:,:),psi_add(:,:),Psi_mod(:,:),Psi_mod_1(:,:),Psi_mod_0(:,:),betafunc(:,:)
complex(dp), allocatable :: A(:,:),B(:,:),C(:,:),C1(:,:),D(:,:),Q(:,:),Uh(:,:),U(:,:),Si(:,:)
complex(dp), allocatable :: Vloc(:,:),HCC(:,:),Deeq_so(:,:,:,:),Deeq(:,:,:),HV(:,:)
complex(dp), allocatable :: HLL(:,:),TLL(:,:),HLLL(:,:),TLLL(:,:),HLLLL(:,:),TLLLL(:,:)
complex(dp), allocatable :: dal(:),dbe(:),work(:),rwork(:),kval(:),id(:,:)

character(256) :: nome

logical :: en_select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

 a0=a_1(1)*bohr

 if(abs(1.0_dp - a_1(1)*bohr/ac1) > 1.0d-6) then
    write(*,*)'Severe warning: ac1 is not the same of the QE simulation'
    write(*,*)a_1(1)*bohr, ac1, 1.0_dp - a_1(1)*bohr/ac1
    stop
 end if
 if(abs(1.0_dp - a_2(2)*bohr/ac2) > 1.0d-6) then
    write(*,*)'Severe warning: ac2 is not the same of the QE simulation'
    write(*,*)a_2(2)*bohr, ac2, 1.0_dp - a_2(2)*bohr/ac2
    stop
 end if
 
 if(abs(1.0_dp - a_3(3)*bohr/ac3) > 1.0d-6) then
    write(*,*)'Severe warning: ac1 is not the same of the QE simulation'
    write(*,*)a_3(3)*bohr, ac3, 1.0_dp - a_3(3)*bohr/ac3
    stop
 end if
 
write(*,*)
write(*,*)' Lattice parameter a0 = ', a0, 'cm'
write(*,*)

read(10,'(a)') comment
write(*,*) comment
read(10,*) b_1(:), b_2(:), b_3(:)
write(*,*)

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
      read(10,'(2e25.15)') Vloc(i,is) ! local term
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

!!$IF ( n_spin /= 4 ) THEN
!!$   allocate(hkl(nn,nn))
!!$   allocate(Deeq(nn,nn,N_typ_at))
!!$ELSE
!!$   allocate(Deeq_so(nn,nn,n_spin,N_typ_at))
!!$END IF

IF ( lspinorb ) THEN
   allocate(Deeq_so(nn,nn,n_spin,N_typ_at))
ELSE
   allocate(hkl(nn,nn))
   allocate(Deeq(nn,nn,N_typ_at))
END IF
   

read(10,'(a)') comment
do j=1,N_typ_at
   read(10,'(a)') comment
   read(10,*) i,N_projs(j)
   write(*,*) 'Nprojs', i,N_projs(j)
  IF(.not. lspinorb) THEN ! IF ( n_spin /= 4 ) THEN
      read(10,*)hkl(1:N_projs(j),1:N_projs(j))
      Deeq(1:N_projs(j),1:N_projs(j),j)=hkl(1:N_projs(j),1:N_projs(j))
   ELSE
      do is=1,N_spin
         read(10,*)comment
         read(10,'(2e25.15)')Deeq_so(1:N_projs(j),1:N_projs(j),is,j)
      end do
   END IF
end do

if(allocated(hkl)) deallocate(hkl)
read(10,'(a)') comment
read(10,*) N_beta
write(*,*) 'N_beta =',N_beta

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
do n2=-1000,1000 !minval(miller(1,:)),maxval(miller(1,:))
   vec(1:3)=n2*b_1(1:3)
   if(t0*vec(1)**2*(2.0_dp*pi/a0)**2 < 4.0_dp*Ecutoff )nrx=nrx+1
end do
nrx=(nrx-1)
write(*,*)'NRX =',nrx
!nrx=maxval(miller(1,:))-minval(miller(1,:))+1
nry=0
do n2=-1000,1000!minval(miller(2,:)),maxval(miller(2,:))
   vec(2:3)=n2*b_2(2:3)
   if(t0*vec(2)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff )nry=nry+1
end do
nry=nry-1
write(*,*)'NRY =',nry
nrz=0
do n3=-1000,1000!minval(miller(3,:)),maxval(miller(3,:))
   vec(2:3)=n3*b_3(2:3)
   if(t0*vec(3)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff )nrz=nrz+1
end do
nrz=nrz-1
write(*,*)'NRZ =',nrz

ndx=nrx
ndy=nry
ndz=nrz

ncf=min((nrx/2-2),8)
allocate(coeff(0:ncf))

call coefficienti (2*ncf,coeff)
write(*,*)
write(*,*)'Discretization order =',ncf
write(*,*)

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

nrx0=nrx
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
    
write(*,*)
write(*,*)'Kyz number =',iyz
write(*,*)

do ikx=1,Nkx
   write(*,*)'ikx =',ikx,'is =',is
   read(10,'(a)') comment
   write(*,'(a)') comment
   if(NSPIN==2) then
!      read(10,*) i
      write(*,*) 'Collinear magnetic system, current_spin = ',is
   end if
   read(10,*) k_vec(:,ikx+(iyz-1)*nkx)
   k_vec(:,ikx+(iyz-1)*nkx)=k_vec(:,ikx+(iyz-1)*nkx)/bohr/(2.0_dp*pi/a0)
   write(*,*) k_vec(:,ikx+(iyz-1)*nkx)

   read(10,'(a)') comment
   read(10,*) N_PW(ikx+(iyz-1)*nkx)
   npw=N_PW(ikx+(iyz-1)*nkx)
   write(*,*)'npw =',npw
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
   write(*,'(a, i4, a)')' Using',N_bands,' KS solutions from the DFT simulation'


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

!! check of the ordering of k_vec
do ikx=1,nkx-1
   if(k_vec(1,ikx+(iyz-1)*nkx) .eq. k_vec(1,ikx+1+(iyz-1)*nkx))then
      write(*,*)'problem with kvec(1) definition or ordering, please check the manual'
      write(*,*)1,k_vec(1,ikx+(iyz-1)*nkx), k_vec(1,ikx+1+(iyz-1)*nkx)
      write(*,*)'Simulation aborted'
      stop
   end if
   if(k_vec(2,ikx+(iyz-1)*nkx) .ne. k_vec(2,ikx+1+(iyz-1)*nkx))then
      write(*,*)'WARNING : pssbl pb  with kvec(2) definition or ordering, please check the manual'
      write(*,*)2,k_vec(2,ikx+(iyz-1)*nkx), k_vec(2,ikx+1+(iyz-1)*nkx)
      write(*,*)'Simulation continues at your own risk'
   end if
   if(k_vec(3,ikx+(iyz-1)*nkx) .ne. k_vec(3,ikx+1+(iyz-1)*nkx))then
      write(*,*)'WARNING : pssbl pb  with kvec(3) definition or ordering, please check the manual'
      write(*,*)3,k_vec(3,ikx+(iyz-1)*nkx), k_vec(3,ikx+1+(iyz-1)*nkx)
      write(*,*)'Simulation continues at your own risk'
   end if
end do



!allocate(Psi_old(nkx*nrx*Ngt*npol,N_bands*nkx))
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
              exp(cmplx(0.0_dp,1.0_dp,kind(dp))*(k_vec(1,ikx+(iyz-1)*nkx)+Gx(m1))*(2.0_dp*pi)*dble(n-1)/dble(nrx))/sqrt(dble(Nkx*nrx))
      end do
   end do
end do
allocate(B(nkx*nrx,nkx*N_bands))
do jp=1,npol
   do jgt=1,Ngt
      call zgemm ('c','n',nkx*nrx,nkx*N_bands,nkx*nrx,alpha,Uh,nkx*nrx,&
           A(1+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt:jgt*nkx*nrx+(jp-1)*nkx*nrx*ngt,1:nkx*N_Bands),nkx*nrx,beta,B(1:nkx*nrx,1:nkx*n_bands),nkx*nrx)
      do n=1,nkx*nrx
         Psi_mod(n+(jgt-1)*nkx*nrx+(jp-1)*nkx*nrx*ngt,1:Nkx*N_bands)=B(n,1:Nkx*N_bands)*sqrt(dble(Nkx))
      end do
   end do
end do
deallocate(A,B)

M=N_bands
NM=nkx*N_bands
NMODES=NM
write(*,*)'NM=',NM
write(*,*)'Total number of KS solutions from the DFT simulation =',NM

if(.not. refine)then
write(*,*)
write(*,*)'*******************************************************************************'
write(*,*)'*******************************************************************************'
write(*,'(a,I3,a)')' WARNING: by using Nomp = ',Nomp,' the required memory is approximatively '
write(*,'(F8.1,a)')     dble(nrx*Ngt)*dble(N_beta*nkx)*16.0d-9+& !betafunc
     dble(Nomp*nkx*nrx*npol)*dble(nkx*nrx*Ngt*npol)*16.0d-9+&  ! HCC
     dble(Nomp*2*Nrx*npol)*dble(nrx*Ngt*npol)*16.0d-9+&  ! HLL,TLL
     dble(Nomp*nrx)*dble(nrx*ngt)*16.0d-9+&  ! Q
     dble(Nomp*2*nkx*nrx)*dble(nkx*nrx)*16.0d-9+&  ! B, U
     dble(2*nm)*dble(nrx*Ngt*npol)*16.0d-9, ' Gb'
write(*,*)
write(*,'(a)')' BE SURE THAT YOUR SYSTEM CAN HANDLE THIS SIMULATION!'
write(*,*)
write(*,*)'*******************************************************************************'
write(*,*)'*******************************************************************************'
write(*,'(a)')' Also note that, by decreasing Nomp, the allocated memory can be reduced up to '
write(*,'(F8.1,a)')    dble(nrx*Ngt)*dble(N_beta*nkx)*16.0d-9+& !betafunc
     dble(nkx*nrx*npol)*dble(nkx*nrx*Ngt*npol)*16.0d-9+&  ! HCC
     dble(2*Nrx*npol)*dble(nrx*Ngt*npol)*16.0d-9+&  ! HLL,TLL
     dble(nrx)*dble(nrx*ngt)*16.0d-9+&  ! Q
     dble(2*nkx*nrx)*dble(nkx*nrx)*16.0d-9+&  ! B, U
     dble(2*nm)*dble(nrx*Ngt*npol)*16.0d-9, ' Gb'
write(*,*)'*******************************************************************************'
write(*,*)'*******************************************************************************'
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
      U(1+(jp-1)*NGt:jp*NGt,j)=exp( cmplx(0.0_dp,-1.0_dp,kind(dp))*KGt(2,1:NGt)*2.0_dp*pi/a0*dble(iy)*Dy+&
           cmplx(0.0_dp,-1.0_dp,kind(dp))*KGt(3,1:NGt)*2.0_dp*pi/a0*dble(iz)*Dz )/sqrt(dble((nry)*(nrz)))
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
   
   call ZGEMM('n','n',NM,(nry)*(nrz),NGt*npol,alpha,B,NM,U,NGt*npol,beta,A,NM)

   do iy=1,nry
      do iz=1,nrz
         tmp=0.0d0
         do i=1,nband_v
            tmp=tmp+(A(i,iy+(iz-1)*(nry))*conjg( A(i,iy+(iz-1)*(nry)) ))
         end do
         write(2000+ix,*)(iy-1)*dy*1e8,(iz-1)*dz*1e8,dble(tmp)
               end do
      write(2000+ix,*)
   end do
   close(2000+ix)
end do
deallocate(A,B,U)

end if

deallocate(Psi_mod)
allocate(PSI_MOD(Nrx*Ngt*npol,NM))
!allocate(Psi_old(nrx*Ngt*npol,N_bands*nkx))
!Psi_mod=0.0_dp

if(.not. refine)then
!write(*,*)
!t1=SECNDS(0.0)
!write(*,*)'Starting the orthonormalization...'
ic=1
!call MGS(nrx0*Ngt*npol,NM,C(1:nrx0*ngt*npol,1:NM))
!psi_old=C
!!call ortonorma(nrx0*Ngt*npol,NM,C,PSI_mod)
!psi_mod=C!psi_old
!C=psi_mod
do jp=1,npol
   do jgt=1,Ngt
      do ix=1,nrx0
         PSI_MOD(ix+(ic-1)*Nrx0+(jgt-1)*nrx+(jp-1)*nrx*ngt,1+(ic-1)*NM:ic*NM)=C(ix+(jgt-1)*Nrx0+(jp-1)*nrx0*ngt,1:NM)      
      end do  
   end do
end do
   
!do j=1,NM
!   do i=1,Nrx*Ngt*npol
!      if(C(i,j)/=C(i,j))then
!         write(*,*)'prob with the MGS',j
!         write(*,*)'stopping the simulaton'
!         stop
!         exit
!      end if
!   end do
!end do

!t2=SECNDS(t1)
!write(*,'(a,F7.3,a)')' Orthonormalization done in ',t2,' s'
!write(*,*)
end if

deallocate(C)
deallocate(ks_fun)



if(ncell==1)then
   !!!ic=1
if(nspin==1) then
   open(unit=13,file=TRIM(outdir)//'PPsi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
else
if(is==1)   open(unit=13,file=TRIM(outdir)//'PPsi_Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
if(is==2)   open(unit=13,file=TRIM(outdir)//'PPsi_Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
end if

if(.not. refine)then
   write(*,*)
   write(*,*)'Writing PPSI files...'
   do j=1,NM
      do i=1,Nrx*Ngt*npol
         write(13,'(2e25.15)')PSI_MOD(i,j)
      end do
   end do
   write(*,*)'done'
else if(refine)then
   write(*,*)
   write(*,*)'reading PPSI files...',Nrx*Ngt*npol
   do j=1,NM
      do i=1,Nrx*Ngt*npol
         read(13,'(2e25.15)')tmp1,tmp2
         PSI_MOD(i,j)=cmplx(tmp1,tmp2,kind(dp))
      end do
   end do
   write(*,*)'done'

end if
close(13)
end if
write(*,*)


allocate(Si(NM,NM))
call ZGEMM('c','n',NM,NM,Nrx*Ngt*npol,alpha,PSI_MOD,Nrx*Ngt*npol,PSI_MOD,Nrx*Ngt*npol,beta,Si,NM)

!do i=1,nm
!   do j=1,nm
!      write(1112,*)i,j,(Si(i,j))
!      write(1113,*)i,j,abs(Si(i,j))
!   end do
!   write(1112,*)
!   write(1113,*)
!end do

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

if(ncell==2)then
   allocate(C1(nrx*Ngt*npol,nm1))
   C1=0.0_dp
   nome=trim(indir1)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat'
   allocate(PSI_MOD_1(Nrx*Ngt*npol,NM1))
   open(unit=13,file=nome,status='unknown')
   do j=1,NM1
      do i=1,Nrx*Ngt*npol
         read(13,'(2e25.15)')tmp1,tmp2
         PSI_MOD_1(i,j)=cmplx(tmp1,tmp2,kind(dp))
      end do
   end do
   close(13)
   nome=trim(indir0)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat'
   allocate(PSI_MOD_0(Nrx*Ngt*npol,NM0))
   open(unit=13,file=nome,status='unknown')
   do j=1,NM0
      do i=1,Nrx*Ngt*npol
         read(13,'(2e25.15)')tmp1,tmp2
         PSI_MOD_0(i,j)=cmplx(tmp1,tmp2,kind(dp))
      end do
   end do
   close(13)
end if

write(*,*)
write(*,*)'Transforming H in the reduced basis (this can be long) ...'
write(*,*)
t1=SECNDS(0.0)

!$omp parallel default(none) private(igt,jgt,ikx,ikb,nt,na,m1,n1,n,nn,l,ivec,&
!$omp i,j,ip,jp,tmp,vec,B,U,HLL,TLL,HCC,Q,D)&
!$omp shared(nomp,is,ic,ncell,ind_kg,nkx,nrx,iyz,npol,nm,nm1,ngt,ncf,N_projs,ityp,N_typ_at,N_spin,&
!$omp N_at,N_beta,ind_cube,ind_miller,miller,miller_2D,min_mill_1,max_mill_1,b_2,b_3,&
!$omp a0,t0,Dx,Dy,Dz,coeff,k_vec,Uh,betafunc,Deeq,Deeq_so,PSI_MOD,PSI_MOD_1,A,C,C1,vloc)

allocate(HCC(nkx*nrx*npol,nkx*nrx*Ngt*npol))
allocate(B(nkx*nrx,nkx*nrx),U(nkx*nrx,nkx*nrx))
allocate(Q(nrx,nrx*ngt))
allocate(HLL(Nrx*npol,nrx*Ngt*npol))
allocate(TLL(Nrx*npol,nrx*Ngt*npol))

!$omp do
do igt=1,Ngt
HCC=0.0_dp

HLL=0.0_dp
TLL=0.0_dp

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


end do !enddo ikx


   
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



if(igt <= (Ngt/nomp) .and. mod(igt,10)==0) write(*,*)floor(dble(igt-1)/dble(Ngt/nomp)*100), '% done'


   vec(2:3)= ( miller_2D(2,igt)*b_2(2:3) + miller_2D(3,igt)*b_3(2:3) + k_vec(2:3,1+(iyz-1)*nkx) )*(2.0_dp*pi/a0) 
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
   if(ncell==2)then
!!! multiplication of the PSI_{i+1}
      call ZGEMM('n','n',nrx,nm1,nrx*Ngt*npol,alpha,TLL(1+(ip-1)*nrx:nrx+(ip-1)*nrx,:),nrx,&
           PSI_MOD_1,nrx*Ngt*npol,beta,C1(1+(igt-1)*nrx+(ip-1)*nrx*ngt:igt*nrx+(ip-1)*nrx*ngt,:),nrx)
   end if
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
write(*,*)'100% done in ',t2,' s'

allocate(HLLL(NM,NM),TLLL(NM,NM))
call ZGEMM('c','n',NM,NM,nrx*Ngt*npol,alpha,PSI_MOD,nrx*Ngt*npol,A,nrx*Ngt*npol,beta,HLLL,NM)
call ZGEMM('c','n',NM,NM,nrx*Ngt*npol,alpha,PSI_MOD,nrx*Ngt*npol,C,nrx*Ngt*npol,beta,TLLL,NM)
deallocate(A)
deallocate(C)

if(ncell==2)then

   allocate(HV(NM0,NM1))
!!! multiplication of the PSI_{i}^\dagger
   call ZGEMM('c','n',NM0,NM1,nrx*Ngt*npol,alpha,PSI_MOD_0,nrx*Ngt*npol,C1,nrx*Ngt*npol,beta,HV,NM0)
   deallocate(C1)
   
!   allocate(A(nm1,nm0))
!   A=transpose(conjg(HV))

!!! saving on the disk the coupling matrix corresponding to the ihet heterostructure 
   open(unit=13,file=TRIM(outdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nhet_'//TRIM(STRINGA(ihet))//'.dat',status='unknown')
   do i=1,nm0
      do j=1,nm1
         write(13,'(2e25.15)')HV(i,j)
      end do
   end do
   close(13)

!   deallocate(A)
   deallocate(HV)
   
   deallocate(psi_mod_0)
   deallocate(psi_mod_1)
   
end if


if(ncell==1)then
   
write(*,*)'Computing the Hamiltonian blocks...'


if(iyz==1)then
   allocate(A(NM,NM))
   allocate(E(NM))
   allocate(B(NM,NM))

   A=HLLL+TLLL+transpose(conjg(TLLL))
   call SUB_DEF_Z0_GEN(1,NM,NM,A,Si,E)!!!   call SUB_DEF_Z(1,NM,NM,A,E,B)
   ref=E(nband_v)
   write(*,*)'TOP VB AT GAMMA',ref

   deallocate(A)
   deallocate(B)
   deallocate(E)
end if

HLLL=HLLL-ref*Si !!! THIS sets the zero energy point at the top of EV(Gamma)

end if
end if

deallocate(Uh)

if(ncell==1)then

   
if(refine)then
   allocate(HLLL(NM,NM),TLLL(NM,NM))
end if

ic=1

open(unit=13,file=TRIM(outdir)//'HH00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
if(.not. refine)then
   do i=1,nm
      do j=1,nm
         write(13,'(2e25.15)')HLLL(i,j)
      end do
   end do
else if(refine)then
   do i=1,nm
      do j=1,nm
         read(13,'(2e25.15)')tmp1,tmp2
         HLLL(i,j)=cmplx(tmp1,tmp2,kind(dp))
      end do
   end do
end if
close(13)



open(unit=13,file=TRIM(outdir)//'HH01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
if(.not. refine)then
do i=1,nm
   do j=1,nm
      write(13,'(2e25.15)')TLLL(i,j)
   end do
end do
else if(refine)then
do i=1,nm
   do j=1,nm
      read(13,'(2e25.15)')tmp1,tmp2
      TLLL(i,j)=cmplx(tmp1,tmp2,kind(dp))
   end do
end do
end if
close(13)


if(gap_corr)then

   write(*,*)
   write(*,*)'WARNING!'
   write(*,*)
   write(*,*)'You are using the gap correction option. Use it with caution: it can result in an unphysical deformation of band dispersion'
   write(*,*)
   
   if(nm <= nband_v)then
      write(*,*)'error in gap_corr: nm <= nband_v'
      stop
   end if
   
   n=40000
   allocate(A(NM,NM),B(NM,NM),C(NM,NM))
   allocate(E(NM))
   allocate(TLLLL(NM,NM),HLLLL(NM,NM))
   HLLLL=0.0_dp
   TLLLL=0.0_dp

   do ikx=1,n+1
      
      if(mod(ikx-1,1000)==0) write(*,*)'ikx',ikx-1,dble(ikx-1)/dble(n)*2.0_dp*pi
      A=HLLL+TLLL*exp(im*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
           transpose(conjg(TLLL))*exp(-im*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
      call SUB_DEF_Z_GEN(1,NM,NM,A,Si,E,B) !!    call SUB_DEF_Z(1,NM,NM,A,E,B)

      A=0.0_dp
      forall ( i = 1 : nband_v ) A(i,i)=E(i)
      forall ( i = nband_v+1:nm) A(i,i)=E(i)+delta_gap
      
   call ZGEMM('n','n',NM,NM,NM,alpha,B,NM,A,NM,beta,C,NM)
   call ZGEMM('n','c',NM,NM,NM,alpha,C,NM,B,NM,beta,A,NM)
   
      HLLLL=HLLLL+A/dble(n+1)
      TLLLL=TLLLL+A/dble(n+1)*exp(-im*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
      
   end do
   
   deallocate(E,B)
   allocate(E(N_bands),B(NM,N_bands))
   allocate(U(NM,NM))
   if(NM .ne. nkx*n_bands)then
      write(*,*)'error',NM,n_bands,nkx
      stop
   end if

   HLLL=HLLLL
   TLLL=TLLLL
   deallocate(TLLLL,HLLLL)
   deallocate(A,B,C)
   deallocate(E)  
   deallocate(U)  
   
end if

n=40
M=min(12*nband_v/10,NM)

allocate(hkl(n+1,M))
allocate(E(M))
allocate(A(NM,NM))
do ikx=1,n+1
   A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
        transpose(conjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
   call SUB_DEF_Z0_GEN(1,M,NM,A,Si,E)  !!!!! call SUB_DEF_Z0(1,M,NM,A,E)
   do k=1,M
      hkl(ikx,k)=e(k)
   end do
end do
deallocate(A)
deallocate(E)
do k=1,M
   ii=k
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

if(nk1 > nkx)then
   write(*,*) 'ERROR: inconsistency of selected Bloch states, nkx is too large'
   write(*,*) nk1, nkx
   stop
end if
if(nkplus > nkx)then
   write(*,*) 'ERROR: inconsistency of selected Bloch states, nkplus is too large'
   write(*,*) nkplus, nkx
   stop
end if

allocate(kx1(nk1))
if(nk1==8)kx1=(/-0.375_dp, -0.25_dp, -0.125_dp, 0.0_dp, 0.125_dp, 0.25_dp, 0.375_dp, 0.5_dp/)
if(nk1==6)kx1=(/-1.0_dp/3.0_dp, -1.0_dp/6.0_dp, 0.0_dp, 1.0_dp/6.0_dp,  1.0_dp/3.0_dp, 0.5_dp/)
if(nk1==4)kx1=(/-0.25_dp, 0.0_dp, 0.25_dp, 0.5_dp/)
if(nk1==2)kx1=(/0.0_dp, 0.5_dp/)
if(nk1==1)kx1=(/0.0_dp/)
nnn=nkplus


if(nkplus == 2 .or. nkplus == 4 .or. nkplus == 6 .or. nkplus == 8)then
   write(*,*)'nkplus=',nkplus
else
   write(*,*)'pb with nkplus: incorrect value', nkplus
   stop
end if


allocate(kx2(nnn))
if(nnn==2)kx2=(/ 0.0, 0.5/)
if(nnn==4)kx2=(/-0.25, 0.0, 0.25, 0.5/)
if(nnn==6)kx1=(/-1.0_dp/3.0_dp, -1.0_dp/6.0_dp, 0.0_dp, 1.0_dp/6.0_dp,  1.0_dp/3.0_dp, 0.5_dp/)
if(nnn==8)kx2=(/-0.375, -0.25, -0.125, 0.0, 0.125, 0.25, 0.375, 0.5/)
nplus=mplus
nn=nplus*nnn
!!!!
mm1=nf-ni+1
if(.not.refine) nkx_add = 0
mmm=nk1*mm1+nn + nkx_add
!!!!
allocate(B(NM,mmm))
allocate(A(NM,NM))
allocate(E(nf-ni+1))
allocate(C(NM,nf-ni+1))
do ikx=1,nk1
   A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*kx1(ikx)*2.0_dp*pi)+transpose(conjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*kx1(ikx)*2.0_dp*pi)
   call SUB_DEF_Z_GEN(ni,nf,NM,A,Si,E,C)   !!!call SUB_DEF_Z(ni,nf,NM,A,E,C)
   B(:,1+(ikx-1)*mm1:ikx*mm1)=C(:,1:mm1)

end do
deallocate(C,E)

if(nplus>0)then
   allocate(E(nplus))
   allocate(C(NM,nplus))
   do ikx=1,nnn
!!!      write(*,*)'adding kx',kx2(ikx)
      A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*kx2(ikx)*2.0_dp*pi)+&
           transpose(conjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*kx2(ikx)*2.0_dp*pi)
      call SUB_DEF_Z_GEN(nf+1,nf+nplus,NM,A,Si,E,C)      
      B(:,nk1*mm1+1+(ikx-1)*nplus:nk1*mm1+ikx*nplus)=C(:,1:nplus)
!!!      write(*,*)nf+1,nf+nplus,nk1*mm1+1+(ikx-1)*nplus,nk1*mm1+ikx*nplus
   end do
   deallocate(C,E)
end if

if(Nkx_add > 0)then
   allocate ( psi_add(nrx*Ngt*npol,nkx_add) )
   allocate(E(1))
   allocate(C(NM,1))
   do ikx=1,Nkx_add
      A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*xk_add(ikx)*2.0_dp*pi)+ &
           transpose(conjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*xk_add(ikx)*2.0_dp*pi)
      call SUB_DEF_Z_GEN(nbnd_add(ikx),nbnd_add(ikx),NM,A,Si,E,C)  
      B(:,nk1*mm1+nn+ikx)=C(:,1)
      write(*,*)'nbnd',nbnd_add(ikx),'additional kx',xk_add(ikx),E(1)
   end do
   deallocate(C,E)
end if
   
deallocate(A)
write(*,*)
write(*,*)' Reduced basis order =',mmm
write(*,*)

!!!$!!!call MGS(NM,mmm,B) !!! Modified Gram-Schmidt to orthonormalize the Psi functions at different kx
!!!$!!allocate(A(NM,mmm))
!!!$!!call ortonorma(NM,mmm,B,A)
!!!$!!B=A
!!!$!!deallocate(A)

NMODES=mmm

ic=1
allocate(C(Nrx*Ngt*npol,NModes))
call ZGEMM('n','n',Nrx*Ngt*npol,NModes,NM,alpha,psi_mod,Nrx*Ngt*npol,B,NM,beta,C,Nrx*Ngt*npol)


deallocate(Si)
allocate(Si(NModes,NModes))
call ZGEMM('c','n',NModes,NModes,nrx*ngt*npol,alpha,C,nrx*ngt*npol,C,nrx*ngt*npol,beta,Si,NModes)

!!!!!!!!!!!!!!!!!!!!

   !allocate(A(nrx*ngt*npol,nmodes))
   !A=0.0_dp

i=0
do ikx=1,nk1
do j=1,nkx
   if(abs(kx1(ikx)-k_vec(1,j+(iyz-1)*nkx))<1.0d-6)then
      i=i+1
   end if
end do
end do
if(i /= nk1)then
   write(*,*)'pb with nk1',i,nk1
   stop
end if

allocate(in_kx(nmodes),ind_bnd(nmodes))
do ikx=1,nk1
   do j=1,nkx
      if(abs(kx1(ikx)-k_vec(1,j+(iyz-1)*nkx))<1.0d-3)then
         write(*,*)'ikx=',j,'kx1=',k_vec(1,j+(iyz-1)*nkx)
         do i=1,mm1
            in_kx(i+(ikx-1)*mm1)=k_vec(1,j+(iyz-1)*nkx)
            ind_bnd(i+(ikx-1)*mm1)=i+ni-1
         end do
      end if
   end do
end do

if(nplus>0)then

i=0
do ikx=1,nnn
do j=1,nkx
   if(abs(kx2(ikx)-k_vec(1,j+(iyz-1)*nkx))<1.0d-3)then
      i=i+1
   end if
end do
end do
if(i /= nkplus)then
   write(*,*)'pb with nkplus'
   stop
end if
do ikx=1,nnn
   do j=1,nkx
      if(abs(kx2(ikx)-k_vec(1,j+(iyz-1)*nkx))<1.0d-3)then
         write(*,*)'ikx=',j,'kx2=',k_vec(1,j+(iyz-1)*nkx)
         do i=1,nplus
            in_kx(nk1*mm1 + i+(ikx-1)*nplus)=k_vec(1,j+(iyz-1)*nkx)
            ind_bnd(nk1*mm1 + i+(ikx-1)*nplus)=mm1+i+ni-1
         end do
      end if
   end do
end do

end if

if(nkx_add > 0)then
   do ikx=1,nkx_add
      in_kx(nk1*mm1+nn+ikx) = xk_add(ikx)
      ind_bnd(nk1*mm1 + nn + ikx) = nbnd_add(ikx)
   end do
end if



do i=1,NModes
   do j=1,NModes
      write(700+iyz,*)i,j,abs(Si(i,j))
   end do
   write(700+iyz,*)
end do
close(700+iyz)

allocate(inds(NModes))

do i=1,nmodes
   do ikx=1,Nkx
      do n=1,N_bands
         if(in_kx(i)==k_vec(1,ikx+(iyz-1)*nkx) .and. ind_bnd(i) == n)then 
            inds(i)=n+(ikx-1)*N_bands
         end if
      end do
   end do
   
   do ikx=1,nkx_add
      if(in_kx(i) == xk_add(ikx) .and.  ind_bnd(i) == nbnd_add(ikx) ) then
         inds(i)=nk1*mm1 + nn + ikx
      end if
   end do
   
      write(*,*)'inds',i,inds(i)
end do



!allocate(inds(NModes))
!jj=0
!do i=1,NModes
!   ii=0
!   do j=1,NModes
!      if ( abs(Si(i,j)) > 0.999 ) ii=ii+1
!   end do
!!   if(ii==1)then
!      jj=jj+1
!      inds(jj)=i
!!   end if
!end do
!write(*,*)'number of selected states after cleaning= ',jj
!allocate(A(NM,jj))
!do i=1,jj
!   A(1:NM,i)=B(1:NM,inds(i))
!end do
!deallocate(B)
!allocate(B(NM,jj))
!B=A
!deallocate(A)
!
!allocate(A(Nrx*Ngt*npol,jj))
!do i=1,jj
!   A(1:Nrx*Ngt*npol,i)=C(1:Nrx*Ngt*npol,inds(i))
!end do
!deallocate(C)
!allocate(C(Nrx*Ngt*npol,jj))
!C=A
!deallocate(A)
!
!deallocate(Si)
!allocate(Si(jj,jj))
!call ZGEMM('c','n',jj,jj,nrx*ngt*npol,alpha,C,nrx*ngt*npol,C,nrx*ngt*npol,beta,Si,jj)
!
!allocate(E(jj))
!do i=1,jj
!   E(i)=in_kx(inds(i))
!end do
!deallocate(in_kx)
!allocate(in_kx(jj))
!in_kx=E
!deallocate(E)
!allocate(E(jj))
!do i=1,jj
!   E(i)=dble(ind_bnd(inds(i)))
!end do
!deallocate(ind_bnd)
!allocate(ind_bnd(jj))
!ind_bnd=int(E)
!deallocate(E)
!nmodes=jj
!mmm=jj
!
!deallocate(inds)
!   !stop
!!!!!!!!!!!!!!!!!!!!!!


  ! if(ncell==1)then
      if( nspin == 1 )then
         open(unit=13,file=TRIM(outdir)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
      else
         if(is==1)open(unit=13,file=TRIM(outdir)//'Psi_Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
         if(is==2)open(unit=13,file=TRIM(outdir)//'Psi_Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
      end if
      do j=1,NModes
         do i=1,Nrx*Ngt*npol
            write(13,'(2e25.15)')PSI_MOD(i,inds(j))
         end do
      end do
      do j=1,NModes
         write(13,*)j,in_kx(j),ind_bnd(j)
      end do
      close(13)
      
      !ic=1
      !if( nspin == 1 )then
      !   open(unit=13,file='Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
      !else
      !   if(is==1) open(unit=13,file='Bloch_up_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
      !   if(is==2) open(unit=13,file='Bloch_dw_nkyz_'//TRIM(STRINGA(iyz))//'_ncell_'//TRIM(STRINGA(imat))//'.dat',status='unknown')    
      !end if
      !do i=1,NM
      !   do j=1,nmodes
      !      write(13,*)B(i,j)
      !   end do
      !end do
      !close(13)
   
   
!   deallocate(A)
   deallocate(C)
   deallocate(ind_bnd,in_kx)
   

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
        transpose(conjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
   call SUB_DEF_Z0_GEN(1,mmm,mmm,A,Si,E) !!!call SUB_DEF_Z0(1,mmm,mmm,A,E)
   do k=1,mmm
      hkl(ikx,k)=e(k)
   end do
end do
deallocate(A)
deallocate(E)
do k=1,mmm
   if( nspin == 1 )then
      open(unit=300,file=TRIM(outdir)//'Edisp_'//TRIM(STRINGA(k))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   else
      if(is==1)  open(unit=300,file=TRIM(outdir)//'Edisp_up_'//TRIM(STRINGA(k))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
      if(is==2)  open(unit=300,file=TRIM(outdir)//'Edisp_dw_'//TRIM(STRINGA(k))//'_'//TRIM(STRINGA(iyz))//'.dat',status='unknown')
   end if
   do ikx=1,n+1
      write(300,*)dble(ikx-1-n/2)/dble(n),hkl(ikx,k)
   end do
   close(300)
end do
deallocate(hkl)


!!!!stop
NM=mmm
NMODES=NM
write(*,*)'NUMBER OF BLOCH STATES=',NMODES

!if(ncell==1)then
if( nspin == 1 )then
   open(unit=13,file=TRIM(outdir)//'H00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
else
   if(is==1) open(unit=13,file=TRIM(outdir)//'H00_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
   if(is==2) open(unit=13,file=TRIM(outdir)//'H00_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
end if
do i=1,nm
   do j=1,nm
      write(13,'(2e25.15)')HLLL(inds(i),inds(j))
   end do
end do
close(13)


!   allocate(A(nm,nm))
!   do i=1,nm
!      do j=1,nm
!         A(i,j)=dconjg(TLLL(inds(j),inds(i)))
!      end do
!   end do

   if( nspin == 1 )then
      open(unit=13,file=TRIM(outdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
   else
      if(is==1) open(unit=13,file=TRIM(outdir)//'H01_up_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
      if(is==2) open(unit=13,file=TRIM(outdir)//'H01_dw_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(imat))//'.dat',status='unknown')
   end if
   do i=1,nm
      do j=1,nm
         write(13,'(2e25.15)')TLLL(inds(i),inds(j))
      end do
   end do
   close(13)
!   deallocate(A)

deallocate(inds)

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

HL(:,:,iyz)=HLL(:,:)
TL(:,:,iyz)=TLL(:,:)

deallocate(B)

write(*,*)iyz,'order of the Hamiltonian blocks =',nmodes

end if   !end if ncell==1

if(allocated(Si)) deallocate(Si)
if(allocated(psi_add))deallocate(psi_add)
if(allocated(HLLL))deallocate(HLLL)
if(allocated(TLLL))deallocate(TLLL)
if(allocated(PSI_MOD)) deallocate(PSI_MOD)
if(allocated(HLL))deallocate(HLL)
if(allocated(TLL))deallocate(TLL)
write(*,*)
write(*,*)'END Kyz number =',iyz
write(*,*)
write(*,*)



end do !end do iyz

end do ! fine end do is
read(10,'(a)') comment

close(10)
write(*,*)' END of DFT INPUT FILE'

IF(N_SPIN/=4)THEN
   deallocate(Deeq)
ELSE
   deallocate(Deeq_so)
END IF



!stop

end subroutine read_QE_output



SUBROUTINE SUB_DEF_Z_GEN(Mi,Mf,ny,A,C,subband,U)
implicit none
integer :: ny,mi,mf
real(dp) :: subband(1:(Mf-Mi+1))
complex(dp) :: A(1:NY,1:NY),B(1:NY,1:NY),C(1:NY,1:NY),U(1:NY,1:(Mf-Mi+1))
integer :: INFO
integer, allocatable :: iwork(:),ifail(:)
complex(dp), allocatable :: work(:)
real(dp), allocatable :: rwork(:)
REAL(DP), EXTERNAL :: DLAMCH

B=C

allocate(WORK(20*ny))
allocate(RWORK(7*ny))
allocate(IWORK(5*ny))
allocate(ifail(ny))


call zhegvx(1,'V','I','U',ny,A,ny,B,ny,0.0,0.0,mi,mf,2*DLAMCH('S'),(Mf-Mi+1),subband,U,ny,WORK,20*ny,RWORK,IWORK,ifail,INFO )

deallocate(ifail)
deallocate(work)
deallocate(rwork)
deallocate(iwork)

if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: SUB_DEF_Z_GEN HAS FAILED. INFO=',INFO
   stop
endif     

END SUBROUTINE SUB_DEF_Z_GEN


SUBROUTINE SUB_DEF_Z0_GEN(Mi,Mf,ny,A,C,subband)
implicit none
integer :: ny,mi,mf
real(dp) :: subband(1:(Mf-Mi+1))
complex(dp) :: A(1:NY,1:NY),B(1:NY,1:NY),C(1:NY,1:NY),U(1:NY,1:(Mf-Mi+1))
integer :: INFO
integer, allocatable :: iwork(:),ifail(:)
complex(dp), allocatable :: work(:)
real(dp), allocatable :: rwork(:)
REAL(DP), EXTERNAL :: DLAMCH

B=C

allocate(WORK(20*ny))
allocate(RWORK(7*ny))
allocate(IWORK(5*ny))
allocate(ifail(ny))

CALL ZHEGVX (1,'N','I','U',ny,A,ny,B,ny,0.0,0.0,mi,mf,2*DLAMCH('S'), mf-mi+1, subband,U,ny,work,20*ny,rwork,iwork,ifail,info)

deallocate(ifail)
deallocate(work)
deallocate(rwork)
deallocate(iwork)

if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: SUB_DEF_Z0_GEN HAS FAILED. INFO=',INFO
   stop
endif     

END SUBROUTINE SUB_DEF_Z0_GEN

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

 
 subroutine MGS(np,nm,Q) !!!! MODIFIED GRAM SCHMIDT ALGORITHM
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
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 subroutine ortonorma(np,nm,PSI,PHI) !!!! 
 implicit none
 
 integer, intent(IN) :: np, nm
 complex(dp), intent(IN)  :: PSI(np,nm)
 complex(dp), intent(OUT) :: PHI(np,nm)
 integer :: i,m
 complex(dp), allocatable :: S(:,:),S_m05(:,:)

 
 allocate(S(NM,NM))
 call ZGEMM('c','n',NM,NM,np,alpha,PSI,np,PSI,np,beta,S,NM)

 allocate(S_m05(NM,NM))
 call A_POWER(-0.5_dp,nm,S,S_m05)

 PHI=0.0_dp
 do i=1,NM
    do m=1,NM
       PHI(1:np,i)=PHI(1:np,i)+S_m05(m,i)*PSI(1:np,m)
    end do
 end do
 deallocate(S,S_m05)
 
 
end subroutine ortonorma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine A_POWER(p,NM,A,B)  ! computes A**(p) of the hermitian matrix A
   implicit none
   
   integer, intent(in)      :: nm
   real(dp), intent(in)     :: p
   complex(dp),intent(in)   :: A (nm,nm) 
   complex(dp),intent(out)  :: B (nm,nm) 

   integer                  :: i
   complex(dp), allocatable :: U(:,:),C(:,:)
   real(dp), allocatable    :: E(:)

   allocate(U(nm,nm),C(nm,nm))
   allocate(E(nm))
   
   B=A
   
   call SUB_DEF_Z(1,NM,nm,B,E,U)

   B=0.0_dp
   do i=1,nm
      B(i,i)=E(i)**(p)
   end do
   
   call zgemm('n','c',nm,nm,nm,alpha,B,nm,U,nm,beta,C,nm)
   call zgemm('n','n',nm,nm,nm,alpha,U,nm,C,nm,beta,B,nm)
   
   deallocate(U,C)
   deallocate(E)
   
 end subroutine A_POWER
 


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
END MODULE Pseudopot_so_gen
