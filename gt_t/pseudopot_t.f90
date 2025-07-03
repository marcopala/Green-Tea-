! Copyright or © or Copr. Marco Pala (February 24, 2022)

! e-mail:  marco.pala@c2n.upsaclay.fr ;  marco.pala@cnrs.fr

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

CONTAINS
  
subroutine read_QE_output

implicit none
  
character(len=80) :: comment
integer(k15), allocatable :: miller_2D(:,:)
integer                   :: i,ikx,iyz,jyz,irr,ik,nu,xx,ii,j,l,k,m,n,mm,nn,ll,nrx0,ncell,m1,ix,jx,iy,iz
integer                   :: nm,nadd,nbnd,ngmax,n2,n3,jgt,igt,ip,im,jj,kk,iq,nksq,i0,i1,j0,j1,k0,k1

real(dp)                  :: a_1(3),a_2(3),a_3(3)
real(dp)                  :: b_1(3),b_2(3),b_3(3)
real(dp)                  :: vec(3), t0, tx, ty,tz, a0
real(dp)                  :: Ecutoff,refec,refev,tmp1,tmp2,tmp3
real(dp),    allocatable  :: E(:), KGt(:,:), Gx(:), bb_ev(:), bb_ec(:), hkl(:,:)
real(dp),    allocatable  :: xk(:,:), kread(:,:)

complex(dp), allocatable  :: el_ph_mat(:,:,:,:)
complex(dp), allocatable  :: A(:,:),B(:,:),C(:,:),D(:,:),U(:,:),Uk(:,:)
complex(dp), allocatable  :: HLL(:,:),TLL(:,:),HLLL(:,:),TLLL(:,:) 
!complex(dp), allocatable  :: dens_z(:,:,:), dens_yz(:,:,:),
complex(dp), allocatable  :: dens_x(:,:,:,:), dens_xyz(:,:,:,:)
complex(dp)               :: tmp, zdotc

real(4) :: t1,t2

external zdotc

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


t0=hbareV*hbareV/2.0_dp/m0 !eV cm**2
a0=ac ! cm

write(*,*)'ryd=',ryd
write(*,*)'t0=',t0,t0*(2.0_dp*pi/a0)**2
write(*,*)'a0 (nm)',a0*1d7

a_1=(/ ac1, 0.0_dp, 0.0_dp /)
a_2=(/ 0.0_dp, ac2, 0.0_dp /)
a_3=(/ 0.0_dp, 0.0_dp, ac3 /)
write(*,*)a_1
write(*,*)a_2
write(*,*)a_3


b_1=(/ 2.0_dp*pi/ac1, 0.0_dp, 0.0_dp /)
b_2=(/ 0.0_dp, 2.0_dp*pi/ac2, 0.0_dp /)
b_3=(/ 0.0_dp, 0.0_dp, 2.0_dp*pi/ac3 /)

write(*,*) b_1
write(*,*) b_2
write(*,*) b_3
write(*,*)
write(*,*)dot_product(a_1,b_1)/(2.0_dp*pi),dot_product(a_1,b_2),dot_product(a_1,b_3)
write(*,*)dot_product(a_2,b_1),dot_product(a_2,b_2)/(2.0_dp*pi),dot_product(a_2,b_3)
write(*,*)dot_product(a_3,b_1),dot_product(a_3,b_2),dot_product(a_3,b_3)/(2.0_dp*pi)
write(*,*)

b_1=b_1/(2.0_dp*pi/a0)!/bohr !units of 2pi/a0
b_2=b_2/(2.0_dp*pi/a0)!/bohr !units of 2pi/a0
b_3=b_3/(2.0_dp*pi/a0)!/bohr !units of 2pi/a0
write(*,*) 'b1', b_1
write(*,*) 'b2', b_2
write(*,*) 'b3', b_3

write(*,*) 'NPOL =',NPOL

write(*,*) 'NKyz points=',Nkyz

Ecutoff=ryd*Ecut

ngmax=1000
nrx=0
do n2=-ngmax,ngmax !minval(miller(1,:)),maxval(miller(1,:))
   vec(1:3)=n2*b_1(1:3)
   if(t0*vec(1)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff )nrx=nrx+1
end do
nrx=(nrx-1)
!nrx=maxval(miller(1,:))-minval(miller(1,:))+1
nry=0
do n2=-ngmax,ngmax !minval(miller(2,:)),maxval(miller(2,:))
   vec(2:3)=n2*b_2(2:3)
   if(t0*vec(2)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff )nry=nry+1
end do
nry=nry-1
write(*,*)'NRY',nry
nrz=0
do n3=-ngmax,ngmax !minval(miller(3,:)),maxval(miller(3,:))
   vec(2:3)=n3*b_3(2:3)
   if(t0*vec(3)**2*(2.0_dp*pi/a0)**2 <= 4.0_dp*Ecutoff )nrz=nrz+1
end do
nrz=nrz-1
write(*,*)'NRZ',nrz


Ngt=0
do n2=-nry/2,nry/2 !-maxval(miller(2,:)),maxval(miller(2,:))
do n3=-nrz/2,nrz/2 !maxval(miller(3,:)),maxval(miller(3,:))
vec(2:3)=n2*b_2(2:3)+n3*b_3(2:3)

if(t0*dot_product(vec(2:3),vec(2:3))*(2.0_dp*pi/a0)**2 <= 1.0_dp*Ecutoff) then ! eV
   Ngt=Ngt+1
end if

end do
end do

write(*,*)'Ngt=',Ngt

allocate(miller_2D(3,Ngt))
miller_2D=0

j=0
do n2=-nry/2,nry/2 !-maxval(miller(2,:)),maxval(miller(2,:))
do n3=-nrz/2,nrz/2 !-maxval(miller(3,:)),maxval(miller(3,:))
vec(2:3)=n2*b_2(2:3)+n3*b_3(2:3) !!!! THIS IMPLIES THAT b_1 is orthonal to the transverse plane

if(t0*dot_product(vec(2:3),vec(2:3))*(2.0_dp*pi/a0)**2 <= 1.0_dp*Ecutoff)then ! eV
!!!if(t0*vec(2)**2*(2.0_dp*pi/a0)**2 <= Ecutoff .and. t0*vec(3)**2*(2.0_dp*pi/a0)**2 <= Ecutoff )then 
   j=j+1
   miller_2D(2,j)=n2
   miller_2D(3,j)=n3

end if

end do
end do


ncell=1
nrx0=nrx/ncell
if(mod(nrx,nrx0)/=0)then
   write(*,*)'nrx0problem'
   stop
end if
write(*,*)'ncell',ncell,ncell/2
write(*,*)'NRX',nrx,'NRY',nry,'NRZ',nrz
write(*,*)'max miller_2D 2',maxval(miller_2D(2,:))
write(*,*)'max miller_2D 3',maxval(miller_2D(3,:))

allocate(Gx(nrx))
Gx=0.0_dp
do i=-nrx/2,nrx/2-1 !maxval(miller(1,:)),maxval(miller(1,:))!-1
   m1=i
   if ( m1 < 0 ) then
     m1 = m1 + nrx + 1
   else
      m1 = m1 + 1
   end if
   Gx(m1)=dble(i)*b_1(1)
end do


Dx=a_1(1)/dble(nrx)
Dy=a_2(2)/dble(nry)
Dz=a_3(3)/dble(nrz)


write(*,*)'a0/Dx',a0/Dx,'Dx',Dx
write(*,*)'a0/Dy',a0/Dy,'Dy',Dy
write(*,*)'a0/Dz',a0/Dz,'Dz',Dz

    nx=nrx-1
    ny=nry-1
    nz=nrz-1

    write(*,*)'nx+1=',nx+1
    write(*,*)'ny+1=',ny+1
    write(*,*)'nz+1=',nz+1

    allocate(Nez(Ndeltaz))
    Nez=0
    do i=1,Ndeltaz
       do iz=1,NRZ
          if(dble(iz) <= dble(i)*dble(NRZ)/dble(Ndeltaz) .and. dble(iz) > dble(i-1)*dble(NRZ)/dble(Ndeltaz)) then
             nez(i)=nez(i)+1
          end if
       end do
    end do

    allocate(Ney(Ndeltay))
    Ney=0
    do i=1,Ndeltay
       do iy=1,NRY
          if(dble(iy) <= dble(i)*dble(NRY)/dble(Ndeltay) .and. dble(iy) > dble(i-1)*dble(NRY)/dble(Ndeltay)) ney(i)=ney(i)+1
       end do
    end do

    allocate(Nex(Ndeltax))
    Nex=0
    do i=1,Ndeltax
       do ix=1,NRX
          if(dble(ix) <= dble(i)*dble(NRX)/dble(Ndeltax) .and. dble(ix) > dble(i-1)*dble(NRX)/dble(Ndeltax)) nex(i)=nex(i)+1
       end do
    end do
    
    if(sum(nex(:)) /= nrx )then
       write(*,*)'pb wth nex and nrx'
       stop
    end if
    
    allocate(kc_min(Nkyz,num_mat))
    allocate(kv_max(Nkyz,num_mat))
    write(*,*)'NPOL=',npol

    write(*,*)'Nex',Nex
    write(*,*)'Ney',Ney
    write(*,*)'Nez',Nez
    

do iyz=1,Nkyz  ! Loops over the transverse k vectors
if(k_selec(iyz))then
   write(*,*)'ikyz',iyz
   ky(iyz-((iyz-1)/nky))=k_vec(2,iyz)
   write(*,*)'iky',iyz-((iyz-1)/nky),ky(iyz-((iyz-1)/nky))    !!!! iyz = iy + (iz-1)*nky
   kz(1+((iyz-1)/nky))=k_vec(3,iyz)
   write(*,*)'ikz',1+((iyz-1)/nky),kz(1+((iyz-1)/nky)) !!!! iyz = iy + (iz-1)*nky
end if
end do



do iyz=1,Nkyz
if(k_selec(iyz))then
   
do im=1,num_mat
   nm=NM_mat(im)

   allocate(ULCBB(iyz,im)%H(Nrx*NGt*npol,NM_mat(im)))
   allocate(PSIBB(iyz,im)%H(Nrx*NGt*npol,NM_mat(im)))
   NM=NM_mat(im)
   allocate(A(nrx*ngt*npol,NM))
   
   if(magnetic)then
   open(unit=13,file=TRIM(inputdir)//'Psi_Bloch_'//TRIM(updw)//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
   write(*,*)'reading ',TRIM(inputdir)//'Psi_Bloch_'//TRIM(updw)//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   else
   open(unit=13,file=TRIM(inputdir)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
   write(*,*)'reading ',TRIM(inputdir)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   endif
   do j=1,NM
      do ip=1,npol
         do jgt=1,Ngt
            do ix=1,nrx   
               read(13,'(2e25.15)')tmp1,tmp2
               A(jgt+(ix-1)*ngt+(ip-1)*Ngt*nrx,j)=dcmplx(tmp1,tmp2)
            end do
         end do
      end do
   end do
  ! if(phonons)then
      allocate( ind_kx( iyz,im)%n(NM) )
      ind_kx( iyz,im)%n=0.0_dp
      allocate( ind_bnd(iyz,im)%i(NM) )
      ind_bnd(iyz,im)%i=0
      
      do j=1,NM
         read(13,*)i,ind_kx(iyz,im)%n(j),ind_bnd(iyz,im)%i(j)
      end do
  ! end if
   close(13)
  

   !!! debugging line
!!$if(allocated(C))deallocate(C)
!!$allocate(C(NM,NM))
!!$C=0.0_dp
!!$call ZGEMM('c','n',NM,NM,nrx*ngt*npol,alpha,A,nrx*ngt*npol,A,nrx*ngt*npol,beta,C,NM)
!!$do i=1,nm
!!$   do j=1,nm
!!$      write(329,*)i,j,abs(C(i,j))
!!$   end do
!!$   write(329,*)
!!$end do
!!$allocate(E(nm))
!!$call sub_def_Z0(1,NM,NM,C,E)
!!$do i=1,nm
!!$   write(330,*)i,e(i)
!!$end do
!!$deallocate(E)
!!$deallocate(C)
   !!end debugging lines
   
   allocate(psipsi(iyz,im)%H(NM_mat(im),NM_mat(im)))
   do i=1,nm
      do j=1,nm
         psipsi(iyz,im)%H(i,j) = dot_product(A(:,j),A(:,i))
      end do
   end do

!!$   do i=1,nm
!!$      do j=1,nm
!!$         write(19,*)i,j,abs(psipsi(iyz,im)%H(i,j))
!!$      end do
!!$       write(19,*)
!!$    end do
!!$    stop

   
   allocate(Si_p05(iyz,im)%H(NM_mat(im),NM_mat(im)))
   allocate(Si_m05(iyz,im)%H(NM_mat(im),NM_mat(im)))
   allocate(Si(iyz,im)%H(NM_mat(im),NM_mat(im)))
   allocate(C(NM,NM))
   call ZGEMM('c','n',NM,NM,nrx*ngt*npol,alpha,A,nrx*ngt*npol,A,nrx*ngt*npol,beta,C,NM)
   Si(iyz,im)%H=C  !(C+transpose(dconjg(C)))/2.0_dp


!do i=1,NM
!   do j=1,NM
!      write(700+iyz,*)i,j,abs(C(i,j))
!   end do
!   write(700+iyz,*)
!end do
!close(700+iyz)
   
   deallocate(C)
   
   call A_POWER( 0.5_dp,nm,Si(iyz,im)%H,Si_p05(iyz,im)%H)
   call A_POWER(-0.5_dp,nm,Si(iyz,im)%H,Si_m05(iyz,im)%H)

   ULCBB(iyz,im)%H=A   !!! this is the matrix containing the URBFs
!!!!!!!!!
   if(phonons)then
      PSIBB(iyz,im)%H=A !!! I am stocking the psi to rebuild the unk
      do ip=1,npol
         do j=1,NM  
            do ix=1,Nrx
               PSIBB(iyz,im)%H(1+(ix-1)*NGt+(ip-1)*Ngt*nrx:ngt+(ix-1)*NGt+(ip-1)*Ngt*nrx,j)=&
                    A(1+(ix-1)*ngt+(ip-1)*Ngt*nrx:ngt+(ix-1)*ngt+(ip-1)*Ngt*nrx,j)&
                    *exp(-dcmplx(0.0_dp,1.0_dp)*ind_kx(iyz,im)%n(j)*2.0_dp*pi*dble(ix-1)/dble(Nrx))
            end do
         end do
      end do
   end if
!!!!!!!!
   
deallocate(A)
end do
end if
end do


do iyz=1,Nkyz
if(k_selec(iyz))then
do im=1,num_mat
   nm=NM_mat(im)
   write(*,*)iyz,'im',im,nm
   allocate(HL(iyz,im)%H(NM_mat(im),NM_mat(im)))
   if(magnetic)then
   open(unit=13,file=TRIM(inputdir)//'H00_'//TRIM(updw)//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
   write(*,*)'reading ',TRIM(inputdir)//'H00_'//TRIM(updw)//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   else
   open(unit=13,file=TRIM(inputdir)//'H00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
   write(*,*)'reading ',TRIM(inputdir)//'H00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'  
   end if
   do i=1,nm
      do j=1,nm
         read(13,'(2e25.15)')tmp1,tmp2
         if(i>=j)            HL(iyz,im)%H(i,j)=dcmplx(tmp1,tmp2)
      end do
   end do
   do i=1,nm
      do j=1,nm
           if(i<j) HL(iyz,im)%H(i,j)=conjg(HL(iyz,im)%H(j,i))
      end do
   end do
   close(13)
   !HL(iyz,im)%H=(HL(iyz,im)%H+transpose(dconjg(HL(iyz,im)%H)))/2.0_dp
   HL(iyz,im)%H = HL(iyz,im)%H + off_set(im)*Si(iyz,im)%H  !!! offset of the material
   write(*,*)im,'offset=',off_set(im)
   allocate(TL(iyz,im)%H(NM_mat(im),NM_mat(im)))   !!! H10
   nm=NM_mat(im)
   if(magnetic)then
      open(unit=13,file=TRIM(inputdir)//'H01_'//TRIM(updw)//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
      write(*,*)'reading ',TRIM(inputdir)//'H01_'//TRIM(updw)//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   else
      open(unit=13,file=TRIM(inputdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
      write(*,*)'reading ',TRIM(inputdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   endif
   allocate(A(nm,nm))
   do i=1,nm
      do j=1,nm
         read(13,'(2e25.15)')tmp1,tmp2
         A(i,j)=dcmplx(tmp1,tmp2)
      end do
   end do
   close(13)
   TL(iyz,im)%H=transpose(dconjg(A))!!! this transforms H01 into H10
   deallocate(A)
end do

if(num_het > 0)then

   l=num_mat
   do im=1,num_het
      l=l+1 
      write(*,*)'het n. =',im,mat_l(im),', NM_mat lft =',NM_mat(mat_l(im))
      write(*,*)'het n. =',im,mat_r(im),', NM_mat rgt =',NM_mat(mat_r(im))
      allocate(A(NM_mat(mat_l(im)),NM_mat(mat_r(im))))
      open(unit=13,file=TRIM(inputdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nhet_'//TRIM(STRINGA(im))//'.dat',status='unknown')
      do i=1,NM_mat(mat_l(im))
         do j=1,NM_mat(mat_r(im))
            read(13,'(2e25.15)')tmp1,tmp2
            A(i,j)=dcmplx(tmp1,tmp2)
         end do
      end do
      close(13)
      allocate(TL(iyz,l)%H(NM_mat(mat_r(im)),NM_mat(mat_l(im)))) !!! H10 for the HJ
      TL(iyz,l)%H=transpose(dconjg(A)) !!!! H01 -> H10
      deallocate(A)
   end do

end if
   

n=60
allocate(bb_ev(n+1),bb_ec(n+1))

do im=1,num_mat
   NM=nm_mat(im)
   M=min(nband_val(im)+10,nm)
   allocate(hkl(n+1,M))
   
   allocate(E(M))
   allocate(A(NM,NM),HLLL(NM,NM),TLLL(NM,NM),C(NM,NM))
   HLLL=HL(iyz,im)%H
   TLLL=TL(iyz,im)%H
   C=Si(iyz,im)%H
   

!do i=1,NM
!   do j=1,NM
!      write(700+iyz,*)i,j,abs(C(i,j))
!   end do
!   write(700+iyz,*)
!end do
!close(700+iyz)

   do ikx=1,n+1
      A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
           transpose(dconjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
      call SUB_DEF_Z0_GEN(1,M,NM,A,C,E)
      do k=1,M
         hkl(ikx,k)=e(k)
      end do
      if( nband_val(im)>0)     bb_ev(ikx)=E(nband_val(im))
      bb_ec(ikx)=E(nband_val(im)+1)
   end do
   deallocate(A,HLLL,TLLL,C)
   deallocate(E)
   do k=1,min(10,M)
      if(nband_val(im)>0)then
         ii=nband_val((im))-1*min(10,M)/2+k
      else
         ii=nband_val((im))+k
      end if

      open(unit=300+k,file=TRIM(outdir)//'Edisp_'//TRIM(STRINGA(ii))//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
      do ikx=1,n+1
         write(300+k,*)dble(ikx-1-n/2)/dble(n),hkl(ikx,ii)
      end do
      close(300+k)
   end do


kc_min(iyz,im) = 0.0_dp
kv_max(iyz,im) = 0.0_dp
refec= 1.0d3
refev=-1.0d3
do ikx=n/2+1,n
   if(bb_ec(ikx) <= refec) then
      refec=bb_ec(ikx)
      kc_min(iyz,im)=dble(ikx-1-n/2)/dble(n)
   end if
   if(bb_ev(ikx) >= refev) then
      refev=bb_ev(ikx)
      kv_max(iyz,im)=dble(ikx-1-n/2)/dble(n)
   end if
end do
write(*,*)iyz,im,'kc_min',kc_min(iyz,im),refec
write(*,*)iyz,im,'kv_max',kv_max(iyz,im),refev
deallocate(hkl)

end do

deallocate(bb_ec,bb_ev)


if(allocated(KGt))deallocate(KGt)
allocate(KGt(4,1*Ngt))
KGt=0.0_dp
do j=1,Ngt
   kgt(2,j)=miller_2D(2,j)*b_2(2)+miller_2D(3,j)*b_3(2) + k_vec(2,iyz)
   kgt(3,j)=miller_2D(2,j)*b_2(3)+miller_2D(3,j)*b_3(3) + k_vec(3,iyz)
   kgt(4,j)=(kgt(2,j)**2+kgt(3,j)**2)
end do

if(.not.allocated(kgt_kyz))then
   allocate(kgt_kyz(4,1:1*Ngt,Nkyz))
   kgt_kyz=0
end if

KGt_kyz(1:4,1:Ngt,iyz)=KGt(1:4,1:1*Ngt)

end if
end do ! endo do iyz


if(.not. onlyT .or. phonons )then

if(.not.allocated(Uk)) allocate(Uk(NGt*npol,(nry)*(nrz)))

do im=1,num_mat

   write(*,*)
   write(*,*)'material n.',im
   write(*,*)
   
   NM=NM_mat(im)

   do iyz=1,Nkyz
      if(k_selec(iyz))then
      
     !!! allocate(U_psi(iyz,im)%K(1:NM_mat(im)*NM_mat(im),1:(Ndeltay+1),1:(Ndeltaz+1),1:Nrx))
      
      
!!!! INTERPOLATION ON THE COARSE GRID
     
      t1=SECNDS(0.0)
      
      do iy=1,NRY
         do iz=1,NRZ
            j=iy+(iz-1)*(NRY)
            Uk(1:NGt,j)=exp(dcmplx(0.0_dp,-1.0_dp)*KGt_kyz(2,1:NGt,iyz)*2.0_dp*pi/a0*dble(iy)*Dy+&
                            dcmplx(0.0_dp,-1.0_dp)*KGt_kyz(3,1:NGt,iyz)*2.0_dp*pi/a0*dble(iz)*Dz )/sqrt(dble((nry)*(nrz)))
         end do
      end do
      if(npol>1)      Uk(1+NGt:npol*NGt,1:NRY*NRZ)=Uk(1:NGt,1:NRY*NRZ)


     ! allocate(D(1:NM_mat(im)*NM_mat(im),1:(Ndeltay+1),1:(Ndeltaz+1),1:Nrx))!proxy of U_psi
      allocate(U(Nrx*NGt*npol,NM_mat(im)))!proxy of ULCBB
      U=ULCBB(iyz,im)%H
if (Jz) then
   allocate(dzpsipsi(iyz,im)%H(nrz,NM_mat(im)*NM_mat(im)))
   dzpsipsi(iyz,im)%H=0.0_dp
!allocate(dxpsipsi(iyz,im)%H(nrz,NM_mat(im)*NM_mat(im)))
!dxpsipsi(iyz,im)%H=0.0_dp
end if
      call omp_set_num_threads(Nomp_PP)
      allocate(dens_xyz(nm,nrx,ndeltay+1,ndeltaz+1))
      dens_xyz=0.0_dp
     
      !allocate(C(NM,(nry)*(nrz)))
      !allocate(D(NM,(nry)*(nrz)))

      !$omp parallel default(none) private(ix,iy,iz,ip,i,j,k,jj,kk,jgt,A,B,tmp1,tmp2) &
      !$omp shared(nm,nrx,nry,nrz,ndeltax,ndeltay,nex,ndeltaz,npol,ngt,xx,ney,nez,Uk,U,dens_xyz,im,iyz,Jz,dzpsipsi)
      allocate(A(NM,(nry)*(nrz)))
      allocate(B(NM,ngt*npol))
      
      !$omp do
      do ix=1,Nrx
         write(*,*) 'ix =',ix
         do ip=1,npol
            do jgt=1,Ngt
               do i=1,NM
                  B(i,jgt+(ip-1)*ngt)=(conjg(U(jgt+(ix-1)*Ngt+(ip-1)*nrx*ngt,i)))
               end do
            end do
         end do
         call ZGEMM('n','n',NM,(nry)*(nrz),NGt*npol,alpha,B,NM,Uk(1:Ngt*npol,1:nry*nrz),NGt*npol,beta,A,NM)

!!$         do iz=1,nrz
!!$            tmp1=0.0_dp
!!$            tmp2=0.0_dp
!!$            do i=1,NM
!!$               if(ind_bnd(iyz,im)%i(i) <= nband_val(im))then
!!$                  tmp1=tmp1+dble(A(i,nry/2+(iz-1)*(nry)))
!!$                  tmp2=tmp2+dimag(A(i,nry/2+(iz-1)*(nry)))
!!$               end if
!!$            end do
!!$            write(2000+100*im+iyz,*)ix,iz,sqrt(tmp1**2+tmp2**2)
!!$         end do
!!$         write(2000+100*im+iyz,*)
!!$         do iz=1,nrz
!!$            tmp1=0.0_dp
!!$            tmp2=0.0_dp
!!$            do i=1,NM
!!$               if(ind_bnd(iyz,im)%i(i) > nband_val(im))then
!!$                  tmp1=tmp1+dble(A(i,nry/2+(iz-1)*(nry)))
!!$                  tmp2=tmp2+dimag(A(i,nry/2+(iz-1)*(nry)))
!!$               end if
!!$            end do
!!$            write(3000+100*im+iyz,*)ix,iz,sqrt(tmp1**2+tmp2**2)
!!$         end do
!!$         write(3000+100*im+iyz,*)
         
      if (Jz) then      
         do iy=1,nry
               do i=1,NM
                  do j=1,NM
                     
                     do iz=1,nrz-1
                        dzpsipsi(iyz,im)%H(iz,i+(j-1)*NM)=dzpsipsi(iyz,im)%H(iz,i+(j-1)*NM) + &
                             (conjg( A(j,iy+(iz-1)*(nry)) )*(A(i,iy+(iz)*(nry)) - A(i,iy+(iz-1)*(nry))) &
                             - conjg( A(j,iy+(iz)*(nry)) - A(j,iy+(iz-1)*(nry)) )*A(i,iy+(iz-1)*(nry)))
                     end do
                     iz=nrz
                     dzpsipsi(iyz,im)%H(iz,i+(j-1)*NM)=dzpsipsi(iyz,im)%H(iz,i+(j-1)*NM) + &
                          (conjg( A(j,iy+(iz-1)*(nry)) )*(A(i,iy+(0)*(nry)) - A(i,iy+(iz-1)*(nry))) &
                          - conjg( A(j,iy+(0)*(nry)) - A(j,iy+(iz-1)*(nry)) )*A(i,iy+(iz-1)*(nry)))
                     
               end do
            end do
         end do
      end if
!         if( ix == 1 ) C = A
!         if( ix == 2 ) D = A
         
           iz=0
           do k=1,Ndeltaz
              do kk=1,Nez(k)
                 iz=iz+1
           iy=0
           do j=1,Ndeltay
              do jj=1,Ney(j)
                 iy=iy+1
                 dens_xyz(1:nm,ix,j,k)=dens_xyz(1:nm,ix,j,k)+A(1:nm,iy+(iz-1)*(nry))/dble(Ney(j)*Nez(k))
              end do
           end do
            end do
         end do 

         do k=1,Ndeltaz+1
            dens_xyz(1:nm,ix,Ndeltay+1,k)=dens_xyz(1:nm,ix,1,k)
         end do
         do j=1,Ndeltay+1
            dens_xyz(1:nm,ix,j,Ndeltaz+1)=dens_xyz(1:nm,ix,j,1)
         end do
         
      end do
      !$omp end do
      
      deallocate(A)
      deallocate(B)
      !$omp end parallel    
      deallocate(U)

    !  do iy=1,nry
    !     do iz=1,nrz
    !        do i=1,NM
    !           do j=1,NM
    !              
    !              dxpsipsi(iyz,im)%H(iz,i+(j-1)*NM)=dxpsipsi(iyz,im)%H(iz,i+(j-1)*NM) + &
    !                   (  conjg(D(j,iy+(iz-1)*(nry)) ) * (D(i,iy+(iz-1)*(nry)) - C(i,iy+(iz-1)*(nry))) &
    !                   -conjg(D(j,iy+(iz-1)*(nry)) - C(j,iy+(iz-1)*(nry)) ) * D(i,iy+(iz-1)*(nry)) )
    !              
    !           end do
    !        end do
    !     end do
    !  end do
    !  deallocate(C,D)

!!$!$omp parallel default(none) private(ix,iy,iz,i,j,ii,xx,k,kk) &
!!$!$omp shared(iyz,im,nm,nrx,nry,nrz,nex,nez,ndeltax,ndeltaz, C, dzpsipsi)
!!$   
!!$!$omp do
!!$do iy=1,nry
!!$         write(*,*) 'iy =',iy
!!$do i=1,NM
!!$   do j=1,NM
!!$ix=0
!!$do xx=1,ndeltax
!!$   do ii=1,nex(xx)
!!$      ix=ix+1
!!$         iz=0
!!$         do k=1,ndeltaz
!!$            do kk=1,nez(k)
!!$               iz=iz+1
!!$               if(iz<=nrz)      dzpsipsi(iyz,im)%M(i+(j-1)*NM,xx,k)=dzpsipsi(iyz,im)%M(i+(j-1)*NM,xx,k) + &
!!$                    conjg(C(j,ix+(iy-1)*nrx+(iz-1)*nrx*nry))*C(i,ix+(iy-1)*nrx+(iz)*nrx*nry)/dble(nex(xx)*nez(k)*nry)
!!$            end do
!!$         end do
!!$      end do
!!$   end do
!!$end do
!!$end do
!!$end do
!!$      !$omp end do
!!$      !$omp end parallel 
!!$
!!$deallocate(C)

allocate(dens_x(NM_mat(im), Ndeltax, Ndeltay+1, Ndeltaz+1))
do i=1,nm
   do iz = 1, ndeltaz+1
      do iy = 1, ndeltay+1
         ix=0
         do k = 1, ndeltax
            dens_x(i,k,iy,iz)=0.0_dp
            do kk=1,Nex(k)
               ix=ix+1
               dens_x(i,k,iy,iz)=dens_x(i,k,iy,iz)+dens_xyz(i,ix,iy,iz)/dble(Nex(k))
            end do
         end do
      end do
   end do
end do
allocate(psi_charge(iyz,im)%K(NM_mat(im)*NM_mat(im), Ndeltax, Ndeltay+1, Ndeltaz+1))
psi_charge(iyz,im)%K=0.0_dp

do i=1,nm
   do j=1,nm
do iz = 1, ndeltaz+1
do iy = 1, ndeltay+1
do ix = 1, ndeltax
      psi_charge(iyz,im)%K(i+(j-1)*NM, ix, iy, iz) =  conjg(dens_x(j,ix,iy,iz)) * dens_x(i,ix,iy,iz)
   end do
!do iz = 1, ndeltaz+1
!do iy = 1, ndeltay+1
!   ix=0
!   do k = 1, ndeltax
!      do kk=1,Nex(k)
!         ix=ix+1
!         psi_charge(iyz,im)%K(i+(j-1)*NM, k, iy, iz) = psi_charge(iyz,im)%K(i+(j-1)*NM, k, iy, iz) + &  
!              ( conjg(dens_xyz(j,ix,iy,iz))*dens_xyz(i,ix,iy,iz) )/dble(Nex(k))
!      end do
!   end do
end do
end do
end do
end do

do iy = ndeltay/2, ndeltay/2
do ix = 1, Ndeltax
do iz = 1, ndeltaz+1
   tmp1=0.0_dp
   tmp2=0.0_dp
   do i=1,NM*NM
      tmp1=tmp1 + dble(psi_charge(iyz,im)%K(i, ix, iy, iz))
      tmp2=tmp2 + dimag(psi_charge(iyz,im)%K(i, ix, iy, iz))
   end do
   write(4000+100*im+iyz,*)ix,iz,sqrt(tmp1**2+tmp2**2)
end do
write(1000+100*im+iyz,*)
end do
end do
write(*,*)'interpolation done'
!!stop
!do i=1,nm
!do j=1,nm
!do iz = 1, ndeltaz
!   do iy = 1, ndeltay
!      do ix = 1, ndeltax
!         dzpsipsi(iyz,im)%H(i+(j-1)*NM,iz)=dzpsipsi(iyz,im)%H(i+(j-1)*NM,iz) + &
!              conjg(dens_x(j,ix,iy,iz))*dens_x(i,ix,iy,iz+1)/dble(ndeltay*ndeltax)
!      end do
!   end do   
!end do
!do iz = 1, ndeltaz
!   do iy = 1, ndeltay
!      do ix = 1, ndeltax-1
!         dxpsipsi(iyz,im)%H(i+(j-1)*NM,iz)=dxpsipsi(iyz,im)%H(i+(j-1)*NM,iz) + &
!              conjg(dens_x(j,ix,iy,iz))*dens_x(i,ix+1,iy,iz)/dble(ndeltay*ndeltax)
!      end do
!      dxpsipsi(iyz,im)%H(i+(j-1)*NM,iz)=dxpsipsi(iyz,im)%H(i+(j-1)*NM,iz) + &
!           conjg(dens_x(j,ndeltax,iy,iz))*dens_x(i,1,iy,iz)/dble(ndeltay*ndeltax)
!   end do   
!end do
!end do
!end do

deallocate(dens_xyz)
deallocate(dens_x)



!!$!$omp parallel default(none) private(ix,iy,iz,ip,i,j,jgt,dens_z,dens_yz,A,B,C) &
!!$!$omp shared(iyz,im,ney,nez,nm,nrx,nry,nrz,ndeltay,ndeltaz,dx,dy,dz,npol,ngt,Uk,U,D)
!!$
!!$
!!$allocate(dens_z(NM*NM,(nry),Ndeltaz+1))
!!$allocate(dens_yz(NM*NM,Ndeltay+1,Ndeltaz+1))
!!$
!!$allocate(A(NM,(nry)*(nrz)))
!!$allocate(B(NM,ngt*npol))
!!$allocate(C(NM*NM,(nry)*(nrz)))
!!$
!!$!$omp do 
!!$do ix=1,Nrx
!!$write(*,*) 'ix =',ix
!!$   do ip=1,npol
!!$      do jgt=1,Ngt
!!$         do i=1,NM
!!$            !!B(i,jgt+(ip-1)*ngt)=(dconjg(ULCBB(iyz,im)%H(jgt+(ix-1)*Ngt+(ip-1)*nrx*ngt,i)))
!!$            B(i,jgt+(ip-1)*ngt)=(dconjg(U(jgt+(ix-1)*Ngt+(ip-1)*nrx*ngt,i)))
!!$         end do
!!$      end do
!!$   end do
!!$   call ZGEMM('n','n',NM,(nry)*(nrz),NGt*npol,alpha,B,NM,Uk(1:Ngt*npol,1:nry*nrz),NGt*npol,beta,A,NM)
!!$
!!$   do i=1,NM
!!$      do j=1,NM
!!$         C(i+(j-1)*NM,:)=dconjg(A(i,:))*A(j,:)
!!$      end do
!!$   end do
!!$
!!$   dens_z=0.0_dp
!!$   dens_yz=0.0_dp
!!$   
!!$   iz=1
!!$   do iy=1,NRY
!!$      dens_z(1:NM*NM,iy,iz)=C(1:NM*NM,iy+(iz-1)*(nry))
!!$   end do
!!$   dens_z(1:NM*NM,1:nry,Ndeltaz+1)=dens_z(1:NM*NM,1:nry,iz)
!!$   iz=0
!!$   do i=1,Ndeltaz-1
!!$      do j=1,Nez(i)
!!$         iz=iz+1
!!$         do iy=1,NRY
!!$            dens_z(1:NM*NM,iy,i+1)=dens_z(1:NM*NM,iy,i+1)+C(1:NM*NM,iy+(iz-1)*(nry))/dble(2*Nez(i))
!!$         end do
!!$      end do
!!$   end do
!!$   iz=Nez(1)
!!$   do i=2,Ndeltaz
!!$      do j=1,Nez(i)
!!$         iz=iz+1
!!$         do iy=1,NRY
!!$            dens_z(1:NM*NM,iy,i)=dens_z(1:NM*NM,iy,i)+C(1:NM*NM,iy+(iz-1)*(nry))/dble(2*Nez(i))
!!$         end do
!!$      end do
!!$   end do
!!$   
!!$   iy=1
!!$   dens_yz(1:NM*NM,iy,1:Ndeltaz+1)=dens_z(1:NM*NM,iy,1:Ndeltaz+1)
!!$   dens_yz(1:NM*NM,Ndeltay+1,1:Ndeltaz+1)=dens_yz(1:NM*NM,iy,1:Ndeltaz+1)
!!$   iy=0
!!$   do i=1,Ndeltay-1
!!$      do j=1,Ney(i)
!!$         iy=iy+1
!!$         dens_yz(1:NM*NM,i+1,:)=dens_yz(1:NM*NM,i+1,:)+dens_z(1:NM*NM,iy,:)/dble(2*Ney(i))
!!$      end do
!!$   end do
!!$   iy=Ney(1)
!!$   do i=2,Ndeltay
!!$      do j=1,Ney(i)
!!$         iy=iy+1
!!$         dens_yz(1:NM*NM,i,:)=dens_yz(1:NM*NM,i,:)+dens_z(1:NM*NM,iy,:)/dble(2*Ney(i))
!!$      end do
!!$   end do
!!$   
!!$   !$omp critical 
!!$   !U_psi(iyz,im)%K(1:NM*NM,1:Ndeltay+1,1:Ndeltaz+1,ix)=dens_yz(1:NM*NM,1:Ndeltay+1,1:Ndeltaz+1)
!!$   D(1:NM*NM,1:Ndeltay+1,1:Ndeltaz+1,ix)=dens_yz(1:NM*NM,1:Ndeltay+1,1:Ndeltaz+1)
!!$   !$omp end critical
!!$   
!!$end do !end do ix
!!$
!!$!$omp end do
!!$
!!$deallocate(A,B,C)
!!$deallocate(dens_z,dens_yz)
!!$!$omp end parallel
!!$
!!$
!!$      
!!$allocate(U_psi(iyz,im)%K(1:NM_mat(im)*NM_mat(im),1:(Ndeltay+1),1:(Ndeltaz+1),1:Nrx))      
!!$U_psi(iyz,im)%K=D
!!$
!!$deallocate(U,D)
!!$
t2=SECNDS(t1)
WRITE(*,*)'Time spent to compute the interpolation (s)',t2

 
if(iyz==1)then
   allocate(E(NM))
   allocate(A(NM,NM),HLL(NM,NM),TLL(NM,NM),C(NM,NM))
   HLL=HL(iyz,im)%H
   TLL=TL(iyz,im)%H
   C=SI(iyz,im)%H

   A=HLL+TLL*exp(cmplx(0.0_dp,1.0_dp)*kc_min(iyz,im)*2.0_dp*pi)+&
        transpose(dconjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*kc_min(iyz,im)*2.0_dp*pi)
   call SUB_DEF_Z0_GEN(1,NM,NM,A,C,E) !!!call SUB_DEF_Z(1,NM,NM,A,E,B)
   ref_ec(im)=E(nband_val(im)+1)
   write(*,*)'ref_ec at the first Kyz = ',nband_val(im)+1,ref_ec(im)
if (nband_val(im)>0)then
   A=HLL+TLL*exp(cmplx(0.0_dp,1.0_dp)*kv_max(iyz,im)*2.0_dp*pi)+&
        transpose(dconjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*kv_max(iyz,im)*2.0_dp*pi)
   call SUB_DEF_Z0_GEN(1,NM,NM,A,C,E) !!!call SUB_DEF_Z(1,NM,NM,A,E,B)
   ref_ev(im)=E(nband_val(im))
end if
   write(*,*)'ref_ev at the first Kyz = ',nband_val(im),ref_ev(im)

   E_GAP=ref_ec(im)-ref_ev(im)
   write(*,*)'E_GAP at the first Kyz = ',e_gap
   
   deallocate(A,HLL,TLL,C)
   deallocate(E)
end if

write(*,*)'End ikyz =',iyz
end if
end do ! end do iyz

do iyz=1,Nkyz
   if(k_selec(iyz))then
      allocate(E(NM))
      allocate(A(NM,NM),HLL(NM,NM),TLL(NM,NM),C(NM,NM))
      HLL=HL(iyz,im)%H
      TLL=TL(iyz,im)%H
      C=SI(iyz,im)%H
      
      A=HLL+TLL*exp(cmplx(0.0_dp,1.0_dp)*kc_min(iyz,im)*2.0_dp*pi)+&
           transpose(dconjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*kc_min(iyz,im)*2.0_dp*pi)
      call SUB_DEF_Z0_GEN(1,NM,NM,A,C,E) 
      ref_ec(im)=E(nband_val(im)+1)
      write(*,*)'ref_ec at the first Kyz =',nband_val(im)+1,ref_ec(im)
      if (nband_val(im)>0)then
         A=HLL+TLL*exp(cmplx(0.0_dp,1.0_dp)*kv_max(iyz,im)*2.0_dp*pi)+&
              transpose(dconjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*kv_max(iyz,im)*2.0_dp*pi)
         call SUB_DEF_Z0_GEN(1,NM,NM,A,C,E) 
         ref_ev(im)=E(nband_val(im))
      end if
      write(*,*)'ref_ev at the first Kyz =',nband_val(im),ref_ev(im)
      
      E_GAP=ref_ec(im)-ref_ev(im)
      write(*,*)'E_GAP at the first Kyz =',e_gap
      
      deallocate(A,HLL,TLL,C)
      deallocate(E)
      exit
   end if
end do

if(phonons)then
write(*,*)
write(*,*)'Computing the el-ph matrix'
t1=SECNDS(0.0)

if(dfpt)then 
!!! this part implements the use of el-ph matrix elements computed by QE within the DFPT. This is a work in progress, its use is for developpers and is not yet completely tested. 

write(*,*)
write(*,*)'Reading the DFPT file ' ,TRIM(input_file_DFPT)
write(*,*)
write(*,*)'WARNING: this is a work in progress, still to be tested' 
write(*,*)

open(unit=13,file=TRIM(input_file_DFPT),status='unknown')
read(13,'(A)') comment !'nbnd='
read(13,'(I)')nbnd
write(*,*)comment, nbnd
read(13,'(A)') comment !'nmodes='
read(13,'(I)')nqmodes
write(*,*)comment, nqmodes
read(13,'(A)') comment !'nqs='
read(13,'(I)')nqs
write(*,*)comment, nqs
nadd=nqs

allocate(x_q(3,nqs))
allocate(omega_q(nqs,nqmodes))
allocate(el_ph_mat(nqs, nqmodes, nbnd, nbnd))
el_ph_mat=0.0_dp


read(13,'(a)') comment ! 'q-points:'
write(*,*) comment
read(13,*) ((x_q(i, j), i = 1, 3), j = 1, nqs)
write(*,*) ((x_q(i, j), i = 1, 3), j = 1, nqs)

x_q(2, 1:nqs)=x_q(2, 1:nqs)/ac1*ac2
x_q(3, 1:nqs)=x_q(3, 1:nqs)/ac1*ac3

!!! dble check
if(mod(nqs,nkyz) /= 0)then
   write(*,*)'problem with the list of q or k vectors'
   write(*,*)nqs,nqx,nkyz
   stop
end if
nqx=nqs/nkyz
write(*,*)
write(*,*)'NQx=',nqx
write(*,*)'NQyz=',nkyz
write(*,*)
!!!

!!! CHECK that the q_yz correspond to the k_yz
l=0
do jyz=1,nkyz
do jx=1,nqx
   iq=jx+(jyz-1)*nqx
!   write(*,*)iq,x_q(2:3,iq)
   do iyz=1,nkyz
      IF ( ABS( x_q(2,iq)-k_vec(2,iyz) ) < 1.0E-3_DP .AND. &
           ABS( x_q(3,iq)-k_vec(3,iyz) ) < 1.0E-3_DP ) THEN
         ind_q( jx, iyz ) = iq
         !write(*,*)iq, jx, iyz
         l=l+1
      END IF
   end do
end do
end do
if( l /= nqs)then
   write(*,*)'inconsistency between x_q and k_vec',l,nqs,nqx*nkyz
   stop
end if

do jyz=1,nkyz
do jx=1,nqx
   
   iq=jx+(jyz-1)*nqx ! it supposes that iq = jx+(jyz-1)*nqx
   
   read(13,'(A)') comment !'nksq='
   read(13,'(I)')nksq
   write(*,*)comment, nksq

   read(13, '(a3,3e20.10)') comment, tmp1, tmp2, tmp3! x_q(1:3, iq)
   write(*,*) comment, tmp1, tmp2, tmp3! x_q(1:3,iq)
   read(13, '(a16,i6)') comment,i
   write(*,*)  comment,i
   if (i /= iq) stop  
   write(*,*)'iq=',i
   read(13,'(a)') comment ! list of k points:
   write(*,*) comment

   if(dot_product(x_q(:,iq),x_q(:,iq))>1.0d-16)then
      allocate(xk(4,2*nksq))
      do nn=1,nksq
         READ(13,'(I4,4e20.10)')i,xk(1:4,2*nn-1)
      end do
      READ(13,'(a)') comment !list of k+q points:
      write(*,*) comment
      do nn=1,nksq
         READ(13,'(I4,4e20.10)')i,xk(1:4,2*nn)
      end do
   else
      allocate(xk(4,nksq))
      READ(13,*)
      do nn=1,nksq
         READ(13,'(I4,4e20.10)')i,xk(1:4,nn)
      end do
   end if

   read(13,'(a)') comment !squared hbar omega (Ry)
   write(*,*) comment
   do ll=1,nqmodes
      read(13,*) i,omega_q(iq,ll) ! ! squared hbar omega (Ry)
      write(*,*) i,omega_q(iq,ll)
      !if (  omega_q(iq,ll) > 0.0_dp ) then
         omega_q(iq,ll)=ryd*dsqrt(abs(omega_q(iq,ll))) ! hbar omega in eV
         write(*,*) 'omega_q',iq,ll,omega_q(iq,ll)
      !end if
   end do

   read(13,'(a)') comment !matrix elements (Ryd)
   write(*,*)comment
   read(13,*)i
   write(*,*)'ios=',i
   if( i/= 0 )then
      write(*,*)'pb wt ios'
      stop
   end if
   READ(13,'(A)')comment !reading k-points.dat
   write(*,*)comment
   
   OPEN( unit = 11, file = TRIM(inputdir)//'kpoints.dat', status = 'unknown')
   READ(11,*) nkread
   WRITE(*,*) nkread
   ALLOCATE ( kread(3,nkread) )
   DO ii = 1, nkread
      READ(11,*)  kread(1:3,ii)
   END DO
   CLOSE(11)
   
   READ(13,*)nkread
   DO ii = 1, nkread
      READ(13,*)kread(1:3,ii)
      write(*,*)kread(1:3,ii)
   END DO


   DO irr = 1, nkread !ciclo sui punti k della base UCRBF (kx,ky,kz)
      READ(13,'(i,a6)')i,comment!,kread(1,irr) !legge kx
      READ(13,'(4X,A3,3e20.10)') ,comment,tmp1,tmp2,tmp3
      write(*,*)comment,tmp1,tmp2,tmp3
      READ(13,'(I2,A5,3e20.10)') nn,comment,tmp1,tmp2,tmp3
      write(*,*)comment,tmp1,tmp2,tmp3
      READ(13,' (A)') comment!  ibnd jbnd mode    |g| [meV] ",/,30("-") )'  
      !write(*,*)comment
      READ(13,' (A)')
      do i = 1,nbnd
         do j = 1,nbnd
            do nu = 1,nqmodes
               READ(13,'(3i5,1e25.15)') ii, jj, ll, &
                    tmp1!, tmp2 !  
   !            write(*,'(3i5,2e25.15)') ii, jj, ll, tmp1, tmp2
               el_ph_mat(iq, ll, ii, jj) = tmp1
            enddo
         enddo
      enddo
   END DO

   deallocate(xk,kread)

end do
end do !end iq
close(13)
write(*,*)'End reading the DFPT file ',TRIM(input_file_DFPT)


do jyz=1,NKyz ! index of q_yz
    if(k_selec(jyz))then
      do jx=1,nqx ! index of qx
                  
         do iyz = 1,NKyz   ! index of k_yz
            if(k_selec(iyz))then
               
               allocate(el_ph_mtrx(iyz,jx,jyz,im)%M(nqmodes,NM,NM))
               el_ph_mtrx(iyz,jx,jyz,im)%M=0.0d0
                              
               jj = ind_kyz( k_vec(2:3,iyz) - k_vec(2:3,jyz) )  ! index of k_yz' = k_yz - q_yz 
               
               if(k_selec(jj))then
                  
               do nu= 1 , nqmodes
               do m = 1 , NM
               do n = 1 , NM      

                  el_ph_mtrx(iyz, jx, jyz, im)%M(nu, n, m) = &
                       el_ph_mat(jx+(jyz-1)*nqx, nu, ind_bnd(iyz,im)%i(n), ind_bnd(jj,im)%i(m))
                  
               end do
               end do
               end do

               end if
               
            end if
            
         end do
      end do
      
   end if
end do

deallocate(el_ph_mat)

!stop
end if
!!! end of the dfpt part



if (.not. dfpt) then
      
   do jyz=1,NKyz ! index of q_yz
      if(k_selec(jyz))then
      do iyz = 1,NKyz   ! index of k_yz
         if(k_selec(iyz))then
            
           ! allocate(el_ph_mtrx(iyz,1,jyz,im)%M(1,NM_mat(im),NM_mat(im)))
           ! el_ph_mtrx(iyz,1,jyz,im)%M=0.0d0

            
            nn=nqx-1
            do jx=1,nqx
               
           allocate(el_ph_mtrx(iyz,jx,jyz,im)%M(1,NM_mat(im),NM_mat(im)))
            el_ph_mtrx(iyz,jx,jyz,im)%M=0.0d0
            
               jj = ind_kyz( k_vec(2:3,iyz) - k_vec(2:3,jyz) )  ! index of k_yz' = k_yz - q_yz 
               
               if(k_selec(jj))then
                  allocate(C(nm,nm))
                  C=0.0_dp
                  ip=1 !!! no polarization is assumed
                                    
                  call omp_set_num_threads(Nrx)
                  !$omp parallel default(none) private(ix,n,m) &
                  !$omp shared(nrx,nm,ngt,im,nband_val,ip,jj,iyz,jx,nn,ind_kx,PSIBB,C)
                  
                  !$omp do 
                  do ix=1,Nrx
                     
                     do m=1,NM
                        do n=1,NM
                           !$omp atomic
                           C(n,m)=C(n,m)+zdotc(ngt,PSIBB(iyz,im)%H(1+(ix-1)*NGt+(ip-1)*Ngt*nrx:ngt+(ix-1)*NGt+(ip-1)*Ngt*nrx,n),1, &
                                PSIBB(jj,im)%H(1+(ix-1)*NGt+(ip-1)*Ngt*nrx:ngt+(ix-1)*NGt+(ip-1)*Ngt*nrx,m),1) * &
                                exp( -cmplx(0.0_dp,1.0_dp) * ( ind_kx(iyz,im)%n(n) - ind_kx(jj,im)%n(m) - dble(2*(jx-1)-dble(nn))/dble(2*nn)  ) * &
                                2.0_dp*pi*dble(ix-1)/dble(Nrx) )
                        end do
                     end do
                     
                  end do
                  !$omp end do
                  !$omp end parallel
                  
                  ll=1
                  do j=1,NM
                     do i=1,NM
                        el_ph_mtrx(iyz,jx,jyz,im)%M(ll,i,j)=C(i,j)
                        !!!el_ph_mtrx(iyz,1,jyz,im)%M(ll,i,j)=el_ph_mtrx(iyz,1,jyz,im)%M(ll,i,j)+C(i,j)!/dble(nqx)
                     end do
                  end do
                  deallocate(C)
                  
               end if

            end do
         end if
      end do
      end if
   end do
   
end if ! if not dfpt

t2=SECNDS(t1)
WRITE(*,*)'Time spent to compute the el-ph matrix (s)',t2
write(*,*)

endif  !end if phonons

!!!stop

end do ! end do im

deallocate(Uk)

end if



do i=1,nky
   if(abs(ky(i)-0.0d0)<1.0d-3 .or. abs(ky(i)-0.5d0)<1.0d-3) deg_ky(i)=1.0_dp
   write(*,*)'Ky',i,ky(i),deg_ky(i)
end do
do j=1,nkz
   if(abs(kz(j)-0.0d0)<1.0d-3 .or. abs(kz(j)-0.5d0)<1.0d-3) deg_kz(j)=1.0_dp
   write(*,*)'Kz',j,kz(j),deg_kz(j)
end do

write(*,*)'Nkyz=',Nkyz
allocate(deg_kyz(Nkyz))
deg_kyz=deg_ky*deg_kz

forall(j=1:nkz, i=1:nky) deg_kyz(i+(j-1)*nky)=deg_ky(i)*deg_kz(j)

do j=1,nkz
   do i=1,nky
      write(*,*)'deg_kyz',i+(j-1)*nky,deg_kyz(i+(j-1)*nky)
   end do
end do


!!!stop

end subroutine read_QE_output




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
 !!!! MODIFIED GRAM-SCHMIDT ALGORITHM

   implicit none
 
   integer, intent(IN) :: np, nm
   complex(dp), intent(INOUT) ::  Q(np,nm)
   integer :: j,k
   complex(dp), allocatable :: X(:,:),R(:,:)
 
   allocate(X(np,nm))  
   allocate(R(nm,nm)) 

   X=Q
   R=0.0_dp
   Q=0.0_dp
   R(1,1) = sqrt(dble(dot_product(X(:,1),X(:,1))))
   Q(:,1) = X(:,1)/R(1,1)
   
   do k=2,nm
      do j=1,k-1
         R(j,k) = dot_product(Q(:,j),X(:,k))
         X(:,k) =  X(:,k) - R(j,k)*Q(:,j) 
      end do
      R(k,k) = sqrt(dble(dot_product(X(:,k),X(:,k))))
      Q(:,k) = X(:,k)/R(k,k) 
   end do
   
   deallocate(R,X)
   
 end subroutine MGS
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
   
   call SUB_DEF_Z(1,NM,NM,B,E,U)

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

 integer function ind_Kyz(v)
   implicit none

   INTEGER  :: iz,iy,l
   real(dp) :: v(2),vv(2)

   ind_kyz=0

   vv=v

   if(abs(v(1))>0.5_dp) vv(1)=abs(abs(v(1))-1.0_dp)
   if(abs(v(2))>0.5_dp) vv(2)=abs(abs(v(2))-1.0_dp)
   
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
      write(*,*)'pb w ind_kyz',ind_kyz
      write(*,*)v(1),v(2)
      write(*,*)vv(1),vv(2)
      stop
   end if
   
 end function ind_Kyz
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SUB_DEF_0(Mi,Mf,ny,A,subband)

  implicit none
  integer     :: ny,mi,mf
  real(dp)    :: subband(1:(Mf-Mi+1))
  complex(dp) :: A(1:NY,1:NY),Uii(1:NY,1:Mf-Mi+1)
  integer     :: INFO
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
  
END SUBROUTINE SUB_DEF_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Pseudopot_so_gen
