MODULE Pseudopot_so_gen

  USE static
  USE indata

  IMPLICIT NONE

  SAVE


CONTAINS
  
subroutine read_QE_output
  implicit none
  
integer(k15), allocatable :: miller_2D(:,:)
integer :: i,ikx,iyz,ii,j,l,k,m,n,nkx,nrx0,ncell,n2,n3,m1,jgt,ix,iy,iz,ip,im
integer :: nm,ngmax

real(dp) :: a_1(3),a_2(3),a_3(3)
real(dp) :: b_1(3),b_2(3),b_3(3)
real(dp) :: vec(3),t0,a0
real(dp) :: Ecutoff,refec,refev
complex(dp) :: tmp
real(4) :: t1,t2
real(dp), allocatable :: E(:), KGt(:,:), Gx(:), bb_ev(:), bb_ec(:), hkl(:,:)
complex(dp), allocatable :: A(:,:),B(:,:),C(:,:),U(:,:)
complex(dp), allocatable :: HLL(:,:),TLL(:,:),HLLL(:,:),TLLL(:,:) 
complex(dp), allocatable :: dens_z(:,:,:), dens_yz(:,:,:)


if( .not. wr_ham ) write(*,*)'READING HAMILTONIANS'

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



if(g_spin==1) npol = 2
if(g_spin==2) npol = 1
write(*,*) 'NPOL =',NPOL


write(*,*) 'NKyz points=',Nkyz


Ecutoff=ryd*Ecut!/3.0_dp

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
do n2=-nry/2,nry/2!-maxval(miller(2,:)),maxval(miller(2,:))
do n3=-nrz/2,nrz/2!maxval(miller(3,:)),maxval(miller(3,:))
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
do n2=-nry/2,nry/2!-maxval(miller(2,:)),maxval(miller(2,:))
do n3=-nrz/2,nrz/2!-maxval(miller(3,:)),maxval(miller(3,:))
vec(2:3)=n2*b_2(2:3)+n3*b_3(2:3) !!!! THIS IMPLIES THAT b_1 is orthonal to the transverse plane

if(t0*dot_product(vec(2:3),vec(2:3))*(2.0_dp*pi/a0)**2 <= 1.0_dp*Ecutoff)then ! eV
!if(t0*vec(2)**2*(2.0_dp*pi/a0)**2 <= Ecutoff .and. t0*vec(3)**2*(2.0_dp*pi/a0)**2 <= Ecutoff )then 
   j=j+1
   miller_2D(2,j)=n2
   miller_2D(3,j)=n3
!   write(21,*)j,miller_2D(2,j),miller_2D(3,j)
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
write(*,*)'NRX',nrx,nrx0,'NRY',nry,'NRZ',nrz
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
  ! write(24,*)m1,gx(m1),int(gx(m1)/b_1(1))
end do


Dx=a_1(1)/dble(nrx)
Dy=a_2(2)/dble(nry)
Dz=a_3(3)/dble(nrz)


write(*,*)'a0/Dx',a0/Dx,Dx
write(*,*)'a0/Dy',a0/Dy,Dy
write(*,*)'a0/Dz',a0/Dz,Dz

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
!             write(*,*)i,nez(i),dble(i)*dble(NRZ)/dble(Ndeltaz)
          end if
       end do
    end do
  !  write(*,*)'nez',nez(1:Ndeltaz)
    do i=1,Ndeltaz
       do j=1,Nez(i)
   !       write(*,*)i,j,j+(i-1)*nez(i)
       end do
    end do
    allocate(Ney(Ndeltay))
    Ney=0
    do i=1,Ndeltay
       do iy=1,NRY
          if(dble(iy) <= dble(i)*dble(NRY)/dble(Ndeltay) .and. dble(iy) > dble(i-1)*dble(NRY)/dble(Ndeltay)) ney(i)=ney(i)+1
       end do
    end do
 !   write(*,*)'ney',ney(1:Ndeltay)
   ! do i=1,Ndeltay
   !    do j=1,Ney(i)
    !      write(*,*)i,j,j+(i-1)*ney(i)
    !   end do
    !end do
    allocate(Nex(Ndeltax))
    Nex=0
    do i=1,Ndeltax
       do ix=1,NRX
          if(dble(ix) <= dble(i)*dble(NRX)/dble(Ndeltax) .and. dble(ix) > dble(i-1)*dble(NRX)/dble(Ndeltax)) nex(i)=nex(i)+1
       end do
    end do
!    write(*,*)'nex',nex(1:Ndeltax)
!    do i=1,Ndeltax
!       do j=1,Nex(i)
!     !     write(*,*)i,j,j+(i-1)*nex(i)
!       end do
!    end do


    allocate(kc_min(Nkyz,num_mat))
    allocate(kv_max(Nkyz,num_mat))
    write(*,*)'NPOL=',npol

!!stop

do iyz=1,Nkyz  ! Loops over the transverse k vectors

!   nn=maxval(N_PW(1+(iyz-1)*nkx:iyz*nkx))
!write(*,*)'size of Hk',dble(nn*npol)**2*nkx*16.0d-9   
!   if(wr_ham == .true.)then
!      allocate(Hk(nn*npol,nn*npol,nkx))
!      Hk=0.0_dp
!   end if
   
write(*,*)'ikyz',iyz

   ky(iyz-((iyz-1)/nky))=k_vec(2,iyz)
   write(*,*)'iky',iyz-((iyz-1)/nky),ky(iyz-((iyz-1)/nky))    !!!! iyz = iy + (iz-1)*nky
   kz(1+((iyz-1)/nky))=k_vec(3,iyz)
   write(*,*)'ikz',1+((iyz-1)/nky),kz(1+((iyz-1)/nky)) !!!! iyz = iy + (iz-1)*nky

end do


do iyz=1,Nkyz
do im=1,num_mat
   nm=NM_mat(im)
   write(*,*)iyz,'im',im,nm
   allocate(HL(iyz,im)%H(NM_mat(im),NM_mat(im)))
   open(unit=13,file=TRIM(inputdir)//'H00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
   write(*,*)'reading ',TRIM(inputdir)//'H00_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   do i=1,nm
      do j=1,nm
         read(13,*)tmp
         HL(iyz,im)%H(i,j)=tmp
      end do
      HL(iyz,im)%H(i,i)=HL(iyz,im)%H(i,i)+off_set(im)
   end do
   close(13)

   allocate(TL(iyz,im)%H(NM_mat(im),NM_mat(im)))
   nm=NM_mat(im)
   write(*,*)iyz,'im',im,nm
   open(unit=13,file=TRIM(inputdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
   write(*,*)'reading ',TRIM(inputdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   do i=1,nm
      do j=1,nm
         read(13,*)tmp
         TL(iyz,im)%H(i,j)=tmp
      end do
   end do
   close(13)
  
end do
if(num_mat>=2)then

   l=num_mat
   do im=1,num_het
      l=l+1
      allocate(TL(iyz,l)%H(NM_mat(mat_1(im)),NM_mat(mat_0(im))))
      write(*,*)iyz,'ihl',l,NM_mat(mat_1(im)),NM_mat(mat_0(im))
      
      open(unit=13,file=TRIM(inputdir)//'H01_nkyz_'//TRIM(STRINGA(iyz))//'_nhet_'//TRIM(STRINGA(im))//'.dat',status='unknown')
      if(htype(im) .eq. 'n')then
         do i=1,NM_mat(mat_1(im))
            do j=1,NM_mat(mat_0(im))
               read(13,*)tmp
               TL(iyz,l)%H(i,j)=tmp
            end do
         end do
      else if(htype(im) .eq. 'c')then
         do j=1,NM_mat(mat_0(im))
            do i=1,NM_mat(mat_1(im))
               read(13,*)tmp
               TL(iyz,l)%H(i,j)=conjg(tmp)
            end do
         end do
      else
         write(*,*)'wrong htype'
         exit
      end if
      close(13)
      
   end do

end if


n=40
allocate(bb_ev(n+1),bb_ec(n+1))


do im=1,num_mat
   M=min(nband_val(im)+10,nm)
   allocate(hkl(n+1,M))
   NM=nm_mat(im)

   allocate(E(M))
   allocate(A(NM,NM),HLLL(NM,NM),TLLL(NM,NM))
   HLLL=HL(iyz,im)%H
   TLLL=TL(iyz,im)%H
   do ikx=1,n+1
      A=HLLL+TLLL*exp(cmplx(0.0_dp,1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)+&
           transpose(dconjg(TLLL))*exp(cmplx(0.0_dp,-1.0_dp)*dble(ikx-1-n/2)/dble(n)*2.0_dp*pi)
      call SUB_DEF_Z0(1,M,NM,A,E)
      do k=1,M
         hkl(ikx,k)=e(k)
      end do
if( nband_val(im)>0)     bb_ev(ikx)=E(nband_val(im))
      bb_ec(ikx)=E(nband_val(im)+1)
   end do
   deallocate(A,HLLL,TLLL)
   deallocate(E)
   do k=1,min(10,M)
      if(nband_val(im)>0)then
         ii=nband_val((im))-1*min(10,M)/2+k
      else
         ii=nband_val((im))+k
      end if
      write(*,*)k,nband_val((im)),ii
      open(unit=300+k,file='Edisp_'//TRIM(STRINGA(ii))//'_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
      do i=1,nm
         do ikx=1,n+1
            write(300+k,*)dble(ikx-1-n/2)/dble(n),hkl(ikx,ii)
         end do
         close(300+k)
      end do
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
write(*,*)iyz,'kc_min',kc_min(iyz,im),refec
write(*,*)iyz,'kv_max',kv_max(iyz,im),refev
deallocate(hkl)

end do

deallocate(bb_ec,bb_ev)

write(*,*)'shape TL', shape(TL)
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

end do ! endo do iyz


if(.not. onlyT)then


do im=1,num_mat
!if(schottky_type(im)==0)then      
do iyz=1,Nkyz

   allocate(ULCBB(iyz,im)%H(Nrx*NGt*npol,NM_mat(im)))
   NM=NM_mat(im)
   write(*,*)Nrx*NGt*npol,NM_mat(im),Nrx*NGt*npol*NM_mat(im)

   open(unit=13,file=TRIM(inputdir)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat',status='unknown')
   write(*,*)'reading ',TRIM(inputdir)//'Psi_Bloch_nkyz_'//TRIM(STRINGA(iyz))//'_nmat_'//TRIM(STRINGA(im))//'.dat'
   do j=1,NM
      do ip=1,npol
         do jgt=1,Ngt
            do ix=1,nrx   
               read(13,*)tmp
               ULCBB(iyz,im)%H(jgt+(ix-1)*ngt+(ip-1)*Ngt*nrx,j)=tmp
            end do
         end do
      end do
   end do
   close(13)

   allocate(form_factor(iyz,im)%F(NM_mat(im),NM_mat(im)))
 
end do
!end if
end do

allocate(U(NGt*npol,(nry)*(nrz)))

do im=1,num_mat
!if(schottky_type(im)==0)then
   
   NM=NM_mat(im)
   write(*,*)im,NM_mat(im),NM_mat(im)


   do iyz=1,Nkyz
      
      allocate(U_psi(iyz,im)%K(1:NM_mat(im)*NM_mat(im),1:(Ndeltay+1),1:(Ndeltaz+1),1:Nrx))
      
      write(*,*)'dimension of U_psi',iyz,im, (dble(NM_mat(im))**2)*dble(Ndeltay+1)*dble(Ndeltaz+1)*dble(nrx)*16.0d-9,'Gb'
      
!!!! INTERPOLATION ON THE COARSE GRID
      
      t1=SECNDS(0.0)
      
      do iy=1,NRY
         do iz=1,NRZ
            j=iy+(iz-1)*(NRY)
            U(1:NGt,j)=exp( cmplx(0.0_dp,-1.0_dp)*KGt_kyz(2,1:NGt,iyz)*2.0_dp*pi/a0*dble(iy)*Dy+&
                            cmplx(0.0_dp,-1.0_dp)*KGt_kyz(3,1:NGt,iyz)*2.0_dp*pi/a0*dble(iz)*Dz )/sqrt(dble((nry)*(nrz)))
         end do
      end do
      if(npol>1)      U(1+NGt:npol*NGt,1:NRY*NRZ)=U(1:NGt,1:NRY*NRZ)

call omp_set_num_threads(Nomp)!call omp_set_num_threads(3)

!$omp parallel default(none) private(ix,iy,iz,ip,i,j,jgt,dens_z,dens_yz,A,B,C) &
!$omp shared(iyz,im,ney,nez,nm,nrx,nry,nrz,ndeltay,ndeltaz,dx,dy,dz,npol,ngt,U,ULCBB,U_psi,FORM_FACTOR)


allocate(dens_z(NM*NM,(nry),Ndeltaz+1))
allocate(dens_yz(NM*NM,Ndeltay+1,Ndeltaz+1))

allocate(A(NM,(nry)*(nrz)))
allocate(B(NM,ngt*npol))
allocate(C(NM*NM,(nry)*(nrz)))

form_factor(iyz,im)%F=0.0d0

 !$omp do 
   do ix=1,Nrx
   write(*,*)'ix',ix

   do ip=1,npol
      do jgt=1,Ngt
         do i=1,NM
            B(i,jgt+(ip-1)*ngt)=(conjg(ULCBB(iyz,im)%H(jgt+(ix-1)*Ngt+(ip-1)*nrx*ngt,i)))
         end do
      end do
   end do
   call ZGEMM('n','n',NM,(nry)*(nrz),NGt*npol,alpha,B,NM,U,NGt*npol,beta,A,NM)

!   do iy=1,NRY
!      do iz=1,NRZ
!         do i=1,NM
!            A(i,iy+(iz-1)*(nry))=A(i,iy+(iz-1)*(nry))*exp(-cmplx(0.0_dp,1.0_dp)*ELCH/HBAR*Hz(ix)*1.0d-4*dx*dy*dble(iy-iz))   !!! bz is in Tesla
!         end do
!      end do
!   end do

   
!!!
!!$   do iy=1,nry
!!$      do iz=1,nrz
!!$         tmp=0.0d0
!!$         do i=1,nm!nband_val(im)
!!$            tmp=tmp+(A(i,iy+(iz-1)*(nry))*dconjg( A(i,iy+(iz-1)*(nry)) ))
!!$         end do
!!$         write(1000+ix,*)(iy-1)*dy*1e7,(iz-1)*dz*1e7,dble(tmp)
!!$               end do
!!$      write(1000+ix,*)
!!$   end do
!!!!
   
   do i=1,NM
      do j=1,NM    
         form_factor(iyz,im)%F(i,j)=form_factor(iyz,im)%F(i,j)+&
              sum(dconjg(A(i,:))*A(j,:)*dconjg(A(j,:))*A(i,:))/(nrx*dx*dy*dz) 
      end do
   end do

   do i=1,NM
      do j=1,NM
         C(i+(j-1)*NM,:)=dconjg(A(i,:))*A(j,:)
      end do
   end do

   dens_z=0.0_dp
   dens_yz=0.0_dp
   
   iz=1
   do iy=1,NRY
      dens_z(1:NM*NM,iy,iz)=C(1:NM*NM,iy+(iz-1)*(nry))
   end do
   dens_z(1:NM*NM,1:nry,Ndeltaz+1)=dens_z(1:NM*NM,1:nry,iz)
   iz=0
   do i=1,Ndeltaz-1
      do j=1,Nez(i)
         iz=iz+1
         do iy=1,NRY
            dens_z(1:NM*NM,iy,i+1)=dens_z(1:NM*NM,iy,i+1)+C(1:NM*NM,iy+(iz-1)*(nry))/dble(2*Nez(i))
         end do
      end do
   end do
   iz=Nez(1)
   do i=2,Ndeltaz
      do j=1,Nez(i)
         iz=iz+1
         do iy=1,NRY
            dens_z(1:NM*NM,iy,i)=dens_z(1:NM*NM,iy,i)+C(1:NM*NM,iy+(iz-1)*(nry))/dble(2*Nez(i))
         end do
      end do
   end do
   
   iy=1
   dens_yz(1:NM*NM,iy,1:Ndeltaz+1)=dens_z(1:NM*NM,iy,1:Ndeltaz+1)
   dens_yz(1:NM*NM,Ndeltay+1,1:Ndeltaz+1)=dens_yz(1:NM*NM,iy,1:Ndeltaz+1)
   iy=0
   do i=1,Ndeltay-1
      do j=1,Ney(i)
         iy=iy+1
         dens_yz(1:NM*NM,i+1,:)=dens_yz(1:NM*NM,i+1,:)+dens_z(1:NM*NM,iy,:)/dble(2*Ney(i))
      end do
   end do
   iy=Ney(1)
   do i=2,Ndeltay
      do j=1,Ney(i)
         iy=iy+1
         dens_yz(1:NM*NM,i,:)=dens_yz(1:NM*NM,i,:)+dens_z(1:NM*NM,iy,:)/dble(2*Ney(i))
      end do
   end do

   U_psi(iyz,im)%K(1:NM*NM,1:Ndeltay+1,1:Ndeltaz+1,ix)=dens_yz(1:NM*NM,1:Ndeltay+1,1:Ndeltaz+1)
   
end do
 !$omp end do

deallocate(A,B,C)
deallocate(dens_z,dens_yz)
!$omp end parallel

!!!! END

       t2=SECNDS(t1)
        WRITE(*,*)iyz,im,'TIME SPENT TO COMPUTE the interpolation POINT (s)',t2

!stop

if(iyz==1)then
   allocate(E(NM))
   allocate(B(NM,NM))
   allocate(A(NM,NM),HLL(NM,NM),TLL(NM,NM))
   HLL=HL(iyz,im)%H
   TLL=TL(iyz,im)%H

   A=HLL+TLL*exp(cmplx(0.0_dp,1.0_dp)*kc_min(iyz,im)*2.0_dp*pi)+&
        transpose(dconjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*kc_min(iyz,im)*2.0_dp*pi)
   call SUB_DEF_Z(1,NM,NM,A,E,B)
   ref_ec(im)=E(nband_val(im)+1)
   write(*,*)'ref_ec at Kyz=0 ',im,ref_ec(im)
if (nband_val(im)>0)then
   A=HLL+TLL*exp(cmplx(0.0_dp,1.0_dp)*kv_max(iyz,im)*2.0_dp*pi)+&
        transpose(dconjg(TLL))*exp(cmplx(0.0_dp,-1.0_dp)*kv_max(iyz,im)*2.0_dp*pi)
   call SUB_DEF_Z(1,NM,NM,A,E,B)
   ref_ev(im)=E(nband_val(im))
end if
   write(*,*)'ref_ev at Kyz=0 ',im,ref_ev(im)

   E_GAP=ref_ec(im)-ref_ev(im)
   write(*,*)'E_GAP at Kyz=0 = ',e_gap
   
   deallocate(A,HLL,TLL)
   deallocate(B)
   deallocate(E)
end if

write(*,*)'fine ikyz',iyz
end do
!deallocate(Upsi)
!end if
end do

deallocate(U)

end if

do i=1,nky
   if(abs(ky(i)-0.0d0)<1.0d-3 .or. abs(ky(i)-0.5d0)<1.0d-3) deg_ky(i)=1.0_dp
   write(*,*)'Ky',i,ky(i),deg_ky(i)
end do
do j=1,nkz
   if(abs(kz(j)-0.0d0)<1.0d-3 .or. abs(kz(j)-0.5d0)<1.0d-3) deg_kz(j)=1.0_dp
   write(*,*)'Kz',j,kz(j),deg_kz(j)
end do

write(*,*)'Kkyz',Nkyz
allocate(deg_kyz(Nkyz))
deg_kyz=deg_ky*deg_kz

forall(j=1:nkz, i=1:nky) deg_kyz(i+(j-1)*nky)=deg_ky(i)*deg_kz(j)

do j=1,nkz
   do i=1,nky
      write(*,*)'deg_kyz',i+(j-1)*nky,deg_kyz(i+(j-1)*nky)
   end do
end do


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
b(3)=1.0_dp !b=[0 0 1 0 0 0 0 0 0 0 0]';
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
