MODULE negf

USE static
USE indata

INTEGER, ALLOCATABLE :: NM(:), flag(:,:)
INTEGER              :: iyz, jyz

CONTAINS

SUBROUTINE negf_mixed(POT3D,charge_n,charge_p,ISDcurrent,IDScurrent,IDScurrentb,ss,gg,ext_iter)

IMPLICIT NONE

REAL(DP), INTENT(IN)     :: POT3D(1:NTOT_X,1:NTOT_Y,1:NTOT_Z)
INTEGER,  INTENT(IN)     :: ss,gg,ext_iter

REAL(DP), INTENT(OUT)    :: ISDcurrent,IDScurrent,IDScurrentb
REAL(DP), INTENT(OUT)    :: charge_n(1:NTOT_X,1:NTOT_Y,1:NTOT_Z)
REAL(DP), INTENT(OUT)    :: charge_p(1:NTOT_X,1:NTOT_Y,1:NTOT_Z)

INTEGER                  :: i, j, l, n, xx, ix, iy, iz, ip, ii, pp, nn
INTEGER                  :: Nop, ref_index, nmax, ee, nee, nsol, SCBA_iter

REAL(DP)                 :: Gcon, en, epsilon, emin_local
REAL(DP), ALLOCATABLE    :: con(:),cone(:),conb(:)

COMPLEX(DP), ALLOCATABLE :: g_lesser_diag_local(:,:,:), g_lesser(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: g_greater_diag_local(:,:,:), g_greater(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: g_r_diag_local(:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_lesser_ph(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_greater_ph(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_lesser_ph_prev(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_greater_ph_prev(:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_r_ph_prev(:,:,:)

REAL(DP), ALLOCATABLE    :: emin_yz(:), emax_yz(:), SCBA_error(:), SCBA_x(:), SCBA_f(:)
REAL(DP)                 :: SCBA_scalar
REAL(DP)                 :: n_bose_g

REAL(DP), ALLOCATABLE    :: degeneracy(:), trans(:,:,:,:), cur(:,:,:,:)
REAL(DP), ALLOCATABLE    :: ldos(:,:,:,:), zdos(:,:,:,:), ndos(:,:,:,:), pdos(:,:,:,:)
REAL(DP)                 :: tr,tre,ttr,sumt,kappax

COMPLEX(DP), allocatable :: Hi(:,:,:,:),pot(:,:)
COMPLEX(dp), allocatable :: dosn(:,:,:,:),dosp(:,:,:,:)
COMPLEX(dp), allocatable :: U(:,:), GBB(:,:), Gcc(:,:), A(:,:), B(:,:)
REAL(DP), allocatable    :: subband(:,:,:), neutr(:)
REAL(dp), allocatable    :: pt(:), CR(:,:), Ry(:), Rz(:), dens_yz(:,:),dens_z(:,:)
REAL(dp), allocatable    :: omega(:,:), derror(:,:), derror_old(:,:)

real(4) :: t1,t2

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  ISDcurrent=0.0_dp
  IDScurrent=0.0_dp
  charge_n=0.0_dp
  charge_p=0.0_dp
  Gcon=0.0_dp

  allocate(RY(NY+1))
  do iy=1,NY+1
     RY(iy)=dble(iy)*Dy
  end do
  allocate(RZ(NZ+1))
  do iz=1,NZ+1
     RZ(iz)=dble(iz)*Dz
  end do

  write(*,*)'NTOT_X=',nTOT_X
  write(*,*)'NTOT_Y=',nTOT_Y
  write(*,*)'NTOT_Z=',nTOT_Z


  write(*,*)'ac',ac
  write(*,*)'dx',dx
  write(*,*)'dy',dy
  write(*,*)'dz',dz
  write(*,*)'nx',nx
  write(*,*)'ny',ny
  write(*,*)'nz',nz

open(unit=10,file='poty_'//TRIM(STRINGA(ext_iter))//'.dat',status='unknown')
do j=Ntot_y/2,NTOT_Y/2
   do i=1,NTOT_X
      do l=1,NTOT_Z
          write(10,*)i,l,POT3D(i,j,l)
       enddo
    write(10,*)
    enddo
 enddo
close(10)

open(unit=10,file='potx.dat',status='unknown')
 do l=1,NTOT_Z
    do j=1,Ntot_y
       do i=NTOT_X/2,NTOT_X/2
          write(10,*)j,l,POT3D(i,j,l)
       enddo
    enddo
    write(10,*)
 enddo
close(10)

open(unit=10,file='potz.dat',status='unknown')
 do l=NTOT_Z/2,NTOT_Z/2
    do j=1,Ntot_y
       do i=1,NTOT_X
          write(10,*)i,j,POT3D(i,j,l)
       enddo
    write(10,*)
    enddo
 enddo
close(10)


allocate(con(NKYZ),cone(NKYZ),conb(NKYZ))
  con=0.0_dp
  cone=0.0_dp
  conb=0.0_dp

  write(*,*)'NKGt=',Ngt*npol
  write(*,*)'NKGt*Nrx=',Ngt*npol*Nrx

  allocate(NM(ncx_d))
  
  do l=1,ncx_d
     NM(l)=nm_mat(imat(l))
  end do


NMAX=maxval(NM)
write(*,*)'NMAX',NMAX

allocate(Hi(NMAX,NMAX,ncx_d,nkyz))
Hi=0.0_dp

allocate(degeneracy(nkyz))
do iyz=1,NKYZ !!! loop over transverse k-vectors
   if(k_selec(iyz))then
   t1=SECNDS(0.0)
   
   degeneracy(iyz)=deg_kyz(iyz)*g_spin 
   write(*,*)iyz,'degeneracy',deg_kyz(iyz),degeneracy(iyz)

   if(allocated(kgt))deallocate(kgt)
   allocate(KGt(4,1*Ngt))
   KGt(1:4,1:1*Ngt)=KGt_kyz(1:4,1:1*Ngt,iyz)


if(.not. onlyT)then
allocate(U(Ngt,(Ny+1)*(Nz+1)))
do iy=1,NY+1
   do iz=1,NZ+1
      j=iy+(iz-1)*(NY+1)
      U(1:Ngt,j)=exp(dcmplx(0.0_dp,-1.0_dp)*KGt_kyz(2,1:Ngt,iyz)*2.0_dp*pi/ac*Ry(iy)+&
           dcmplx(0.0_dp,-1.0_dp)*KGt_kyz(3,1:Ngt,iyz)*2.0_dp*pi/ac*Rz(iz))/sqrt(dble((Ny+1)*(Nz+1)))
   end do
end do

call omp_set_num_threads(Nomp) !!! this sets the environment variable

!$omp parallel default(none) private(xx,ix,iy,iz,i,j,n,ip,pt,Gcc,GBB,A,B,CR,t1,t2) &
!$omp shared(NCX_D,Nrx,Ndeltax,Ndeltay,Ndeltaz,Ny,Nz,NMAX,NM,NGT,npol,iyz,imat,ac,&
!$omp ULCBB,pot3D,KGt_kyz,Ry,Rz,Hi,U,deltay,deltaz,dy,dz,to2_lft,to2_bot )

!allocate(pot(Ngt,Ngt))
allocate(Gcc(Ngt,Ngt))
allocate(pt((Ny+1)*(Nz+1)))
allocate(A(Ngt,(Ny+1)*(Nz+1)))
allocate(CR(2,2))

  !$omp do
do xx=1,Ncx_D

do ip=1,npol
!   Gcc=0.0_dp
!   do ix=1,Ndeltax
!      pt=0.0_dp
!      do iy=1,Ny+1  
!      do iz=1,Nz+1  
!         CR(1,1)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3),  to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3))
!         CR(2,1)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3+1),to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3))
!         CR(1,2)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3),  to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3+1))
!         CR(2,2)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3+1),to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3+1))
!         call surf_interp(iy,iz,CR,pt(iy+(iz-1)*(NY+1)))
!      end do
!      end do
!      forall (i = 1:Ngt, j = 1:(Ny+1)*(Nz+1) ) A(i,j) = U(i,j) * pt(j)
!      
!      call zgemm('n','c',Ngt,ngt,(Ny+1)*(Nz+1),alpha,A,ngt,U(1:ngt,1:(Ny+1)*(Nz+1)),ngt,beta,pot,ngt)
!      
!      Gcc(1:ngt,1:Ngt)=Gcc(1:Ngt,1:Ngt)+pot(1:Ngt,1:Ngt)/dble(Ndeltax)
!      
!   end do
   
   !!!! TRANSFORMING INTO THE Bloch states BASIS

   allocate(GBB(Ngt,NM(xx)))
   allocate(B(NM(xx),NM(xx)))
   do n=1,Nrx
      
      if( mod(n,ceiling(dble(Nrx)/dble(Ndeltax))) == 1)then
         ix=ceiling(dble(n)/dble(Nrx)*dble(Ndeltax))
         if(ix > Ndeltax)then
            write(*,*)'error in the interpolation',ix,Ndeltax
            stop
         end if
         
      pt=0.0_dp
      do iy=1,Ny+1
      do iz=1,Nz+1
         CR(1,1)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3),  to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3))
         CR(2,1)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3+1),to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3))
         CR(1,2)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3),  to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3+1))
         CR(2,2)=pot3D(ix+(xx-1)*Ndeltax, to2_lft+ceiling(dble(iy-1)*dy/deltay+1.0d-3+1),to2_bot+ceiling(dble(iz-1)*dz/deltaz+1.0d-3+1))
         call surf_interp(iy,iz,CR,pt(iy+(iz-1)*(NY+1)))
      end do
      end do

      forall (i = 1:Ngt, j = 1:(Ny+1)*(Nz+1) ) A(i,j) = U(i,j) * pt(j)

      call zgemm('n','c',Ngt,ngt,(Ny+1)*(Nz+1),alpha,A,ngt,U(1:ngt,1:(Ny+1)*(Nz+1)),ngt,beta,Gcc,ngt)

      end if

      
      call ZGEMM('n','n',ngt,NM(xx),ngt,alpha,Gcc(1:ngt,1:ngt),ngt,&
           ULCBB(iyz,imat(xx))%H(1+(n-1)*Ngt+(ip-1)*ngt*nrx:n*Ngt+(ip-1)*ngt*nrx,1:NM(xx)),Ngt,beta,GBB(1:ngt,1:NM(xx)),Ngt)      
      call ZGEMM('c','n',NM(xx),NM(xx),Ngt,alpha,ULCBB(iyz,imat(xx))%H(1+(n-1)*Ngt+(ip-1)*ngt*nrx:n*Ngt+(ip-1)*ngt*nrx,1:NM(xx)),&
           Ngt,GBB(1:ngt,1:NM(xx)),Ngt,beta,B(1:NM(xx),1:NM(xx)),NM(xx))

      Hi(1:NM(xx),1:NM(xx),xx,iyz)=Hi(1:NM(xx),1:NM(xx),xx,iyz)+B(1:NM(xx),1:NM(xx))
      
   end do
   deallocate(GBB,B)
end do
Hi(1:NM(xx),1:NM(xx),xx,iyz)=(Hi(1:NM(xx),1:NM(xx),xx,iyz)+transpose(conjg(Hi(1:NM(xx),1:NM(xx),xx,iyz))))/2.0_dp

end do
!$omp end do

deallocate(pt)
!deallocate(pot)
deallocate(Gcc)
deallocate(A)
deallocate(CR)

!$omp end parallel


if(allocated(U))deallocate(U)
end if

t2=SECNDS(t1)
if(.not.onlyT) write(*,*)'potential transformed in ',t2,'s'

end if
end do !fine loop kyz

nsol=nsolv+nsolc
allocate(subband(nsol,ncx_d,nkyz))
allocate(emin_yz(NKYZ),emax_yz(NKYZ))

do iyz=1,NKYZ
if(k_selec(iyz))then
call omp_set_num_threads(NCX_D) 
!$omp parallel default(none) private(xx,ref_index,kappax,A) &
!$omp shared(NCX_D,iyz,NM,imat,ihet,nband_val,Hi,HL,TL,nsolv,nsolc,subband,kc_min,kv_max,chtype)
!$omp do
do xx=1,ncx_d
   ref_index=nband_val(imat(xx))!!!+off_k_nvb(iyz)
   
   allocate(A(NM(xx),NM(xx)))
   
   kappax=0.0d0

   Hi(1:NM(xx),1:NM(xx),xx,iyz)=Hi(1:NM(xx),1:NM(xx),xx,iyz)+HL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx))!HL(1:NM(xx),1:NM(xx),iyz,im(xx))
   if(nband_val(imat(xx))>0)then
      kappax=kv_max(iyz,imat(xx))
      A(1:NM(xx),1:NM(xx))=Hi(1:NM(xx),1:NM(xx),xx,iyz)+&
           (TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx)))*exp(cmplx(0.0_dp,1.0_dp)*kappax*2.0_dp*pi)+&
           transpose(dconjg(TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx))))*exp(cmplx(0.0_dp,-1.0_dp)*kappax*2.0_dp*pi)
      call sub_def_0(ref_index-nsolv+1,ref_index,NM(xx),A,subband(1:nsolv,xx,iyz))
   end if
   kappax=kc_min(iyz,imat(xx))
   A(1:NM(xx),1:NM(xx))=Hi(1:NM(xx),1:NM(xx),xx,iyz)+&
        (TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx)))*exp(cmplx(0.0_dp,1.0_dp)*kappax*2.0_dp*pi)+&
        transpose(dconjg(TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx))))*exp(cmplx(0.0_dp,-1.0_dp)*kappax*2.0_dp*pi)
   call sub_def_0(ref_index+1,ref_index+nsolc,NM(xx),A,subband(nsolv+1:nsolv+nsolc,xx,iyz))
   deallocate(A)
end do
  !$omp end do
  !$omp end parallel


if(nsolv==0 .and. chtype /= 'n')then
   write(*,*)'pb w nsolv or chtype'
end if

!! devo definire il punto di neutralita'
if(.not.allocated(neutr))allocate(neutr(ncx_d))
do xx=1,ncx_d
   if(nsolv>0)then
      neutr(xx)=(subband(nsolv,xx,iyz)+subband(nsolv+1,xx,iyz))/2.0_dp
   else
      neutr(xx)=subband(nsolv+1,xx,iyz)-100.0d-3
   end if
   write(199,*)xx,neutr(xx)
enddo

  epsilon=5*Nop_g*Eop!(BOLTZ*TEMP)

  SELECT CASE (chtype)
  CASE ('n')
        emin_yz(iyz)=min(min(mud,mus),MINVAL(subband(nsolv+1,:,iyz)))-epsilon
        emax_yz(iyz)=max(max(mus,mud),MAXVAL(subband(nsolv+1,:,iyz)))+NKT*(BOLTZ*TEMP)
        
   CASE ('p')
        emin_yz(iyz)=min(min(mud,mus),MINVAL(subband(nsolv,:,iyz)))-NKT*(BOLTZ*TEMP)
        emax_yz(iyz)=max(max(mus,mud),MAXVAL(subband(nsolv,:,iyz)))+epsilon

  CASE ('t')
     if(mud.ge.mus)then

        emin_yz(iyz)=mus-NKT*(BOLTZ*TEMP)!MIN(mud,MINVAL(subband(nsolv+1,:)))-NKT*(BOLTZ*TEMP)
        emax_yz(iyz)=mud+NKT*(BOLTZ*TEMP)!MAX(mus,MAXVAL(subband(nsolv,:)))+NKT*(BOLTZ*TEMP)
        if(MINVAL(subband(nsolv+1,:,iyz)).gt.mus+NKT*(BOLTZ*TEMP))emax=mud
        if(MAXVAL(subband(nsolv,:,iyz)).lt.mud-NKT*(BOLTZ*TEMP))  emin=mus

        if(MINVAL(subband(nsolv+1,:,iyz)).gt.mud+1*(BOLTZ*TEMP) .and. &
             MAXVAL(subband(nsolv,:,iyz)).lt.mus-1*(BOLTZ*TEMP)) emin_yz(iyz)=emax_yz(iyz)
   else
        emin_yz(iyz)=MIN(mud,MINVAL(subband(nsolv+1,:,iyz)))-NKT*(BOLTZ*TEMP)
        emax_yz(iyz)=MAX(mus,MAXVAL(subband(nsolv,:,iyz)))+NKT*(BOLTZ*TEMP)
        if(MINVAL(subband(nsolv+1,:,iyz)).gt.mus+NKT*(BOLTZ*TEMP))emax=mud
        if(MAXVAL(subband(nsolv,:,iyz)).lt.mud-NKT*(BOLTZ*TEMP))  emin=mus
        if(MINVAL(subband(nsolv+1,:,iyz)).gt.mus+NKT*(BOLTZ*TEMP) .and. &
             MAXVAL(subband(nsolv,:,iyz)).lt.mud-NKT*(BOLTZ*TEMP)) emin_yz(iyz)=emax_yz(iyz)
     end if
     !if(emin.gt.emax) emin_yz(iyz)=emax_yz(iyz)
  END SELECT

  do j=1,nsol
     do xx=1,ncx_d
        write(400+j,*)xx,subband(j,xx,iyz)
     enddo
  enddo

  end if
end do  ! fine loop kyz

  emax=maxval(emax_yz(:))
  emin=minval(emin_yz(:))

  if(onlyT)then
     emin=min(emin,min(mus,mud))-NKT*(BOLTZ*TEMP)
     emax=max(emax,max(mus,mud))+NKT*(BOLTZ*TEMP)
  end if
  
  deallocate(emin_yz,emax_yz)

  Nop=FLOOR((emax-emin)/Eop+0.5_dp)

  write(*,*)'GLOBAL ENERGY RANGE'
  write(*,*)'mus=',mus,'mud=',mud
  write(*,*)'Emax=',emax
  write(*,*)'Emin=',emin
  write(*,*)'Nop=',Nop

  allocate(flag(Nop,NKYZ))
  
!!!!!!!!!!!!do iyz=1,NKYZ

  IF(Nop.eq.0)THEN
     write(*,*) 'NOP null, stopping...'
     stop
  END IF
  
  allocate(dosn(1:NMAX,1:NMAX,1:NcX_d,1:NKYZ))
  allocate(dosp(1:NMAX,1:NMAX,1:NcX_d,1:NKYZ))

  dosn=0.0_dp
  dosp=0.0_dp

  ALLOCATE(cur(1:Nsub,1:Nop,1:Ncx_d-1,1:NKYZ))
  ALLOCATE(ldos(1:Nsub,1:Nop,1:Ncx_d,1:NKYZ))
  ALLOCATE(zdos(1:Nsub,1:Nop,1:Ncx_d,1:NKYZ))
  ALLOCATE(ndos(1:Nsub,1:Nop,1:Ncx_d,1:NKYZ))
  ALLOCATE(pdos(1:Nsub,1:Nop,1:Ncx_d,1:NKYZ))
  ALLOCATE(trans(1:Nsub,1:Nop,4,1:NKYZ))


  ! ==========  Start of parallel resolution ===========================
  !write(*,*) 'starting the parallel loop'
  
  call omp_set_num_threads(min(Nomp,Nop)) !!! this sets the environment variable OMP_NUM_THREADS 

  ALLOCATE(sigma_greater_ph(1:Nop,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(sigma_lesser_ph(1:Nop,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(sigma_greater_ph_prev(1:Nop,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(sigma_lesser_ph_prev(1:Nop,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(sigma_r_ph_prev(1:Nop,1:NMAX,1:Ncx_d))
  ALLOCATE(SCBA_error(1:NMAX),SCBA_x(1:NMAX),SCBA_f(1:NMAX))

  do ee=1,Nsub

     t1=SECNDS(0.0)

     emin_local=en_global(ee) 
     flag=0

if(phonons)then 

     sigma_greater_ph_prev=0.0_dp
     sigma_lesser_ph_prev=0.0_dp
     sigma_r_ph_prev=0.0_dp
     

     SCBA_scalar=1.0_dp
     SCBA_iter=0  

     ALLOCATE(derror(Nop,Ncx_D),derror_old(Nop,Ncx_D),omega(Nop,Ncx_D))
     derror_old=0.0_dp
     omega=0.0_dp
     ALLOCATE(g_lesser(Nop,NMAX,Ncx_d,NKYZ),g_greater(Nop,NMAX,Ncx_d,NKYZ))
     DO WHILE((SCBA_scalar.gt.SCBA_tolerance).and.(SCBA_iter.lt.SCBA_max_iter))
        SCBA_iter=SCBA_iter+1
        
        g_lesser = 0.0d0
        g_greater= 0.0d0

        do iyz=1,NKYZ !loop over kyz    
if(k_selec(iyz))then
           !$omp parallel default(none)  private(jyz,nee,xx,pp,nn,EN,tr,tre,g_lesser_diag_local,g_greater_diag_local,g_r_diag_local) &
           !$omp shared(ee,iyz,Nop,nm,imat,ncx_d,ncy,ncz,nmax,nkyz,emin,emin_local,Eop,mus,mud,hi,&
           !$omp sigma_lesser_ph_prev,sigma_greater_ph_prev,sigma_r_ph_prev,cur,ldos,ndos,pdos,zdos,chtype,neutr,&
           !$omp degeneracy,g_spin,w,temp,trans,g_lesser,g_greater,form_factor,flag,k_selec)

        ALLOCATE(g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
        ALLOCATE(g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
        ALLOCATE(g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
    
       !$omp do
        do nee=1,Nop
           EN=emin+emin_local+dble(nee-1)*Eop
           
           sigma_r_ph_prev(nee,1:NMAX,1:NCX_D)=DCMPLX(0.0_dp*DBLE(sigma_r_ph_prev(nee,1:NMAX,1:NCX_D)),&
                0.5_dp*((aimag(sigma_greater_ph_prev(nee,1:NMAX,1:NCX_D,iyz))) - (aimag(sigma_lesser_ph_prev(nee,1:NMAX,1:NCX_D,iyz)))))
     
           call RGF(NMAX,flag(nee,iyz),EN,mus,mud,Hi(1:NMAX,1:NMAX,1:Ncx_d,iyz),&
                sigma_lesser_ph_prev(nee,1:NMAX,1:Ncx_d,iyz),sigma_greater_ph_prev(nee,1:NMAX,1:Ncx_d,iyz),sigma_r_ph_prev(nee,1:NMAX,1:Ncx_d),&
                g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d),&
                tr,tre,cur(ee,nee,1:Ncx_d-1,iyz))
           
           do xx=1,ncx_d
              do nn=1,nm(xx)

                 do jyz = 1,NKYZ
              if(k_selec(jyz))then
                 g_lesser (nee,1:nm(xx),xx,jyz)=g_lesser (nee,1:nm(xx),xx,jyz)+&
                      degeneracy(iyz)/g_spin/dble(NCY*NCZ)*g_lesser_diag_local(nn,nn,xx) * form_factor(jyz,iyz,imat(xx))%F(1:nm(xx),nn)
                 
                 g_greater(nee,1:nm(xx),xx,jyz)=g_greater(nee,1:nm(xx),xx,jyz)+&
                      degeneracy(iyz)/g_spin/dble(NCY*NCZ)*g_greater_diag_local(nn,nn,xx)* form_factor(jyz,iyz,imat(xx))%F(1:nm(xx),nn)
              end if
              end do
              
              end do
           end do
           
        enddo  ! end do nee
        !$omp end do


  DEALLOCATE(g_lesser_diag_local)
  DEALLOCATE(g_greater_diag_local)
  DEALLOCATE(g_r_diag_local)

  !$omp end parallel

  end if
end do !End of the loop over kyz
     


     !  self-energies calculation
  sigma_lesser_ph=0.0d0
  sigma_greater_ph=0.0d0
  !sigma_r_ph=0.0d0

  do jyz = 1,NKYZ
     if(k_selec(jyz))then
  !$omp parallel default(none)  private(nee,ii,pp,l) shared(jyz,Nop,Nop_g,sigma_lesser_ph,sigma_greater_ph, &
  !$omp g_lesser,g_greater,Dop_g,n_bose_g,Dac,ncx_d,nm,eop,temp)

  !$omp do
  DO nee=1,Nop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! g-type optical phonon !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     n_bose_g=1.0_dp/(EXP((Nop_g*Eop)/(BOLTZ*TEMP))-1.0_dp)
     
     IF(nee.le.l*Nop_g)THEN
        IF(nee.le.Nop-l*Nop_g)THEN
           ! E+Eop_g
           
           DO ii=1,Ncx_d
              DO pp=1,NM(ii)
                 
                 sigma_lesser_ph(nee,pp,ii,jyz)=sigma_lesser_ph(nee,pp,ii,jyz)+&
                      Dop_g*(n_bose_g+1.0d0)*g_lesser(nee+l*Nop_g,pp,ii,jyz)
                 
                 sigma_greater_ph(nee,pp,ii,jyz)=sigma_greater_ph(nee,pp,ii,jyz)+&
                      Dop_g*(n_bose_g)*(g_greater(nee+l*Nop_g,pp,ii,jyz))
                 
              END DO
           END DO
           
        END IF
           
     ELSE
        
        IF(nee.le.Nop-l*Nop_g)THEN
           ! E+Eop_g and E-Eop_g
           
           DO ii=1,Ncx_d
              DO pp=1,NM(ii)
      
                 sigma_lesser_ph(nee,pp,ii,jyz)= sigma_lesser_ph(nee,pp,ii,jyz)+&
                      Dop_g*(n_bose_g+1.0d0)*g_lesser(nee+l*Nop_g,pp,ii,jyz)+&   
                      Dop_g*(n_bose_g)*g_lesser(nee-l*Nop_g,pp,ii,jyz) 
      
                 sigma_greater_ph(nee,pp,ii,jyz)=sigma_greater_ph(nee,pp,ii,jyz)+&
                      Dop_g*(n_bose_g)*(g_greater(nee+l*Nop_g,pp,ii,jyz))+&
                      Dop_g*(n_bose_g+1.0d0)*(g_greater(nee-l*Nop_g,pp,ii,jyz))
                 
              END DO
           END DO
           
        ELSE 
           ! E-Eop_g
           DO ii=1,Ncx_d
              DO pp=1,NM(ii)       
      
                 sigma_lesser_ph(nee,pp,ii,jyz)=sigma_lesser_ph(nee,pp,ii,jyz)+&
                      Dop_g*(n_bose_g)*g_lesser(nee-l*Nop_g,pp,ii,jyz)  
      
                 sigma_greater_ph(nee,pp,ii,jyz)=sigma_greater_ph(nee,pp,ii,jyz)+&
                      Dop_g*(n_bose_g+1.0d0)*(g_greater(nee-l*Nop_g,pp,ii,jyz))
                 
              END DO
           END DO
           
        END IF
     END IF
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! Acoustic !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
     DO ii=1,Ncx_d
        DO pp=1,NM(ii)       
      
           sigma_lesser_ph (nee,pp,ii,jyz)=sigma_lesser_ph (nee,pp,ii,jyz) + &
                Dac*g_lesser (nee,pp,ii,jyz) 
      
           sigma_greater_ph(nee,pp,ii,jyz)=sigma_greater_ph(nee,pp,ii,jyz) + &
                Dac*g_greater(nee,pp,ii,jyz) 

        END DO
     END DO

        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!IMAGINARY PART!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     DO ii=1,Ncx_d
        DO pp=1,NM(ii)
        
           sigma_lesser_ph(nee,pp,ii,jyz)=DCMPLX(0.0_dp*dble(sigma_lesser_ph(nee,pp,ii,jyz)),abs(aimag(sigma_lesser_ph(nee,pp,ii,jyz))))
        
           sigma_greater_ph(nee,pp,ii,jyz)=DCMPLX(0.0_dp*dble(sigma_greater_ph(nee,pp,ii,jyz)),-abs(aimag(sigma_greater_ph(nee,pp,ii,jyz))))
        
        end do
     end do
     
  END DO
  !$omp end do  ! end of self-energies calculation
  !$omp end parallel
end if
end do
  
  if(scba_scalar > 1.0d-3)then
     omega=0.0_dp
  else
     omega=1.0_dp-scba_alpha
  end if
  
  SCBA_f=0.0_dp
  SCBA_x=0.0_dp
  do nee=1,Nop
     DO ii=1,Ncx_D
        do jyz=1,NKYZ
        if(k_selec(jyz))then
        SCBA_f(1:NM(ii))=SCBA_f(1:NM(ii))+((dimag(sigma_lesser_ph(nee,1:NM(ii),ii,jyz)-sigma_lesser_ph_prev(nee,1:NM(ii),ii,jyz)))**2 +&
          (dimag(sigma_greater_ph(nee,1:NM(ii),ii,jyz)-sigma_greater_ph_prev(nee,1:NM(ii),ii,jyz)))**2   )
          
        SCBA_x(1:NM(ii))=SCBA_x(1:NM(ii))+((dimag(sigma_lesser_ph(nee,1:NM(ii),ii,jyz)))**2 + (dimag(sigma_greater_ph(nee,1:NM(ii),ii,jyz)))**2)
        
        derror(nee,ii)=     abs(sum(dimag(sigma_lesser_ph(nee,1:NM(ii),ii,jyz)-sigma_lesser_ph_prev(nee,1:NM(ii),ii,jyz))))/dble(NM(ii)) 
        end if
        end do
     END DO
  end do
  derror_old=derror

  SCBA_error(:)=SCBA_f(:)/(SCBA_x(:)+1.0d-16)
  SCBA_scalar=0.0_dp !errore
  DO i=1,NMAX
!!!$omp atomic
     SCBA_scalar=SCBA_scalar+sqrt(SCBA_error(i))/dble(NMAX)
  END DO
  

!do nee=1,Nop
!EN=emin+emin_local+dble(nee-1)*Eop
!do xx=1,ncx_d
!   write(656+ee*1000,*)(xx-1)*ac1*1.0d7,EN, sum(aimag( sigma_lesser_ph(nee,1:NMAX,xx)))
!   write(657+ee*1000,*)(xx-1)*ac1*1.0d7,EN,-sum(aimag(sigma_greater_ph(nee,1:NMAX,xx)))
!end do
!write(656+ee*1000,*)
!write(657+ee*1000,*)
!end do
!close(656+ee*1000)
!close(657+ee*1000)
     


  DO nee=1, Nop
     DO ii=1,Ncx_D
        DO jyz=1,NKYZ
           sigma_lesser_ph_prev(nee,1:NM(ii),ii,jyz)=sigma_lesser_ph_prev(nee,1:NM(ii),ii,jyz) + &
                (1.0d0-omega(nee,ii))*(sigma_lesser_ph(nee,1:NM(ii),ii,jyz)-sigma_lesser_ph_prev(nee,1:NM(ii),ii,jyz))     
           sigma_greater_ph_prev(nee,1:NM(ii),ii,jyz)=sigma_greater_ph_prev(nee,1:NM(ii),ii,jyz) + &
                (1.0d0-omega(nee,ii))*(sigma_greater_ph(nee,1:NM(ii),ii,jyz)-sigma_greater_ph_prev(nee,1:NM(ii),ii,jyz))
        END DO
     END DO
  END DO
  
  write(*,*)ee,'scba_iter',SCBA_iter,SCBA_scalar


  END DO!Fine autoconsistenza di Born
!stop

  DEALLOCATE(g_lesser,g_greater)
  DEALLOCATE(omega,derror,derror_old)

end if !fine endif

  do iyz=1,NKYZ !loop over kyz
if(k_selec(iyz))then
     !$omp parallel default(none)  private(nee,xx,EN,tr,tre,g_lesser_diag_local,g_greater_diag_local,g_r_diag_local) &
     !$omp shared(ee,iyz,Nop,ncx_d,nkyz,nmax,emin,emin_local,Eop,mus,mud,hi,sigma_lesser_ph_prev,sigma_greater_ph_prev,sigma_r_ph_prev, &
     !$omp cur,ldos,ndos,pdos,zdos,chtype,neutr,degeneracy,w,temp,con,cone,conb,trans,dosn,dosp,flag,phonons)

  ALLOCATE(g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
  ALLOCATE(g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
  ALLOCATE(g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d))

        
  !$omp do
  do nee=1,Nop
     EN=emin+emin_local+dble(nee-1)*Eop
         
     if(phonons)then
        call RGF(NMAX,flag(nee,iyz),EN,mus,mud,Hi(1:NMAX,1:NMAX,1:Ncx_d,iyz),&
             0.0_dp*sigma_r_ph_prev(nee,1:NMAX,1:Ncx_d),0.0_dp*sigma_r_ph_prev(nee,1:NMAX,1:Ncx_d),0.0_dp*sigma_r_ph_prev(nee,1:NMAX,1:Ncx_d),&
             g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d),&
             tr,tre,cur(ee,nee,1:Ncx_d-1,iyz))
        conb(iyz) =conb(iyz) +degeneracy(iyz)*w(ee)*(tr)
     end if
     
     sigma_r_ph_prev(nee,1:NMAX,1:NCX_D)=DCMPLX(0.0_dp*DBLE(sigma_r_ph_prev(nee,1:NMAX,1:NCX_D)),&
          0.5_dp*((aimag(sigma_greater_ph_prev(nee,1:NMAX,1:NCX_D,iyz))) - (aimag(sigma_lesser_ph_prev(nee,1:NMAX,1:NCX_D,iyz))))) 

     call RGF(NMAX,flag(nee,iyz),EN,mus,mud,Hi(1:NMAX,1:NMAX,1:Ncx_d,iyz),&
          sigma_lesser_ph_prev(nee,1:NMAX,1:Ncx_d,iyz),sigma_greater_ph_prev(nee,1:NMAX,1:Ncx_d,iyz),sigma_r_ph_prev(nee,1:NMAX,1:Ncx_d),&
          g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d),&
          tr,tre,cur(ee,nee,1:Ncx_d-1,iyz))
              
     do xx=1,ncx_d
        ldos(ee,nee,xx,iyz)=-traccia(dimag(g_r_diag_local(1:NMAX,1:NMAX,xx)))
        ndos(ee,nee,xx,iyz)=traccia(dimag(g_lesser_diag_local(1:NMAX,1:NMAX,xx)))
        pdos(ee,nee,xx,iyz)=-traccia(dimag(g_greater_diag_local(1:NMAX,1:NMAX,xx)))
        zdos(ee,nee,xx,iyz)= traccia(dimag( g_greater_diag_local(1:NMAX,1:NMAX,xx) - g_lesser_diag_local(1:NMAX,1:NMAX,xx) &
             -g_r_diag_local(1:NMAX,1:NMAX,xx)+transpose(dconjg(g_r_diag_local(1:NMAX,1:NMAX,xx)))))
     end do
     trans(ee,nee,1,iyz)=tr
     trans(ee,nee,2,iyz)=tre
     trans(ee,nee,3,iyz)=tr /(1.0_dp/(1.0_dp+exp((EN-mus)/(BOLTZ*TEMP)))-1.0_dp/(1.0_dp+exp((EN-mud)/(BOLTZ*TEMP))))
     trans(ee,nee,4,iyz)=tre/(1.0_dp/(1.0_dp+exp((EN-mud)/(BOLTZ*TEMP)))-1.0_dp/(1.0_dp+exp((EN-mus)/(BOLTZ*TEMP))))
     SELECT CASE (chtype)
     CASE ('t') 
        do xx=1,ncx_d    
           if(EN .ge. neutr(xx))then 
              dosn(1:NMAX,1:NMAX,xx,iyz)=dosn(1:NMAX,1:NMAX,xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_lesser_diag_local(1:NMAX,1:NMAX,xx))/(2.0_dp*pi)
           else if(EN .lt. neutr(xx))then ! 
              dosp(1:NMAX,1:NMAX,xx,iyz)=dosp(1:NMAX,1:NMAX,xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_greater_diag_local(1:NMAX,1:NMAX,xx))/(2.0_dp*pi)
           end if
        enddo
     CASE ('n')           
        do xx=1,ncx_d
           if(EN .ge. neutr(xx))then !if(EN.ge.subband(nsolv+1,xx)-epsilon)then ! if(EN.gt.neutr(xx))then
              dosn(1:NMAX,1:NMAX,xx,iyz)=dosn(1:NMAX,1:NMAX,xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_lesser_diag_local(1:NMAX,1:NMAX,xx))/(2.0_dp*pi)
           end if
        enddo
     CASE ('p')
        do xx=1,ncx_d
           if(EN .lt. neutr(xx))then !if(EN.le.subband(nsolv,xx)+epsilon)then ! if(EN.lt.neutr(xx))then
              dosp(1:NMAX,1:NMAX,xx,iyz)=dosp(1:NMAX,1:NMAX,xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_greater_diag_local(1:NMAX,1:NMAX,xx))/(2.0_dp*pi)
           end if
        end do
     END SELECT
     
     con(iyz) = con(iyz)  + degeneracy(iyz)*w(ee)*(tr)
  
     cone(iyz)= cone(iyz) + degeneracy(iyz)*w(ee)*(tre)
     
  enddo  ! Nop
  !$omp end do
     


  DEALLOCATE(g_lesser_diag_local)
  DEALLOCATE(g_greater_diag_local)
  DEALLOCATE(g_r_diag_local)

  !$omp end parallel
     
  write(*,*)ee,con(iyz),cone(iyz)
  end if
  end do !end loop kyz
     
enddo     ! Nsub
  

  deallocate(sigma_greater_ph_prev)
  deallocate(sigma_lesser_ph_prev)
  deallocate(sigma_r_ph_prev)
  deallocate(sigma_greater_ph)
  deallocate(sigma_lesser_ph)
  deallocate(SCBA_error,SCBA_x)
  deallocate(flag)

write(*,*)'Writing files'

  do iyz=1,NKYZ  
if(k_selec(iyz))then
  ! ==========  End of parallel resolution ===========================
if(iyz < 10)then

if(phonons)then
   OPEN(UNIT=10,FILE='scat_pdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=20,FILE='scat_ndens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=30,FILE='scat_LDOS_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE='scat_Jdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE='scat_Jx_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=50,FILE='scat_Jspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=51,FILE='scat_Tspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=61,FILE='scat_subv_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=62,FILE='scat_subc_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
else
   OPEN(UNIT=10,FILE='bal_pdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=20,FILE='bal_ndens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=30,FILE='bal_LDOS_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=60,FILE='bal_ZDOS_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE='bal_Jdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE='bal_Jx_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=50,FILE='bal_Jspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=51,FILE='bal_Tspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=61,FILE='bal_subv_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=62,FILE='bal_subc_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
end if



  do xx=1,ncx_d
     if(nsolv > 0) write(61,*)(xx-1)*ac1*1.0d7,subband(nsolv,xx,iyz)
     write(62,*)(xx-1)*ac1*1.0d7,subband(nsolv+1,xx,iyz)
   end do
   
do nee=1,Nop
   do ee=1,Nsub
     emin_local=en_global(ee) 
        EN=emin+emin_local+dble(nee-1)*Eop
        do xx=1,ncx_d
           write(10,*)(xx-1)*ac1*1.0d7,EN,pdos(ee,nee,xx,iyz)
           write(20,*)(xx-1)*ac1*1.0d7,EN,ndos(ee,nee,xx,iyz)
           write(30,*)(xx-1)*ac1*1.0d7,EN,ldos(ee,nee,xx,iyz)
           write(60,*)(xx-1)*ac1*1.0d7,EN,zdos(ee,nee,xx,iyz)
        end do
        do xx=1,ncx_d-1
           write(40,*)(xx-1)*ac1*1.0d7,EN,cur(ee,nee,xx,iyz)
        end do
        write(10,*)
        write(20,*)
        write(30,*)
        write(60,*)
        write(40,*)
        write(50,*)EN,trans(ee,nee,1,iyz),trans(ee,nee,2,iyz)
        write(51,*)EN,trans(ee,nee,3,iyz),trans(ee,nee,4,iyz)
     end do
  end do
  do xx=1,ncx_d-1
     ttr=0.0_dp
     do ee=1,Nsub
        do nee=1,Nop
           emin_local=en_global(ee) 
           EN=emin+emin_local+dble(nee-1)*Eop
           ttr=ttr+w(ee)*cur(ee,nee,xx,iyz)*Eop
        end do
     end do
     write(41,*)(xx-1)*ac1*1.0d7,ttr
  end do
end if


    CLOSE (UNIT=10)
    CLOSE (UNIT=20)
    CLOSE (UNIT=60)
    CLOSE (UNIT=30)
    CLOSE (UNIT=40)
    CLOSE (UNIT=41)
    CLOSE (UNIT=50)
    CLOSE (UNIT=51)
    CLOSE (UNIT=61)
    CLOSE (UNIT=62)

!!!stop

if(.not.onlyT)then    
   
   Write(*,*) 'Transforming the carrier density',', ikyz',iyz

   t1=SECNDS(0.0)
   
   call omp_set_num_threads(Nomp)
  !$omp parallel default(none) private(xx,n,i,j,ix,iy,iz,dens_yz,dens_z) &
  !$omp shared(chtype,iyz,imat,NCX_D,Nrx,Ndeltax,Ndeltay,Ndeltaz,NM,Ny,Nz,ac,Ry,Rz,dosn,dosp, &
  !$omp U_PSI,DX,DY,DZ,charge_n,charge_p,Ney,Nez,NTOT_Y,NTOT_Z,to2_lft,to2_bot)

   allocate(dens_yz(Ndeltay+1,Ndeltaz+1))

  !$omp do
   do xx=1,NCX_D
     !write(*,*)'xx',xx,imat(xx)
      do ix=1,Nrx
         if((chtype .eq. 'n') .or. (chtype .eq. 't'))then
            dens_yz=0.0_dp
            do i=1,NM(xx)
               do j=1,NM(xx)
                  dens_yz(1:Ndeltay+1,1:Ndeltaz+1)=dens_yz(1:Ndeltay+1,1:Ndeltaz+1)+&
                       dimag(dosn(i,j,xx,iyz)*U_psi(iyz,imat(xx))%K(i+(j-1)*NM(xx),1:Ndeltay+1,1:Ndeltaz+1,ix))
               end do
            end do
            do i=1,Ndeltay+1
               do j=1,Ndeltaz+1
                  if(dens_yz(i,j)<0.0_dp)dens_yz(i,j)=0.0_dp
               end do
            end do
            charge_n(1+(xx-1)*Ndeltax,1+to2_lft:Ndeltay+1+to2_lft ,1+to2_bot:Ndeltaz+1+to2_bot )=&
                 charge_n(1+(xx-1)*Ndeltax,1+to2_lft:Ndeltay+1+to2_lft ,1+to2_bot:Ndeltaz+1+to2_bot )+&
                 dens_yz(1:Ndeltay+1,1:Ndeltaz+1)
         end if
         if((chtype .eq. 'p') .or. (chtype .eq. 't'))then
            dens_yz=0.0_dp
            do i=1,NM(xx)
               do j=1,NM(xx)
                  dens_yz(1:Ndeltay+1,1:Ndeltaz+1)=dens_yz(1:Ndeltay+1,1:Ndeltaz+1)-&
                       dimag(dosp(i,j,xx,iyz)*U_psi(iyz,imat(xx))%K(i+(j-1)*NM(xx),1:Ndeltay+1,1:Ndeltaz+1,ix))
               end do
            end do
            do i=1,Ndeltay+1
               do j=1,Ndeltaz+1
                  if(dens_yz(i,j)<0.0_dp)dens_yz(i,j)=0.0_dp
               end do
            end do
            charge_p(1+(xx-1)*Ndeltax,1+to2_lft:Ndeltay+1+to2_lft ,1+to2_bot:Ndeltaz+1+to2_bot )=&
                 charge_p(1+(xx-1)*Ndeltax,1+to2_lft:Ndeltay+1+to2_lft ,1+to2_bot:Ndeltaz+1+to2_bot )+&
                 dens_yz(1:Ndeltay+1,1:Ndeltaz+1)
         end if
         
      end do
      
      do j=2,Ndeltax
         charge_n(j+(xx-1)*Ndeltax,1:NTOT_y,1:NTOT_Z)=charge_n(1+(xx-1)*Ndeltax,1:NTOT_y,1:NTOT_Z)
         charge_p(j+(xx-1)*Ndeltax,1:NTOT_y,1:NTOT_Z)=charge_p(1+(xx-1)*Ndeltax,1:NTOT_y,1:NTOT_Z)
      end do
      
   end do
!$omp end do

   deallocate(dens_yz)

  !$omp end parallel

end if

end if
end do ! end do iyz

if(phonons)then
   OPEN(UNIT=30,FILE='scat_LDOS_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE='scat_Jdens_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE='scat_Jx_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
else
   OPEN(UNIT=30,FILE='bal_LDOS_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE='bal_Jdens_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE='bal_Jx_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
end if


do nee=1,Nop
   do ee=1,Nsub
     emin_local=en_global(ee) !local minimum
     EN=emin+emin_local+dble(nee-1)*Eop
     do xx=1,ncx_d
     sumt=0.0_dp
     do iyz=1,nkyz
if(k_selec(iyz))   sumt=sumt+ldos(ee,nee,xx,iyz)
     end do
     write(30,*)(xx-1)*ac1*1.0d7,EN,sumt
     end do
     do xx=1,ncx_d-1
     sumt=0.0_dp
     do iyz=1,nkyz
if(k_selec(iyz))     sumt=sumt+cur(ee,nee,xx,iyz)
     end do        
     write(40,*)(xx-1)*ac1*1.0d7,EN,sumt
     end do
     write(30,*)
     write(40,*)
  end do
end do

do xx=1,ncx_d-1
   ttr=0.0_dp
   do iyz=1,NKYZ
if(k_selec(iyz))then
      do ee=1,Nsub
         do nee=1,Nop
            emin_local=en_global(ee) 
            EN=emin+emin_local+dble(nee-1)*Eop
            ttr=ttr+w(ee)*cur(ee,nee,xx,iyz)*Eop
         end do
       end do
end if
    end do
    write(41,*)(xx-1)*ac1*1.0d7,ttr
 end do
 close(30)
 close(40)
 close(41)

  deallocate(trans)
  deallocate(ldos)
  deallocate(zdos)
  deallocate(ndos)
  deallocate(pdos)
  deallocate(cur)
  deallocate(dosn)
  deallocate(dosp)

deallocate(subband)
deallocate(Hi)
     

  charge_n=charge_n/(DX*NCY*DY*NCZ*DZ)/(dble(Nrx)/dble(Ndeltax))
  charge_p=charge_p/(DX*NCY*DY*NCZ*DZ)/(dble(Nrx)/dble(Ndeltax))

charge_n(NTOT_X,1:NTOT_Y,1:NTOT_Z)=charge_n(NTOT_X-1,1:NTOT_Y,1:NTOT_Z)
charge_p(NTOT_X,1:NTOT_Y,1:NTOT_Z)=charge_p(NTOT_X-1,1:NTOT_Y,1:NTOT_Z)

open(unit=10,file='caricx_'//TRIM(STRINGA(ext_iter))//'.dat',status='unknown')
if(chtype == 'n')then
do xx=1,1!Ncx_d
   do ix=1,1!Ndeltax
      do iz=1,NTOT_Z
         do iy=1,NTOT_Y
            write(10,*)(iy-1)*deltay*1d7,(iz-1)*deltaz*1d7,charge_n(ix+(xx-1)*Ndeltax,iy,iz)
         end do
         write(10,*)
      end do
   end do
end do
else
do xx=1,1!Ncx_d
   do ix=1,1!Ndeltax
      do iz=1,NTOT_Z
         do iy=1,NTOT_Y
            write(10,*)(iy-1)*deltay*1d7,(iz-1)*deltaz*1d7,charge_p(ix+(xx-1)*Ndeltax,iy,iz)
         end do
         write(10,*)
      end do
   end do
end do
end if
close(10)

open(unit=10,file='caricy_'//TRIM(STRINGA(ext_iter))//'.dat',status='unknown')
if(chtype == 'n')then
do xx=1,Ncx_d
   do ix=1,Ndeltax
      do iz=1,NTOT_Z
         do iy=NTOT_Y/2+1,NTOT_Y/2+1
            write(10,*)ix+(xx-1)*Ndeltax,iz,charge_n(ix+(xx-1)*Ndeltax,iy,iz)
         end do
      end do
      write(10,*)
   end do
end do
else
do xx=1,Ncx_d
   do ix=1,Ndeltax
      do iz=1,NTOT_Z
         do iy=NTOT_Y/2,NTOT_Y/2
            write(10,*)ix+(xx-1)*Ndeltax,iz,charge_p(ix+(xx-1)*Ndeltax,iy,iz)
         end do
      end do
      write(10,*)
   end do
end do
end if
close(10)

deallocate(KGt)
deallocate(NM)
        t2=SECNDS(t1)
        WRITE(*,*)'TIME SPENT TO COMPUTE THE CHARGE (s)',t2
        WRITE(*,*)     

  IDScurrent=sum(con(1:NKYZ))/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
  ISDcurrent=sum(cone(1:NKYZ))/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
  if(phonons)then
      IDScurrentb=sum(conb(1:NKYZ))/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
   else
      IDScurrentb=IDScurrent
   end if
  
  write(*,*)'IDScurrent =',IDScurrent,ISDcurrent

  Gcon=sum(con(1:NKYZ))/abs(mud-mus)

  write(*,*)'Gcon=',Gcon,sum(cone(1:nkyz))/abs(mud-mus)
  write(*,*)

  deallocate(con,cone,conb)

!!!!stop
END SUBROUTINE negf_mixed


subroutine RGF(nmax,ff,E,mul,mur,Hii,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur) 

  implicit none
  
  INTEGER     :: i,j,l,nmax,ff
  REAL(dp)    :: E,mul,mur,tr,tre
  REAL(dp)    :: cur(Ncx_d-1)
  COMPLEX(dp) :: sigma_lesser_ph(nmax,Ncx_d),sigma_greater_ph(nmax,Ncx_d),sigma_r_ph(nmax,Ncx_d)
  COMPLEX(dp) :: sig(nmax,nmax),sigmal(nmax,nmax),sigmar(nmax,nmax)
  COMPLEX(dp) :: Gn(nmax,nmax),Gp(nmax,nmax),G00(nmax,nmax),GN0(nmax,nmax)
  COMPLEX(dp) :: Gl(nmax,nmax,Ncx_d),Gln(nmax,nmax,Ncx_d),Glp(nmax,nmax,Ncx_d)
  COMPLEX(dp) :: ldos(nmax,nmax,ncx_d),ndens(nmax,nmax,ncx_d),pdens(nmax,nmax,ncx_d)
  COMPLEX(dp) :: Hii(nmax,nmax,ncx_d),H00(nmax,nmax),H10(nmax,nmax)
  COMPLEX(dp) :: A(nmax,nmax),B(nmax,nmax),C(nmax,nmax),Id(nmax,nmax)
  COMPLEX(dp) :: z

if(ff == 0)then

  z=dcmplx(E,0.0d-6)

  Id=0.0_dp
  do i=1,nmax
     Id(i,i)=dcmplx(1.0_dp,0.0_dp)
  enddo

  Gln=0.0_dp
  Glp=0.0_dp
  Gl=0.0_dp
  ldos=0.0_dp
  ndens=0.0_dp
  pdens=0.0_dp
  cur=0.0_dp


! self energy on the left contact
  l=1
  H00(1:nm(l),1:nm(l))=Hii(1:nm(l),1:nm(l),l)
  do  i=1,nm(l)
     H00(i,i)=H00(i,i)+sigma_r_ph(i,l)
  enddo

  H10(1:nm(l),1:nm(l))=TL(iyz,ihet(l))%H(1:NM(l),1:NM(l))
  call sancho(nm(l),E,H00(1:nm(l),1:nm(l)),transpose(dconjg(H10(1:nm(l),1:nm(l)))),G00(1:nm(l),1:nm(l))) 
  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           ff=l
           write(*,*)'NaN warning! Pb w lft snch. Eliminating E=',E
           G00=0.0_dp
           exit
        end if
     end do
  end do
  if( traccia(dimag(- G00(:,:))) < 0.0d0 )then 
     ff=l
     write(*,*)'pb w snch l',E,traccia(dimag(- G00(:,:)))
  end if

  
  H10(1:nm(1),1:nm(1))=TL(iyz,ihet(l))%H(1:NM(1),1:NM(1))  !!!! H_{1,0}
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,H10(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,A(1:nm(l),1:nm(l)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),H10(1:nm(l),1:nm(l)),nm(l),beta,sigmal(1:nm(l),1:nm(l)),nm(l))  
  A(1:nm(l),1:nm(l))=z*Id(1:nm(l),1:nm(l))-H00(1:nm(l),1:nm(l))-sigmal(1:nm(l),1:nm(l))
                  
  call invert(A(1:nm(l),1:nm(l)),nm(l),nm(l))
  Gl(1:nm(l),1:nm(l),l)=A(1:nm(l),1:nm(l))
 
  sig(1:nm(l),1:nm(l))=-(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm((E-mul)/(BOLTZ*TEMP)) 
  do  i=1,nm(l)
     sig(i,i)=sig(i,i)+sigma_lesser_ph(i,l)
  end do
  
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  Gln(1:nm(l),1:nm(l),l)=C(1:nm(l),1:nm(l))

  sig(1:nm(l),1:nm(l))=(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm(-(E-mul)/(BOLTZ*TEMP))
  do  i=1,nm(l)
     sig(i,i)=sig(i,i)+sigma_greater_ph(i,l)
  end do

  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  Glp(1:nm(l),1:nm(l),l)=C(1:nm(l),1:nm(l))

  Do l=2,Ncx_d-1
     !write(*,*)'l=',l
     H00(1:nm(l),1:nm(l))=Hii(1:nm(l),1:nm(l),l)
     do i=1,nm(l)
        H00(i,i)=H00(i,i)+sigma_r_ph(i,l)
     end do
     H10(1:nm(l),1:nm(l-1))=TL(iyz,ihet(l))%H(1:NM(l),1:NM(l-1)) !!!! H_{l,l-1}
     call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Gl(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
     call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
     A(1:nm(l),1:nm(l))=z*Id(1:nm(l),1:nm(l))-H00(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))              
     call invert(A(1:nm(l),1:nm(l)),nm(l),nm(l))
     Gl(1:nm(l),1:nm(l),l)=A(1:nm(l),1:nm(l))

     sig(1:nm(l-1),1:nm(l-1))=Gln(1:nm(l-1),1:nm(l-1),l-1)
     call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),sig(1:nm(l-1),1:nm(l-1)),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
     call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
     do i=1,nm(l)
        C(i,i)=C(i,i)+sigma_lesser_ph(i,l)
     end do

     call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),C(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l)) 
     call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,Gn(1:nm(l),1:nm(l)),nm(l))

     Gln(1:nm(l),1:nm(l),l)=Gn(1:nm(l),1:nm(l))

     sig(1:nm(l-1),1:nm(l-1))=Glp(1:nm(l-1),1:nm(l-1),l-1)
     call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),sig(1:nm(l-1),1:nm(l-1)),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l))
     call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
     do i=1,nm(l)
        C(i,i)=C(i,i)+sigma_greater_ph(i,l)
     end do

     call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),C(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
     call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,Gp(1:nm(l),1:nm(l)),nm(l))
     Glp(1:nm(l),1:nm(l),l)=Gp(1:nm(l),1:nm(l))

  enddo

! self energy on the right contact
  l=Ncx_d

  H00(1:nm(l),1:nm(l))=Hii(1:nm(l),1:nm(l),l)
  do i=1,nm(l)
     H00(i,i)=H00(i,i)+sigma_r_ph(i,l)
  end do

  H10(1:nm(l),1:nm(l))=TL(iyz,ihet(l+1))%H(1:NM(l),1:NM(l))    !!! H(N+1,N)
  call sancho(nm(l),E,H00(1:nm(l),1:nm(l)),H10(1:nm(l),1:nm(l)),G00(1:nm(l),1:nm(l)))
  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           ff=l
           write(*,*)'NaN warning! Pb w rht snch. Eliminating E=',E
           G00=0.0_dp
           exit
        end if
     end do
  end do
  if( traccia(dimag(- G00(:,:))) < 0.0d0 )then 
     ff=l
     write(*,*)'pb w snch r',E,traccia(dimag(- G00(:,:)))
  end if

  call zgemm('c','n',nm(l),nm(l),nm(l),alpha,H10(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,A(1:nm(l),1:nm(l)),nm(l)) 
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),H10(1:nm(l),1:nm(l)),nm(l),beta,sigmar(1:nm(l),1:nm(l)),nm(l)) 

  H10(1:nm(l),1:nm(l-1))=TL(iyz,ihet(l))%H(1:NM(l),1:NM(l-1))  !!! H(N,N-1)
  call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Gl(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))

  G00(1:nm(l),1:nm(l))=z*Id(1:nm(l),1:nm(l))-H00(1:nm(l),1:nm(l))-sigmar(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))   
  call invert(G00(1:nm(l),1:nm(l)),nm(l),nm(l))

  call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Gln(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  sig(1:nm(l),1:nm(l))=-(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm((E-mur)/(BOLTZ*TEMP))+C(1:nm(l),1:nm(l))
  do i=1,nm(l)
     sig(i,i)=sig(i,i)+sigma_lesser_ph(i,l)
  end do

  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,G00(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,Gn(1:nm(l),1:nm(l)),nm(l)) 

  ndens(1:nm(l),1:nm(l),l)=Gn(1:nm(l),1:nm(l))

  call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Glp(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  sig(1:nm(l),1:nm(l))=(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm(-(E-mur)/(BOLTZ*TEMP))+C(1:nm(l),1:nm(l))
  do i=1,nm(l)
     sig(i,i)=sig(i,i)+sigma_greater_ph(i,l)
  end do

  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,G00(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,Gp(1:nm(l),1:nm(l)),nm(l))

  pdens(1:nm(l),1:nm(l),l)=Gp(1:nm(l),1:nm(l))
  ldos(1:nm(l),1:nm(l),l)=G00(1:nm(l),1:nm(l))
 
  A(1:nm(l),1:nm(l))=-(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm((E-mur)/(BOLTZ*TEMP))
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gp(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  A(1:nm(l),1:nm(l))=(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm(-(E-mur)/(BOLTZ*TEMP))
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gn(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  tr=-traccia(dble(B(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))))

!-------------------------
   
     if( abs(traccia(dimag( pdens(:,:,l) - ndens(:,:,l) - 2.0_dp*ldos(:,:,l)  ))) > 1.0d-1 )then
!        write(*,*)'pssbl pb w E =',E
        !ff=l  
     end if

  
  do l=Ncx_d-1,1,-1

     H10(1:nm(l+1),1:nm(l))=TL(iyz,ihet(l+1))%H(1:NM(l+1),1:NM(l))  !!! H(l+1,l)
     
     call zgemm('n','c',nm(l+1),nm(l),nm(l),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),Gl(1:nm(l),1:nm(l),l),nm(l),beta, B(1:nm(l+1),1:nm(l)), nm(l+1)) 
     if(chtype == 'p')then
        call zgemm('n','n',nm(l+1),nm(l),nm(l+1),alpha,Gp(1:nm(l+1),1:nm(l+1)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta, C(1:nm(l+1),1:nm(l)), nm(l+1))
     else
        call zgemm('n','n',nm(l+1),nm(l),nm(l+1),alpha,Gn(1:nm(l+1),1:nm(l+1)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta, C(1:nm(l+1),1:nm(l)), nm(l+1))
     end if
     if(chtype == 'p')then
        call zgemm('n','n',nm(l+1),nm(l),nm(l),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),Glp(1:nm(l),1:nm(l),l),nm(l),beta, B(1:nm(l+1),1:nm(l)), nm(l+1))
     else
        call zgemm('n','n',nm(l+1),nm(l),nm(l),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),Gln(1:nm(l),1:nm(l),l),nm(l),beta, B(1:nm(l+1),1:nm(l)), nm(l+1))
     end if
     call zgemm('n','n',nm(l+1),nm(l),nm(l+1),alpha,G00(1:nm(l+1),1:nm(l+1)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta, A(1:nm(l+1),1:nm(l)), nm(l+1))
     B(1:nm(l+1),1:nm(l))=C(1:nm(l+1),1:nm(l))+A(1:nm(l+1),1:nm(l))
     call zgemm('c','n',nm(l),nm(l),nm(l+1),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta,A(1:nm(l),1:nm(l)),nm(l))      !!! G<_i+1,i
     cur(l)=2.0_dp*traccia(dble(A(1:nm(l),1:nm(l))))


     !-------------------------

     call zgemm('n','c',nm(l),nm(l+1),nm(l),alpha,Gl(1:nm(l),1:nm(l),l),nm(l),H10(1:nm(l+1),1:nm(l)),nm(l+1),beta,B(1:nm(l),1:nm(l+1)),nm(l)) 
     call zgemm('n','n',nm(l),nm(l+1),nm(l+1),alpha,B(1:nm(l),1:nm(l+1)),nm(l),G00(1:nm(l+1),1:nm(l+1)),nm(l+1),beta,GN0(1:nm(l),1:nm(l+1)),nm(l))      !!! G_i,i+1

     call zgemm('n','n',nm(l),nm(l),nm(l+1),alpha,GN0(1:nm(l),1:nm(l+1)),nm(l),H10(1:nm(l+1),1:nm(l)),nm(l+1),beta,B(1:nm(l),1:nm(l)),nm(l))
     call zgemm('n','n',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),Gl(1:nm(l),1:nm(l),l),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
     
     G00(1:nm(l),1:nm(l)) = Gl(1:nm(l),1:nm(l),l) + C(1:nm(l),1:nm(l))                                       !!! G_i,i

!-------------------------
     
     call zgemm('n','c',nm(l),nm(l+1),nm(l),alpha,Gl(1:nm(l),1:nm(l),l),nm(l),H10(1:nm(l+1),1:nm(l)),nm(l+1),beta,A(1:nm(l),1:nm(l+1)),nm(l))  
     call zgemm('n','n',nm(l),nm(l+1),nm(l+1),alpha,A(1:nm(l),1:nm(l+1)),nm(l),Gn(1:nm(l+1),1:nm(l+1)),nm(l+1),beta,C(1:nm(l),1:nm(l+1)),nm(l))     

     call zgemm('n','n',nm(l),nm(l),nm(l+1),alpha,C(1:nm(l),1:nm(l+1)),nm(l),H10(1:nm(l+1),1:nm(l)),nm(l+1),beta,B(1:nm(l),1:nm(l)),nm(l))
     call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),Gl(1:nm(l),1:nm(l),l),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))

     Gn(1:nm(l),1:nm(l)) = Gln(1:nm(l),1:nm(l),l) + C(1:nm(l),1:nm(l)) 

     call zgemm('n','n',nm(l),nm(l),nm(l+1),alpha,GN0(1:nm(l),1:nm(l+1)),nm(l),H10(1:nm(l+1),1:nm(l)),nm(l+1),beta,sig(1:nm(l),1:nm(l)),nm(l)) 
     call zgemm('n','n',nm(l),nm(l),nm(l),alpha,sig(1:nm(l),1:nm(l)),nm(l),Gln(1:nm(l),1:nm(l),l),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))     
     
     Gn(1:nm(l),1:nm(l)) = Gn(1:nm(l),1:nm(l))+C(1:nm(l),1:nm(l))-transpose(dconjg(C(1:nm(l),1:nm(l)))) !!!! I am using  (G<)^+ = -G< 	 !!! G<_i,i
!-------------------------
     
     call zgemm('n','n',nm(l),nm(l+1),nm(l+1),alpha,A(1:nm(l),1:nm(l+1)),nm(l),Gp(1:nm(l+1),1:nm(l+1)),nm(l+1),beta,C(1:nm(l),1:nm(l+1)),nm(l))     
     call zgemm('n','n',nm(l),nm(l),nm(l+1),alpha,C(1:nm(l),1:nm(l+1)),nm(l),H10(1:nm(l+1),1:nm(l)),nm(l+1),beta,B(1:nm(l),1:nm(l)),nm(l))
     call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),Gl(1:nm(l),1:nm(l),l),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))

     Gp(1:nm(l),1:nm(l)) = Glp(1:nm(l),1:nm(l),l) + C(1:nm(l),1:nm(l))  
     call zgemm('n','n',nm(l),nm(l),nm(l),alpha,sig(1:nm(l),1:nm(l)),nm(l),Glp(1:nm(l),1:nm(l),l),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))     
     
     Gp(1:nm(l),1:nm(l)) = Gp(1:nm(l),1:nm(l))+C(1:nm(l),1:nm(l))-transpose(dconjg(C(1:nm(l),1:nm(l))))   !!! G<_i,i
!-------------------------    
     ndens(1:nm(l),1:nm(l),l)=Gn(1:nm(l),1:nm(l))
     pdens(1:nm(l),1:nm(l),l)=Gp(1:nm(l),1:nm(l))
     ldos(1:nm(l),1:nm(l),l)=G00(1:nm(l),1:nm(l))
  
     
     if( abs(traccia(dimag( pdens(:,:,l) - ndens(:,:,l) - 2.0_dp*ldos(:,:,l)  ))) > 1.0d-1 )then
!        write(*,*)'pssbl pb w E =',E
        !ff=l
     end if

  enddo
  l=1
  A(1:nm(l),1:nm(l))=-(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm((E-mul)/(BOLTZ*TEMP))
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gp(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  A(1:nm(l),1:nm(l))=(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm(-(E-mul)/(BOLTZ*TEMP))
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gn(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  tre=-traccia(dble(B(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))))

 if ( ff /= 0 ) then
     write(*,*)'pssbl pb w E =',E,ff
     tr=0.0_dp
     tre=0.0_dp
     ndens=0.0_dp
     pdens=0.0_dp
     ldos=0.0_dp
     cur=0.0_dp
  end if
elseif ( ff /= 0 ) then
     write(*,*)'ignoring E =',E
     tr=0.0_dp
     tre=0.0_dp
     ndens=0.0_dp
     pdens=0.0_dp
     ldos=0.0_dp
     cur=0.0_dp
  end if

end subroutine RGF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Function dferm(a)
    real(dp) a,dferm
    dferm=-1.0_dp/(1.0_dp+DExp(-a))
  End Function dferm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Sancho-Rubio
subroutine sancho(nm,E,H00,H10,G00)
  implicit none
  integer, intent(IN) :: nm
  integer :: i
  integer :: nmax=500
  COMPLEX(DP) :: z
  real(dp), intent(IN) :: E
  REAL(DP) :: TOL=1.0D-15, error

  COMPLEX(DP), INTENT(IN) :: H00(nm,nm), H10(nm,nm)
  COMPLEX(DP), INTENT(OUT) :: G00(nm,nm)

  COMPLEX(DP), ALLOCATABLE :: A(:,:), B(:,:), C(:,:)
  COMPLEX(DP), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)

  Allocate( H_BB(nm,nm) )
  Allocate( H_SS(nm,nm) )
  Allocate( H_01(nm,nm) )
  Allocate( H_10(nm,nm) )
  Allocate( Id(nm,nm) )
  Allocate( A(nm,nm) )
  Allocate( B(nm,nm) )
  Allocate( C(nm,nm) )
  
  z = dcmplx(E,eta)
  
  Id=0.0_dp
  do i=1,nm
     Id(i,i)=1.0_dp
  enddo

  H_BB = 0.5_dp*(H00+transpose(dconjg(H00))) 
  H_10 = H10
  H_01 = TRANSPOSE( DCONJG( H10 ) )
  H_SS = H_BB


  do i = 1,nmax

     A = z*Id - H_BB
     call invert(A,nm,nm)
     call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,A,nm,beta,B,nm)
     call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,A,nm,beta,C,nm)
     call zgemm('n','n',nm,nm,nm,alpha,C,nm,H_10,nm,beta,A,nm)
     H_SS = H_SS + A
     H_BB = H_BB + A
     call zgemm('n','n',nm,nm,nm,alpha,B,nm,H_01,nm,beta,A,nm)
     H_BB = H_BB + A
     call zgemm('n','n',nm,nm,nm,alpha,C,nm,H_01,nm,beta,A,nm)
     H_01 = A
     call zgemm('n','n',nm,nm,nm,alpha,B,nm,H_10,nm,beta,A,nm)
     H_10 = A

     error = abs(sum (sum(H_10, dim=2), dim=1) ) + abs(sum (sum(H_01, dim=2), dim=1) )

     IF ( abs(error) .lt. TOL ) THEN
        EXIT
     END IF
     
  enddo

  if (i==nmax)then
     write(*,*)'pb nmax',E,error
     G00=0.0_dp
  else
     !!!! G00 = dcmplx(E,eta)*Id - H_SS
     G00 = E*Id - H_SS
     call invert(G00,nm,nm)
  end if


  Deallocate( A )
  Deallocate( B )
  Deallocate( C )

  Deallocate( H_BB )
  Deallocate( H_SS )
  Deallocate( H_01 )
  Deallocate( H_10 )
  Deallocate( Id )

end subroutine sancho



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Old Sancho-Rubio
subroutine oldsancho(nm,E,H00,H10,G00)

  integer, intent(IN) :: nm
  integer i,nmax
  COMPLEX(DP) :: z
  real(dp), intent(IN) :: E
  REAL(DP) :: TOL=1.0D-12, error

  COMPLEX(DP), INTENT(IN) :: H00(nm,nm), H10(nm,nm)
  COMPLEX(DP), INTENT(OUT) :: G00(nm,nm)

  COMPLEX(DP), ALLOCATABLE :: A(:,:), B(:,:), C(:,:)
  COMPLEX(DP), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)

  Allocate( H_BB(nm,nm) )
  Allocate( H_SS(nm,nm) )
  Allocate( H_01(nm,nm) )
  Allocate( H_10(nm,nm) )
  Allocate( Id(nm,nm) )
  Allocate( A(nm,nm) )
  Allocate( B(nm,nm) )
  Allocate( C(nm,nm) )

  nmax=100!500


  z = E+dcmplx(0.0_dp,eta)

  Id=0.0_dp
  do i=1,nm
     Id(i,i)=1.0_dp
  enddo

  H_BB = H00
  H_10 = H10
  H_01 = TRANSPOSE( DCONJG( H10 ) )
  H_SS = H00

  do i = 1, nmax

     A = z*Id - H_BB
     call invert(A,nm,nm)
     call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_10,nm,beta,B,nm)
     call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm)
     H_SS = H_SS + C
     H_BB = H_BB + C
     call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,C,nm)
     call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_01,nm,beta,B,nm)
     call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,A,nm)
     H_BB = H_BB + A
     H_10 = C
     call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm)
     H_01 = C
     error=traccia(abs(C))
    ! write(*,*)i,error
     IF ( abs(error) .lt. TOL ) THEN
!       write(*,*) 'SR: Exited, abs(error)=',abs(error)
       EXIT
     END IF

     IF (i == nmax) THEN
        write(*,*) 'SEVERE warning: nmax reached in sancho!!!',error,E
        EXIT
     END IF

  enddo

  if (i==nmax)then
     G00=0.0_dp
  else
     G00 = E*Id - H_SS
     call invert(G00,nm,nm)
  end if
!  GBB = z*Id - H_BB
!  call invert(GBB,nm,nm)

  Deallocate( A )
  Deallocate( B )
  Deallocate( C )

  Deallocate( H_BB )
  Deallocate( H_SS )
  Deallocate( H_01 )
  Deallocate( H_10 )
  Deallocate( Id )

end subroutine oldsancho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function traccia(A)
  real(dp) :: A(:,:)
  real(dp) :: traccia
  integer :: i,j

  traccia=dble(sum( a(:,:), &
       reshape( (/((i==j,i=1,size(a,1)),j=1,size(a,2))/), &
       (/size(a,1),size(a,2)/) ) ) ) 

end function traccia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Real(dp) Function ferm(a)
  real(dp) a
  ferm = (1.0_dp+DExp(a))**(-1.0_dp)
End Function ferm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine invert(A,n,np)

  implicit none
      
  integer :: info,n,np
      
  integer, dimension(np) :: ipiv
      
  complex(dp) :: A(n,n)    
  complex(dp) :: work(np*np)   

  call zgetrf(n,n,A,n,ipiv,info)
  call zgetri(n,A,n,ipiv,work,np*np,info)

end subroutine invert


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

SUBROUTINE surf_interp(iy,iz,mat,interp)
  implicit none

  integer,intent(in)   :: iy, iz
  real(dp),intent(in)  :: mat(2,2)
  real(dp),intent(out) :: interp 
  real(dp)             ::  vecl(2), vecr(2)
  
  vecl(1) = ceiling(dble(iz-1)*dz/deltaz+1.0d-3+1)*deltaz - (iz-1)*dz
  vecl(2) = - ceiling(dble(iz-1)*dz/deltaz+1.0d-3)*deltaz + (iz-1)*dz
  
  call dgemv('n', 2, 2, dalpha, mat, 2, vecl, 1, dbeta, vecr, 1)
  
  vecl(1) = ceiling(dble(iy-1)*dy/deltay+1.0d-3+1)*deltay - (iy-1)*dy
  vecl(2) = - ceiling(dble(iy-1)*dy/deltay+1.0d-3)*deltay + (iy-1)*dy

  interp=(vecl(1)*vecr(1)+vecl(2)*vecr(2))/deltaz/deltay
 
END SUBROUTINE surf_interp
         

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE negf


