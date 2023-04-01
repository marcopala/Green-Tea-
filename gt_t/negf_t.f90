! Copyright or © or Copr. Marco Pala (February 24, 2022)

! e-mail: marco.pala@c2n.upsaclay.fr ;   marco.pala@c2n.upsaclay.fr

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

! In this respect, the user's atention is drawn to the risks associated
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

MODULE negf
  
USE static  
USE indata
use mpi

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

INTEGER                  :: i, j, l, n, ix, iy, iz, ip, ii, xx, ll, pp, nn, nde
INTEGER                  :: Nop,ref_index, nmax, ee, nee, nsol, SCBA_iter, comp_ch(NCX_D)

REAL(DP)                 :: Gcon, en, epsilon, emin_local
REAL(DP), ALLOCATABLE    :: con(:),cone(:),conb(:)

COMPLEX(DP), ALLOCATABLE :: g_lesser_diag_local(:,:,:), g_lesser(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: g_greater_diag_local(:,:,:), g_greater(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: g_r_diag_local(:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_lesser_ph(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_greater_ph(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_lesser_ph_prev(:,:,:,:,:)
COMPLEX(DP), ALLOCATABLE :: sigma_greater_ph_prev(:,:,:,:,:)
!COMPLEX(DP), ALLOCATABLE :: sigma_r_ph_prev(:,:,:,:)

REAL(DP), ALLOCATABLE    :: emin_yz(:), emax_yz(:), SCBA_error(:), SCBA_x(:), SCBA_f(:)
REAL(DP)                 :: SCBA_scalar
REAL(DP)                 :: n_bose_g

REAL(DP), ALLOCATABLE    :: degeneracy(:), trans(:,:,:,:), cur(:,:,:,:)
REAL(DP), ALLOCATABLE    :: ldos(:,:,:,:), zdos(:,:,:,:), ndos(:,:,:,:), pdos(:,:,:,:), R_in(:,:,:,:), R_out(:,:,:,:)
REAL(DP)                 :: tr,tre,ttr,sumt,kappax,tmp,n_2D
!COMPLEX(DP)              :: tmp


COMPLEX(DP), allocatable :: Hi(:,:,:,:),pot(:,:)
COMPLEX(dp), allocatable :: dosn(:,:,:,:),dosp(:,:,:,:)
COMPLEX(dp), allocatable :: U(:,:), GBB(:,:), Gcc(:,:), A(:,:), B(:,:), C(:,:),D(:,:),CC(:,:,:),DD(:,:,:)
REAL(DP), allocatable    :: subband(:,:,:), neutr(:,:)
REAL(dp), allocatable    :: pt(:), CR(:,:), Ry(:), Rz(:), dens_yz(:,:),dens_z(:,:)
REAL(dp), allocatable    :: omega(:,:), derror(:,:,:)

real(4) :: t1,t2,t3,t4

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

open(unit=10,file=TRIM(outdir)//'poty_'//TRIM(STRINGA(ext_iter))//'.dat',status='unknown')
do j=Ntot_y/2,NTOT_Y/2
   do i=1,NTOT_X
      do l=1,NTOT_Z
          write(10,*)i,l,POT3D(i,j,l)
       enddo
    write(10,*)
    enddo
 enddo
close(10)

open(unit=10,file=TRIM(outdir)//'potx.dat',status='unknown')
 do l=1,NTOT_Z
    do j=1,Ntot_y
       do i=NTOT_X/2,NTOT_X/2
          write(10,*)j,l,POT3D(i,j,l)
       enddo
    enddo
    write(10,*)
 enddo
close(10)

open(unit=10,file=TRIM(outdir)//'potz.dat',status='unknown')
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

comp_ch=0
do xx=1,NCX_D
   if( schottky_source .and. xx <= source_len/Ndeltax ) comp_ch(xx)=11
   if( schottky_drain  .and. xx > (source_len+2*spacer+gate_len)/Ndeltax ) comp_ch(xx)=12
end do

allocate(Hi(NMAX,NMAX,ncx_d,nkyz))
Hi=0.0_dp

allocate(degeneracy(nkyz))
do iyz=1,NKYZ !!! loop over transverse k-vectors
   if(k_selec(iyz))then
   
   degeneracy(iyz)=deg_kyz(iyz)*g_spin
   write(*,*) 
   write(*,*)iyz,'degeneracy',deg_kyz(iyz),degeneracy(iyz)
   write(*,*)
   if(allocated(kgt))deallocate(kgt)
   allocate(KGt(4,1*Ngt))
   KGt(1:4,1:1*Ngt)=KGt_kyz(1:4,1:1*Ngt,iyz)


if(.not. onlyT .or. in_pot)then
      
write(*,*)'Transforming the potential'
   t1=SECNDS(0.0)
allocate(U(Ngt,(Ny+1)*(Nz+1)))
do iy=1,NY+1
   do iz=1,NZ+1
      j=iy+(iz-1)*(NY+1)
      U(1:Ngt,j)=exp(cmplx(0.0_dp,-1.0_dp)*KGt_kyz(2,1:Ngt,iyz)*2.0_dp*pi/ac*Ry(iy)+&
           cmplx(0.0_dp,-1.0_dp)*KGt_kyz(3,1:Ngt,iyz)*2.0_dp*pi/ac*Rz(iz))/sqrt(dble((Ny+1)*(Nz+1)))
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

   
   !!!! TRANSFORMING POT INTO THE Bloch state BASIS

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

Hi(1:NM(xx),1:NM(xx),xx,iyz)=(Hi(1:NM(xx),1:NM(xx),xx,iyz)+transpose(dconjg(Hi(1:NM(xx),1:NM(xx),xx,iyz))))/2.0_dp

end do
!$omp end do nowait

deallocate(pt)
!deallocate(pot)
deallocate(Gcc)
deallocate(A)
deallocate(CR)

!$omp end parallel

   if(allocated(U))deallocate(U)

   t2=SECNDS(t1)
   write(*,*)'Potential transformed in ',t2,'s'
end if


end if
end do !fine loop kyz

nsol=nsolv+nsolc
allocate(subband(nsol,ncx_d,nkyz))
allocate(emin_yz(NKYZ),emax_yz(NKYZ))
emin_yz=0.0_dp
emax_yz=0.0_dp


!! neutrality point for the semiconductors
if(.not.allocated(neutr))allocate(neutr(ncx_d,nkyz))

do iyz=1,NKYZ
   if(k_selec(iyz))then
      
call omp_set_num_threads(NCX_D) 

!$omp parallel default(none) private(xx,ref_index,kappax,A,B) &
!$omp shared(NCX_D,iyz,NM,nrx,ngt,npol,imat,ihet,nband_val,ULCBB,Hi,Si,HL,TL,nsolv,nsolc,subband,kc_min,kv_max,chtype)
!$omp do
do xx=1,ncx_d
   ref_index=nband_val(imat(xx))!!!+off_k_nvb(iyz)
   
   allocate(B(NM(xx),NM(xx)))
   B(1:NM(xx),1:NM(xx))=Si(iyz,imat(xx))%H(1:NM(xx),1:NM(xx))
      
   kappax=0.0d0

   allocate(A(NM(xx),NM(xx)))
   Hi(1:NM(xx),1:NM(xx),xx,iyz)=Hi(1:NM(xx),1:NM(xx),xx,iyz)+HL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx))
   if(nband_val(imat(xx))>0)then
      kappax=kv_max(iyz,imat(xx))
      A(1:NM(xx),1:NM(xx))=Hi(1:NM(xx),1:NM(xx),xx,iyz)+&
           (TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx)))*exp(cmplx(0.0_dp,1.0_dp)*kappax*2.0_dp*pi)+&
           transpose(dconjg(TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx))))*exp(cmplx(0.0_dp,-1.0_dp)*kappax*2.0_dp*pi)
      call SUB_DEF_Z0_GEN(ref_index-nsolv+1,ref_index,NM(xx),A,B, subband(1:nsolv,xx,iyz))
   end if
   kappax=kc_min(iyz,imat(xx))
   A(1:NM(xx),1:NM(xx))=Hi(1:NM(xx),1:NM(xx),xx,iyz)+&
        (TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx)))*exp(cmplx(0.0_dp,1.0_dp)*kappax*2.0_dp*pi)+&
        transpose(dconjg(TL(iyz,imat(xx))%H(1:NM(xx),1:NM(xx))))*exp(cmplx(0.0_dp,-1.0_dp)*kappax*2.0_dp*pi)
   call SUB_DEF_Z0_GEN(ref_index+1,ref_index+nsolc,NM(xx),A, B, subband(nsolv+1:nsolv+nsolc,xx,iyz))

   deallocate(A,B)
end do
!$omp end do nowait
!$omp end parallel

if(nsolv==0 .and. chtype /= 'n')then
   write(*,*)'pb w nsolv or chtype'
end if


do xx=1,ncx_d
   if(nsolv>0)then
      neutr(xx,iyz)=(subband(nsolv,xx,iyz)+subband(nsolv+1,xx,iyz))/2.0_dp
   else
      neutr(xx,iyz)=subband(nsolv+1,xx,iyz)-100.0d-3
   end if
enddo

  epsilon=max(5*Nop_g*Eop,5*(BOLTZ*TEMP))

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

  write(*,*)
  write(*,*)iyz,'emin =',emin_yz(iyz),'emax =',emax_yz(iyz)
  
end if
end do  ! end of loop over kyz


  emax=maxval(emax_yz(:))
  emin=minval(emin_yz(:))

  if(onlyT .and. .not. in_pot)then
     emin=min(emin,min(mus,mud))-NKT*(BOLTZ*TEMP)
     emax=max(emax,max(mus,mud))+NKT*(BOLTZ*TEMP)
  end if
  
  deallocate(emin_yz,emax_yz)

  Nop=FLOOR((emax-emin)/Eop+0.5_dp)

  write(*,*)
  write(*,*)'GLOBAL ENERGY RANGE'
  write(*,*)'Emin=',emin,'Emax=',emax
  write(*,*)'Nop=',Nop
  write(*,*)
  write(*,*)'mus=',mus
  write(*,*)'mud=',mud
  write(*,*)
  
  allocate(flag(Nop,NKYZ))
  

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
  if(phonons)allocate(R_in(1:Nsub,1:Nop,1:Ncx_d,1:NKYZ),R_out(1:Nsub,1:Nop,1:Ncx_d,1:NKYZ))

  ! ==========  Start of parallel resolution ===========================
  !write(*,*) 'starting the parallel loop'
  
  call omp_set_num_threads(min(Nomp,Nop)) !!! this sets the environment variable OMP_NUM_THREADS 

  ALLOCATE(sigma_greater_ph(1:Nop,1:NMAX,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(sigma_lesser_ph(1:Nop,1:NMAX,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(sigma_greater_ph_prev(1:Nop,1:NMAX,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(sigma_lesser_ph_prev(1:Nop,1:NMAX,1:NMAX,1:Ncx_d,1:nkyz))
  ALLOCATE(SCBA_error(1:NMAX),SCBA_x(1:NMAX),SCBA_f(1:NMAX))

  do ee=1,Nsub

     emin_local=en_global(ee) 
     flag=0

if(phonons)then 

     sigma_greater_ph_prev=0.0_dp
     sigma_lesser_ph_prev=0.0_dp
     

     SCBA_scalar=1.0_dp
     SCBA_iter=0  

     ALLOCATE(derror(Nop,Ncx_D,nkyz),omega(Nop,Ncx_D))
     derror=1.0_dp
     omega=0.0_dp
     ALLOCATE(g_lesser(Nop,NMAX,NMAX,ncx_d,NKYZ),g_greater(Nop,NMAX,NMAX,Ncx_d,NKYZ))

        g_lesser  = 0.0d0
        g_greater = 0.0d0
        
     DO WHILE((SCBA_scalar.gt.SCBA_tolerance).and.(SCBA_iter.lt.SCBA_max_iter))
   
        t1=SECNDS(0.0)
        SCBA_iter=SCBA_iter+1
        
        
        do iyz=1,NKYZ !loop over ( kyz' = kyz - qyz ) ! this is the idex of k_yz' = k_yz - q_yz 
           if(k_selec(iyz))then
              
            
              !$omp parallel default(none)  private(jyz,ii,ll,ix,nee,xx,pp,nn,EN,tr,tre,g_lesser_diag_local,g_greater_diag_local,g_r_diag_local,A,B,C,D) &
              !$omp shared(scba_iter,ee,iyz,Nop,nm,imat,ncx_d,ncy,ncz,nmax,nkyz,nkx,nqmodes,k_vec,emin,emin_local,Eop,mus,mud,hi,&
              !$omp sigma_lesser_ph_prev,sigma_greater_ph_prev,cur,ldos,ndos,pdos,zdos,chtype,neutr,omega_q,ind_q,&
              !$omp degeneracy,g_spin,w,temp,trans,g_lesser,g_greater,el_ph_mtrx,flag,k_selec,dfpt,Si,derror,scba_tolerance)

        ALLOCATE(g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
        ALLOCATE(g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
        ALLOCATE(g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
        allocate(A(NMAX,NMAX),B(NMAX,NMAX),C(NMAX,NMAX),D(NMAX,NMAX))
    
        !$omp do
        do nee=1,Nop
        
        if(maxval(abs(derror(nee,:,iyz))) > scba_tolerance)then
           EN=emin+emin_local+dble(nee-1)*Eop
     
           call RGF(scba_iter,NMAX,flag(nee,iyz),EN,mus,mud,Hi(1:NMAX,1:NMAX,1:Ncx_d,iyz),&
                sigma_lesser_ph_prev(nee,1:nmax,1:NMAX,1:Ncx_d,iyz),sigma_greater_ph_prev(nee,1:nmax,1:NMAX,1:Ncx_d,iyz),&
                g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d),&
                tr,tre,cur(ee,nee,1:Ncx_d-1,iyz))

           do xx=1,ncx_d
              g_lesser (nee,1:nm(xx),1:nm(xx),xx,iyz) = g_lesser_diag_local (1:nm(xx),1:nm(xx),xx)
              g_greater(nee,1:nm(xx),1:nm(xx),xx,iyz) = g_greater_diag_local(1:nm(xx),1:nm(xx),xx)
           end do
           
    !!!! DEBUGGING       FILES
!!$           DO xx=1,Ncx_d
!!$              write(3200+scba_iter+10*(ee-1),*)xx,En, traccia(aimag(sigma_lesser_ph_prev(nee,1:nm(xx),1:nm(xx),xx,iyz)))
!!$              write(3300+scba_iter+10*(ee-1),*)xx,En, traccia(-aimag(sigma_greater_ph_prev(nee,1:nm(xx),1:nm(xx),xx,iyz)))
!!$              write(3400+scba_iter+10*(ee-1),*)xx,En, traccia(aimag(g_lesser(nee,1:nm(xx),1:nm(xx),xx,iyz)))
!!$              write(3500+scba_iter+10*(ee-1),*)xx,En, traccia(-aimag(g_greater(nee,1:nm(xx),1:nm(xx),xx,iyz)))
!!$              write(3600+scba_iter+10*(ee-1),*)xx,En, traccia(-aimag(g_r_diag_local(1:nm(xx),1:nm(xx),xx)))
!!$           end do
!!$          write(3200+scba_iter+10*(ee-1),*)
!!$          write(3300+scba_iter+10*(ee-1),*)
!!$          write(3400+scba_iter+10*(ee-1),*)
!!$          write(3500+scba_iter+10*(ee-1),*)
!!$          write(3600+scba_iter+10*(ee-1),*)

        end if
     end do  ! end do nee
     !$omp end do nowait
     DEALLOCATE(g_lesser_diag_local)
     DEALLOCATE(g_greater_diag_local)
     DEALLOCATE(g_r_diag_local)
     DEALLOCATE(A,B,C,D)
     !$omp end parallel
        
    !!!! DEBUGGING       FILES
!!$      close(3200+scba_iter+10*(ee-1))
!!$      close(3300+scba_iter+10*(ee-1))
!!$      close(3400+scba_iter+10*(ee-1))
!!$      close(3500+scba_iter+10*(ee-1))
!!$      close(3600+scba_iter+10*(ee-1))
      
     end if
  end do !End of the loop over ( kyz + qyz )
  t2=SECNDS(t1)
  write(*,*)'time spent to compute G',t2

  !!!!!  self-energies calculation
  sigma_lesser_ph=0.0d0
  sigma_greater_ph=0.0d0
  
  t1=SECNDS(0.0)

  if (dfpt) then
        
     !$omp parallel default(none)  private(NdE,ll,ix,ii,iyz,jyz,nee,xx,A,B,CC,DD) shared(Nop,nmax,nkx,nkyz,g_spin,ncx_d,ncy,ncz,nm, &
     !$omp sigma_lesser_ph,sigma_greater_ph,g_lesser,g_greater,SCBA_iter,Eop,temp,omega_q,ind_q,nqmodes,k_selec,dfpt,el_ph_mtrx,k_vec,imat,degeneracy)
     allocate(A(NMAX,NMAX),B(NMAX,NMAX),CC(NMAX,NMAX,Ncx_d),DD(NMAX,NMAX,Ncx_d))     
     DO nee=1,Nop
        do jyz = 1,NKYZ   ! index of k_yz 
           if(k_selec(jyz))then
              do iyz =1,nkyz ! index of k_yz'
                 if(k_selec(iyz))then
                    ii = ind_kyz(  k_vec(2:3,jyz) - k_vec(2:3,iyz) ) ! this is the idex of q_yz = k_yz - k_yz'
                    do ll=1,nqmodes
                       do ix=1,nkx
                          NdE=ceiling(omega_q(ind_q(ix,ii),ll)/Eop)
                          !$omp do
                          DO xx=1,Ncx_d
                             A(1:nm(xx),1:nm(xx))=el_ph_mtrx(jyz,ix,ii,imat(xx))%M(ll,1:nm(xx),1:nm(xx))
                             IF(nee <= Nde)THEN
                                IF(nee <= Nop-Nde)THEN
                                   ! E+Eop
                                   CC(1:nm(xx),1:nm(xx),xx)=g_lesser(nee+NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP)) + 1.0_dp )
                                   DD(1:nm(xx),1:nm(xx),xx)=g_greater(nee+NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP)) ) 
                                END IF
                             ELSE
                                IF(nee <= Nop-Nde)THEN
                                   ! E+Eop and E-Eop
                                   CC(1:nm(xx),1:nm(xx),xx)=g_lesser(nee+NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP)) + 1.0_dp ) + &
                                        g_lesser(nee-NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP))  )
                                   DD(1:nm(xx),1:nm(xx),xx)=g_greater(nee+NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP)) ) + &
                                        g_greater(nee-NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP)) + 1.0_dp )
                                ELSE 
                                   ! E-Eop
                                   CC(1:nm(xx),1:nm(xx),xx)=g_lesser(nee-NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP)) )
                                   DD(1:nm(xx),1:nm(xx),xx)=g_greater(nee-NdE,1:nm(xx),1:nm(xx),xx,iyz)*( bose(NdE*Eop/(BOLTZ*TEMP)) + 1.0_dp ) 
                                END IF
                             END IF
                   call zgemm('n','n',nm(xx),nm(xx),nm(xx),alpha,A(1:nm(xx),1:nm(xx)),nm(xx),CC(1:nm(xx),1:nm(xx),xx),nm(xx),beta,B(1:nm(xx),1:nm(xx)),nm(xx)) 
                   call zgemm('n','c',nm(xx),nm(xx),nm(xx),alpha,B(1:nm(xx),1:nm(xx)),nm(xx),A(1:nm(xx),1:nm(xx)),nm(xx),beta,CC(1:nm(xx),1:nm(xx),xx),nm(xx))
                   call zgemm('n','n',nm(xx),nm(xx),nm(xx),alpha,A(1:nm(xx),1:nm(xx)),nm(xx),DD(1:nm(xx),1:nm(xx),xx),nm(xx),beta,B(1:nm(xx),1:nm(xx)),nm(xx)) 
                   call zgemm('n','c',nm(xx),nm(xx),nm(xx),alpha,B(1:nm(xx),1:nm(xx)),nm(xx),A(1:nm(xx),1:nm(xx)),nm(xx),beta,DD(1:nm(xx),1:nm(xx),xx),nm(xx))
                             sigma_lesser_ph(nee,1:nm(xx),1:nm(xx),xx,jyz) =sigma_lesser_ph(nee,1:nm(xx),1:nm(xx),xx,jyz) +&
                                  degeneracy(iyz)/g_spin/dble(NCY*NCZ*nkx) * CC(1:nm(xx),1:nm(xx),xx)
                             sigma_greater_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)=sigma_greater_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)+&
                                  degeneracy(iyz)/g_spin/dble(NCY*NCZ*nkx) * DD(1:nm(xx),1:nm(xx),xx)
                          END DO
                          !$omp end do nowait
                                                    
                       end do
                    end do
                 end if
              end do
           end if
        end do
  
     END DO

  
     DEALLOCATE(A,B,CC,DD)

     !$omp end parallel
  end if
                    
  if (.not. dfpt) then


     !$omp parallel default(none)  private(ii,iyz,jyz,nee,xx,A,B,CC,DD) shared(Nop,nmax,nkyz,g_spin,ncx_d,ncy,ncz,nm,Nop_g,sigma_lesser_ph,sigma_greater_ph,sigma_lesser_ph_prev,sigma_greater_ph_prev, &
     !$omp g_lesser,g_greater,SCBA_iter,Dop_g,n_bose_g,Dac,Eop,temp,k_selec,dfpt,el_ph_mtrx,k_vec,imat,degeneracy,derror,scba_tolerance)
     
     allocate(A(NMAX,NMAX),B(NMAX,NMAX),CC(NMAX,NMAX,Ncx_d),DD(NMAX,NMAX,Ncx_d))
    
     DO nee=1,Nop
        do jyz = 1,NKYZ   ! index of k_yz 
           if(k_selec(jyz))then

        if(maxval(abs(derror(nee,:,jyz))) > scba_tolerance)then
          
              do iyz = 1,NKYZ ! this is the idex of k_yz' 
                 if(k_selec(iyz))then
                    ii = ind_kyz( k_vec(2:3,jyz) - k_vec(2:3,iyz) ) ! this is the idex of q_yz = k_yz - k_yz'
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! acoustic + optical phonons !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
                    n_bose_g=1.0_dp/(EXP((Nop_g*Eop)/(BOLTZ*TEMP)) - 1.0_dp)
              
                    !$omp do
                    do xx=1,ncx_d
                       
                    IF(nee.le.Nop_g)THEN
                       IF(nee.le.Nop-Nop_g)THEN
                          ! E+Eop_g
                       
                             A(1:nm(xx),1:nm(xx))=el_ph_mtrx(jyz,1,ii,imat(xx))%M(1,1:nm(xx),1:nm(xx))
                             CC(1:nm(xx),1:nm(xx),xx)=Dop_g*g_lesser (nee+Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*( n_bose_g + 1.0_dp ) +&
                                  Dac*g_lesser (nee,1:nm(xx),1:nm(xx),xx,iyz)

                             DD(1:nm(xx),1:nm(xx),xx)=Dop_g*g_greater (nee+Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*( n_bose_g ) +&
                                  Dac*g_greater (nee,1:nm(xx),1:nm(xx),xx,iyz)
                       
                          END IF
                       
                    ELSE
                       
                       IF(nee.le.Nop-Nop_g)THEN
                          ! E+Eop_g and E-Eop_g
                                                     
                             A(1:nm(xx),1:nm(xx))=el_ph_mtrx(jyz,1,ii,imat(xx))%M(1,1:nm(xx),1:nm(xx))
                             CC(1:nm(xx),1:nm(xx),xx)=Dop_g*g_lesser (nee+Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*( n_bose_g + 1.0_dp ) +&
                                  Dop_g*g_lesser (nee-Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*  ( n_bose_g )  + &
                                  Dac*g_lesser (nee,1:nm(xx),1:nm(xx),xx,iyz)
                                                    
                             DD(1:nm(xx),1:nm(xx),xx)=Dop_g*g_greater (nee+Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*( n_bose_g ) +&
                                  Dop_g*g_greater (nee-Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*  ( n_bose_g +1.0_dp )  + &
                                  Dac*g_greater (nee,1:nm(xx),1:nm(xx),xx,iyz)
                                                                           
                       ELSE 
                          ! E-Eop_g
                             A(1:nm(xx),1:nm(xx))=el_ph_mtrx(jyz,1,ii,imat(xx))%M(1,1:nm(xx),1:nm(xx))
                             CC(1:nm(xx),1:nm(xx),xx)=Dop_g*g_lesser (nee-Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*  ( n_bose_g )  + &
                                  Dac*g_lesser (nee,1:nm(xx),1:nm(xx),xx,iyz)
                            
                             DD(1:nm(xx),1:nm(xx),xx)=Dop_g*g_greater (nee-Nop_g,1:nm(xx),1:nm(xx),xx,iyz)*  ( n_bose_g +1.0_dp )  + &
                                  Dac*g_greater (nee,1:nm(xx),1:nm(xx),xx,iyz)

                          END IF
                       END IF
                       
                    call zgemm('n','n',nm(xx),nm(xx),nm(xx),alpha,A(1:nm(xx),1:nm(xx)),nm(xx),CC(1:nm(xx),1:nm(xx),xx),nm(xx),beta,B(1:nm(xx),1:nm(xx)),nm(xx)) 
                    call zgemm('n','c',nm(xx),nm(xx),nm(xx),alpha,B(1:nm(xx),1:nm(xx)),nm(xx),A(1:nm(xx),1:nm(xx)),nm(xx),beta,CC(1:nm(xx),1:nm(xx),xx),nm(xx))                                
                    call zgemm('n','n',nm(xx),nm(xx),nm(xx),alpha,A(1:nm(xx),1:nm(xx)),nm(xx),DD(1:nm(xx),1:nm(xx),xx),nm(xx),beta,B(1:nm(xx),1:nm(xx)),nm(xx)) 
                    call zgemm('n','c',nm(xx),nm(xx),nm(xx),alpha,B(1:nm(xx),1:nm(xx)),nm(xx),A(1:nm(xx),1:nm(xx)),nm(xx),beta,DD(1:nm(xx),1:nm(xx),xx),nm(xx))
                    
                    sigma_lesser_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)=sigma_lesser_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)+&
                         degeneracy(iyz)/g_spin/dble(NCY*NCZ)*&
                         CC(1:nm(xx),1:nm(xx),xx)
                    sigma_greater_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)=sigma_greater_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)+&
                         degeneracy(iyz)/g_spin/dble(NCY*NCZ)*&
                         DD(1:nm(xx),1:nm(xx),xx)
                             
                 end do
                 !$omp end do nowait
                 
              end if
           end do

           else
        sigma_lesser_ph(nee,1:nmax,1:nmax,1:ncx_d,jyz)  = sigma_lesser_ph_prev(nee,1:nmax,1:nmax,1:ncx_d,jyz)
        sigma_greater_ph(nee,1:nmax,1:nmax,1:ncx_d,jyz) = sigma_greater_ph_prev(nee,1:nmax,1:nmax,1:ncx_d,jyz)
     end if
  end if
end do
     
end DO
  
     DEALLOCATE(A,B,CC,DD)

    !$omp end parallel
  end if

t2=SECNDS(t1)
WRITE(*,*)'TIME SPENT TO COMPUTE SIGMA_LESSER_PH (s)',t2
        

  if(scba_scalar > 1.0d-3)then
     omega=0.0_dp
  else
     omega=1.0_dp-scba_alpha
  end if
  
  SCBA_f=0.0_dp
  SCBA_x=0.0_dp
  derror=0.0_dp
  !$omp parallel default(none)  private(nee,xx,i,jyz) shared(Nop,nkyz,ncx_d,k_selec,NM,SCBA_f,SCBA_x,derror, &
  !$omp sigma_lesser_ph,sigma_greater_ph,sigma_lesser_ph_prev,sigma_greater_ph_prev)

  !$omp do
  do nee=1,Nop
     do jyz=1,NKYZ
        if(k_selec(jyz))then
           DO xx=1,Ncx_D
              do i=1,NM(xx)
                 SCBA_f(i)=SCBA_f(i)+((dimag(sigma_lesser_ph(nee,i,i,xx,jyz)-sigma_lesser_ph_prev(nee,i,i,xx,jyz)))**2 +&
                      (dimag(sigma_greater_ph(nee,i,i,xx,jyz)-sigma_greater_ph_prev(nee,i,i,xx,jyz)))**2   )
                 
                 SCBA_x(i)=SCBA_x(i)+((dimag(sigma_lesser_ph(nee,i,i,xx,jyz)))**2 + (dimag(sigma_greater_ph(nee,i,i,xx,jyz)))**2)
              end do
              derror(nee,xx,jyz)= derror(nee,xx,jyz)+abs(traccia(dimag(sigma_lesser_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)-sigma_lesser_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz)))-traccia(dimag(sigma_greater_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)-sigma_greater_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz))))
           end do
        end if
     END DO
  end do
  !$omp end do  nowait
  !$omp end parallel
  
  SCBA_error(:)=SCBA_f(:)/(SCBA_x(:)+1.0d-16)
  SCBA_scalar=sum(sqrt(SCBA_error(:)))/dble(NMAX)


  !$omp parallel default(none)  private(nee,xx,jyz) shared(Nop,nkyz,ncx_d,k_selec,NM,omega,&
  !$omp sigma_lesser_ph,sigma_greater_ph,sigma_lesser_ph_prev,sigma_greater_ph_prev)
  !$omp do
  DO xx=1,Ncx_D
     DO jyz=1,NKYZ
        if(k_selec(jyz))then
           do nee=1,Nop
           sigma_lesser_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz)=sigma_lesser_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz) + &
                (1.0d0-omega(nee,xx))*(sigma_lesser_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)-sigma_lesser_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz))     
           sigma_greater_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz)=sigma_greater_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz) + &
                (1.0d0-omega(nee,xx))*(sigma_greater_ph(nee,1:nm(xx),1:nm(xx),xx,jyz)-sigma_greater_ph_prev(nee,1:nm(xx),1:nm(xx),xx,jyz))
        end do
        end if
     END DO
  END DO
  !$omp end do nowait
  !$omp end parallel

  
  write(*,*)ee,'scba_iter',SCBA_iter,SCBA_scalar
  
END DO ! END of the SCBA di Born
!!!stop

  DEALLOCATE(g_lesser,g_greater)
  DEALLOCATE(omega,derror)

end if ! endif phonons

do iyz=1,NKYZ !loop over kyz
   if(k_selec(iyz))then
   !$omp parallel default(none) &
   !$omp private(nee,xx,i,j,ix,iy,iz,EN,tr,tre,tmp,g_lesser_diag_local,g_greater_diag_local,g_r_diag_local,A) &
   !$omp shared(scba_iter,ee,iyz,Nop,ncx_d,nkyz,nmax,nm,nrx,nrz,imat,emin,emin_local,Eop,mus,mud,hi,&
   !$omp AJ,BJ,CJ,sigma_lesser_ph_prev,sigma_greater_ph_prev, &
   !$omp cur,ldos,ndos,pdos,zdos,chtype,neutr,degeneracy,w,temp,con,cone,conb,trans,dosn,dosp, &
   !$omp flag,phonons,R_in,R_out)

  ALLOCATE(g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
  ALLOCATE(g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
  ALLOCATE(g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d))
  ALLOCATE(A(1:nmax,1:nmax))
        
  !$omp do 
  do nee=1,Nop
     EN=emin+emin_local+dble(nee-1)*Eop
         
     if(phonons)then
        !!! elastic calculation with the same band profile of the inelastic calculation
        call RGF(scba_iter,NMAX,flag(nee,iyz),EN,mus,mud,Hi(1:NMAX,1:NMAX,1:Ncx_d,iyz),&
             0.0_dp*sigma_lesser_ph_prev(nee,1:NMAX,1:NMAX,1:Ncx_d,iyz),0.0_dp*sigma_greater_ph_prev(nee,1:NMAX,1:NMAX,1:Ncx_d,iyz),&
             g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d),&
             tr,tre,cur(ee,nee,1:Ncx_d-1,iyz))
        conb(iyz) =conb(iyz) +degeneracy(iyz)*w(ee)*(tr)
     end if
     
     call RGF(scba_iter,NMAX,flag(nee,iyz),EN,mus,mud,Hi(1:NMAX,1:NMAX,1:Ncx_d,iyz),&
          sigma_lesser_ph_prev(nee,1:NMAX,1:NMAX,1:Ncx_d,iyz),sigma_greater_ph_prev(nee,1:NMAX,1:NMAX,1:Ncx_d,iyz),&
          g_lesser_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_greater_diag_local(1:NMAX,1:NMAX,1:Ncx_d),g_r_diag_local(1:NMAX,1:NMAX,1:Ncx_d),&
          tr,tre,cur(ee,nee,1:Ncx_d-1,iyz))
     
     do xx=1,ncx_d
        ldos(ee,nee,xx,iyz)=-traccia(dimag(g_r_diag_local(1:NM(xx),1:NM(xx),xx)))
        ndos(ee,nee,xx,iyz)=traccia(dimag(g_lesser_diag_local(1:NM(xx),1:NM(xx),xx)))
        pdos(ee,nee,xx,iyz)=-traccia(dimag(g_greater_diag_local(1:NM(xx),1:NM(xx),xx)))
        zdos(ee,nee,xx,iyz)= traccia(dimag( g_greater_diag_local(1:NM(xx),1:NM(xx),xx) - g_lesser_diag_local(1:NM(xx),1:NM(xx),xx) &
             -g_r_diag_local(1:NM(xx),1:NM(xx),xx)+transpose(dconjg(g_r_diag_local(1:NM(xx),1:NM(xx),xx)))))
     end do
     if(phonons)then  
        do xx=1,ncx_d
           call zgemm('n','n',nm(xx),nm(xx),nm(xx),alpha,sigma_lesser_ph_prev(nee,1:NM(xx),1:NM(xx),xx,iyz),nm(xx),g_greater_diag_local(1:NM(xx),1:NM(xx),xx),nm(xx),beta,A(1:nm(xx),1:nm(xx)),nm(xx))
           R_in (ee,nee,xx,iyz)=traccia(dble(A(1:nm(xx),1:nm(xx))))
           call zgemm('n','n',nm(xx),nm(xx),nm(xx),alpha,sigma_greater_ph_prev(nee,1:NM(xx),1:NM(xx),xx,iyz),nm(xx),g_lesser_diag_local(1:NM(xx),1:NM(xx),xx),nm(xx),beta,A(1:nm(xx),1:nm(xx)),nm(xx))
           R_out(ee,nee,xx,iyz)=traccia(dble(A(1:nm(xx),1:nm(xx))))        
        end do
     end if
     trans(ee,nee,1,iyz)=tr
     trans(ee,nee,2,iyz)=tre
     trans(ee,nee,3,iyz)=tr /(1.0_dp/(1.0_dp+exp((EN-mus)/(BOLTZ*TEMP)))-1.0_dp/(1.0_dp+exp((EN-mud)/(BOLTZ*TEMP))))
     trans(ee,nee,4,iyz)=tre/(1.0_dp/(1.0_dp+exp((EN-mud)/(BOLTZ*TEMP)))-1.0_dp/(1.0_dp+exp((EN-mus)/(BOLTZ*TEMP))))
     SELECT CASE (chtype)
     CASE ('t') 
        do xx=1,ncx_d    
           if(EN .ge. neutr(xx,iyz))then 

              dosn(1:NM(xx),1:NM(xx),xx,iyz)=dosn(1:NM(xx),1:NM(xx),xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_lesser_diag_local(1:NM(xx),1:NM(xx),xx))/(2.0_dp*pi)

           else if(EN .lt. neutr(xx,iyz))then 

              dosp(1:NM(xx),1:NM(xx),xx,iyz)=dosp(1:NM(xx),1:NM(xx),xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_greater_diag_local(1:NM(xx),1:NM(xx),xx))/(2.0_dp*pi)

           end if
        enddo
     CASE ('n')           
        do xx=1,ncx_d
           if(EN .ge. neutr(xx,iyz))then 

              dosn(1:NM(xx),1:NM(xx),xx,iyz)=dosn(1:NM(xx),1:NM(xx),xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_lesser_diag_local(1:NM(xx),1:NM(xx),xx))/(2.0_dp*pi)

           end if
        enddo
     CASE ('p')
        do xx=1,ncx_d
           if(EN .lt. neutr(xx,iyz))then

              dosp(1:NM(xx),1:NM(xx),xx,iyz)=dosp(1:NM(xx),1:NM(xx),xx,iyz)+&
                   degeneracy(iyz)*w(ee)*(g_greater_diag_local(1:NM(xx),1:NM(xx),xx))/(2.0_dp*pi)

           end if
        end do
     END SELECT
     
     con(iyz) = con(iyz)  + degeneracy(iyz)*w(ee)*(tr)
  
     cone(iyz)= cone(iyz) + degeneracy(iyz)*w(ee)*(tre)
     
  enddo  ! Nop
  !$omp end do nowait
     
  DEALLOCATE(g_lesser_diag_local)
  DEALLOCATE(g_greater_diag_local)
  DEALLOCATE(g_r_diag_local)
  DEALLOCATE(A)

  !$omp end parallel

  do nee= 2,Nop-1
     if(flag(nee,iyz)/=0)then
        ldos(ee,nee,1:ncx_d,iyz)=( ldos(ee,nee-1,1:ncx_d,iyz) + ldos(ee,nee+1,1:ncx_d,iyz) )/2.0_dp
        ndos(ee,nee,1:ncx_d,iyz)=( ndos(ee,nee-1,1:ncx_d,iyz) + ndos(ee,nee+1,1:ncx_d,iyz) )/2.0_dp
        pdos(ee,nee,1:ncx_d,iyz)=( pdos(ee,nee-1,1:ncx_d,iyz) + pdos(ee,nee+1,1:ncx_d,iyz) )/2.0_dp
        cur(ee,nee,1:Ncx_d-1,iyz)=( cur(ee,nee-1,1:Ncx_d-1,iyz) + cur(ee,nee+1,1:Ncx_d-1,iyz) )/2.0_dp
        trans(ee,nee,1:4,iyz)=( trans(ee,nee-1,1:4,iyz) + trans(ee,nee+1,1:4,iyz) )/2.0_dp
     end if
  end do
  
  write(*,*)ee,con(iyz),cone(iyz)
  end if
  end do !end loop kyz
     
enddo     ! Nsub


  deallocate(sigma_greater_ph_prev)
  deallocate(sigma_lesser_ph_prev)
  deallocate(sigma_greater_ph)
  deallocate(sigma_lesser_ph)
  deallocate(SCBA_error,SCBA_x)
  deallocate(flag)

write(*,*)
write(*,*)'Writing output files ...'
write(*,*)

do iyz=1,NKYZ  
if(k_selec(iyz))then
  ! ==========  End of parallel resolution ===========================
if(iyz <= 10)then

if(phonons)then
   OPEN(UNIT=10,FILE=TRIM(outdir)//'scat_pdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=20,FILE=TRIM(outdir)//'scat_ndens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=30,FILE=TRIM(outdir)//'scat_LDOS_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=60,FILE=TRIM(outdir)//'scat_ZDOS_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE=TRIM(outdir)//'scat_Jdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE=TRIM(outdir)//'scat_Jx_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=50,FILE=TRIM(outdir)//'scat_Jspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=51,FILE=TRIM(outdir)//'scat_Tspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=61,FILE=TRIM(outdir)//'scat_subv_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=62,FILE=TRIM(outdir)//'scat_subc_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
else
   OPEN(UNIT=10,FILE=TRIM(outdir)//'bal_pdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=20,FILE=TRIM(outdir)//'bal_ndens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=30,FILE=TRIM(outdir)//'bal_LDOS_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=60,FILE=TRIM(outdir)//'bal_ZDOS_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE=TRIM(outdir)//'bal_Jdens_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE=TRIM(outdir)//'bal_Jx_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=50,FILE=TRIM(outdir)//'bal_Jspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=51,FILE=TRIM(outdir)//'bal_Tspectrum_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=61,FILE=TRIM(outdir)//'bal_subv_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=62,FILE=TRIM(outdir)//'bal_subc_convergence_kyz_'//TRIM(STRINGA(iyz))//'_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
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
         write(40,*)(xx-1)*ac1*1.0d7,EN,degeneracy(iyz)*cur(ee,nee,xx,iyz)
      end do
      write(10,*)
      write(20,*)
      write(30,*)
      write(60,*)
      write(40,*)
      write(50,*)EN,degeneracy(iyz)*trans(ee,nee,1,iyz),degeneracy(iyz)*trans(ee,nee,2,iyz)
      write(51,*)EN,degeneracy(iyz)*trans(ee,nee,3,iyz),degeneracy(iyz)*trans(ee,nee,4,iyz)
   end do
end do
do xx=1,ncx_d-1
   ttr=0.0_dp
   do ee=1,Nsub
      do nee=1,Nop
         ttr=ttr+degeneracy(iyz)*w(ee)*cur(ee,nee,xx,iyz)
      end do
   end do
   write(41,*)(xx-1)*ac1*1.0d7,ttr/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
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
   
   Write(*,*) 'Transforming the carrier density',', ikyz =',iyz

   t1=SECNDS(0.0)
   
   call omp_set_num_threads(Nomp)
  !$omp parallel default(none) private(xx,n,i,j,ix,iy,iz,dens_yz,dens_z) &
  !$omp shared(chtype,iyz,imat,NCX_D,Nrx,Ndeltax,Ndeltay,Ndeltaz,NM,Ny,Nz,ac,Ry,Rz,dosn,dosp, &
  !$omp U_PSI,DX,DY,DZ,charge_n,charge_p,Ney,Nez,NTOT_Y,NTOT_Z,to2_lft,to2_bot,comp_ch)

   allocate(dens_yz(Ndeltay+1,Ndeltaz+1))

  !$omp do
   do xx=1,NCX_D
   if(comp_ch(xx)==0)then
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
   end if
   end do
   !$omp end do nowait
   
   deallocate(dens_yz)
   
   !$omp end parallel
   
end if

end if
end do ! end do iyz

if(phonons)then
   OPEN(UNIT=30,FILE=TRIM(outdir)//'scat_LDOS_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE=TRIM(outdir)//'scat_Jdens_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE=TRIM(outdir)//'scat_Jx_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=42,FILE=TRIM(outdir)//'scat_Junbalance_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=43,FILE=TRIM(outdir)//'scat_Qdens_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=44,FILE=TRIM(outdir)//'scat_Rin_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=45,FILE=TRIM(outdir)//'scat_Rout_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
else
   OPEN(UNIT=30,FILE=TRIM(outdir)//'bal_LDOS_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=40,FILE=TRIM(outdir)//'bal_Jdens_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
   OPEN(UNIT=41,FILE=TRIM(outdir)//'bal_Jx_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
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
            if(k_selec(iyz))     sumt=sumt+cur(ee,nee,xx,iyz)*degeneracy(iyz)
         end do
         write(40,*)(xx-1)*ac1*1.0d7,EN,sumt/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
      end do
      
      do  xx=1,ncx_d-2
         sumt=0.0_dp
         do iyz=1,nkyz
            if(k_selec(iyz)) sumt=sumt+(cur(ee,nee,xx+1,iyz)-cur(ee,nee,xx,iyz))*degeneracy(iyz)
         end do
         write(42,*)(xx)*ac1*1.0d7,EN,sumt/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
      end do
      
      do  xx=1,ncx_d-2
         sumt=0.0_dp
         do iyz=1,nkyz
            if(k_selec(iyz)) sumt=sumt+(cur(ee,nee,xx+1,iyz)-cur(ee,nee,xx,iyz))*degeneracy(iyz)
         end do
         write(43,*)(xx)*ac1*1.0d7,EN,-EN*sumt/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
      end do
      
      write(30,*)
      write(40,*)
      write(42,*)
      write(43,*)
   end do
end do

if(phonons)then
   do nee=1,Nop
      do ee=1,Nsub
         emin_local=en_global(ee) !local minimum
         EN=emin+emin_local+dble(nee-1)*Eop
         do xx=1,ncx_d
            sumt=0.0_dp
            do iyz=1,nkyz
               if(k_selec(iyz))   sumt=sumt+R_in(ee,nee,xx,iyz)*degeneracy(iyz)
            end do
            write(44,*)(xx-1)*ac1*1.0d7,EN,sumt
         end do
         do xx=1,ncx_d
            sumt=0.0_dp
            do iyz=1,nkyz
               if(k_selec(iyz))   sumt=sumt+R_out(ee,nee,xx,iyz)*degeneracy(iyz)
            end do
            write(45,*)(xx-1)*ac1*1.0d7,EN,sumt
         end do
         
      write(44,*)
      write(45,*)
   end do
end do
close(44)
close(45)
end if
      
do xx=1,ncx_d-1
   ttr=0.0_dp
   do iyz=1,NKYZ
      if(k_selec(iyz))then
         do ee=1,Nsub
            do nee=1,Nop
               ttr=ttr+degeneracy(iyz)*w(ee)*cur(ee,nee,xx,iyz)
            end do
         end do
      end if
   end do
   write(41,*)(xx-1)*ac1*1.0d7,ttr/dble(NCY*NCZ)/(hbar*2.0_dp*pi)*ELCH**2
end do
close(30)
close(40)
close(41)
close(42)
close(43)

deallocate(trans)
deallocate(ldos)
deallocate(zdos)
deallocate(ndos)
deallocate(pdos)
deallocate(cur)
deallocate(dosn)
deallocate(dosp)
if(phonons)deallocate(R_in,R_out)

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
OPEN(UNIT=11,FILE=TRIM(outdir)//'charge_n2D_convergence_vd_'//TRIM(STRINGA(ss))//'_vg_'//TRIM(STRINGA(gg))//'.dat',STATUS='replace')
if(chtype == 'n')then
do xx=1,Ncx_d
   do ix=1,Ndeltax
      n_2D=0.0_dp
      do iz=1,NTOT_Z
         do iy=NTOT_Y/2+1,NTOT_Y/2+1
            write(10,*)ix+(xx-1)*Ndeltax,iz,charge_n(ix+(xx-1)*Ndeltax,iy,iz)
         end do
         do iy=1,NTOT_Y
            n_2D=n_2D+charge_n(ix+(xx-1)*Ndeltax,iy,iz)*deltaz/NTOT_Y
         end do
      end do
      write(11,*)ix+(xx-1)*Ndeltax,n_2D
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
close(11)

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
  
  write(*,*)
  write(*,*)'IDScurrent =',IDScurrent,ISDcurrent
  write(*,*)
  write(*,*)'Bal_IDScurrent=', IDScurrentb
  
  Gcon=sum(con(1:NKYZ))/abs(mud-mus)

  write(*,*)'Gcon=',Gcon,sum(cone(1:nkyz))/abs(mud-mus)
  write(*,*)

  deallocate(con,cone,conb)

!!!!stop
END SUBROUTINE negf_mixed


subroutine RGF(it,nmax,ff,E,mul,mur,Hii,sigma_lesser_ph,sigma_greater_ph,ndens,pdens,ldos,tr,tre,cur) 

  implicit none
  
  INTEGER     :: it,i,j,l,nmax,ff
  REAL(dp)    :: E,mul,mur,tr,tre,tmp
  REAL(dp)    :: cur(Ncx_d-1)
  COMPLEX(dp) :: Hii(nmax,nmax,ncx_d),sigma_lesser_ph(nmax,nmax,Ncx_d),sigma_greater_ph(nmax,nmax,Ncx_d)
  COMPLEX(dp), allocatable :: sig(:,:),sigmal(:,:),sigmar(:,:),sigma_r_ph(:,:)
  COMPLEX(dp), allocatable :: Gn(:,:),Gp(:,:),G00(:,:),GN0(:,:)
  COMPLEX(dp), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:)
  COMPLEX(dp) :: ldos(nmax,nmax,ncx_d),ndens(nmax,nmax,ncx_d),pdens(nmax,nmax,ncx_d)
  COMPLEX(dp), allocatable :: H00(:,:),H10(:,:)
  COMPLEX(dp), allocatable :: A(:,:),B(:,:),C(:,:),D(:,:),Id(:,:)
  COMPLEX(dp) :: z
  real(qp)    :: zeta
  logical     :: flag_nan

!!!  seuil=1.0d0
  
if(ff == 0)then
   allocate( Gl(nmax,nmax,Ncx_d), Gln(nmax,nmax,Ncx_d), Glp(nmax,nmax,Ncx_d) )
   allocate( Gn(nmax,nmax), Gp(nmax,nmax), G00(nmax,nmax), GN0(nmax,nmax) )
   allocate( sig(nmax,nmax), sigmal(nmax,nmax), sigmar(nmax,nmax), sigma_r_ph(nmax,nmax) )
   allocate( H00(nmax,nmax), H10(nmax,nmax), A(nmax,nmax), B(nmax,nmax), C(nmax,nmax), Id(nmax,nmax) )
   
  z=cmplx(E,1.0d-6,kind=dp)

  Gln=0.0_dp
  Glp=0.0_dp
  Gl=0.0_dp
  ldos=0.0_dp
  ndens=0.0_dp
  pdens=0.0_dp
  cur=0.0_dp

  do l=1,ncx_d
     sigma_lesser_ph(1:nm(l),1:nm(l),l)=(sigma_lesser_ph(1:nm(l),1:nm(l),l) - &
          transpose(dconjg(sigma_lesser_ph(1:nm(l),1:nm(l),l))))/2.0_dp
     sigma_greater_ph(1:nm(l),1:nm(l),l)=(sigma_greater_ph(1:nm(l),1:nm(l),l) - &
          transpose(dconjg(sigma_greater_ph(1:nm(l),1:nm(l),l))))/2.0_dp
  end do
! self energy on the left contact
  l=1
  Id(1:nm(l),1:nm(l))=Si(iyz,imat(l))%H(1:NM(l),1:NM(l))
  
  sigma_r_ph(1:nm(l),1:nm(l))=(sigma_greater_ph(1:nm(l),1:nm(l),l) - sigma_lesser_ph(1:nm(l),1:nm(l),l))/2.0_dp

  H00(1:nm(l),1:nm(l))=Hii(1:nm(l),1:nm(l),l)+sigma_r_ph(1:nm(l),1:nm(l))
  
  H10(1:nm(l),1:nm(l)) = TL(iyz,ihet(l))%H(1:NM(l),1:NM(l)) !!!! H_{1,0}
  !write(*,*)'ls',E
  call sancho(nm(l),E,H00(1:nm(l),1:nm(l)),transpose(dconjg(H10(1:nm(l),1:nm(l)))),id(1:nm(l),1:nm(l)),G00(1:nm(l),1:nm(l))) 
  !write(*,*)'ls',E,-traccia(dimag(G00(1:nm(l),1:nm(l))))
  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           flag_nan=.true.
           exit
        end if
     end do
  end do
  
  if(traccia(dimag(G00(1:nm(l),1:nm(l)))) <= 0.0_dp)flag_nan=.true.
  if(flag_nan)call oldsancho(nm(l),E,H00(1:nm(l),1:nm(l)),transpose(dconjg(H10(1:nm(l),1:nm(l)))),id(1:nm(l),1:nm(l)),G00(1:nm(l),1:nm(l)))
  flag_nan=.false.
  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           ff=l
           write(*,*)'NaN warning! Pb w lft snch. Eliminating E =',E
           exit
        end if
     end do
  end do
  
  !  if( traccia(dimag(- G00(1:nm(l),1:nm(l)))) <= 0.0d0 )then
!     G00=0.0_dp
!     ff=l
!     !write(*,*)'pb w snch l, E=',E,traccia(dimag(- G00(1:NM(l),1:NM(l))))
!  end if
     
  H10(1:nm(1),1:nm(1))=TL(iyz,ihet(l))%H(1:NM(1),1:NM(1))  !!!! H_{1,0}
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,H10(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,A(1:nm(l),1:nm(l)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),H10(1:nm(l),1:nm(l)),nm(l),beta,sigmal(1:nm(l),1:nm(l)),nm(l))  
  A(1:nm(l),1:nm(l))=z*Id(1:nm(l),1:nm(l))-H00(1:nm(l),1:nm(l))-sigmal(1:nm(l),1:nm(l))
  call invert(A(1:nm(l),1:nm(l)),nm(l),nm(l))
  Gl(1:nm(l),1:nm(l),l)=A(1:nm(l),1:nm(l))

  zeta=(E-mul)/(BOLTZ*TEMP)
  sig(1:nm(l),1:nm(l))=-(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm(zeta) 
  sig(1:nm(l),1:nm(l))=sig(1:nm(l),1:nm(l))+sigma_lesser_ph(1:nm(l),1:nm(l),l)
  
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  !C(1:nm(l),1:nm(l))=(C(1:nm(l),1:nm(l))-transpose(dconjg(C(1:nm(l),1:nm(l)))))/2.0_dp
  Gln(1:nm(l),1:nm(l),l)=C(1:nm(l),1:nm(l))

  sig(1:nm(l),1:nm(l))=(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm(-zeta)
  sig(1:nm(l),1:nm(l))=sig(1:nm(l),1:nm(l))+sigma_greater_ph(1:nm(l),1:nm(l),l)

  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  !C(1:nm(l),1:nm(l))=(C(1:nm(l),1:nm(l))-transpose(dconjg(C(1:nm(l),1:nm(l)))))/2.0_dp
  Glp(1:nm(l),1:nm(l),l)=C(1:nm(l),1:nm(l))

  !Gl(1:nm(l),1:nm(l),l)=(Gl(1:nm(l),1:nm(l),l)+transpose(conjg(Gl(1:nm(l),1:nm(l),l))))/2.0_dp + &
  !     (Glp(1:nm(l),1:nm(l),l)-Gln(1:nm(l),1:nm(l),l))/2.0_dp

  
  Do l=2,Ncx_d-1
 
     Id(1:nm(l),1:nm(l))=Si(iyz,imat(l))%H(1:NM(l),1:NM(l))
     sigma_r_ph(1:nm(l),1:nm(l))=(sigma_greater_ph(1:nm(l),1:nm(l),l) - sigma_lesser_ph(1:nm(l),1:nm(l),l))/2.0_dp
     
     H00(1:nm(l),1:nm(l))=Hii(1:nm(l),1:nm(l),l)+sigma_r_ph(1:nm(l),1:nm(l))
            
     H10(1:nm(l),1:nm(l-1))=TL(iyz,ihet(l))%H(1:NM(l),1:NM(l-1)) !!!! H_{l,l-1}
   
     call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Gl(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
     call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
     A(1:nm(l),1:nm(l))=z*Id(1:nm(l),1:nm(l))-H00(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))              
     call invert(A(1:nm(l),1:nm(l)),nm(l),nm(l))
     Gl(1:nm(l),1:nm(l),l)=A(1:nm(l),1:nm(l))     

     sig(1:nm(l-1),1:nm(l-1))=Gln(1:nm(l-1),1:nm(l-1),l-1)
     call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),sig(1:nm(l-1),1:nm(l-1)),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
     call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
     C(1:nm(l),1:nm(l))=C(1:nm(l),1:nm(l))+sigma_lesser_ph(1:nm(l),1:nm(l),l)

     call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),C(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l)) 
     call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,Gn(1:nm(l),1:nm(l)),nm(l))
     !Gn(1:nm(l),1:nm(l))=(Gn(1:nm(l),1:nm(l))-transpose(dconjg(Gn(1:nm(l),1:nm(l)))))/2.0_dp
     Gln(1:nm(l),1:nm(l),l)=Gn(1:nm(l),1:nm(l))

     sig(1:nm(l-1),1:nm(l-1))=Glp(1:nm(l-1),1:nm(l-1),l-1)
     call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),sig(1:nm(l-1),1:nm(l-1)),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l))
     call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))

     C(1:nm(l),1:nm(l))=C(1:nm(l),1:nm(l))+sigma_greater_ph(1:nm(l),1:nm(l),l)

     call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),C(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
     call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),A(1:nm(l),1:nm(l)),nm(l),beta,Gp(1:nm(l),1:nm(l)),nm(l))
     !Gp(1:nm(l),1:nm(l))=(Gp(1:nm(l),1:nm(l))-transpose(dconjg(Gp(1:nm(l),1:nm(l)))))/2.0_dp
     Glp(1:nm(l),1:nm(l),l)=Gp(1:nm(l),1:nm(l))

     !Gl(1:nm(l),1:nm(l),l)=(Gl(1:nm(l),1:nm(l),l)+transpose(conjg(Gl(1:nm(l),1:nm(l),l))))/2.0_dp + &
     !     (Glp(1:nm(l),1:nm(l),l)-Gln(1:nm(l),1:nm(l),l))/2.0_dp
     
  enddo

  
! self energy on the right contact
  l=Ncx_d

  Id(1:nm(l),1:nm(l))=Si(iyz,imat(l))%H(1:NM(l),1:NM(l))
 
  sigma_r_ph(1:nm(l),1:nm(l))=(sigma_greater_ph(1:nm(l),1:nm(l),l) - sigma_lesser_ph(1:nm(l),1:nm(l),l))/2.0_dp
  !if( traccia(dimag(- sigma_r_ph(1:nm(l),1:nm(l)))) < 0.0d0 ) sigma_r_ph=0.0_dp

  H00(1:nm(l),1:nm(l))=Hii(1:nm(l),1:nm(l),l)+sigma_r_ph(1:nm(l),1:nm(l))
  
  H10(1:nm(l),1:nm(l))=TL(iyz,ihet(l+1))%H(1:NM(l),1:NM(l))   !!! H(N+1,N)
 ! write(*,*)'rs',E
  call sancho(nm(l),E,H00(1:nm(l),1:nm(l)),H10(1:nm(l),1:nm(l)),id(1:nm(l),1:nm(l)),G00(1:nm(l),1:nm(l)))
 ! write(*,*)'rs',E,-traccia(dimag(G00(1:nm(l),1:nm(l))))
  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           !ff=l
           !write(*,*)'NaN warning! Pb w lft snch. Eliminating E =',E
           !write(*,*)'NaN warning! Pb w lft snch. E =',E,'trying again'
           flag_nan=.true.
           exit
        end if
     end do
  end do
  if(traccia(dimag(G00(1:nm(l),1:nm(l)))) <= 0.0_dp)flag_nan=.true.
  if(flag_nan)call oldsancho(nm(l),E,H00(1:nm(l),1:nm(l)),H10(1:nm(l),1:nm(l)),id(1:nm(l),1:nm(l)),G00(1:nm(l),1:nm(l)))
  flag_nan=.false.
  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           ff=l
           write(*,*)'NaN warning! Pb w rgt snch. Eliminating E =',E
           exit
        end if
     end do
  end do
  
  !if( traccia(dimag(- G00(1:nm(l),1:nm(l)))) <= 0.0d0 )then 
  !   G00=0.0_dp
  !   ff=l
  !   !write(*,*)'pb w snch r, E=',E,traccia(dimag(- G00(1:nm(l),1:nm(l))))
  !end if

  call zgemm('c','n',nm(l),nm(l),nm(l),alpha,H10(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,A(1:nm(l),1:nm(l)),nm(l)) 
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),H10(1:nm(l),1:nm(l)),nm(l),beta,sigmar(1:nm(l),1:nm(l)),nm(l)) 

  H10(1:nm(l),1:nm(l-1))=TL(iyz,ihet(l))%H(1:NM(l),1:NM(l-1)) !!! H(N,N-1)
  call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Gl(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))

  G00(1:nm(l),1:nm(l))=z*Id(1:nm(l),1:nm(l))-H00(1:nm(l),1:nm(l))-sigmar(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))   
  call invert(G00(1:nm(l),1:nm(l)),nm(l),nm(l))

  call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Gln(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))

  zeta=(E-mur)/(BOLTZ*TEMP)
  sig(1:nm(l),1:nm(l))=-(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm(zeta)+C(1:nm(l),1:nm(l))

  sig(1:nm(l),1:nm(l))=sig(1:nm(l),1:nm(l))+sigma_lesser_ph(1:nm(l),1:nm(l),l)

  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,G00(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l)) 
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,Gn(1:nm(l),1:nm(l)),nm(l)) 

  !Gn(1:nm(l),1:nm(l))=(Gn(1:nm(l),1:nm(l))-transpose(dconjg(Gn(1:nm(l),1:nm(l)))))/2.0_dp
  ndens(1:nm(l),1:nm(l),l)=Gn(1:nm(l),1:nm(l))

  call zgemm('n','n',nm(l),nm(l-1),nm(l-1),alpha,H10(1:nm(l),1:nm(l-1)),nm(l),Glp(1:nm(l-1),1:nm(l-1),l-1),nm(l-1),beta,B(1:nm(l),1:nm(l-1)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l-1),alpha,B(1:nm(l),1:nm(l-1)),nm(l),H10(1:nm(l),1:nm(l-1)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  
  sig(1:nm(l),1:nm(l))=(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm(-zeta)+C(1:nm(l),1:nm(l))
  sig(1:nm(l),1:nm(l))=sig(1:nm(l),1:nm(l))+sigma_greater_ph(1:nm(l),1:nm(l),l)

  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,G00(1:nm(l),1:nm(l)),nm(l),sig(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  call zgemm('n','c',nm(l),nm(l),nm(l),alpha,B(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,Gp(1:nm(l),1:nm(l)),nm(l))

  !Gp(1:nm(l),1:nm(l))=(Gp(1:nm(l),1:nm(l))-transpose(dconjg(Gp(1:nm(l),1:nm(l)))))/2.0_dp
  pdens(1:nm(l),1:nm(l),l)=Gp(1:nm(l),1:nm(l))

  ldos(1:nm(l),1:nm(l),l) = G00(1:nm(l),1:nm(l))
!  ldos(1:nm(l),1:nm(l),l)=(G00(1:nm(l),1:nm(l))+transpose(conjg(G00(1:nm(l),1:nm(l)))))/2.0_dp + &
!       (Gp(1:nm(l),1:nm(l))-Gn(1:nm(l),1:nm(l)))/2.0_dp
  
  A(1:nm(l),1:nm(l))=-(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm(zeta)
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gp(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  A(1:nm(l),1:nm(l))=(sigmar(1:nm(l),1:nm(l))-transpose(dconjg(sigmar(1:nm(l),1:nm(l)))))*ferm(-zeta)
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gn(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  tr=-traccia(dble(B(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))))

!-------------------------
  tmp=abs(traccia(dimag( pdens(1:nm(l),1:nm(l),l) - ndens(1:nm(l),1:nm(l),l)&
       - ldos(1:nm(l),1:nm(l),l) + transpose(dconjg(ldos(1:nm(l),1:nm(l),l))) )))
  if( tmp  > seuil )then
     !!!write(*,*)l,E,tmp
     ff=l  
  end if
  
  
  do l=Ncx_d-1,1,-1

     H10(1:nm(l+1),1:nm(l))=TL(iyz,ihet(l+1))%H(1:NM(l+1),1:NM(l)) !!! H(l+1,l)
     
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

     !Gn(1:nm(l),1:nm(l))=(Gn(1:nm(l),1:nm(l))-transpose(dconjg(Gn(1:nm(l),1:nm(l)))))/2.0_dp
     ndens(1:nm(l),1:nm(l),l)=Gn(1:nm(l),1:nm(l))
     !Gp(1:nm(l),1:nm(l))=(Gp(1:nm(l),1:nm(l))-transpose(dconjg(Gp(1:nm(l),1:nm(l)))))/2.0_dp
     pdens(1:nm(l),1:nm(l),l)=Gp(1:nm(l),1:nm(l))
     ldos(1:nm(l),1:nm(l),l) = G00(1:nm(l),1:nm(l))
     !ldos(1:nm(l),1:nm(l),l)=(G00(1:nm(l),1:nm(l))+transpose(conjg(G00(1:nm(l),1:nm(l)))))/2.0_dp + &
     !  (Gp(1:nm(l),1:nm(l))-Gn(1:nm(l),1:nm(l)))/2.0_dp

     tmp=abs(traccia(dimag( pdens(1:nm(l),1:nm(l),l) - ndens(1:nm(l),1:nm(l),l) &
          - ldos(1:nm(l),1:nm(l),l) + transpose(dconjg(ldos(1:nm(l),1:nm(l),l))) )))
     if( tmp  > seuil  )then
        !!!write(*,*)l,E,tmp
        ff=l  
     end if
     
  end do
  
  l=1
  zeta=(E-mul)/(BOLTZ*TEMP)
  A(1:nm(l),1:nm(l))=-(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm(zeta)
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gp(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  A(1:nm(l),1:nm(l))=(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm(-zeta)
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gn(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  tre=-traccia(dble(B(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))))
  
  if ( ff /= 0 ) then
     write(*,*)'ignoring E =',E
     tr=0.0_dp
     tre=0.0_dp
     ndens=0.0_dp
     pdens=0.0_dp
     ldos=0.0_dp
     cur=0.0_dp
  end if

  deallocate( Gl, Gln, Glp )
  deallocate( Gn, Gp, G00, GN0 )
  deallocate( sig, sigmal, sigmar, sigma_r_ph )
  deallocate( H00, H10, A , B, C, Id )
  
elseif ( ff /= 0 ) then
   !!!write(*,*)'ignoring E =',E
   tr=0.0_dp
   tre=0.0_dp
   ndens=0.0_dp
   pdens=0.0_dp
   ldos=0.0_dp
   cur=0.0_dp
end if

end subroutine RGF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sancho-Rubio
subroutine sancho(nm,E,H00,H10,S,G00)
  implicit none
  integer, intent(IN) :: nm
  integer :: i
  integer :: nmax = 100
  COMPLEX(DP) :: z
  real(dp), intent(IN) :: E
  REAL(DP) :: TOL=1.0D-15, error

  COMPLEX(DP), INTENT(IN) :: H00(nm,nm), H10(nm,nm), S(nm,nm)
  COMPLEX(DP), INTENT(OUT) :: G00(nm,nm)

  COMPLEX(DP), ALLOCATABLE :: A(:,:), B(:,:), C(:,:)
  COMPLEX(DP), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:)

  Allocate( H_BB(nm,nm) )
  Allocate( H_SS(nm,nm) )
  Allocate( H_01(nm,nm) )
  Allocate( H_10(nm,nm) )
  Allocate( A(nm,nm) )
  Allocate( B(nm,nm) )
  Allocate( C(nm,nm) )
  
  z = cmplx(E,eta,kind=dp)
  
  H_BB = H00  
  H_10 = H10
  H_01 = TRANSPOSE( DCONJG( H10 ) )
  H_SS = H_BB

  do i = 1,nmax
     
     A = z*S - H_BB
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

     error = abs( sum(sum(H_10, dim=2), dim=1) ) + abs(sum(sum(H_01, dim=2), dim=1) )
!     write(*,*)i,E,error
     IF ( abs(error) .lt. TOL ) THEN
        EXIT
     END IF
     
  enddo

  if (i == nmax)then
     write(*,*)'pb nmax',E,error
     G00=0.0_dp
  else
     G00 = E*S - H_SS
     call invert(G00,nm,nm)
  end if


  Deallocate( A )
  Deallocate( B )
  Deallocate( C )
  Deallocate( H_BB )
  Deallocate( H_SS )
  Deallocate( H_01 )
  Deallocate( H_10 )

end subroutine sancho


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Old Sancho-Rubio
subroutine oldsancho(nm,E,H00,H10,Id,G00)

  integer, intent(IN) :: nm
  integer i,nmax
  COMPLEX(DP) :: z
  real(dp), intent(IN) :: E
  REAL(DP) :: TOL=1.0D-12, error

  COMPLEX(DP), INTENT(IN) :: H00(nm,nm), H10(nm,nm), Id(nm,nm)
  COMPLEX(DP), INTENT(OUT) :: G00(nm,nm)

  COMPLEX(DP), ALLOCATABLE :: A(:,:), B(:,:), C(:,:)
  COMPLEX(DP), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:)

  Allocate( H_BB(nm,nm) )
  Allocate( H_SS(nm,nm) )
  Allocate( H_01(nm,nm) )
  Allocate( H_10(nm,nm) )
  Allocate( A(nm,nm) )
  Allocate( B(nm,nm) )
  Allocate( C(nm,nm) )

  nmax=100

  z = E+cmplx(0.0_dp,1.0d-3,kind=dp)


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
     error=traccia(abs(C(1:nm,1:nm)))
    ! write(*,*)i,error
     IF ( abs(error) .lt. TOL ) THEN
!       write(*,*) 'SR: Exited, abs(error)=',abs(error)
       EXIT
     END IF

     IF (i == nmax) THEN
        write(*,*) 'warning: nmax reached in sancho!!!',error,E
        write(*,*) 'try to reduce eta'
        EXIT
     END IF

  enddo

  if (i==nmax)then
     G00=0.0_dp
  else
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

end subroutine oldsancho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function traccia(A)
  real(dp) :: A(:,:)
  real(dp) :: traccia
  integer :: i,j

  traccia=dble(sum( a(:,:), &
       reshape( (/((i==j,i=1,size(a,1)),j=1,size(a,2))/), &
       (/size(a,1),size(a,2)/) ) ) ) 

end function traccia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function bose(a)
  real(dp) a,bose
  bose=1.0_dp/(EXP(a) - 1.0_dp)
end function bose

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Real(qp) Function ferm(a)
  real(qp) a
  ferm = 1.0_qp/(1.0_qp + Exp(a))
End Function ferm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine invert(A,n,np)

  implicit none
      
  integer :: info,n,np
      
  integer, dimension(np) :: ipiv
      
  complex(dp) :: A(n,n)    
  complex(dp) :: work(np*np)   

  call zgetrf(n,n,A,n,ipiv,info)
  call zgetri(n,A,n,ipiv,work,np*np,info)

end subroutine invert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer function ind_Kyz(v)
   implicit none

   INTEGER  :: iz,iy,l
   real(dp) :: v(2),vv(2)

   ind_kyz=0

   vv=v

   if(abs(V(1))>0.5_dp) vv(1)=abs(abs(v(1))-1.0_dp)
   if(abs(V(2))>0.5_dp) vv(2)=abs(abs(v(2))-1.0_dp)
   
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
      write(*,*)'pb w ind_kyz'
      stop
   end if
   
 end function ind_Kyz
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE negf


