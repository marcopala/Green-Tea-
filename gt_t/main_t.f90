! Copyright or © or Copr. Marco Pala (February 24, 2022)

! e-mail: marco.pala@c2n.upsaclay.fr ;   marco.pala@cnrs.fr

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

program main


use static
use pseudopot_so_gen
use indata

use poisson
use negf

implicit none 

INTEGER               :: gg, ss, bb, ii, jj, transport_iter, MAXITER
REAL(DP)              :: vg, bz, transport_error!, car, carp
LOGICAL               :: init, sweepbz

REAL(DP), ALLOCATABLE :: POT_3D(:)
REAL(DP), ALLOCATABLE :: EC_3D(:)
REAL(DP), ALLOCATABLE :: EV_3D(:)
REAL(DP), ALLOCATABLE :: outer_pot3D_old(:)
REAL(DP), ALLOCATABLE :: pot3D(:,:,:)
REAL(DP), ALLOCATABLE :: plotpot(:,:,:)
REAL(DP), ALLOCATABLE :: rho_3D_n(:)
REAL(DP), ALLOCATABLE :: rho_3D_p(:)
REAL(DP), ALLOCATABLE :: rho_n(:,:,:)
REAL(DP), ALLOCATABLE :: rho_p(:,:,:)
REAL(DP), ALLOCATABLE :: drho_3D_n(:)
REAL(DP), ALLOCATABLE :: drho_3D_p(:)
REAL(DP), ALLOCATABLE :: Fn(:)
REAL(DP), ALLOCATABLE :: Fp(:)
REAL(DP), ALLOCATABLE :: Hz(:) ! vertical magnetic field
REAL(DP)              :: conductance, conductanceb, ISDcurrent, IDScurrent, IDScurrentb

real(4) :: t1,t2

gg=0

conductance=0.0d0
conductanceb=0.0d0
IDScurrent=0.0d0
ISDcurrent=0.0d0
IDScurrentb=0.0d0

call indata_readinput()
call indata_grid()

call indata_structure_init()
call indata_build_doping(dop_vec,coul,source_dop_val,drain_dop_val,channel_dop_val,&
     NUMN_3D,LWORK_3D) 

call read_QE_output()
!stop


IF (VGMIN.eq.VGMAX) THEN
   NUMVG=0
ELSE
   NUMVG=ceiling((VGMAX-VGMIN)/DELTAVG)  
END IF 
IF (VDMIN.eq.VDMAX) THEN
   NUMVD=0
ELSE
   NUMVD=ceiling((VDMAX-VDMIN)/DELTAVD)  
END IF
IF (BZMIN.eq.BZMAX) THEN
   NUMBZ=0
ELSE
   NUMBZ=ceiling((BZMAX-BZMIN)/DELTABZ)  
END IF


!DO bb=0, NUMBZ
  DO ss=0, NUMVD
     open(23,file=TRIM(outdir)//'bal_current_vg_'//TRIM(STRINGA(ss))//'.dat',status='replace')
     close(23)
  END DO
  DO gg=0, NUMVG
     open(23,file=TRIM(outdir)//'bal_current_vd_'//TRIM(STRINGA(gg))//'.dat',status='replace')
     close(23)
  END DO
!END DO

  
if(phonons)then
!DO bb=0, NUMBZ
   DO ss=0, NUMVD
      open(23,file=TRIM(outdir)//'scat_current_vg_'//TRIM(STRINGA(ss))//'.dat',status='replace')
      close(23)
   END DO
   DO gg=0, NUMVG
      open(23,file=TRIM(outdir)//'scat_current_vd_'//TRIM(STRINGA(gg))//'.dat',status='replace')
      close(23)
   END DO
!END DO
end if

open(23, file = TRIM(outdir)//'convergence_LOG.dat',status='replace')
close(23)


   ALLOCATE(POT_3D(0:LWORK_3D-1))
   ALLOCATE(outer_pot3D_old(0:LWORK_3D-1))
   ALLOCATE(EC_3D(0:LWORK_3D-1))
   ALLOCATE(EV_3D(0:LWORK_3D-1))
   ALLOCATE(pot3D(1:NTOT_X,1:NTOT_Y,1:NTOT_Z))
   ALLOCATE(plotpot(1:NTOT_X,1:NTOT_Y,1:NTOT_Z))
   ALLOCATE(rho_3D_n(0:LWORK_3D-1))
   ALLOCATE(rho_3D_p(0:LWORK_3D-1))
   ALLOCATE(rho_n(1:NTOT_X,1:NTOT_Y,1:NTOT_Z))
   ALLOCATE(rho_p(1:NTOT_X,1:NTOT_Y,1:NTOT_Z))
   ALLOCATE(drho_3D_n(0:LWORK_3D-1))
   ALLOCATE(drho_3D_p(0:LWORK_3D-1))
   ALLOCATE(Fn(0:LWORK_3D-1))
   ALLOCATE(Fp(0:LWORK_3D-1))
   ALLOCATE(Hz(1:NTOT_X))

   Fn=3.0d0
   Fp=-3.0d0
   Hz=0.0d0
   rho_n=0.0d0
   rho_p=0.0d0
   mus=0.0d0
   mud=0.0d0
   bz=0.0d0

   outer_pot3D_old=0.0d0
   EC_3D=0.0d0
   EV_3D=0.0d0
   rho_3D_n=0.0d0
   rho_n=0.0d0
   rho_3D_p=0.0d0
   rho_p=0.0d0

!  write(*,*)'potelectr=',potelectr

DO gg=0, NUMVG
   DO ss=0, NUMVD
      sweepbz=.TRUE.
      DO bb=0, NUMBZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Stimulus definition!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         mus=0.0d0
         mud=-(VDMIN+(1.0d0*ss)*DELTAVD)
         vg=VGMIN+(1.0d0*gg)*DELTAVG
         bz=BZMIN+(1.0d0*bb)*DELTABZ
         potelectr=vg
         potelectr11=-mus
         potelectr12=-mud
        

         DO ii=1, NTOT_X
            if( ii <= NTOT_X ) Hz(ii) = bz - bz*(ii-NTOT_X+drain_len)/dble(drain_len) ! drain
            if(ii<=source_len+gate_len+2*spacer) Hz(ii) = bz ! channel
            if(ii<=source_len) Hz(ii) = bz*(ii-1)/dble(source_len) ! source
         END DO

         if( bz /= 0.0_dp )then
            write(*,*)'BE CAREFUL: the magnetic field is oriented along the z-axis and works only for non-periodic systems!'
            write(*,*) 'MAGNETIC FIELD NOT YET implemented'
            stop
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)'///////////////////////////////////////////////////////////////////////'
  write(*,*)'//VG=',vg,'VD=',-mud,'BZ=',bz,'//'
  write(*,*)'///////////////////////////////////////////////////////////////////////'


  open(23, file = TRIM(outdir)//'convergence_LOG.dat',status='old',position='append')
  write(23,*) ':::::::::::::::::::::Bias definition::::::::::::::::::::::::'
  write(23,*) '::VG=',vg,'VD=',-mud,'BZ=',bz,'::'
  write(23,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
  close(23)


  transport_error=1.0d0
  transport_iter=0
  init=.TRUE.
  maxiter=max_iter_outer

  IF(sweepbz)THEN
     if(ss.eq.0 .and. gg.eq.0 .and. bb.eq.0)then
        write(*,*)'////////////////////////////////////////////////////////////////'
        if(.not. in_pot)then
           write(*,*)'Poisson init guess'
           call pot_init(POT_3D,whichkind_3D,map_3D,ntot_x,ntot_y,ntot_z,lwork_3D)
           write(*,*)'POT INIT DONE',mus,mud
        else
           open(10,file='Init_Potential.dat',status='unknown')
           do ii=0,LWORK_3D-1
              read(10,*)POT_3D(ii)
           end do
           close(10)
        end if
     end if
  END IF


t1=SECNDS(0.0)

DO WHILE ((transport_error.ge.ERROR_OUTER).and.(transport_iter.le.MAXITER))
  transport_iter=transport_iter+1
  write(*,*)'ITERATION n.', transport_iter


  EC_3D=0.0d0
  EV_3D=0.0d0
  rho_3D_n=0.0d0
  rho_n=0.0d0
  drho_3D_n=0.0d0
  rho_3D_p=0.0d0
  rho_p=0.0d0
  drho_3D_p=0.0d0

  CALL shift_potential(EC_3D,-POT_3D, 3.0d0,E_GAP,whichkind_3D,map_3D,list_3D_ord,&
       whichkind_3D_ord,epsilon_3D,NTOT_X,NTOT_Y,NTOT_Z,NUMEL_3D,LWORK_3D)
  CALL shift_potential(EV_3D,-POT_3D,-3.0d0,0.0d0,whichkind_3D,map_3D,list_3D_ord,&
       whichkind_3D_ord,epsilon_3D,NTOT_X,NTOT_Y,NTOT_Z,NUMEL_3D,LWORK_3D)

  write(*,*)'e_gap',E_GAP
  write(*,*)'diel_sc',diel_sc

  CALL map_potential(pot3D,-POT_3D,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

  write(*,*)
  write(*,*)' negf starts'
  call negf_mixed(POT3D,rho_n,rho_p,ISDcurrent,IDScurrent,IDScurrentb,ss,gg,transport_iter)
  write(*,*)

  !if(onlyT)then
  !   write(*,*)
  !   write(*,*)
  !   write(*,*) ' Simulation ended. &
  !        Only transmission and current in the flat-band configuration were computed. '
  !   write(*,*)
  !   write(*,*)
  !   write(*,*) ' This is an open source software distributed under the CeCILL license.'
  !   write(*,*) ' Please cite  '
  !   write(*,*) ' M. G. Pala, P. Giannozzi, and D. Esseni, Phys. Rev. B 102, 045410 (2020)'
  !   write(*,*) ' DOI: https://doi.org/10.1103/PhysRevB.102.045410'
  !   write(*,*)
  !   write(*,*)
  !   stop
  !end if

  
  OPEN(UNIT=22,FILE='charge_y'//TRIM(STRINGA(transport_iter))//'.dat',STATUS='UNKNOWN')
  OPEN(UNIT=23,FILE='charge_x'//TRIM(STRINGA(transport_iter))//'.dat',STATUS='UNKNOWN')
  OPEN(UNIT=24,FILE='charge_z'//TRIM(STRINGA(transport_iter))//'.dat',STATUS='UNKNOWN')
  DO ii=1, NTOT_X
     DO jj=1, NTOT_Z
        write(22,*)rho_n(ii,NTOT_Y/2,jj),rho_p(ii,NTOT_Y/2,jj)
     END DO
     write(22,*)
  END DO
  DO ii=1, NTOT_Y
     DO jj=1, NTOT_Z
        write(23,*)rho_n(NTOT_X/2,ii,jj),rho_p(NTOT_X/2,ii,jj)
     END DO
     write(23,*)
  END DO
  DO ii=1, NTOT_Y
     DO jj=1, NTOT_X
        write(24,*)rho_n(jj,ii,NTOT_Z/2),rho_p(jj,ii,NTOT_Z/2)
     END DO
     write(24,*)
  END DO
  CLOSE(22)
  CLOSE(23)
  CLOSE(24)


  write(*,*)'////////////////////////////////////////////////////////'
  write(*,*)'                        NEGF SOLVED                     '
  write(*,*)'////////////////////////////////////////////////////////'
  write(*,*)'////////////////////////////////////////////////////////'
  write(*,*)'                     SOLVING  POISSON                   '
  write(*,*)'////////////////////////////////////////////////////////'

  

  CALL map_charge(rho_3D_n,rho_n,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D) 
  CALL map_charge(rho_3D_p,rho_p,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D) 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     
     CALL map_potential(plotpot,EC_3D,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='ec_y.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='ec_x.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)

     CALL map_potential(plotpot,EV_3D,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='ev_y.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='ev_x.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)
     CALL map_potential(plotpot,POT_3D,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='potpoisson_y.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='potpoisson_x.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)


     write(*,*)'IMREF CALCULATION...'
     CALL poisson_imref_n(Fn, EC_3D, rho_3D_n, LWORK_3D)
     CALL poisson_imref_p(Fp, EV_3D, rho_3D_p, LWORK_3D)
     write(*,*)'...IMREF DONE'

     
     CALL map_potential(plotpot,Fn,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='fn_y.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='fn_x.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)
     
     CALL map_potential(plotpot,Fp,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='fp_y.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='fp_x.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)
     
     CALL map_potential(plotpot,dop_vec,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='dop_y.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='dop_x.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)
     
     CALL map_potential(plotpot,rho_3D_n,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='rho_n.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='rho_n.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)
     CALL map_potential(plotpot,rho_3D_p,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='rho_p.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='rho_p.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)

     
     CALL map_potential(plotpot,drho_3D_n,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='drho_n.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='drho_n.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)
     CALL map_potential(plotpot,drho_3D_p,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

     OPEN(UNIT=22,FILE='drho_p.dat',STATUS='UNKNOWN')
     OPEN(UNIT=23,FILE='drho_p.dat',STATUS='UNKNOWN')
     DO ii=1, NTOT_X
        DO jj=1, NTOT_Z
           write(22,*)plotpot(ii,NTOT_Y/2,jj)
        END DO
        write(22,*)
     END DO
     DO ii=1, NTOT_Y
        DO jj=1, NTOT_Z
           write(23,*)plotpot(NTOT_X/2,ii,jj)
        END DO
        write(23,*)
     END DO
     CLOSE(22)
     CLOSE(23)
     


  
  CALL poisson_nonlin_selfconsistent(POT_3D,EC_3D,EV_3D,rho_3D_n,drho_3D_n,Fn,rho_3D_p,drho_3D_p,Fp,dop_vec,coul,LWORK_3D,nnz_3D)

  write(*,*)'////////////////////////////////////////////////////////'
  write(*,*)'                      POISSON SOLVED                    '
  write(*,*)'////////////////////////////////////////////////////////'


  transport_error=0.d0
  DO ii=0,LWORK_3D-1
     transport_error=transport_error+(POT_3D(ii)-outer_pot3D_old(ii))**2
  END DO
  transport_error=sqrt(transport_error/LWORK_3D)

  write(*,*)'AlphaPot=',alphapot

  POT_3D(:)=outer_pot3D_OLD(:) + alphapot*(POT_3D(:)-outer_pot3D_OLD(:))
  outer_pot3D_OLD(:)=POT_3D(:)


  write(*,*) 'Error of the self-consistent cycle = ',transport_error
  open(23, file = TRIM(outdir)//'convergence_LOG.dat',status='old',position='append')
  write(23,*) 'Iter=', transport_iter, 'Err=', transport_error
  close(23)
  
  CALL map_potential(plotpot,rho_3D_n,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)

  OPEN(UNIT=22,FILE='chargepoisson_y.dat',STATUS='UNKNOWN')
  OPEN(UNIT=23,FILE='chargepoisson_x.dat',STATUS='UNKNOWN')
  DO ii=1, NTOT_X
     DO jj=1, NTOT_Z
        write(22,*)plotpot(ii,NTOT_Y/2,jj)
     END DO
     write(22,*)
  END DO
  DO ii=1, NTOT_Y
     DO jj=1, NTOT_Z
        write(23,*)plotpot(NTOT_X/2,ii,jj)
     END DO
     write(23,*)
  END DO
  CLOSE(22)
  CLOSE(23)
  CALL map_potential(plotpot,rho_3D_p,whichkind_3D,map_3D,NTOT_X,NTOT_Y,NTOT_Z,LWORK_3D)
  OPEN(UNIT=22,FILE='chargepoissop_y.dat',STATUS='UNKNOWN')
  OPEN(UNIT=23,FILE='chargepoissop_x.dat',STATUS='UNKNOWN')
  DO ii=1, NTOT_X
     DO jj=1, NTOT_Z
        write(22,*)plotpot(ii,NTOT_Y/2,jj)
     END DO
     write(22,*)
  END DO
  DO ii=1, NTOT_Y
     DO jj=1, NTOT_Z
        write(23,*)plotpot(NTOT_X/2,ii,jj)
     END DO
     write(23,*)
  END DO
  CLOSE(22)
  CLOSE(23)

open(10,file='Last_Potential-vd'//TRIM(STRINGA(ss))//'-vg'//TRIM(STRINGA(gg))//'.dat',status='replace')
do ii=0,LWORK_3D-1
   write(10,*)POT_3D(ii)
end do
close(10)

END DO

t2=SECNDS(t1)
write(*,*)t2,'SECs SPENT FOR THE BIAS POINT ','VG =',vg,'VD =',-mud

IF(NUMBZ.gt.0)sweepbz=.FALSE.

open(23, file = TRIM(outdir)//'convergence_LOG.dat',status='old',position='append')
IF(transport_iter.ge.MAX_ITER_OUTER)THEN
   write(23,*) '::::!!!!MAX ITER NUMBER REACHED!!!!::::'
END IF
write(23,*) 'Convergence achieved in', transport_iter, 'iterations'
write(23,*) 'with error =', transport_error
write(23,*) 'Current =',IDScurrent,ISDcurrent
write(23,*) 'Balistic current =',IDScurrentb
conductance  = abs( IDScurrent  )/( abs(mud-mus) )
conductanceb = abs( IDScurrentb )/( abs(mud-mus) )
write(23,*) 'Conductance=',conductance
close(23)


open(23,file=TRIM(outdir)//'bal_current_vg_'//TRIM(STRINGA(ss))//'.dat',status='old',POSITION='append')
write(23,*) vg, ISDcurrent, IDScurrentB
close(23)
open(23,file=TRIM(outdir)//'bal_current_vd_'//TRIM(STRINGA(gg))//'.dat',status='old',POSITION='append')
write(23,*) -mud, ISDcurrent, IDScurrentB
close(23)



if(phonons)then
open(23,file=TRIM(outdir)//'scat_current_vg_'//TRIM(STRINGA(ss))//'.dat',status='old',POSITION='append')
write(23,*) vg, ISDcurrent, IDScurrent
close(23)
open(23,file=TRIM(outdir)//'scat_current_vd_'//TRIM(STRINGA(gg))//'.dat',status='old',POSITION='append')
write(23,*) -mud, ISDcurrent, IDScurrent
close(23)
end if

END DO
END DO
END DO

   DEALLOCATE(POT_3D)
   DEALLOCATE(outer_pot3D_old)
   DEALLOCATE(EC_3D)
   DEALLOCATE(EV_3D)
   DEALLOCATE(pot3D)
   DEALLOCATE(plotpot)
   DEALLOCATE(rho_3D_n)
   DEALLOCATE(rho_3D_p)
   DEALLOCATE(rho_n)
   DEALLOCATE(rho_p)
   DEALLOCATE(drho_3D_n)
   DEALLOCATE(drho_3D_p)
   DEALLOCATE(Fn)
   DEALLOCATE(Fp)
   DEALLOCATE(Hz)

write(*,*)
write(*,*)
write(*,*) ' This is an open source software distributed under the CeCILL license.' 
write(*,*) ' Please cite  '
write(*,*) ' M. G. Pala, P. Giannozzi, and D. Esseni, Phys. Rev. B 102, 045410 (2020)'
write(*,*) ' DOI: https://doi.org/10.1103/PhysRevB.102.045410'
write(*,*)
write(*,*)


stop


end program main
