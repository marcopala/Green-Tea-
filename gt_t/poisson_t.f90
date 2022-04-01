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

MODULE poisson

  USE static
  USE indata
  USE fermi

  IMPLICIT NONE 
  SAVE

  CONTAINS

SUBROUTINE poisson_nonlin_selfconsistent(pot3D,EC3D,EV3D,outer_rho,outer_drho,Fn,outer_rho_p,outer_drho_p,Fp,doping,coul,LWORK_3D,nnn_3D)
   

    IMPLICIT NONE

    REAL(DP), INTENT(OUT)   :: POT3D(0:LWORK_3D-1) 
    REAL(DP), ALLOCATABLE   :: plotpot(:,:,:)

    REAL(DP), INTENT(INOUT) :: EC3D(0:LWORK_3D-1) 
    REAL(DP), INTENT(INOUT) :: EV3D(0:LWORK_3D-1)    
    REAL(DP), INTENT(IN)    :: Fn(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_rho(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_drho(0:LWORK_3D-1)   
    REAL(DP), INTENT(IN)    :: Fp(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_rho_p(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_drho_p(0:LWORK_3D-1)
    REAL(DP), INTENT(IN)    :: doping(0:LWORK_3D-1)
    INTEGER, INTENT(INOUT)  :: coul(0:LWORK_3D-1)
    INTEGER, INTENT(IN)     :: LWORK_3D
    INTEGER, INTENT(IN)     :: nnn_3D

!!!!!!! Variables for the linear solver !!!!!!!!!!!!!!

    REAL(DP), ALLOCATABLE :: laplacian(:)
    INTEGER,  ALLOCATABLE :: IA(:)
    REAL(DP), ALLOCATABLE :: laplYSMP(:)
    REAL(DP), ALLOCATABLE :: rhs(:)
    REAL(DP), ALLOCATABLE :: valueslap(:) 
    REAL(DP), ALLOCATABLE :: values(:) 
    INTEGER,  ALLOCATABLE :: JA(:)       

!!!!!!! Variables for the PARDISO solver!!!!!!!!!!!!!!!

    INTEGER               :: RCI_request
    INTEGER               :: ipar(1:128)
    REAL(DP)              :: dpar(1:128)
    REAL(DP), ALLOCATABLE :: tmp(:,:)
    INTEGER               :: itercount    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER               :: poisson_iter
    REAL(DP)              :: error

    REAL(DP), ALLOCATABLE :: inner_potold(:)
    REAL(DP), ALLOCATABLE :: rho(:)
    REAL(DP), ALLOCATABLE :: drho(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL(DP)              :: A(1:8,1:8)  
    REAL(DP)              :: XX(1:8), YY(1:8), ZZ(1:8)
    REAL(DP)              :: xi,yi,zi,xj,yj,zj
    INTEGER, ALLOCATABLE  :: IPIV(:)
    INTEGER               :: ii,jj,nel,icc,count
    INTEGER               :: prev_plane, next_plane, out_plane

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*)'STARTING THE ALLOCATION'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ALLOCATE(laplacian(0:LWORK_3D-1))
    ALLOCATE(IA(0:LWORK_3D))
    ALLOCATE(laplYSMP(0:LWORK_3D-1))
    ALLOCATE(rhs(0:LWORK_3D-1))
    ALLOCATE(valueslap(0:nnn_3D-1))
    ALLOCATE(values(0:nnn_3D-1))
    ALLOCATE(JA(0:nnn_3D-1))
    ALLOCATE(tmp(1:LWORK_3D,1:4))
    ALLOCATE(inner_potold(0:LWORK_3D-1))
    ALLOCATE(rho(0:LWORK_3D-1))  
    ALLOCATE(drho(0:LWORK_3D-1))  
    ALLOCATE(IPIV(1:LWORK_3D))  
    ALLOCATE(plotpot(1:NTOT_X,1:NTOT_Y,1:NTOT_Z))

    IA(:)=0 
    valueslap(:)=0.0_dp
    JA(:)=0
    values(:)=0.0_dp
    laplYSMP(:)=0.0_dp
    rhs(:)=0.0_dp
    IPIV(:)=0
    laplacian(:)=0.0_dp
    inner_potold(:)=0.0_dp   
    rho(:)=0.0_dp
    drho(:)=0.0_dp
    A(:,:)=0.0_dp
    XX(:)=0.0_dp
    YY(:)=0.0_dp
    ZZ(:)=0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*)'POISSON ALLOCATION SUCCESSFUL!!',size(POT3D)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    DO ii=1, LWORK_3D-1
    IA(ii) = IA(ii-1) + (connect_3D(ii-1,7)+1)
    END DO

    IA(LWORK_3D)=nnn_3D   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Laplacian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO nel=0, NUMEL_3D-1 
      
    DO jj=1,8 
        XX(jj)=coord_3D_ord(1,list_3D_ord(jj,nel))
        YY(jj)=coord_3D_ord(2,list_3D_ord(jj,nel))
        ZZ(jj)=coord_3D_ord(3,list_3D_ord(jj,nel))
    END DO 

!Flux surfaces
    A(:,:)=0.0_dp
    A(1,1)=deltax*deltay*deltaz/8.0_dp 
    A(2,2)=A(1,1)  
    A(3,3)=A(1,1)  
    A(4,4)=A(1,1)  
    A(5,5)=A(1,1)  
    A(6,6)=A(1,1)  
    A(7,7)=A(1,1)  
    A(8,8)=A(1,1)  

    DO jj=1,8
       xj=coord_3D_ord(1,list_3D_ord(jj,nel))
       yj=coord_3D_ord(2,list_3D_ord(jj,nel))
       zj=coord_3D_ord(3,list_3D_ord(jj,nel))   
    DO ii=jj+1,8
       xi=coord_3D_ord(1,list_3D_ord(ii,nel))
       yi=coord_3D_ord(2,list_3D_ord(ii,nel))
       zi=coord_3D_ord(3,list_3D_ord(ii,nel))
    IF((xi.ne.xj).and.(yi.eq.yj).and.(zi.eq.zj))A(jj,ii)=deltay*deltaz/(4.0_dp*deltax)
    IF((xi.eq.xj).and.(yi.ne.yj).and.(zi.eq.zj))A(jj,ii)=deltax*deltaz/(4.0_dp*deltay)
    IF((xi.eq.xj).and.(yi.eq.yj).and.(zi.ne.zj))A(jj,ii)=deltax*deltay/(4.0_dp*deltaz)
    A(ii,jj)=A(jj,ii)
    END DO
    END DO

    DO jj=0, 7 !Loop on all the element nodes

    IF((whichkind_3D_ord(list_3D_ord(jj+1,nel)).eq.0).or.&
       (whichkind_3D_ord(list_3D_ord(jj+1,nel)).eq.-1) )THEN !Matrix entry

    next_plane =mod(jj+1,4)+4*floor(jj/4.0_dp)+1
    prev_plane =mod(jj+3,4)+4*floor(jj/4.0_dp)+1
    out_plane=mod(jj+4,8)+1   

    JA(IA(list_3D_ord(jj+1,nel)))=list_3D_ord(jj+1,nel)

    valueslap(IA(list_3D_ord(jj+1,nel)))=valueslap(IA(list_3D_ord(jj+1,nel))) + epsilon_3D(nel)*&
    (A(jj+1,next_plane)+A(jj+1,prev_plane)+A(jj+1,out_plane))

    DO icc=1,connect_3D(list_3D_ord(jj+1,nel),7)
              
       IF ( (XX(next_plane)==coord_3D_ord(1,connect_3D(list_3D_ord(jj+1,nel),icc))) .and. &
            (YY(next_plane)==coord_3D_ord(2,connect_3D(list_3D_ord(jj+1,nel),icc))) .and. &
            (ZZ(next_plane)==coord_3D_ord(3,connect_3D(list_3D_ord(jj+1,nel),icc))) )THEN
          
          JA(IA(list_3D_ord(jj+1,nel))+icc)=connect_3D(list_3D_ord(jj+1,nel),icc)
          valueslap(IA(list_3D_ord(jj+1,nel))+icc)=valueslap(IA(list_3D_ord(jj+1,nel))+icc)&
               -epsilon_3D(nel)*(A(jj+1,next_plane))
       END IF

       IF ( (XX(prev_plane)==coord_3D_ord(1,connect_3D(list_3D_ord(jj+1,nel),icc))) .and. &
            (YY(prev_plane)==coord_3D_ord(2,connect_3D(list_3D_ord(jj+1,nel),icc))) .and. &
            (ZZ(prev_plane)==coord_3D_ord(3,connect_3D(list_3D_ord(jj+1,nel),icc))) )THEN
          
          JA(IA(list_3D_ord(jj+1,nel))+icc)=connect_3D(list_3D_ord(jj+1,nel),icc)
          valueslap(IA(list_3D_ord(jj+1,nel))+icc)=valueslap(IA(list_3D_ord(jj+1,nel))+icc)&
               -epsilon_3D(nel)*(A(jj+1,prev_plane))
          
       END IF

       IF ( (XX(out_plane)==coord_3D_ord(1,connect_3D(list_3D_ord(jj+1,nel),icc))) .and. &
            (YY(out_plane)==coord_3D_ord(2,connect_3D(list_3D_ord(jj+1,nel),icc))) .and. &
            (ZZ(out_plane)==coord_3D_ord(3,connect_3D(list_3D_ord(jj+1,nel),icc))) )THEN
          
          JA(IA(list_3D_ord(jj+1,nel))+icc)=connect_3D(list_3D_ord(jj+1,nel),icc)
          valueslap(IA(list_3D_ord(jj+1,nel))+icc)=valueslap(IA(list_3D_ord(jj+1,nel))+icc)&
               -epsilon_3D(nel)*(A(jj+1,out_plane))
          
       END IF
       
           END DO
           
        END IF
        
     END DO
  END DO
  
!////////////////////////////////////////////////////////////////////
!//////////////////// END INITIAL MATRIX ////////////////////////////
!////////////////////////////////////////////////////////////////////

  error=1.0_dp
  poisson_iter=0
  inner_potold(:)=POT3D(:)

  DO WHILE ((error.ge.ERROR_INNER).and.(poisson_iter.le.MAX_ITER_INNER))

  poisson_iter=poisson_iter+1

  CALL poisson_charge_from_imref_n(outer_rho, EC3D, Fn, LWORK_3D)
  CALL poisson_deriv_from_imref_n(outer_drho, EC3D, Fn, LWORK_3D)
  CALL poisson_charge_from_imref_p(outer_rho_p, EV3D, Fp, LWORK_3D)
  CALL poisson_deriv_from_imref_p(outer_drho_p, EV3D, Fp, LWORK_3D)

  rho(:)=-ELCH*(outer_rho(:)-outer_rho_p(:))
  drho(:)=-ELCH*(outer_drho(:)-outer_drho_p(:))

  values(:)=valueslap(:)
  CALL multYSMP(nnn_3D,LWORK_3D,valueslap,JA,IA,inner_potold,laplYSMP)
  rhs(:)=-laplYSMP(:)


!////////////////////////////////////////////////////////////////////
!////////////// MATRICE INIZIALE DI SOLUZIONE ///////////////////////
!////////////////////////////////////////////////////////////////////


    DO nel=0, NUMEL_3D-1 !ciclo su tutti gli elementi
  
    DO jj=1,8 
        XX(jj)=coord_3D_ord(1,list_3D_ord(jj,nel))
        YY(jj)=coord_3D_ord(2,list_3D_ord(jj,nel))
        ZZ(jj)=coord_3D_ord(3,list_3D_ord(jj,nel))
    END DO
 
   
!Flux surfaces
    A(:,:)=0.0_dp

    A(1,1)=deltax*deltay*deltaz/8.0_dp !Subvolume of the box inherent to a node
    A(2,2)=A(1,1)  
    A(3,3)=A(1,1)  
    A(4,4)=A(1,1)  
    A(5,5)=A(1,1)  
    A(6,6)=A(1,1)  
    A(7,7)=A(1,1)  
    A(8,8)=A(1,1)  

    DO jj=1,8
             xj=coord_3D_ord(1,list_3D_ord(jj,nel))
             yj=coord_3D_ord(2,list_3D_ord(jj,nel))
             zj=coord_3D_ord(3,list_3D_ord(jj,nel))   
    DO ii=jj+1,8
             xi=coord_3D_ord(1,list_3D_ord(ii,nel))
             yi=coord_3D_ord(2,list_3D_ord(ii,nel))
             zi=coord_3D_ord(3,list_3D_ord(ii,nel))
    IF((xi.ne.xj).and.(yi.eq.yj).and.(zi.eq.zj))A(jj,ii)=deltay*deltaz/(4.0_dp*deltax)
    IF((xi.eq.xj).and.(yi.ne.yj).and.(zi.eq.zj))A(jj,ii)=deltax*deltaz/(4.0_dp*deltay)
    IF((xi.eq.xj).and.(yi.eq.yj).and.(zi.ne.zj))A(jj,ii)=deltax*deltay/(4.0_dp*deltaz)
    A(ii,jj)=A(jj,ii)
    END DO
    END DO
 
    DO jj=0, 7


    IF((whichkind_3D_ord(list_3D_ord(jj+1,nel)).eq.0).or. &
       (whichkind_3D_ord(list_3D_ord(jj+1,nel)).eq.-1) )THEN !matrix node
   
    next_plane =mod(jj+1,4)+4*floor(jj/4.0_dp)+1
    prev_plane =mod(jj+3,4)+4*floor(jj/4.0_dp)+1
    out_plane=mod(jj+4,8)+1 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        values(IA(list_3D_ord(jj+1,nel)))=values(IA(list_3D_ord(jj+1,nel)))&
                                        -drho(list_3D_ord(jj+1,nel))*A(jj+1,jj+1)            

        rhs(list_3D_ord(jj+1,nel))=rhs(list_3D_ord(jj+1,nel))+&
                             rho(list_3D_ord(jj+1,nel))*A(jj+1,jj+1)+&
                             doping(list_3D_ord(jj+1,nel))*A(jj+1,jj+1)

        !Charged defect
        IF(coul(list_3D_ord(jj+1,nel)).eq.1)THEN
        rhs(list_3D_ord(jj+1,nel))=rhs(list_3D_ord(jj+1,nel))+&
                             doping(list_3D_ord(jj+1,nel))*A(jj+1,jj+1)
        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF(whichkind_3D_ord(list_3D_ord(prev_plane,nel)).eq.1) THEN !
           rhs(list_3D_ord(jj+1,nel))=rhs(list_3D_ord(jj+1,nel))+&
                epsilon_3D(nel)*A(jj+1,prev_plane)*&
                potelectr 
        END IF
        
        IF(whichkind_3D_ord(list_3D_ord(next_plane,nel)).eq.1) THEN !
           
           rhs(list_3D_ord(jj+1,nel))=rhs(list_3D_ord(jj+1,nel))+&
                epsilon_3D(nel)*A(jj+1,next_plane)*&
                potelectr
           
        END IF
        
        IF(whichkind_3D_ord(list_3D_ord(out_plane,nel)).eq.1) THEN !
           rhs(list_3D_ord(jj+1,nel))=rhs(list_3D_ord(jj+1,nel))+&
                epsilon_3D(nel)*A(jj+1,out_plane)*&
                potelectr 
        END IF
        
     END IF
     
  END DO
  
END DO
   

  POT3D(:)=0.0_dp

  ipar=0
  dpar=0.0_dp
  tmp=0.0_dp
  RCI_request=0
  count=0

  CALL dcg_init&
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, &
       ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))

  IF(RCI_request.ne.0)THEN
     write(*,*)'dcg_init error =',RCI_request
  STOP
  END IF

  ipar(5)=5000
  ipar(9)=1
  ipar(10)=0
  dpar(1)=1.0d-5

  CALL dcg_check&
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, &
       ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))

  IF(RCI_request.ne.0)THEN
     write(*,*)'dcg_check error=',RCI_request
  STOP
  END IF

  CALL dcg(LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, &
       ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))

  DO WHILE(RCI_request.eq.1)

  count=count+1

  CALL multYSMP(nnn_3D,LWORK_3D,values,JA,IA,tmp(1:LWORK_3D,1),tmp(1:LWORK_3D,2))

  CALL dcg&
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, &
       ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))
 
  END DO

  IF(RCI_request.eq.0)THEN

  CALL dcg_get&
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, &
       ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4), itercount)  

  ELSE

  write(*,*)'dcg_solver error=',RCI_request
  write(*,*)'iterations =', count
  STOP

  END IF


  error=0._dp
  DO ii=0,LWORK_3D-1
     error=error+POT3D(ii)*POT3D(ii)
  END DO
  error=sqrt(error/dble(LWORK_3D))
  write(*,*) 'Error = ',error


  POT3D(:)=inner_potold(:)+POT3D(:)
 
  inner_potold(:)=POT3D(:) 
 
  CALL shift_potential(EC3D,-POT3D,3.0_dp,E_GAP,whichkind_3D,map_3D,list_3D_ord,&
       whichkind_3D_ord,epsilon_3D,NTOT_X,NTOT_Y,NTOT_Z,NUMEL_3D,LWORK_3D)
  CALL shift_potential(EV3D,-POT3D,-3.0_dp,0.0_dp,whichkind_3D,map_3D,list_3D_ord,&
       whichkind_3D_ord,epsilon_3D,NTOT_X,NTOT_Y,NTOT_Z,NUMEL_3D,LWORK_3D)
 
 
END DO !END DO WHILE


    DEALLOCATE(laplacian)
    DEALLOCATE(IA)
    DEALLOCATE(laplYSMP)
    DEALLOCATE(rhs)
    DEALLOCATE(valueslap)
    DEALLOCATE(values)
    DEALLOCATE(JA)
    DEALLOCATE(tmp)
    DEALLOCATE(inner_potold)
    DEALLOCATE(rho)
    DEALLOCATE(drho)
    deallocate(ipiv)
    deallocate(plotpot)


END SUBROUTINE poisson_nonlin_selfconsistent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine multYSMP(nz,dim,A,JA,IA,B,PR)
  implicit none
  
  integer,  intent(in) :: nz, dim 
  REAL(DP), intent(in) :: A(0:nz-1),B(0:dim-1)
  integer,  intent(in) :: JA(0:nz-1),IA(0:dim)
  REAL(DP),intent(out) :: PR(0:dim-1)

  integer :: ii,jj

  PR(:)=0.0_dp
  do ii=0,dim-1
     do jj=IA(ii), (IA(ii+1)-1)
        PR(ii)=PR(ii)+A(jj)*B(JA(jj))
     end do
  end do
  
end subroutine multYSMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE pot_init(POT3D,which,map,nx,ny,nz,lwork)

IMPLICIT NONE

REAL(DP),INTENT(OUT) :: POT3D(0:lwork-1)

INTEGER,  INTENT(IN) :: which(0:nx*ny*nz-1)
INTEGER,  INTENT(IN) :: map(0:nx*ny*nz-1)
INTEGER,  INTENT(IN) :: nx
INTEGER,  INTENT(IN) :: ny
INTEGER,  INTENT(IN) :: nz
INTEGER,  INTENT(IN) :: lwork
INTEGER              :: ii, x_index

 POT3D=0.0
 DO ii=0, nx*ny*nz-1
 x_index=ii/(ny*nz)

IF(which(ii).le.0)THEN

if(chtype.eq.'n')then

   IF(x_index.le.source_len)THEN
      POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))-mus
   ELSEIF(x_index .gt. source_len+2*spacer+gate_len .and. &
        x_index .lt. source_len+2*spacer+gate_len+drain_len)THEN
      POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))-mud
   ELSEIF(x_index .gt. source_len .and. &
      x_index .le. source_len+2*spacer+gate_len)THEN
      POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))-0.1_dp 
   END IF
   IF(x_index.eq.source_len+2*spacer+gate_len+drain_len)THEN
      POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax))-mud 
   END IF
   
elseif(chtype.eq.'p')then
      
   IF(x_index.le.source_len)THEN
      POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))-mus
   ELSEIF(x_index.gt.source_len+2*spacer+gate_len .and. &
        x_index .lt. source_len+2*spacer+gate_len+drain_len)THEN
      POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))-mud
   ELSEIF(x_index.gt.source_len.and.x_index.le.source_len+2*spacer+gate_len)THEN
      POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))+0.1_dp
   END IF
   IF(x_index.eq.source_len+2*spacer+gate_len+drain_len)THEN
      POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax))-mud 
   END IF

elseif(chtype.eq.'t')then
   
   IF(x_index.le.source_len)THEN
      POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))-mus
   ELSEIF(x_index.gt.source_len+2*spacer+gate_len .and. &
        x_index .lt. source_len+2*spacer+gate_len+drain_len)THEN
      POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))-mud
   ELSEIF(x_index.gt.source_len.and.x_index.le.source_len+2*spacer+gate_len)THEN
      POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))
   END IF
   IF(x_index.eq.source_len+2*spacer+gate_len+drain_len)THEN
      POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax))-mud 
   END IF

end if
   
END IF
END DO

END SUBROUTINE pot_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     
  SUBROUTINE poisson_imref_n(Fn, Ec, rho, lwork)

  IMPLICIT NONE

  REAL(DP), INTENT(INOUT) :: Fn(0:lwork-1)      
  REAL(DP), INTENT(IN)    :: Ec(0:lwork-1) 
  REAL(DP), INTENT(IN)    :: rho(0:lwork-1)
  INTEGER, INTENT(IN)     :: lwork

  REAL(DP) :: N_c, N3D, f, df, TOLERANCE, delta, fn_old
  INTEGER :: ii, nn  
  
  Fn=3.0_dp

  N_c=6.2d17
  N3D=N_c*(TEMP)**(1.5) 
  
  TOLERANCE=1.0d-15
  nn=0

  DO ii=0, lwork-1 

  f=rho(ii)-N3D*FDP0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP))
  df=N3D*(FDM0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)
  delta=-f/df
  fn_old=Fn(ii)
  if(rho(ii).eq.0.0_DP)then
  Fn(ii)=100.0_dp
  else
  nn=0
  DO WHILE((ABS(delta).gt.TOLERANCE).and.(nn.lt.400))

  Fn(ii)=Fn(ii)-delta
  f=rho(ii)-N3D*FDP0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP))
  df=N3D*(FDM0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP)))/(BOLTZ*TEMP)
  delta=-f/df
  
  nn=nn+1
  if(nn.gt.500)then
     write(*,*)'problem in imref',nn
     stop
  endif
  END DO
  
  end if

  IF(ABS(Fn(ii)).gt.100.0_dp)Fn(ii)=3.0_dp
  
END DO

END SUBROUTINE poisson_imref_n

  SUBROUTINE poisson_imref_p(Fp, Ev, rho, lwork)

  IMPLICIT NONE

  REAL(DP), INTENT(INOUT) :: Fp(0:lwork-1)     
  REAL(DP), INTENT(IN)    :: Ev(0:lwork-1) 
  REAL(DP), INTENT(IN)    :: rho(0:lwork-1)
  INTEGER, INTENT(IN)     :: lwork

  REAL(DP) :: N_v, N3D, f, df, TOLERANCE, delta, fp_old
  INTEGER :: ii, nn  
  
  Fp=-3.0_dp

  N_v=6.5d19
  N3D=N_v*(TEMP)**(1.5) 
  
  TOLERANCE=1.0d-15
  nn=0

  DO ii=0, lwork-1
  f=rho(ii)-N3D*FDP0P5(-(Fp(ii)-Ev(ii))/(BOLTZ*TEMP))
  df=N3D*(FDM0P5(-(Fp(ii)-Ev(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)
  delta=f/df
  fp_old=Fp(ii)
  if(rho(ii).lt.1.0_DP)then
  Fp(ii)=100.0_dp
  else
  nn=0 
  DO WHILE((ABS(delta).gt.TOLERANCE).and.(nn.lt.200))
     Fp(ii)=Fp(ii)-delta
     f=rho(ii)-N3D*FDP0P5(-(Fp(ii)-Ev(ii))/(BOLTZ*TEMP))
     df=N3D*(FDM0P5(-(Fp(ii)-Ev(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)
     delta=f/df
     nn=nn+1
  END DO
  
  end if

  if(nn.ge.200)Fp(ii)=fp_old
  
  END DO

END SUBROUTINE poisson_imref_p


SUBROUTINE poisson_charge_from_imref_n(rho, Ec, Fn, lwork)

  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: rho(0:lwork-1)
  REAL(DP), INTENT(IN)  :: Ec(0:lwork-1)
  REAL(DP), INTENT(IN)  :: Fn(0:lwork-1)
  INTEGER,  INTENT(IN)  :: lwork
  REAL(DP) :: N_c, N3D
  INTEGER :: ii


  N_c=6.2d17
  N3D=N_c*(TEMP)**(1.5) 

  rho=0.0_dp
  DO ii=0, lwork-1
    IF(Fn(ii).ne.100.0_dp)rho(ii)=N3D*FDP0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP))   
  END DO
 
END SUBROUTINE poisson_charge_from_imref_n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE poisson_charge_from_imref_p(rho, Ev, Fp, lwork)

  IMPLICIT NONE
  
  REAL(DP), INTENT(OUT) :: rho(0:lwork-1)   
  REAL(DP), INTENT(IN)  :: Ev(0:lwork-1)
  REAL(DP), INTENT(IN)  :: Fp(0:lwork-1)
  INTEGER,  INTENT(IN)  :: lwork
    
  REAL(DP) :: N_v, N3D
  INTEGER :: ii
    
  N_v=6.5d19
  N3D=N_v*(TEMP)**(1.5) 

  rho=0.0_dp
    
  DO ii=0, lwork-1
     IF(Fp(ii).ne.100.0_dp)rho(ii)=N3D*FDP0P5(-(Fp(ii)-Ev(ii))/(BOLTZ*TEMP))   
  END DO
        
END SUBROUTINE poisson_charge_from_imref_p

SUBROUTINE poisson_deriv_from_imref_n(drho, Ec, Fn, lwork)

  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: drho(0:lwork-1)
  REAL(DP), INTENT(IN)  :: Ec(0:lwork-1)
  REAL(DP), INTENT(IN)  :: Fn(0:lwork-1)
  INTEGER,  INTENT(IN)  :: lwork
  REAL(DP) :: N_c, N3D
  INTEGER :: ii
  
  N_c=6.2d17
  N3D=N_c*(TEMP)**(1.5) 

  drho=0.0_dp

  DO ii=0, lwork-1
     IF(Fn(ii).ne.100.0_dp)THEN 
        drho(ii)=N3D*( FDM0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)        
     END IF
  END DO
  
END SUBROUTINE poisson_deriv_from_imref_n

SUBROUTINE poisson_deriv_from_imref_p(drho, Ev, Fp, lwork)
    
    IMPLICIT NONE
    
    REAL(DP), INTENT(OUT) :: drho(0:lwork-1)  
    REAL(DP), INTENT(IN)  :: Ev(0:lwork-1)
    REAL(DP), INTENT(IN)  :: Fp(0:lwork-1)
    INTEGER,  INTENT(IN)  :: lwork
    
    REAL(DP) :: N_v, N3D
    INTEGER :: ii
     
    N_v=6.5d19
    N3D=N_v*(TEMP)**(1.5) 
   
    drho=0.0_dp
    DO ii=0, lwork-1
       IF(Fp(ii).ne.100.0_dp)THEN
        drho(ii)=-N3D*( FDM0P5(-(Fp(ii)-Ev(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)
       END IF
    END DO
    
END SUBROUTINE poisson_deriv_from_imref_p


END MODULE poisson
