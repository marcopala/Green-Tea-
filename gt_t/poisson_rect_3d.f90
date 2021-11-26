MODULE poisson

  USE static
  USE indata
  USE fermi

  IMPLICIT NONE 
  SAVE

  CONTAINS

SUBROUTINE poisson_nonlin_selfconsistent(pot3D,EC3D,EV3D,outer_rho,outer_drho,Fn,outer_rho_p,outer_drho_p,Fp,doping,coul,LWORK_3D,nnz_3D)
   

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: POT3D(0:LWORK_3D-1) 
    REAL(DP), ALLOCATABLE :: plotpot(:,:,:)

    REAL(DP), INTENT(INOUT) :: EC3D(0:LWORK_3D-1) 
    REAL(DP), INTENT(INOUT) :: EV3D(0:LWORK_3D-1)    
    REAL(DP), INTENT(IN)    :: Fn(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_rho(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_drho(0:LWORK_3D-1)   
    REAL(DP), INTENT(IN)    :: Fp(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_rho_p(0:LWORK_3D-1)
    REAL(DP), INTENT(INOUT) :: outer_drho_p(0:LWORK_3D-1)
    REAL(DP), INTENT(IN) :: doping(0:LWORK_3D-1)
    INTEGER, INTENT(INOUT) :: coul(0:LWORK_3D-1)
    INTEGER, INTENT(IN) :: LWORK_3D
    INTEGER, INTENT(IN) :: nnz_3D

!!!!!!! Variables for the linear solver !!!!!!!!!!!!!!

    REAL(DP), ALLOCATABLE :: laplacian(:)
    INTEGER, ALLOCATABLE          :: IA(:)
    REAL(DP), ALLOCATABLE :: laplYSMP(:)
    REAL(DP), ALLOCATABLE :: rhs(:)
    REAL(DP), ALLOCATABLE :: valueslap(:) 
    REAL(DP), ALLOCATABLE :: values(:) 
    INTEGER, ALLOCATABLE          :: JA(:)       

!!!!!!! Variables for the PARDISO solver!!!!!!!!!!!!!!!

    INTEGER                       :: RCI_request
    !INTEGER,ALLOCATABLE           :: ipar(:)
    !REAL(DP),ALLOCATABLE           :: dpar(:)
    INTEGER          :: ipar(1:128)
    REAL(DP)          :: dpar(1:128)
    REAL(DP), ALLOCATABLE :: tmp(:,:)
    INTEGER                       :: itercount    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER                       :: poisson_iter
    REAL(DP)              :: error

    REAL(DP), ALLOCATABLE :: inner_potold(:)
    REAL(DP), ALLOCATABLE :: rho(:)
    REAL(DP), ALLOCATABLE :: drho(:)

!    REAL(DP), ALLOCATABLE :: rhovec(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!number of nodes/element =8!!!!!!!!!!!!!!!!!!
!!!!nuber of fluxes for every node in an element = 3!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !REAL(DP),ALLOCATABLE  :: A(:,:)
    REAL(DP)               :: A(1:8,1:8)  
    REAL(DP)               :: XX(1:8), YY(1:8), ZZ(1:8)
    REAL(DP) :: xi,yi,zi,xj,yj,zj
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER :: ii,jj,nel,cc,count
    INTEGER :: prec_piano, succ_piano, fuori_piano

    LOGICAL :: first


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)'START ALLOCATION'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !ALLOCATE(ipar(1:128))
    !ALLOCATE(dpar(1:128))
    ALLOCATE(laplacian(0:LWORK_3D-1))
    ALLOCATE(IA(0:LWORK_3D))
    ALLOCATE(laplYSMP(0:LWORK_3D-1))
    ALLOCATE(rhs(0:LWORK_3D-1))
    ALLOCATE(valueslap(0:nnz_3D-1))
    ALLOCATE(values(0:nnz_3D-1))
    ALLOCATE(JA(0:nnz_3D-1))
    ALLOCATE(tmp(1:LWORK_3D,1:4))
    ALLOCATE(inner_potold(0:LWORK_3D-1))
    ALLOCATE(rho(0:LWORK_3D-1))  
    ALLOCATE(drho(0:LWORK_3D-1))  
    ALLOCATE(IPIV(1:LWORK_3D))  
!    ALLOCATE(rhovec(1:NTOT_X,1:NTOT_Y,1:NTOT_Z))
    ALLOCATE(plotpot(1:NTOT_X,1:NTOT_Y,1:NTOT_Z))
    !ALLOCATE(A(1:8,1:8))
    !ALLOCATE(XX(1:8))
    !ALLOCATE(YY(1:8))
    !ALLOCATE(ZZ(1:8))
    IA(:)=0 
    valueslap(:)=0.0_dp
    JA(:)=0
    values(:)=0.0_dp
    laplYSMP(:)=0.0_dp
    rhs(:)=0._dp
    IPIV(:)=0
    laplacian(:)=0.0_dp
    inner_potold(:)=0.0_dp   
    rho(:)=0.0_dp
    drho(:)=0.0_dp
    A(:,:)=0.0_dp
    XX(:)=0.0_dp
    YY(:)=0.0_dp
    ZZ(:)=0.0_dp

  first=.TRUE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)'POISSON ALLOCATION SUCCESSFUL!!',size(POT3D)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    DO ii=1, LWORK_3D-1

    IA(ii) = IA(ii-1) + (connect_3D(ii-1,7)+1)

    END DO

    IA(LWORK_3D)=nnz_3D   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!Calcolo della matrice del Laplaciano!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    DO nel=0, NUMEL_3D-1 !ciclo su tutti gli elementi
      
    DO jj=1,8 
        XX(jj)=coord_3D_ord(1,lista_3D_ord(jj,nel))
        YY(jj)=coord_3D_ord(2,lista_3D_ord(jj,nel))
        ZZ(jj)=coord_3D_ord(3,lista_3D_ord(jj,nel))
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
             xj=coord_3D_ord(1,lista_3D_ord(jj,nel))
             yj=coord_3D_ord(2,lista_3D_ord(jj,nel))
             zj=coord_3D_ord(3,lista_3D_ord(jj,nel))   
    DO ii=jj+1,8
             xi=coord_3D_ord(1,lista_3D_ord(ii,nel))
             yi=coord_3D_ord(2,lista_3D_ord(ii,nel))
             zi=coord_3D_ord(3,lista_3D_ord(ii,nel))
    IF((xi.ne.xj).and.(yi.eq.yj).and.(zi.eq.zj))A(jj,ii)=deltay*deltaz/(4.0_dp*deltax)
    IF((xi.eq.xj).and.(yi.ne.yj).and.(zi.eq.zj))A(jj,ii)=deltax*deltaz/(4.0_dp*deltay)
    IF((xi.eq.xj).and.(yi.eq.yj).and.(zi.ne.zj))A(jj,ii)=deltax*deltay/(4.0_dp*deltaz)
    A(ii,jj)=A(jj,ii)
    END DO
    END DO

    DO jj=0, 7 !Loop on all the element nodes

    IF((whichkind_3D_ord(lista_3D_ord(jj+1,nel)).eq.0).or.&
       (whichkind_3D_ord(lista_3D_ord(jj+1,nel)).eq.-1) )THEN !Matrix entry

    succ_piano =mod(jj+1,4)+4*floor(jj/4.0_dp)+1
    prec_piano =mod(jj+3,4)+4*floor(jj/4.0_dp)+1
    fuori_piano=mod(jj+4,8)+1   


    JA(IA(lista_3D_ord(jj+1,nel)))=lista_3D_ord(jj+1,nel)

    valueslap(IA(lista_3D_ord(jj+1,nel)))=valueslap(IA(lista_3D_ord(jj+1,nel))) + epsilon_3D(nel)*&
    (A(jj+1,succ_piano)+A(jj+1,prec_piano)+A(jj+1,fuori_piano))


            DO cc=1,connect_3D(lista_3D_ord(jj+1,nel),7)


              !vicino sul piano (successivo)
              IF ( (XX(succ_piano)==coord_3D_ord(1,connect_3D(lista_3D_ord(jj+1,nel),cc))) .and. &
                   (YY(succ_piano)==coord_3D_ord(2,connect_3D(lista_3D_ord(jj+1,nel),cc))) .and. &
                   (ZZ(succ_piano)==coord_3D_ord(3,connect_3D(lista_3D_ord(jj+1,nel),cc))) )THEN
                 
                JA(IA(lista_3D_ord(jj+1,nel))+cc)=connect_3D(lista_3D_ord(jj+1,nel),cc)
                valueslap(IA(lista_3D_ord(jj+1,nel))+cc)=valueslap(IA(lista_3D_ord(jj+1,nel))+cc)&
                     -epsilon_3D(nel)*(A(jj+1,succ_piano))
              END IF
              !vicino sul piano (precedente)
              IF ( (XX(prec_piano)==coord_3D_ord(1,connect_3D(lista_3D_ord(jj+1,nel),cc))) .and. &
                   (YY(prec_piano)==coord_3D_ord(2,connect_3D(lista_3D_ord(jj+1,nel),cc))) .and. &
                   (ZZ(prec_piano)==coord_3D_ord(3,connect_3D(lista_3D_ord(jj+1,nel),cc))) )THEN
                 
                JA(IA(lista_3D_ord(jj+1,nel))+cc)=connect_3D(lista_3D_ord(jj+1,nel),cc)
                valueslap(IA(lista_3D_ord(jj+1,nel))+cc)=valueslap(IA(lista_3D_ord(jj+1,nel))+cc)&
                     -epsilon_3D(nel)*(A(jj+1,prec_piano))

              END IF
              !vicino fuori piano
              IF ( (XX(fuori_piano)==coord_3D_ord(1,connect_3D(lista_3D_ord(jj+1,nel),cc))) .and. &
                   (YY(fuori_piano)==coord_3D_ord(2,connect_3D(lista_3D_ord(jj+1,nel),cc))) .and. &
                   (ZZ(fuori_piano)==coord_3D_ord(3,connect_3D(lista_3D_ord(jj+1,nel),cc))) )THEN
                 
                JA(IA(lista_3D_ord(jj+1,nel))+cc)=connect_3D(lista_3D_ord(jj+1,nel),cc)
                valueslap(IA(lista_3D_ord(jj+1,nel))+cc)=valueslap(IA(lista_3D_ord(jj+1,nel))+cc)&
                     -epsilon_3D(nel)*(A(jj+1,fuori_piano))

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

  !write(*,*)'ITERAZIONE CONSISTENZA INTERNA:', poisson_iter


  CALL poisson_charge_from_imref_n(outer_rho, EC3D, Fn, LWORK_3D)
  CALL poisson_deriv_from_imref_n(outer_drho, EC3D, Fn, LWORK_3D)
  CALL poisson_charge_from_imref_p(outer_rho_p, EV3D, Fp, LWORK_3D)
  CALL poisson_deriv_from_imref_p(outer_drho_p, EV3D, Fp, LWORK_3D)

  rho(:)=-ELCH*(outer_rho(:)-outer_rho_p(:))
  drho(:)=-ELCH*(outer_drho(:)-outer_drho_p(:))



  values(:)=valueslap(:)
  CALL prodottoYSMP(nnz_3D,LWORK_3D,valueslap,JA,IA,inner_potold,laplYSMP)
  rhs(:)=-laplYSMP(:)


!////////////////////////////////////////////////////////////////////
!////////////// MATRICE INIZIALE DI SOLUZIONE ///////////////////////
!////////////////////////////////////////////////////////////////////



    DO nel=0, NUMEL_3D-1 !ciclo su tutti gli elementi
  
    DO jj=1,8 
        XX(jj)=coord_3D_ord(1,lista_3D_ord(jj,nel))
        YY(jj)=coord_3D_ord(2,lista_3D_ord(jj,nel))
        ZZ(jj)=coord_3D_ord(3,lista_3D_ord(jj,nel))
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
             xj=coord_3D_ord(1,lista_3D_ord(jj,nel))
             yj=coord_3D_ord(2,lista_3D_ord(jj,nel))
             zj=coord_3D_ord(3,lista_3D_ord(jj,nel))   
    DO ii=jj+1,8
             xi=coord_3D_ord(1,lista_3D_ord(ii,nel))
             yi=coord_3D_ord(2,lista_3D_ord(ii,nel))
             zi=coord_3D_ord(3,lista_3D_ord(ii,nel))
    IF((xi.ne.xj).and.(yi.eq.yj).and.(zi.eq.zj))A(jj,ii)=deltay*deltaz/(4.0_dp*deltax)
    IF((xi.eq.xj).and.(yi.ne.yj).and.(zi.eq.zj))A(jj,ii)=deltax*deltaz/(4.0_dp*deltay)
    IF((xi.eq.xj).and.(yi.eq.yj).and.(zi.ne.zj))A(jj,ii)=deltax*deltay/(4.0_dp*deltaz)
    A(ii,jj)=A(jj,ii)
    END DO
    END DO
 
    DO jj=0, 7


    IF((whichkind_3D_ord(lista_3D_ord(jj+1,nel)).eq.0).or. &
       (whichkind_3D_ord(lista_3D_ord(jj+1,nel)).eq.-1) )THEN !nodo di matrice

   
    succ_piano =mod(jj+1,4)+4*floor(jj/4.0_dp)+1
    prec_piano =mod(jj+3,4)+4*floor(jj/4.0_dp)+1
    fuori_piano=mod(jj+4,8)+1 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!Variante aalocando la carica per porzioni
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        values(IA(lista_3D_ord(jj+1,nel)))=values(IA(lista_3D_ord(jj+1,nel)))&
                                        -drho(lista_3D_ord(jj+1,nel))*A(jj+1,jj+1)            

        rhs(lista_3D_ord(jj+1,nel))=rhs(lista_3D_ord(jj+1,nel))+&
                             rho(lista_3D_ord(jj+1,nel))*A(jj+1,jj+1)+&
                             doping(lista_3D_ord(jj+1,nel))*A(jj+1,jj+1)

        !Difetto puntuale

        IF(coul(lista_3D_ord(jj+1,nel)).eq.1)THEN
        rhs(lista_3D_ord(jj+1,nel))=rhs(lista_3D_ord(jj+1,nel))+&
                             doping(lista_3D_ord(jj+1,nel))*A(jj+1,jj+1)
        !write(*,*)'k'
        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           IF(whichkind_3D_ord(lista_3D_ord(prec_piano,nel)).eq.1) THEN !prec nodo sul piano
              rhs(lista_3D_ord(jj+1,nel))=rhs(lista_3D_ord(jj+1,nel))+&
                   epsilon_3D(nel)*A(jj+1,prec_piano)*&
                   potelectr!(whichkind_3D_ord(lista_3D_ord(prec_piano,nel)))
           END IF
           
           IF(whichkind_3D_ord(lista_3D_ord(succ_piano,nel)).eq.1) THEN !!suc nodo sul piano

              rhs(lista_3D_ord(jj+1,nel))=rhs(lista_3D_ord(jj+1,nel))+&
                   epsilon_3D(nel)*A(jj+1,succ_piano)*&
                   potelectr

           END IF

           IF(whichkind_3D_ord(lista_3D_ord(fuori_piano,nel)).eq.1) THEN !nodo fuori piano
             rhs(lista_3D_ord(jj+1,nel))=rhs(lista_3D_ord(jj+1,nel))+&
                   epsilon_3D(nel)*A(jj+1,fuori_piano)*&
                   potelectr!(whichkind_3D_ord(lista_3D_ord(fuori_piano,nel)))
           END IF

           !IF(whichkind_3D_ord(lista_3D_ord(fuori_piano,nel)).eq.11) THEN !nodo fuori piano
           !  rhs(lista_3D_ord(jj+1,nel))=rhs(lista_3D_ord(jj+1,nel))+&
           !        epsilon_3D(nel)*A(jj+1,fuori_piano)*&
           !        potschottky!(whichkind_3D_ord(lista_3D_ord(fuori_piano,nel)))
           !END IF

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
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))

  IF(RCI_request.ne.0)THEN
     write(*,*)'dcg_init error =',RCI_request
  STOP
  END IF

  !Imposto ipar in modo tale che il solutore
  !gestisca in modo autonomo la convergenza sul residuo
  ipar(5)=5000
  ipar(9)=1
  ipar(10)=0
  !Modifica il valore della tolleranza sul residuo
  !Il valore di default e' 1.0e-06
  dpar(1)=1.0d-5

  CALL dcg_check&
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))

  IF(RCI_request.ne.0)THEN
     write(*,*)'dcg_check error=',RCI_request
  STOP
  END IF

  CALL dcg(LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))

  DO WHILE(RCI_request.eq.1)

  count=count+1

  CALL prodottoYSMP(nnz_3d,LWORK_3D,values,JA,IA,tmp(1:LWORK_3D,1),tmp(1:LWORK_3D,2))
!!$  CALL mkl_dcsrgemv&
!!$       ('N',LWORK,values(0:nnz-1),colptr(0:LWORK),rowind(0:nnz-1),tmp(1:LWORK,1),tmp(1:LWORK,2))

  CALL dcg&
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4))
 
  END DO

  IF(RCI_request.eq.0)THEN

  CALL dcg_get&
       (LWORK_3D, POT3D(0:LWORK_3D-1), rhs(0:LWORK_3D-1), RCI_request, ipar(1:128), dpar(1:128), tmp(1:LWORK_3D,1:4), itercount)  
 
!!$  write(*,*)'System successfully solved'
!!$  write(*,*)'Number of iteraitons', itercount

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
  write(*,*) 'Errore = ',error



  POT3D(:)=inner_potold(:)+POT3D(:)
 
  inner_potold(:)=POT3D(:) 
 
  CALL shift_potential(EC3D,-POT3D,3.0_dp,E_GAP,whichkind_3D,map_3D,lista_3D_ord,whichkind_3D_ord,epsilon_3D,NTOT_X,NTOT_Y,NTOT_Z,NUMEL_3D,LWORK_3D)
  CALL shift_potential(EV3D,-POT3D,-3.0_dp,0.0_dp,whichkind_3D,map_3D,lista_3D_ord,whichkind_3D_ord,epsilon_3D,NTOT_X,NTOT_Y,NTOT_Z,NUMEL_3D,LWORK_3D)
 
 
  END DO!Chiusura DO WHILE




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



subroutine prodottoYSMP(nnz,dim,A,JA,IA,B,PR)
  implicit none
  
  integer, intent(in) :: nnz, dim !dimensione del vettore B
  REAL(DP), intent(in) :: A(0:nnz-1),B(0:dim-1)
  integer, intent(in) :: JA(0:nnz-1),IA(0:dim)
  REAL(DP), intent(out) :: PR(0:dim-1)

  integer :: ii,jj

  PR(:)=0.0_dp
  do ii=0,dim-1
     do jj=IA(ii), (IA(ii+1)-1)
        PR(ii)=PR(ii)+A(jj)*B(JA(jj))
     end do
  end do
end subroutine prodottoYSMP 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE pot_init(POT3D,which,map,nx,ny,nz,lwork)

IMPLICIT NONE

REAL(DP), INTENT(OUT) :: POT3D(0:lwork-1)

INTEGER, INTENT(IN) :: which(0:nx*ny*nz-1)
INTEGER, INTENT(IN) :: map(0:nx*ny*nz-1)
INTEGER, INTENT(IN) :: nx
INTEGER, INTENT(IN) :: ny
INTEGER, INTENT(IN) :: nz
INTEGER, INTENT(IN) :: lwork
INTEGER :: ii, x_index

 POT3D=0.0
 DO ii=0, nx*ny*nz-1
 x_index=ii/(ny*nz)

IF(which(ii).le.0)THEN

if(chtype.eq.'n')then
 IF(x_index.le.source_len)THEN
    POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))-mus
 ELSEIF(x_index.gt.source_len+2*spacer+gate_len .and. x_index .lt. source_len+2*spacer+gate_len+drain_len)THEN
    POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))-mud
 ELSEIF(x_index.gt.source_len.and.x_index.le.source_len+2*spacer+gate_len)THEN
    POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))-0.1_dp 
 END IF
 IF(x_index.eq.source_len+2*spacer+gate_len+drain_len)THEN
    POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax))+0.0*off_set(imat(x_index/Ndeltax))-mud 
 END IF
elseif(chtype.eq.'p')then

 IF(x_index.le.source_len)THEN
    POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))-mus
 ELSEIF(x_index.gt.source_len+2*spacer+gate_len .and. x_index .lt. source_len+2*spacer+gate_len+drain_len)THEN
    POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))-mud
 ELSEIF(x_index.gt.source_len.and.x_index.le.source_len+2*spacer+gate_len)THEN
    POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))+0.1_dp
 END IF
 IF(x_index.eq.source_len+2*spacer+gate_len+drain_len)THEN
    POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax))+0.0*off_set(imat(x_index/Ndeltax))-mud 
 END IF
elseif(chtype.eq.'t')then

 IF(x_index.le.source_len)THEN
    POT3D(map(ii))=ref_ev(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))-mus
 ELSEIF(x_index.gt.source_len+2*spacer+gate_len .and. x_index .lt. source_len+2*spacer+gate_len+drain_len)THEN
    POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))-mud
 ELSEIF(x_index.gt.source_len.and.x_index.le.source_len+2*spacer+gate_len)THEN
    POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax+1))+0.0*off_set(imat(x_index/Ndeltax+1))
 END IF
 IF(x_index.eq.source_len+2*spacer+gate_len+drain_len)THEN
    POT3D(map(ii))=ref_ec(imat(x_index/Ndeltax))+0.0*off_set(imat(x_index/Ndeltax))-mud 
 END IF
end if

END IF
 END DO


END SUBROUTINE pot_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     

  SUBROUTINE poisson_imref_n(Fn, Ec, rho, lwork)

  IMPLICIT NONE

  REAL(DP), INTENT(INOUT) :: Fn(0:lwork-1)      !Pseudo-potential

  REAL(DP), INTENT(IN)  :: Ec(0:lwork-1) !Potential
  REAL(DP), INTENT(IN)  :: rho(0:lwork-1) !charge
  INTEGER, INTENT(IN)  :: lwork

  REAL(DP) :: N_c, N3D, f, df, TOLERANCE, delta, fn_old
  !REAL(DP), ALLOCATABLE :: Ec(:)
  INTEGER :: ii, nn  

  !ALLOCATE(Ec(0:lwork-1))

  !Ec=-pot+E_GAP/2.0_dp
  
  Fn=3._dp!2.0_dp

!  dd=1.0d-6 
  N_c=6.2d17!6.2d17
  N3D=N_c*(TEMP)**(1.5) 
  write(*,*)'poisson_imref N3D=',N3D

  TOLERANCE=1.0d-15
  nn=0

  DO ii=0, lwork-1 

  f=rho(ii)-N3D*FDP0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP))
!!$  df=-N3D*(FDP0P5((Fn(ii)+dd-Ec(ii))/(BOLTZ*TEMP)) - FDP0P5((Fn(ii)-dd-Ec(ii))/(BOLTZ*TEMP)))/(2.0_dp*dd)
  df=N3D*(FDM0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)
  delta=-f/df
  fn_old=Fn(ii)
  if(rho(ii).eq.0.0_DP)then
  Fn(ii)=100.0_dp
  else
  nn=0
  DO WHILE((ABS(delta).gt.TOLERANCE).and.(nn.lt.400))
  !DO WHILE((ABS(delta).gt.TOLERANCE))
  Fn(ii)=Fn(ii)-delta
  f=rho(ii)-N3D*FDP0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP))
!!$  df=-N3D*(FDP0P5((Fn(ii)+dd-Ec(ii))/(BOLTZ*TEMP)) - FDP0P5((Fn(ii)-dd-Ec(ii))/(BOLTZ*TEMP)))/(2.0_dp*dd)
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

  REAL(DP), INTENT(INOUT) :: Fp(0:lwork-1)      !Pseudo-potential

  REAL(DP), INTENT(IN)  :: Ev(0:lwork-1) !Potential
  REAL(DP), INTENT(IN)  :: rho(0:lwork-1) !charge
  INTEGER, INTENT(IN)  :: lwork

  REAL(DP) :: N_v, N3D, f, df, TOLERANCE, delta, fp_old
  INTEGER :: ii, nn  


  
  Fp=-3._dp!2.0_dp

!  dd=1.0d-6
  N_v=6.5d17!3.5d17
  N3D=N_v*(TEMP)**(1.5) 
  write(*,*)'poisson_imref N3D=',N3D

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
  INTEGER,          INTENT(IN)  :: lwork

  REAL(DP) :: N_c, N3D
  INTEGER :: ii


  N_c=6.2d17!6.2d17
  N3D=N_c*(TEMP)**(1.5) 

  rho=0.0_dp

  DO ii=0, lwork-1
    IF(Fn(ii).ne.100.0_dp)rho(ii)=N3D*FDP0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP))   
  END DO
 
  !DEALLOCATE(Ec)

END SUBROUTINE poisson_charge_from_imref_n

SUBROUTINE poisson_charge_from_imref_p(rho, Ev, Fp, lwork)

  IMPLICIT NONE
  
  REAL(DP), INTENT(OUT) :: rho(0:lwork-1)
    
  REAL(DP), INTENT(IN)  :: Ev(0:lwork-1)
  REAL(DP), INTENT(IN)  :: Fp(0:lwork-1)
  INTEGER,          INTENT(IN)  :: lwork
    
  REAL(DP) :: N_v, N3D
  INTEGER :: ii
    
  N_v=6.5d17!3.5d17
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
  INTEGER,          INTENT(IN)  :: lwork

  REAL(DP) :: N_c, N3D
  INTEGER :: ii
  !REAL(DP), ALLOCATABLE :: Ec(:)
  !REAL(DP) :: xx_mass, yy_mass, zz_mass
  
  N_c=6.2d17!6.2d17
  N3D=N_c*(TEMP)**(1.5) 

  drho=0.0_dp
!  dd=1.0d-6

  DO ii=0, lwork-1
     IF(Fn(ii).ne.100.0_dp)THEN 
!!$  drho(ii)=N3D*( FDP0P5((Fn(ii)+dd-Ec(ii))/(BOLTZ*TEMP)) &
!!$                 - FDP0P5((Fn(ii)-dd-Ec(ii))/(BOLTZ*TEMP)) )/(2.0_dp*dd)
        drho(ii)=N3D*( FDM0P5((Fn(ii)-Ec(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)
        
     END IF
  END DO
  
  !DEALLOCATE(Ec)

END SUBROUTINE poisson_deriv_from_imref_n

SUBROUTINE poisson_deriv_from_imref_p(drho, Ev, Fp, lwork)
    
    IMPLICIT NONE
    
    REAL(DP), INTENT(OUT) :: drho(0:lwork-1)
    
    REAL(DP), INTENT(IN)  :: Ev(0:lwork-1)
    REAL(DP), INTENT(IN)  :: Fp(0:lwork-1)
    INTEGER,          INTENT(IN)  :: lwork
    
    REAL(DP) :: N_v, N3D
    INTEGER :: ii
     
    N_v=6.5d17!3.5d17
    N3D=N_v*(TEMP)**(1.5) 
   
    drho=0.0_dp
    !dd=1.0d-6
    
    DO ii=0, lwork-1
       IF(Fp(ii).ne.100.0_dp)THEN
!!          drho(ii)=N3D*( FDP0P5(-(Fp(ii)+dd-Ev(ii))/(BOLTZ*TEMP)) - FDP0P5(-(Fp(ii)-dd-Ev(ii))/(BOLTZ*TEMP)) )/(2.0_dp*dd)
        drho(ii)=-N3D*( FDM0P5(-(Fp(ii)-Ev(ii))/(BOLTZ*TEMP)) )/(BOLTZ*TEMP)
       END IF
    END DO
    
END SUBROUTINE poisson_deriv_from_imref_p




END MODULE poisson
