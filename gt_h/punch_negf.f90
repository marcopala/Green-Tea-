!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE punch( what )
  !----------------------------------------------------------------------------
  !! This routine is called at the end of the run to save on a file
  !! the information needed for further processing (phonon etc.).
  !
  !! * what = 'all' : write xml data file and, if io_level > -1, charge 
  !!                  density and wavefunctions. For final data.
  !! * what = 'config' : write xml data file and charge density; also,
  !!                     for nks=1, wavefunctions in plain binary format
  !!                     (see why in comments below). For intermediate 
  !!                     or incomplete results
  !! * what = 'config-only' : write xml data file only
  !! * what = 'config-init' : write xml data file only excluding final results
  !!                          (for dry run, can be called at early stages).
  !
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunpun, iunwfc, nwordwfc, diropn,      &
                                   tmp_dir, prefix, restart_dir, xmlfile, &
                                   create_directory
  USE control_flags,        ONLY : io_level, lscf, lxdm
  USE klist,                ONLY : nks
  USE io_files,             ONLY : psfile, pseudo_dir
  USE wrappers,             ONLY : f_copy
  USE spin_orb,             ONLY : lforcet
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : nspin
  USE ions_base,            ONLY : nsp
  USE pw_restart_new,       ONLY : pw_write_schema, pw_write_binaries
  USE qexsd_module,         ONLY : qexsd_reset_steps
  USE io_rho_xml,           ONLY : write_scf
  USE a2F,                  ONLY : la2F, a2Fsave
  USE wavefunctions,        ONLY : evc
  USE xdm_module,           ONLY : write_xdmdat
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  !! see main comment
  !
  ! ... local variables
  !
  LOGICAL :: exst, only_init, wf_collect
  CHARACTER(LEN=320) :: cp_source, cp_dest
  INTEGER            :: cp_status, nt
  !
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing output data file ",A)' ) &
      TRIM ( restart_dir ( ) )
  iunpun = 4
  !
  ! ...New-style I/O with xml schema and (optionally) hdf5 binaries
  !
  ! ... create the main restart directory (if needed)
  !
  CALL create_directory( restart_dir( ) )
  !
  ! ... wf_collect keeps track whether wfcs are written in portable format
  !
  wf_collect = ( TRIM(what) == 'all' )
  only_init  = ( TRIM(what) == 'config-init' )
  CALL pw_write_schema( only_init, wf_collect )
  !
  IF ( TRIM(what) == 'all' .AND. io_level < 0 ) RETURN
  !
  ! ... charge density - also writes rho%ns if lda+U and rho%bec if PAW
  ! ... do not overwrite the scf charge density with a non-scf one
  ! ... (except in the 'force theorem' calculation of MAE where the
  ! ...  charge density differs from the one read from disk)
  !
  IF (TRIM(what) == 'all' .OR. TRIM(what) == 'config' ) THEN
     IF ( lscf .OR. lforcet ) CALL write_scf( rho, nspin )
  ENDIF
  !
  IF (TRIM(what) == 'all') THEN 
     !
     ! ... copy xml file one level up (FIXME: why?),
     ! ... copy pseudopotential files into the .save directory
     !
     IF (ionode) THEN
        !
        cp_source = xmlfile ( )
        cp_dest   = TRIM(tmp_dir)//TRIM(prefix)//'.xml'
        cp_status = f_copy(cp_source, cp_dest)
        !
        DO nt = 1, nsp
           cp_source = TRIM(pseudo_dir)//psfile(nt)
           cp_dest   = TRIM(restart_dir ( ) ) //psfile(nt)
           IF ( TRIM(cp_source) /= TRIM(cp_dest) ) &
                cp_status = f_copy(cp_source, cp_dest)
        ENDDO
        !
        ! write XDM dispersion data (coefficients and vdw radii) to xdm.dat
        IF (lxdm) THEN
           CALL write_xdmdat()
        ENDIF
     ENDIF
     !
     ! ... wavefunctions in "collected" format - also G- and k+G-vectors
     !
     CALL pw_write_binaries( )

     ! ... if allocated, deallocate variables containing info on ionic steps 
     ! 
     CALL qexsd_reset_steps()
     !
  ELSEIF ( TRIM(what) == 'config' .AND.  nks == 1 ) THEN
     !
     ! ... here we are stopping an incomplete calculations - wavefunctions are 
     ! ... stored in buffers and saved when buffers are closed. For 1 k-point 
     ! ... however there is no buffer: wavefunctions must be saved to file here
     !
     IF (io_level < 1) CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
     CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, 1 )
     IF (io_level < 1) CLOSE ( UNIT=iunwfc, STATUS='keep' )
     CALL infomsg('punch','wavefunctions written to file')
     !
  ENDIF
  !
  ! ... FIXME: for electron-phonon calculations - data should be read from xml file!
  ! 
  IF ( la2F ) CALL a2Fsave()
  !


  CALL simple_output ( )
!!!CALL simple_diag (  )


  RETURN
  !
END SUBROUTINE punch




!
! Copyright (C) 2019-2020 Paolo Giannozzi
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE simple_output (  )
  !----------------------------------------------------------------------------
  !
  ! Not-so-smart but easy-to-read output file for simple cases (e.g. Si)
  !
  USE kinds,        ONLY : dp
  USE io_global,    ONLY : stdout
  USE io_files,     ONLY : iunpun, iunwfc, tmp_dir, prefix, postfix, nwordwfc
  USE cell_base,    ONLY : at, bg, alat, tpiba
  USE ions_base,    ONLy : nat, nsp, ityp, atm, tau
  USE gvect,        ONLY : ngm, mill
  USE scf,          ONLY : vrs
  USE fft_rho,      ONLY : rho_r2g
  USE lsda_mod,     ONLY : nspin
  USE klist,        ONLY : nks, xk, ngk, igk_k
  USE uspp,         ONLY : okvan, nkb, vkb, dvan, dvan_so
  USE uspp_param,   ONLY : nh
  USE fft_base,     ONLY : dffts
  USE buffers,      ONLY : get_buffer
  USE spin_orb,     ONLY : lspinorb, domag
  USE wvfct,        ONLY : nbnd, et
  USE wavefunctions,ONLY : evc
  !

  IMPLICIT NONE
  !
  integer,parameter :: k15 = selected_int_kind(15)
  COMPLEX(dp), ALLOCATABLE :: vaux(:,:)
  CHARACTER(LEN=256) :: fileout
  INTEGER(k15) :: ig, ikb
  INTEGER is, ik, ibnd, na, nt
  !
  fileout = TRIM(tmp_dir) // TRIM( prefix ) // postfix // 'output.dat'
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing simple output data file ",A)' ) &
       fileout
  iunpun = 4
  OPEN ( UNIT = iunpun, FORM = 'formatted', STATUS = 'unknown', &
       FILE = fileout )
  WRITE(iunpun,'("# Primitive lattice vectors a_1, a_2, a_3 (a.u.)")')
  WRITE(iunpun,*) alat*at(:,1), alat*at(:,2), alat*at(:,3)
  WRITE(iunpun,'("# Reciprocal lattice vectors b_1, b_2, b_3 (a.u.)")')
  WRITE(iunpun,*) tpiba*bg(:,1), tpiba*bg(:,2), tpiba*bg(:,3)
  WRITE(iunpun,'("# Number of types of atom")')
  WRITE(iunpun,*) nsp
  WRITE(iunpun,'("# Number of atoms")')
  WRITE(iunpun,*) nat
  WRITE(iunpun,'("# Atomic species and positions (x, y, z, in a.u.)")')
  DO na =1, nat
     nt = ityp(na)
     WRITE(iunpun,*) ityp(na)
     WRITE(iunpun,'(a,3e25.15)') atm(nt), alat*tau(:,na)

  END DO
  WRITE(iunpun,'("# number of G-vectors")')
  WRITE(iunpun,*) ngm, shape(mill)
  WRITE(iunpun,'("# Miller indices: G=i_1*b_1+i_2*b_2+i_3*b_3")')
  do ig=1,ngm
    WRITE(iunpun,'(3i14)') mill(:,ig)
!     WRITE(iunpun,*) mill(1,ig),mill(2,ig),mill(3,ig)
  end do
  WRITE(iunpun,'("# number of spin components")')
  WRITE(iunpun,*) nspin
  IF ( nspin == 4 ) THEN
     WRITE(iunpun,'("# magnetic, spin-orbit?")')
     WRITE(iunpun,*) domag, lspinorb
  END IF
  ALLOCATE (vaux(ngm,nspin))
  CALL rho_r2g (dffts, vrs, vaux )
  WRITE(iunpun,'("# Local potential V(G) (one column per spin component)")')
  DO is=1,nspin
     WRITE(iunpun,'("# spin component n.",i4)')
     do ig=1,ngm
        WRITE(iunpun,'(2e25.15)') vaux(ig,is)
     end do
  END DO
  DEALLOCATE (vaux)

  WRITE(iunpun,'("# US PP?")')
  WRITE(iunpun,*) okvan

  WRITE(iunpun,'("# Nomber of projectors")')
  WRITE(iunpun,*) nh(1:nsp)

  WRITE(iunpun,'("# Nonlocal PP coefficients Dlm")')
  DO nt = 1, nsp
     WRITE(iunpun,'("# Atom  type, number of projectors")')
     WRITE(iunpun,*) nt, nh(nt)
     IF ( nspin /= 4 ) THEN
        WRITE(iunpun,*) dvan(1:nh(nt),1:nh(nt),nt)
     ELSE
        DO is = 1, nspin
           WRITE(iunpun,'("# spin component n.",i4)') is
           WRITE(iunpun,*) dvan_so(1:nh(nt),1:nh(nt),is,nt)
        END DO
     END IF
  END DO
  WRITE(iunpun,'("# number of beta functions")')
  WRITE(iunpun,*) nkb

  WRITE(iunpun,'("# number of k-points")')
  WRITE(iunpun,*) nks

  WRITE(iunpun,'("# number of NPW for k")')
  WRITE(iunpun,*) ngk(1:nks)


  DO ik=1,nks
     WRITE(iunpun,'("# k-point n.",i4)') ik
     write(*,*) 'ik',ik,ngk(ik)
     !write(*,*) 'shape vkb',shape(vkb)
     !write(*,*) 'shape evc',shape(evc)
     WRITE(iunpun,*) tpiba*xk(:,ik)
     WRITE(iunpun,'("# number of plane waves")')
     WRITE(iunpun,*) ngk(ik)
     WRITE(iunpun,'("# index of k+G: (k+G)_i = k + G_index(i)")')
     do ig=1,ngk(ik)
        WRITE(iunpun,'(i8)') igk_k(ig,ik)
     end do
     CALL init_us_2 (ngk(ik), igk_k(1,ik), xk (1, ik), vkb)
     DO ikb=1,nkb
        WRITE(iunpun,'("# beta function n.",i4)') ikb
        WRITE(iunpun,'(2e25.15)') vkb(1:ngk(ik),ikb)
     END DO
     CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
     WRITE(iunpun,'("# number of bands")')
     WRITE(iunpun,*) nbnd
     DO ibnd=1,nbnd
        WRITE(iunpun,'("# band n.",i4)') ibnd
        WRITE(iunpun,'("# eigenvalue (eV):",e25.15)') et(ibnd,ik)*13.6058
        IF ( nspin /= 4 ) THEN
           WRITE(iunpun,'(2e25.15)') evc(1:ngk(ik),ibnd)
        ELSE
           WRITE(iunpun,'(2e25.15)') evc(1:ngk(ik),ibnd)
           WRITE(iunpun,'(2e25.15)') evc(1+maxval(ngk):ngk(ik)+maxval(ngk),ibnd)
        END IF
     END DO
  END DO
  WRITE(iunpun,'("# end of file")')
  !
  CLOSE (unit=iunpun, STATUS='keep')
  !
END SUBROUTINE simple_output
!----------------------------------------------------------------------------
!
SUBROUTINE simple_diag (  )
  !----------------------------------------------------------------------------
  !
  ! Check: read the simple data file, re-diagonalize the matrix
  !
  USE kinds,        ONLY : dp
  USE io_global,    ONLY : stdout
  USE io_files,     ONLY : iunpun, iunwfc, tmp_dir, prefix, postfix, nwordwfc
  USE cell_base,    ONLY : at, bg, alat, tpiba
  USE ions_base,    ONLy : nat, nsp, ityp, atm, tau
  USE gvect,        ONLY : ngm, mill
  USE scf,          ONLY : vrs
  USE fft_rho,      ONLY : rho_r2g
  USE lsda_mod,     ONLY : nspin
  USE spin_orb,     ONLY : lspinorb, domag
  USE klist,        ONLY : nks, xk, ngk, igk_k
  USE uspp,         ONLY : okvan, nkb, vkb, dvan, dvan_so
  USE uspp_param,   ONLY : nh
  USE fft_base,     ONLY : dfftp
  USE fft_rho,      ONLY : rho_r2g
  USE buffers,      ONLY : get_buffer
  USE wvfct,        ONLY : nbnd, et
  USE wavefunctions,ONLY : evc
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=80) :: line
  COMPLEX(dp), ALLOCATABLE :: vaux(:,:), h(:,:), v(:,:)
  CHARACTER(LEN=256) :: fileout
  INTEGER :: npw, ig, is, ik, ikb, ibnd, na, nt, nt_, i, j, ii, jj, ij, &
       ih, jh, i1, i2, i3, ipol
  INTEGER, ALLOCATABLE :: limm(:,:,:)
  REAL(dp) :: g(3)
  REAL(dp), ALLOCATABLE :: et_(:)
  LOGICAL :: debug = .FALSE.
  COMPLEX(dp) :: auu
  !
  fileout = TRIM(tmp_dir) // TRIM( prefix ) // postfix // 'output.dat'
  WRITE( UNIT = stdout, FMT = '(/,5X,"Checking simple output data file ",A)' ) &
       fileout
  iunpun = 4
!write(*,*)'ciao1'
  OPEN ( UNIT = iunpun, FORM = 'formatted', STATUS = 'old', FILE = fileout )
  READ(iunpun,'(a)')  line
  READ(iunpun,*) at(:,1), at(:,2), at(:,3)
  IF ( debug) WRITE(stdout,*) trim(line), at(:,1), at(:,2), at(:,3)
  READ(iunpun,'(a)')  line
  READ(iunpun,*) bg(:,1), bg(:,2), bg(:,3)
  IF ( debug) WRITE(stdout,*) trim(line), bg(:,1), bg(:,2), bg(:,3)
  READ(iunpun,'(a)')  line
  READ(iunpun,*) nsp
  IF ( debug) WRITE(stdout,*) trim(line), nsp
  READ(iunpun,'(a)')  line
  READ(iunpun,*) nat
  IF ( debug) WRITE(stdout,*) trim(line), nat
  READ(iunpun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  DO na =1, nat
     READ(iunpun,*) nt
     READ(iunpun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
  END DO

!write(*,*)'ciao2'

  READ(iunpun,'(a)')  line
  READ(iunpun,*) ngm
  IF ( debug) WRITE(stdout,*) trim(line), ngm
  READ(iunpun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  READ(iunpun,'(3i8)') (mill(:,ig), ig=1,ngm)
  ! if i1=mill(1,ig), i2=mill(2,ig), i3=mill(3,ig), then:
  ! limm(i1,i2,i3) = ig
  i1 = MAXVAL( ABS(mill(1,:)) )
  i2 = MAXVAL( ABS(mill(2,:)) )
  i3 = MAXVAL( ABS(mill(3,:)) )
  ALLOCATE (limm(-i1:i1,-i2:i2,-i3:i3))
  limm = 0
  DO ig=1,ngm
     limm( mill(1,ig), mill(2,ig), mill(3,ig) ) = ig
  ENDDO
  !
  READ(iunpun,'(a)')  line
  READ(iunpun,*) nspin
  IF ( debug) WRITE(stdout,*) trim(line), nspin
  IF ( nspin <= 2 ) THEN
     npol = 1
  ELSE
     READ(iunpun,'(a)')  line
     READ(iunpun,*) domag, lspinorb
     IF ( domag .OR. .NOT.lspinorb ) &
          CALL errore ('simple_diag','spin-orbit with no magnetization only',1)
     npol = 2
  ENDIF
!write(*,*)'ciao3'

  ALLOCATE (vaux(ngm,nspin))
  READ(iunpun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  DO is=1,nspin
     READ(iunpun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
     READ(iunpun,*) (vaux(ig,is), ig=1,ngm)
  END DO
  READ(iunpun,'(a)')  line
  READ(iunpun,*)  okvan

  READ(iunpun,'(a)')  line
  READ(iunpun,*) nh(1:nsp)

!write(*,*)'ciao3'

  IF ( okvan ) CALL errore ('simple_diag','US PP not implemented',1)
  READ(iunpun,'(a)')  line
  IF ( debug) WRITE(stdout,'(a)') line
  DO nt = 1, nsp
     READ(iunpun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
     READ(iunpun,*) nt_, nh(nt)
     IF ( nspin /= 4 ) THEN
        READ(iunpun,*) dvan(1:nh(nt),1:nh(nt),nt)
     ELSE
        DO is=1,nspin
           READ(iunpun,'(a)')  line
           IF ( debug) WRITE(stdout,'(a)') line
           READ(iunpun,*) dvan_so(1:nh(nt),1:nh(nt),is,nt)
        END DO
     END IF
  END DO
     do j=1,Nsp
      do is=1,nspin
 write(700+j,*) Dvan_so(1:Nh(j),1:Nh(j),is,j)
 write(700+j,*)
end do
end do
!stop
  READ(iunpun,'(a)')  line
  READ(iunpun,*)  nkb
  IF ( debug) WRITE(stdout,*) trim(line), nkb
  !
  READ(iunpun,'(a)')  line
  READ(iunpun,*) nks

  READ(iunpun,'(a)')  line
  READ(iunpun,*) ngk(1:nks)

  IF ( debug) WRITE(stdout,*) trim(line), nks
  DO ik=1,nks
     READ(iunpun,'(a)')  line
     READ(iunpun,*) xk(:,ik)
     IF ( debug) WRITE(stdout,*) trim(line), xk(:,ik)
     READ(iunpun,'(a)')  line
     READ(iunpun,*) ngk(ik)
     IF ( debug) WRITE(stdout,*) trim(line), ngk(ik)
     READ(iunpun,'(a)')  line
     IF ( debug) WRITE(stdout,'(a)') line
     READ(iunpun,'(i8)') (igk_k(ig,ik), ig=1,ngk(ik))
     DO ikb=1,nkb
        READ(iunpun,'(a)')  line
        IF ( debug) WRITE(stdout,'(a)') line
        READ(iunpun,'(2e25.15)') vkb(1:ngk(ik),ikb)
     END DO
     CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
     READ(iunpun,'(a)')  line
     READ(iunpun,*)  nbnd
     IF ( debug) WRITE(stdout,*) trim(line), nbnd
     DO ibnd=1,nbnd
        READ(iunpun,'(a)')  line
        IF ( debug) WRITE(stdout,'(a)') line
        READ(iunpun,'(a)')  line
        IF ( debug) WRITE(stdout,'(a)') line
        READ(iunpun,'(2e25.15)') evc(:,ibnd)
     END DO
     WRITE(stdout,'("# Data read for k-point #",i4,", diagonalizing...")') ik
     npw = ngk(ik)
     ALLOCATE ( h(npol*npw,npol*npw), v(npol*npw,npol*npw), et_(npol*npw) )
     h(:,:) = (0.0_dp,0.0_dp)
     !write(*,*)'pre-pre-end',ik
     DO j=1,npw
        !
        ! kinetic energy
        !
        g(:) = mill(1,igk_k(j,ik))*bg(:,1) + &
               mill(2,igk_k(j,ik))*bg(:,2) + &
               mill(3,igk_k(j,ik))*bg(:,3)
        h(j,j)= ( xk(1,ik)+g(1) )**2 + &
                ( xk(2,ik)+g(2) )**2 + &
                ( xk(3,ik)+g(3) )**2
        IF ( npol == 2 ) h(npw+j,npw+j) = h(j,j)
        write(2000,*)j,abs(h(j,j))*13.6056929972990
        !
        ! nonlocal potential
        !
        ikb=0
        DO nt=1,nsp
           DO na=1,nat
              IF ( nt == ityp(na) ) THEN
                 DO ih=1,nh(nt)
                    DO jh=1,nh(nt)
                       IF ( nspin /= 4 ) THEN
                          DO i=j,npw
                             h(i,j) = h(i,j) + vkb(i,ikb+ih) * &
                                  dvan(ih,jh,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                       ELSE
                          DO i=j,npw
                             h(i,j) = h(i,j) + vkb(i,ikb+ih) * &
                                  dvan_so(ih,jh,1,nt) *DCONJG(vkb(j,ikb+jh))
                             h(i+npw,j+npw) = h(i+npw,j+npw) + vkb(i,ikb+ih)* &
                                  dvan_so(ih,jh,4,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                          DO i=1,npw
                             h(i,j+npw) = h(i,j+npw) + vkb(i,ikb+ih) * &
                                  dvan_so(ih,jh,2,nt) *DCONJG(vkb(j,ikb+jh))
                             h(i+npw,j) = h(i+npw,j) + vkb(i,ikb+ih) * &
                                  dvan_so(ih,jh,3,nt) *DCONJG(vkb(j,ikb+jh))
                          END DO
                       END IF
                    END DO
                 END DO
                 ikb = ikb + nh(nt)
              END IF
           END DO
        END DO
        DO i=1,10!npw
              write(3100,*)i,j,abs(h(i,j))*13.6056929972990
              write(4100,*)i,j,abs(h(i,j+npw))*13.6056929972990
              write(5100,*)i,j,abs(h(i+npw,j))*13.6056929972990
              write(6100,*)i,j,abs(h(i+npw,j+npw))*13.6056929972990
           END DO
           write(3100,*)
           write(4100,*)
           write(5100,*)
           write(6100,*)

           !
        ! local potential
        !
        DO i=j,npw
           i1 = mill(1,igk_k(i,ik)) - mill(1,igk_k(j,ik))
           i2 = mill(2,igk_k(i,ik)) - mill(2,igk_k(j,ik))
           i3 = mill(3,igk_k(i,ik)) - mill(3,igk_k(j,ik))
           IF ( ABS(i1) > SIZE(limm,1) .OR. &
                ABS(i2) > SIZE(limm,2) .OR. &
                ABS(i3) > SIZE(limm,3) ) &
                CALL errore ('simple_diag','internal error (1)',i)
           ij = limm ( i1,i2,i3 )
           IF ( ij <= 0 .OR. ij > ngm ) &
                CALL errore ('simple_diag','internal error (2)',i)
           DO ipol = 1, npol
              ii = (ipol-1)*npw + i
              jj = (ipol-1)*npw + j
              h(ii,jj) = h(ii,jj) + vaux(ij,1)
              IF ( i > j ) h(jj,ii) =DCONJG(h(ii,jj))
           END DO
        END DO
     END DO

!write(*,*)'pre-end'

     DO j=1,npw
        DO i=1,10!npw
           write(310+ik-1,*)i,j,abs(h(i,j))*13.6056929972990
           write(410+ik-1,*)i,j,abs(h(i,j+npw))*13.6056929972990
           write(510+ik-1,*)i,j,abs(h(i+npw,j))*13.6056929972990
           write(610+ik-1,*)i,j,abs(h(i+npw,j+npw))*13.6056929972990

        END DO
        write(310+ik-1,*)
        write(410+ik-1,*)
        write(510+ik-1,*)
        write(610+ik-1,*)
        END DO

        do ih=1,5!nbnd
           do is=1,5!nbnd
              auu=0.0
              DO j=1,npw*2
                 do i=1,npw*2
                    auu=auu+h(i,j)*evc(j,is)*dconjg(evc(i,ih))
                 end do
              end DO
              write(1000+ik,*)ih,is,abs(auu)*13.6056929972990
              if(ih==is)write(1200+ik,*)ih,abs(auu)*13.6056929972990
           end do
           write(1000+ik,*)
        end DO

             ! write(*,*)'end'


     !
     CALL cdiagh ( npol*npw, h, npol*npw, et_, v)
     do ibnd=1,nbnd
        write(1100+ik-1,*)ibnd,et_(ibnd)*13.6056929972990
     end do
     WRITE(stdout,'(4f12.6)') (et_(ibnd)*13.6056929972990, ibnd=1,nbnd)
     DEALLOCATE ( et_, v, h )
  END DO
  DEALLOCATE (limm, vaux)
  READ(iunpun,'(a)')  line
  WRITE(stdout,'(a)') line
  !
  CLOSE (unit=iunpun, STATUS='keep')
  !
END SUBROUTINE simple_diag






