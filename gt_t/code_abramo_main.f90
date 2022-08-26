program main
  implicit none


  integer,     parameter :: dp = selected_real_kind(15,307)
  complex(dp), parameter :: alpha=cmplx(1.0_dp,0.0_dp, kind=dp)
  complex(dp), parameter :: beta=cmplx(0.0_dp,0.0_dp, kind=dp)
  
  integer :: i,j,l,ff,nmax,ncx_d
  REAL(dp)    :: E,mul,mur,tr,tre,eta
  REAL(dp), allocatable    :: cur(:)

  complex(dp), allocatable :: sigma_lesser_ph(:,:),sigma_greater_ph(:,:),sigma_r_ph(:,:)
  COMPLEX(dp), allocatable :: lcur(:,:,:),ldos(:,:,:),ndens(:,:,:),pdens(:,:,:)
  COMPLEX(dp), allocatable :: Hii(:,:,:),H10(:,:)


  open(unit=14,file='inoutput_rgf.dat',status='unknown')
  read(14,*) 
  read(14,*) eta
  read(14,*) 
  read(14,*) ncx_d
  read(14,*) 
  read(14,*) nmax

  allocate(Hii(nmax,nmax,ncx_d))
  allocate(H10(nmax,nmax))
  allocate(sigma_lesser_ph(nmax,Ncx_d),sigma_greater_ph(nmax,Ncx_d),sigma_r_ph(nmax,Ncx_d))
  allocate(lcur(nmax,nmax,ncx_d-1),ldos(nmax,nmax,ncx_d),ndens(nmax,nmax,ncx_d),pdens(nmax,nmax,ncx_d))
  allocate(cur(Ncx_d-1))
  
  
  read(14,*) 
  read(14,*) ff
  read(14,*) 
  read(14,*) E
  read(14,*) 
  read(14,*) mul 
  read(14,*) 
  read(14,*) mur
  read(14,*) 
  do i=1,nmax
     do j=1,nmax
        do l=1,ncx_d
           read(14,*) Hii(i,j,l)
        end do
     end do
  end do
  read(14,*) 
  do i=1,nmax
     do j=1,nmax
        read(14,*) H10(i,j) !TL(iyz,ihet(1))%H(i,j)
     end do
  end do
  read(14,*) 
  do j=1,nmax
     do l=1,ncx_d
        read(14,*) sigma_lesser_ph(j,l)
     end do
  end do
  read(14,*) 
  do j=1,nmax
     do l=1,ncx_d
        read(14,*) sigma_greater_ph(j,l)
     end do
  end do
  read(14,*) 
  do j=1,nmax
     do l=1,ncx_d
        read(14,*) sigma_r_ph(j,l)
     end do
  end do
  read(14,*) 
  do i=1,nmax
     do j=1,nmax
        do l=1,ncx_d
           read(14,*) ndens(i,j,l)
        end do
     end do
  end do
  read(14,*) 
  do i=1,nmax
     do j=1,nmax
        do l=1,ncx_d
           read(14,*) pdens(i,j,l)
        end do
     end do
  end do
  read(14,*) 
  do i=1,nmax
     do j=1,nmax
        do l=1,ncx_d
  read(14,*) ldos(i,j,l)
        end do
     end do
  end do
  read(14,*) 
  read(14,*) tr
  read(14,*) 
  read(14,*) tre
  read(14,*) 
  do l=1,ncx_d-1
     read(14,*) cur(l)
  end do
  read(14,*) 
  do i=1,nmax
     do j=1,nmax
        do l=1,ncx_d-1
           read(14,*) lcur(i,j,l)
        end do
     end do
  end do
  close(14)


  write(*,*)'begin'
  
  call RGF(ncx_d,nmax,ff,E,mul,mur,H10,Hii,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,lcur) 

  write(*,*)'end'

  
end program main


subroutine RGF(ncx_d,nmax,ff,E,mul,mur,H10,Hii,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,lcur) 

  implicit none
  
  integer,     parameter :: dp = selected_real_kind(15,307)
  complex(dp), parameter :: alpha=cmplx(1.0_dp,0.0_dp, kind=dp)
  complex(dp), parameter :: beta=cmplx(0.0_dp,0.0_dp, kind=dp)
  REAL(DP),    PARAMETER :: BOLTZ=8.61734d-05  !eV K-1
  REAL(DP),    PARAMETER :: TEMP=300.0  !K

  INTEGER     :: i,j,l,nmax,ff,ncx_d,nm(ncx_d)
  REAL(dp)    :: E,mul,mur,tr,tre,traccia,ferm
  REAL(dp)    :: cur(Ncx_d-1)
  COMPLEX(dp) :: sigma_lesser_ph(nmax,Ncx_d),sigma_greater_ph(nmax,Ncx_d),sigma_r_ph(nmax,Ncx_d)
  COMPLEX(dp) :: sig(nmax,nmax),sigmal(nmax,nmax),sigmar(nmax,nmax)
  COMPLEX(dp) :: Gn(nmax,nmax),Gp(nmax,nmax),G00(nmax,nmax),GN0(nmax,nmax)
  COMPLEX(dp) :: Gl(nmax,nmax,Ncx_d),Gln(nmax,nmax,Ncx_d),Glp(nmax,nmax,Ncx_d)
  COMPLEX(dp) :: lcur(nmax,nmax,ncx_d-1),ldos(nmax,nmax,ncx_d),ndens(nmax,nmax,ncx_d),pdens(nmax,nmax,ncx_d)
  COMPLEX(dp) :: Hii(nmax,nmax,ncx_d),H00(nmax,nmax),H10(nmax,nmax)
  COMPLEX(dp) :: A(nmax,nmax),B(nmax,nmax),C(nmax,nmax),Id(nmax,nmax)
  COMPLEX(dp) :: z
  
  nm=nmax
  
  
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

  !H10(1:nm(l),1:nm(l)) = TL(iyz,ihet(l))%H(1:NM(l),1:NM(l)) !!!! H_{1,0}
  call sancho(nm(l),E,H00(1:nm(l),1:nm(l)),transpose(dconjg(H10(1:nm(l),1:nm(l)))),G00(1:nm(l),1:nm(l))) 
  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           ff=l
           write(*,*)'NaN warning! Pb w lft snch. Eliminating E =',E
           G00=0.0_dp
           exit
        end if
     end do
  end do
  if( traccia(dimag(- G00(:,:))) < 0.0d0 )then 
     ff=l
     write(*,*)'pb w snch l',E,traccia(dimag(- G00(:,:)))
  end if
  
  !H10(1:nm(1),1:nm(1))=TL(iyz,ihet(l))%H(1:NM(1),1:NM(1))  !!!! H_{1,0}
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
     !H10(1:nm(l),1:nm(l-1))=TL(iyz,ihet(l))%H(1:NM(l),1:NM(l-1)) !!!! H_{l,l-1}
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

!  H10(1:nm(l),1:nm(l))=TL(iyz,ihet(l+1))%H(1:NM(l),1:NM(l))   !!! H(N+1,N)
  call sancho(nm(l),E,H00(1:nm(l),1:nm(l)),H10(1:nm(l),1:nm(l)),G00(1:nm(l),1:nm(l)))

  do i=1,nm(l)
     do j=1,nm(l)
        if ( G00(i,j) /= G00(i,j) )then
           ff=l
           write(*,*)'NaN warning! Pb w rht snch. Eliminating E =',E
           G00=0.0_dp
           exit
        end if
     end do
  end do
  
  if( traccia(aimag(- G00(1:nm(l),1:nm(l)))) < 0.0d0 )then 
     ff=l
     write(*,*)'pb w snch r',E,traccia(aimag(- G00(1:nm(l),1:nm(l))))
  end if

  call zgemm('c','n',nm(l),nm(l),nm(l),alpha,H10(1:nm(l),1:nm(l)),nm(l),G00(1:nm(l),1:nm(l)),nm(l),beta,A(1:nm(l),1:nm(l)),nm(l)) 
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),H10(1:nm(l),1:nm(l)),nm(l),beta,sigmar(1:nm(l),1:nm(l)),nm(l)) 

!  H10(1:nm(l),1:nm(l-1))=TL(iyz,ihet(l))%H(1:NM(l),1:NM(l-1)) !!! H(N,N-1)
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
   
     if( abs(traccia(dimag( pdens(:,:,l) - ndens(:,:,l) - 2.0_dp*ldos(:,:,l)  ))) > 1.0d-4 )then
        ff=l  
     end if

    do l=Ncx_d-1,1,-1

!     H10(1:nm(l+1),1:nm(l))=TL(iyz,ihet(l+1))%H(1:NM(l+1),1:NM(l)) !!! H(l+1,l)
     
     call zgemm('n','c',nm(l+1),nm(l),nm(l),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),Gl(1:nm(l),1:nm(l),l),nm(l),beta, B(1:nm(l+1),1:nm(l)), nm(l+1)) 
     !if(chtype == 'p')then
     !   call zgemm('n','n',nm(l+1),nm(l),nm(l+1),alpha,Gp(1:nm(l+1),1:nm(l+1)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta, C(1:nm(l+1),1:nm(l)), nm(l+1))
     !else
        call zgemm('n','n',nm(l+1),nm(l),nm(l+1),alpha,Gn(1:nm(l+1),1:nm(l+1)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta, C(1:nm(l+1),1:nm(l)), nm(l+1))
     !end if
     !if(chtype == 'p')then
     !   call zgemm('n','n',nm(l+1),nm(l),nm(l),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),Glp(1:nm(l),1:nm(l),l),nm(l),beta, B(1:nm(l+1),1:nm(l)), nm(l+1))
     !else
        call zgemm('n','n',nm(l+1),nm(l),nm(l),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),Gln(1:nm(l),1:nm(l),l),nm(l),beta, B(1:nm(l+1),1:nm(l)), nm(l+1))
     !end if
     call zgemm('n','n',nm(l+1),nm(l),nm(l+1),alpha,G00(1:nm(l+1),1:nm(l+1)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta, A(1:nm(l+1),1:nm(l)), nm(l+1))
     B(1:nm(l+1),1:nm(l))=C(1:nm(l+1),1:nm(l))+A(1:nm(l+1),1:nm(l))
     call zgemm('c','n',nm(l),nm(l),nm(l+1),alpha,H10(1:nm(l+1),1:nm(l)),nm(l+1),B(1:nm(l+1),1:nm(l)),nm(l+1),beta,A(1:nm(l),1:nm(l)),nm(l))      !!! G<_i+1,i
     cur(l)=2.0_dp*traccia(dble(A(1:nm(l),1:nm(l))))
     lcur(1:nm(l),1:nm(l),l)=2.0_dp*(A(1:nm(l),1:nm(l)))
       

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
  
     
     if( abs(traccia(dimag( pdens(:,:,l) - ndens(:,:,l) - 2.0_dp*ldos(:,:,l)  ))) > 1.0d-4 )then
        ff=l
     end if

  enddo
  
  l=1
  A(1:nm(l),1:nm(l))=-(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm((E-mul)/(BOLTZ*TEMP))
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gp(1:nm(l),1:nm(l)),nm(l),beta,B(1:nm(l),1:nm(l)),nm(l))
  A(1:nm(l),1:nm(l))=(sigmal(1:nm(l),1:nm(l))-transpose(dconjg(sigmal(1:nm(l),1:nm(l)))))*ferm(-(E-mul)/(BOLTZ*TEMP))
  call zgemm('n','n',nm(l),nm(l),nm(l),alpha,A(1:nm(l),1:nm(l)),nm(l),Gn(1:nm(l),1:nm(l)),nm(l),beta,C(1:nm(l),1:nm(l)),nm(l))
  tre=-traccia(dble(B(1:nm(l),1:nm(l))-C(1:nm(l),1:nm(l))))

 if ( ff /= 0 ) then
     tr=0.0_dp
     tre=0.0_dp
     ndens=0.0_dp
     pdens=0.0_dp
     ldos=0.0_dp
     cur=0.0_dp
  end if
elseif ( ff /= 0 ) then
     !write(*,*)'ignoring E =',E
     tr=0.0_dp
     tre=0.0_dp
     ndens=0.0_dp
     pdens=0.0_dp
     ldos=0.0_dp
     cur=0.0_dp
  end if

end subroutine RGF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function traccia(A)
  integer,     parameter :: dp = selected_real_kind(15,307)
  real(dp) :: A(:,:)
  real(dp) :: traccia
  integer :: i,j

  traccia=0.0_dp
  do i=1,size(a,1)
     traccia=traccia+A(i,i)
  end do
!  traccia=dble(sum( a(:,:), &
!       reshape( (/((i==j,i=1,size(a,1)),j=1,size(a,2))/), &
!       (/size(a,1),size(a,2)/) ) ) ) 

end function traccia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function ferm(a)

  integer,     parameter :: dp = selected_real_kind(15,307)

  real(dp) a
  ferm = (1.0_dp+dexp(a))**(-1.0_dp)
End Function ferm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Sancho-Rubio
subroutine sancho(nm,E,H00,H10,G00)
  implicit none

  integer,     parameter :: dp = selected_real_kind(15,307)

  complex(dp), parameter :: alpha=cmplx(1.0_dp,0.0_dp, kind=dp)
  complex(dp), parameter :: beta=cmplx(0.0_dp,0.0_dp, kind=dp)

  integer, intent(IN) :: nm
  integer :: i
  integer :: nmax=500
  COMPLEX(DP) :: z
  real(dp), intent(IN) :: E
  REAL(DP) :: TOL=1.0D-25, error

  COMPLEX(DP), INTENT(IN) :: H00(nm,nm), H10(nm,nm)
  COMPLEX(DP), INTENT(OUT) :: G00(nm,nm)

  COMPLEX(DP), ALLOCATABLE :: A(:,:), B(:,:), C(:,:)
  COMPLEX(DP), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)
  
  REAL(DP),    PARAMETER :: eta=1.0d-6



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


subroutine invert(A,n,np)
  implicit none
      
  integer,     parameter :: dp = selected_real_kind(15,307)

  integer :: info,n,np      
  integer, dimension(np) :: ipiv      
  complex(dp) :: A(n,n)    
  complex(dp) :: work(np*np)   

  call zgetrf(n,n,A,n,ipiv,info)
  call zgetri(n,A,n,ipiv,work,np*np,info)

end subroutine invert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

