module static


integer,parameter :: k15 = selected_int_kind(15)
integer, parameter :: dp = selected_real_kind(15,307)
REAL(DP), PARAMETER :: DIEL_0=8.854187817d-14  ! F/cm
REAL(DP), PARAMETER :: DIEL_METAL=0.0_dp
real(dp), parameter :: pi=abs(acos(-1.0_dp))
real(dp), parameter :: hbar=1.05457180d-34 ! J*s
real(dp), parameter :: ELCH=1.60217662d-19 ! C
real(dp), parameter :: me=9.10938356d-31  ! Kg
REAL(DP), PARAMETER :: BOLTZ=8.61734d-05  !eV K-1
real(dp), parameter :: bohr=0.529177208590d-8 !0.529177d-08   ! cm
real(dp), parameter :: hbareV=hbar/ELCH    ! eV*s
real(dp), parameter :: m0=me/elch*1.0d-04    ! Kg/C times 1E-4 to account for cm**2
real(dp), parameter :: ryd=m0/2.0_dp*(ELCH**2/diel_0/4.0_dp/pi/hbar)**2

real(dp), parameter :: dalpha=1.0_dp
real(dp), parameter :: dbeta=0.0_dp
complex(dp), parameter :: alpha=cmplx(1.0_dp,0.0_dp, kind=dp)
complex(dp), parameter :: beta=cmplx(0.0_dp,0.0_dp, kind=dp)
complex(dp), parameter :: im=cmplx(0.0_dp,1.0_dp, kind=dp)

end module static
