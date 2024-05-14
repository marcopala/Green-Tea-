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

module static


  ! The following units are used:
  ! eV for energies
  ! cm for lengths 
  ! A for currents
  ! C for charges

  integer,     parameter :: k15 = selected_int_kind(15)
  integer,     parameter :: dp = selected_real_kind(15,307)
  integer,     parameter :: qp = selected_real_kind(33, 4931)
  REAL(DP),    PARAMETER :: DIEL_0=8.854187817d-14  ! F/cm
  REAL(DP),    PARAMETER :: DIEL_METAL=0.0_dp
  real(dp),    parameter :: pi=abs(acos(-1.0_dp))
  real(dp),    parameter :: hbar=1.05457180d-34 ! J*s
  real(dp),    parameter :: ELCH=1.60217662d-19 ! C
  real(dp),    parameter :: me=9.10938356d-31  ! Kg
  REAL(DP),    PARAMETER :: BOLTZ=8.61734d-05  !eV K-1
  real(dp),    parameter :: bohr=0.5291772109030d-08    ! cm
  real(dp),    parameter :: hbareV=hbar/ELCH    ! eV*s
  
  real(dp),    parameter :: m0=me/elch*1.0d-04    ! Kg/C times 1E-4 to account for cm**2 when used to discretize the kinetic operator in real space
  real(dp),    parameter :: ryd=m0/2.0_dp*(ELCH**2/diel_0/4.0_dp/pi/hbar)**2
  real(dp),    parameter :: dalpha=1.0_dp
  real(dp),    parameter :: dbeta=0.0_dp
  complex(dp), parameter :: alpha=dcmplx(1.0_dp,0.0_dp)
  complex(dp), parameter :: beta=dcmplx(0.0_dp,0.0_dp)
  complex(dp), parameter :: im=dcmplx(0.0_dp,1.0_dp)

end module static

