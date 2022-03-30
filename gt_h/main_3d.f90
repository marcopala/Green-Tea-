! Copyright or © or Copr. Marco Pala (February 24, 2022)

! e-mail: marco.pala@cnrs.fr

! This software is a computer program whose purpose is 
! to perform self-consistent simulations of nanosystems with a full ab initio approach
! by using the density functional theory and the non-equilibrium Green's function method.

! This software is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
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
! knowledge of the CeCILL-B license and that you accept its terms.

program main


use static
use pseudopot_so_gen
use indata

implicit none 


call indata_readinput()

call indata_grid()

call read_QE_output()


write(*,*)
write(*,*)
write(*,*) ' This is an open source software distributed under the CECILL-B license "http://www.cecill.info".'
write(*,*) ' Please consider citing  '
write(*,*) ' M. G. Pala, P. Giannozzi, and D. Esseni, Phys. Rev. B 102, 045410 (2020)'
write(*,*) ' DOI: https://doi.org/10.1103/PhysRevB.102.045410'
write(*,*)
write(*,*)


stop


end program main
