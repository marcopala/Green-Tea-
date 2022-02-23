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
write(*,*) ' This is an open source software distributed under the CECILL-B license.'
write(*,*) ' Please consider citing  '
write(*,*) ' M. G. Pala, P. Giannozzi, and D. Esseni, Phys. Rev. B 102, 045410 (2020)'
write(*,*) ' DOI: https://doi.org/10.1103/PhysRevB.102.045410'
write(*,*)
write(*,*)


stop


end program main
