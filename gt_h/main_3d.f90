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
write(*,*) ' This program is an open-source code; please cite  '
write(*,*) ' M. G. Pala, P. Giannozzi, and D. Esseni, Phys. Rev. B 102, 045410 (2020)'
write(*,*)
write(*,*)

stop


end program main
