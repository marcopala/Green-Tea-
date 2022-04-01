####
# This Makefile is adapted to work with either ifort or gfortran.
####

#FC90 = gfortran                                                                        # uncomment this line if you use gfortran
FC90 = ifort                                                                          # uncomment this line if you use ifort

# COMPILATION OPTIONS

#F90_FLAGS =  -O3 -cpp -fopenmp -ffree-line-length-300 -J$(MODDIR)                      # uncomment this line if you use gfortran
F90_FLAGS = -O3 -traceback -qopenmp -fpp -I$(MODDIR) -module $(MODDIR)                           # uncomment this line if you use ifort
#F90_FLAGS =  -check all -traceback -warn noexternal -fstack-protector -assume protect_parens -implicitnone -qopenmp -fpp -I$(MODDIR) -module $(MODDIR) 
# Library options in general
MKL_IFORT=/opt/intel/oneapi/mkl/latest/lib/intel64  
LIBS=-L$(MKL_IFORT) -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64 -lpthread
LIBS=-L$(MKL_IFORT) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64 -lpthread


# module directory
MODDIR = compiled

# Source directory
SRCDIR = .

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets
all : gt_t.x
gt_t.x : main_t.f90 static_t.o indata_t.o pseudopot_t.o fermi.o poisson_t.o negf_t.o 
	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) 

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules
%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) 
