#FC90 = gfortran                                                                        # uncomment this line if you use gfortran
FC90 = ifort                                                                          # uncomment this line if you use ifort

# COMPILATION OPTIONS

#F90_FLAGS =  -O3 -cpp -fopenmp -ffree-line-length-300 -J$(MODDIR)                      # uncomment this line if you use gfortran
F90_FLAGS =  -O3 -traceback -qopenmp -fpp -I$(MODDIR) -module $(MODDIR)                           # uncomment this line if you use ifort
#F90_FLAGS =  -traceback -check bounds -qopenmp -fpp -I$(MODDIR) -module $(MODDIR)

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
vpath %.mod $(MODDIR)

# Targets
all : gt_h.x
gt_h.x : main_3d.f90 static_3d.o indata_3d.o pseudopot_so_gen.o 
	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) 

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules
%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) 
