.SUFFIXES: .F .F90 .f90 .o

OBJS = \
	gfs_tke_edmf_compat.o \
	tridi_mod.o           \
	mfpbltq_mod.o         \
	mfscuq_mod.o          \
	canopy_utils_mod.o    \
	sfc_diff.o            \
	iw3jdn.o              \
	sfc_nst_combined.o    \
	satmedmfvdifq.o

all: $(OBJS)
	ar -ru libphys_gfs.a $(OBJS)

physics_gfs_lib: $(OBJS)
	ar -ru ../libphys.a $(OBJS)

gfs_tke_edmf_compat.o: \
	../mpas_atmphys_constants.o

tridi_mod.o: \
	gfs_tke_edmf_compat.o

mfpbltq_mod.o: \
	gfs_tke_edmf_compat.o

mfscuq_mod.o: \
	gfs_tke_edmf_compat.o

canopy_utils_mod.o: \
	gfs_tke_edmf_compat.o

sfc_diff.o: \
	gfs_tke_edmf_compat.o

sfc_nst_combined.o: \
	gfs_tke_edmf_compat.o

satmedmfvdifq.o: \
	gfs_tke_edmf_compat.o \
	tridi_mod.o \
	mfpbltq_mod.o \
	mfscuq_mod.o \
	canopy_utils_mod.o

clean:
	$(RM) *.o *.mod *.F *.F90 *.f90 *.i libphys_gfs.a

.F.o:
	$(RM) $@ $*.mod
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(FFLAGS) -c $*.F \
	$(CPPINCLUDES) $(FCINCLUDES) \
	-I. -I.. -I../physics_mmm -I../physics_wrf \
	-I../../../framework -I../../../external/esmf_time_f90


.f90.o:
	$(RM) $@ $*.mod
	$(FC) $(FFLAGS) -c $*.f90 \
	$(FCINCLUDES) \
	-I. -I.. -I../physics_mmm -I../physics_wrf \
	-I../../../framework -I../../../external/esmf_time_f90
