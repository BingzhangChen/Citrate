.SUFFIXES: .o .i .f90 .f
.f90.o:
	$(CFT) -c $(FFLAGS) $*.f90 -o $*.o
.f.o:
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

MAIN = DRAM
#Source file
SRCS = dsymv.f lsame.f xerbla.f Readcsv.f90 tridiagonal.f90\
	cholesky.f90 gridinterp.f90 time_interp.f90\
	mtfort90.f90 gammadev.f90 gasdev.f90 multiGauss.f90\
	params.f90 bio_MOD.f90 Monod.f90 NPZD_N2.f90 nif_growth.f90\
 	Advection_center.f90 Diff_center.f90 WAPR.f90 pnorm.f90\
	EFT_phygrowth.f90 EFT_2sp.f90 EFT_size.f90 NPZD_Fix.f90\
	NPZDPhy_size.f90 choose_model.f90\
	NPPZDD.f90 NPZD_2sp.f90 GeiderPhy_size.f90 Geider_disc.f90\
	Geider_simple.f90 EFT_simple.f90 EFT_disc.f90 EFT_PPZDD.f90\
	Geider_Droop.f90\
	FlexEFT_cont.f90 MOD_1D.f90 NPZD_cont.f90 CalcYSSE.f90\
	matmuls.f90 Citrate3.f90 NP_closure.f90\
	syminv.f90 Interface_MOD.f90 transform.f90 UpdateCVM.f90\
	SUB_MOD.f90 Main.f90

OBJS =$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
#Executable file
$(MAIN): $(OBJS)
	$(CFT) $(FFLAGS) $(LFLAGS) -o $(MAIN) $(OBJS)

