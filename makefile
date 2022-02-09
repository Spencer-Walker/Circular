include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
FFLAGS := $(shell pkg-config --cflags fgsl)
FLIBS := $(shell pkg-config --libs fgsl)


forpy_mod: forpy_mod.o 
	-${FLINKER} -c forpy_mod.F90

stark_basisf90: stark_basisf90.o chkopts
	-${FLINKER} ${FFLAGS} ${FLIBS} -o stark_basisf90 stark_basisf90.o forpy_mod.o ${SLEPC_EPS_LIB} `python3-config --ldflags --embed`
	-${RM} stark_basisf90.o

atomic_basisf90: atomic_basisf90.o chkopts
	-${FLINKER} ${FFLAGS} ${FLIBS} -o atomic_basisf90 atomic_basisf90.o forpy_mod.o ${SLEPC_EPS_LIB} `python3-config --ldflags --embed`
	-${RM} atomic_basisf90.o

overlapsf90: overlapsf90.o chkopts
	-${FLINKER} ${FFLAGS} ${FLIBS} -o overlapsf90 overlapsf90.o forpy_mod.o ${SLEPC_EPS_LIB} -fno-lto `python3-config --ldflags --embed`
	-${RM} overlapsf90.o
#------------------------------------------------------------------------------------
DATAPATH = ${SLEPC_DIR}/share/slepc/datafiles/matrices
