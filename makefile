include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
FFLAGS := $(shell pkg-config --cflags fgsl)
FLIBS := $(shell pkg-config --libs fgsl)

basisf90: basisf90.o chkopts
	-${FLINKER} ${FFLAGS} ${FLIBS} -o basisf90 basisf90.o ${SLEPC_EPS_LIB}
	-${RM} basisf90.o

#------------------------------------------------------------------------------------
DATAPATH = ${SLEPC_DIR}/share/slepc/datafiles/matrices
