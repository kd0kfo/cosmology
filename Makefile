COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -g -fno-implicit-templates
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -D__USE_BOINC__
#COMPILER = c++ -fpermissive -w -D__USE_BOINC__ -D__VERBOSE__ -static
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -static
SUFFIX = 

LIBDNSTD_PATH =/home/dcoss/tmp/blah/libdnstd/install
LIBMYGL_PATH = /home/dcoss/tmp/blah/libmygl/install
FFTW_PATH = $(HOME)/opt/fftw
MYINCLUDES = -I $(FFTW_PATH)/include/ -I $(LIBMYGL_PATH)/include/  -I $(LIBDNSTD_PATH)/include/

all: ray_trace_ellipse flatten utilities mass_to_shear mycosmo makecluster

ray_trace_ellipse: ray_trace_ellipse.o DStackinstantiations.o DArrayinstantiations.o 
	${COMPILER} ${MYINCLUDES}  $^ -L${LIBMYGL_PATH}/lib -lmygl ${LIBDNSTD_PATH}/lib/libdnstd.a  -o ray_trace_ellipse${SUFFIX}

flatten: flattenmain.cpp flatten.o  DStackinstantiations.o 
	${COMPILER} ${MYINCLUDES} $^ -L${LIBMYGL_PATH}/lib -lmygl -L${LIBDNSTD_PATH}/lib -ldnstd -o flatten${SUFFIX}

utilities: utilitiesmain.cpp utilities.cpp Functions.o flatten.o  DStackinstantiations.o DArrayinstantiations.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ -L${LIBMYGL_PATH}/lib -lmygl -L${LIBDNSTD_PATH}/lib -L$(FFTW_PATH)/lib -ldnstd  -lfftw3 -lm -o utilities${SUFFIX}

mass_to_shear: mass_to_shear.cpp Functions.o utilities.o flatten.o  DStackinstantiations.o DArrayinstantiations.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYFFT_PATH}/libmyfft.a ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a  -o mass_to_shear${SUFFIX}

mycosmo:  Functions.o utilities.o Command.o CommandWords.o Parser.cpp  flatten.o  DArrayinstantiations.o  DStackinstantiations.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a ${LIBMYFFT_PATH}/libmyfft.a -o mycosmo${SUFFIX} -L $(FFTW_PATH)/lib -lfftw3 -lm

makecluster: makecluster.cpp create_cluster.o makecluster 
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a -o makecluster${SUFFIX}

nfwshear: nfwshear.cpp 
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a -o nfwshear${SUFFIX}

clean: 
	rm -f *.o makecluster${SUFFIX}  mycosmo${SUFFIX} mass_to_shear${SUFFIX} utilities${SUFFIX} flatten${SUFFIX} ray_trace_ellipse${SUFFIX}

.cpp.o: %.h
	${COMPILER} ${MYINCLUDES} -c $< -L$(FFTW_PATH)/include/

