COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -g -O0
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -D__USE_BOINC__
#COMPILER = c++ -fpermissive -w -D__USE_BOINC__ -D__VERBOSE__ -static
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -static
SUFFIX = 

LIBDNSTD_PATH =/home/dcoss/tmp/blah/libdnstd/install
LIBMYGL_PATH = /home/dcoss/tmp/blah/libmygl/install
FFTW_PATH = $(HOME)/opt/fftw
MYINCLUDES = -I $(FFTW_PATH)/include/ -I $(LIBMYGL_PATH)/include/  -I $(LIBDNSTD_PATH)/include/

all: ray_trace_ellipse flatten utilities mass_to_shear physcalc makecluster

ray_trace_ellipse: ray_trace_ellipse.o DStackinstantiations.o DArrayinstantiations.o 
	${COMPILER} ${MYINCLUDES}  $^ -L${LIBMYGL_PATH}/lib -lmygl ${LIBDNSTD_PATH}/lib/libdnstd.a  -o ray_trace_ellipse${SUFFIX}

flatten: flattenmain.cpp flatten.o  DStackinstantiations.o 
	${COMPILER} ${MYINCLUDES} $^ -L${LIBMYGL_PATH}/lib -lmygl -L${LIBDNSTD_PATH}/lib -ldnstd -o flatten${SUFFIX}

utilities: utilitiesmain.cpp utilities.cpp Functions.o flatten.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ -ldnstd  -lfftw3 -lm -L${LIBMYGL_PATH}/lib -lmygl -L${LIBDNSTD_PATH}/lib -L$(FFTW_PATH)/lib -o utilities${SUFFIX}

mass_to_shear: mass_to_shear.cpp Functions.o utilities.o flatten.o  DStackinstantiations.o DArrayinstantiations.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYFFT_PATH}/libmyfft.a ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a  -o mass_to_shear${SUFFIX}

physcalc: physcalc.tab.c physcalc.yy.c
	${COMPILER} ${MYINCLUDES} -o physcalc physcalc.tab.c physcalc.yy.c  -lmygl -ldnstd -lm -L${LIBDNSTD_PATH}/lib -L${LIBMYGL_PATH}/lib  -lfftw3 -L${FFTW_PATH}/lib

makecluster: makecluster.cpp create_cluster.o makecluster 
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a -o makecluster${SUFFIX}

nfwshear: nfwshear.cpp 
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a -o nfwshear${SUFFIX}

clean: 
	rm -f *.o makecluster${SUFFIX}  physcalc${SUFFIX} mass_to_shear${SUFFIX} utilities${SUFFIX} flatten${SUFFIX} ray_trace_ellipse${SUFFIX} *.tab.* *.yy.c

.cpp.o: %.h
	${COMPILER} ${MYINCLUDES} -c $< -L$(FFTW_PATH)/include/

%.tab.c %.tab.h: %.y
	bison -d $<

%.yy.c %.yy.h: %.l
	flex -o$@ $<