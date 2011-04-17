COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -g -O0
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -D__USE_BOINC__
#COMPILER = c++ -fpermissive -w -D__USE_BOINC__ -D__VERBOSE__ -static
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -static
SUFFIX = 

LIBDNSTD_PATH =/home/dcoss/opt/libdnstd/install
LIBMYGL_PATH = /home/dcoss/opt/libmygl/install
FFTW_PATH = $(HOME)/opt/fftw
MYINCLUDES = -I $(FFTW_PATH)/include/ -I $(LIBMYGL_PATH)/include/  -I $(LIBDNSTD_PATH)/include/

all: ray_trace_ellipse flatten  physcalc makecluster

ray_trace_ellipse: ray_trace_ellipse.o 
	${COMPILER} ${MYINCLUDES}  $^ -L${LIBMYGL_PATH}/lib -lmygl ${LIBDNSTD_PATH}/lib/libdnstd.a  -o ray_trace_ellipse${SUFFIX}

flatten: flattenmain.cpp flatten.o  
	${COMPILER} ${MYINCLUDES} $^ -L${LIBMYGL_PATH}/lib -lmygl -L${LIBDNSTD_PATH}/lib -ldnstd -o flatten${SUFFIX}

physcalc: physcalc.tab.c physcalc.yy.c
	${COMPILER} ${MYINCLUDES} -o physcalc physcalc.tab.c physcalc.yy.c  -lmygl -ldnstd -lm -L${LIBDNSTD_PATH}/lib -L${LIBMYGL_PATH}/lib  -lfftw3 -L${FFTW_PATH}/lib

makecluster: makecluster.cpp create_cluster.o makecluster 
	${COMPILER} ${MYINCLUDES} $^ -I${LIBMYGL_PATH}/include/ -I${LIBDNSTD_PATH}/include/ -o makecluster${SUFFIX} -lmygl -ldnstd -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib

clean: 
	rm -f *.o makecluster${SUFFIX}  physcalc${SUFFIX} mass_to_shear${SUFFIX} utilities${SUFFIX} flatten${SUFFIX} ray_trace_ellipse${SUFFIX} *.tab.* *.yy.c

.cpp.o: %.h
	${COMPILER} ${MYINCLUDES} -c $< -L$(FFTW_PATH)/include/

%.tab.c %.tab.h: %.y
	bison -d $<

%.yy.c %.yy.h: %.l
	flex -o$@ $<
