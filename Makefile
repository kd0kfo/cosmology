COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -g
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -D__DEBUG__ -static -D__USE_BOINC__
#COMPILER = c++ -fpermissive -w -D__USE_BOINC__ -D__VERBOSE__ -static
#COMPILER = c++ -fpermissive -w -D__VERBOSE__ -static
SUFFIX = 

LIBDNSTD_PATH = /home/david/opt/libdnstd
LIBMYGL_PATH = /home/david/opt/libmygl
LIBMYFFT_PATH = $(LIBDNSTD_PATH)
MYINCLUDES = -I../libmygl/src/ -I../libdnstd/src/ -I../fft/ -I../libmygl/src/EasyBMP/  -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/

all: ray_trace_ellipse flatten utilities mass_to_shear mycosmo makecluster

ray_trace_ellipse: ray_trace_ellipse.o planeinstantiations.o DStackinstantiations.o DArrayinstantiations.o 
	${COMPILER} ${MYINCLUDES}  $^ -L${LIBMYGL_PATH}/lib -lmygl ${LIBDNSTD_PATH}/lib/libdnstd.a  -o ray_trace_ellipse${SUFFIX}

flatten: flattenmain.cpp flatten.o planeinstantiations.o DStackinstantiations.o 
	${COMPILER} ${MYINCLUDES} $^ -L${LIBMYGL_PATH}/lib -lmygl -L${LIBDNSTD_PATH}/lib -ldnstd -o flatten${SUFFIX}

utilities: utilitiesmain.cpp utilities.cpp Functions.o flatten.o planeinstantiations.o DStackinstantiations.o DArrayinstantiations.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ -L${LIBMYGL_PATH}/lib -lmygl -L${LIBDNSTD_PATH}/lib -ldnstd  -lfftw3 -lm -o utilities${SUFFIX}

mass_to_shear: mass_to_shear.cpp Functions.o utilities.o flatten.o planeinstantiations.o DStackinstantiations.o DArrayinstantiations.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYFFT_PATH}/libmyfft.a ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a  -o mass_to_shear${SUFFIX}

mycosmo:  Functions.o utilities.o Command.o CommandWords.o Parser.cpp  flatten.o planeinstantiations.o DArrayinstantiations.o  DStackinstantiations.o  Rainbow.o
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a ${LIBMYFFT_PATH}/libmyfft.a -o mycosmo${SUFFIX} -lfftw3 -lm

makecluster: makecluster.cpp create_cluster.o makeclusterplaneinstantiations.o 
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a -o makecluster${SUFFIX}

nfwshear: nfwshear.cpp 
	${COMPILER} ${MYINCLUDES} $^ ${LIBMYGL_PATH}/libmygl.a ${LIBDNSTD_PATH}/libdnstd.a -o nfwshear${SUFFIX}

clean: 
	rm -f *.o makecluster${SUFFIX}  mycosmo${SUFFIX} mass_to_shear${SUFFIX} utilities${SUFFIX} flatten${SUFFIX} ray_trace_ellipse${SUFFIX}

.cpp.o: %.h
	${COMPILER} ${MYINCLUDES} -c $< -lfftw3 -lm

