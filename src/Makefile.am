LIBDNSTD_PATH =/home/dcoss/tmp/blah/libdnstd/install
LIBMYGL_PATH = /home/dcoss/tmp/blah/libmygl/install
FFTW_PATH = $(HOME)/opt/fftw
MYINCLUDES = -I $(FFTW_PATH)/include/ -I $(LIBMYGL_PATH)/include/  -I $(LIBDNSTD_PATH)/include/

bin_PROGRAMS: ray_trace_ellipse flatten  physcalc makecluster

ray_trace_ellipse_SOURCES = ray_trace_ellipse.cpp
ray_trace_ellipse_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/
ray_trace_ellipse_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib
ray_trace_ellipse_LDADD = -lmygl -ldnstd


flatten_SOURCES =  flattenmain.cpp
flatten_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/
flatten_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib
flatten_LDADD = -lmygl -ldnstd

physcalc_SOURCES =  physcalc.tab.c physcalc.yy.c
physcalc_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/
physcalc_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib
physcalc_LDADD = -lmygl -ldnstd -lfftw3 

makecluster_SOURCES =  makecluster.cpp create_cluster.cpp
makecluster_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/
makecluster_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib
