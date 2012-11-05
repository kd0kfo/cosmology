#LIBDNSTD_PATH =/home/dcoss/tmp/blah/libdnstd/install
#LIBMYGL_PATH = /home/dcoss/tmp/blah/libmygl/install
#FFTW_PATH = $(HOME)/opt/fftw
#MYINCLUDES = -I $(FFTW_PATH)/include/ -I $(LIBMYGL_PATH)/include/  -I $(LIBDNSTD_PATH)/include/

bin_PROGRAMS = ray_trace_ellipse flatten  physcalc makecluster

ray_trace_ellipse_SOURCES = ray_trace_ellipse.cpp
ray_trace_ellipse_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/ $(NETCDF_CPPFLAGS)
ray_trace_ellipse_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib $(NETCDF_LDFLAGS)
ray_trace_ellipse_LDADD = -lmygl -ldnstd $(NETCDF_LDADD)


flatten_SOURCES =  flattenmain.cpp flatten.cpp
flatten_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/
flatten_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib
flatten_LDADD = -lmygl -ldnstd

AM_YFLAGS = -d
physcalc_SOURCES =  physcalc.yacc.ypp physcalc.cpp physcalc.lex.lpp functions.cpp
physcalc_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/
physcalc_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib
physcalc_LDADD = -lmygl -ldnstd 

makecluster_SOURCES =  makecluster.cpp create_cluster.cpp
makecluster_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/
makecluster_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib
makecluster_LDADD = -lmygl -ldnstd 


if BUILD_WITH_MPI
ray_trace_ellipse_SOURCES +=  mpi_utils.cpp
ray_trace_ellipse_CPPFLAGS += -I${MPI_PATH}/include
ray_trace_ellipse_LDFLAGS += -L${MPI_PATH}/lib -lmpi
endif

if BUILD_WITH_FFT
physcalc_CPPFLAGS += -I${FFTW3_PATH}/include
physcalc_LDFLAGS += -L${FFTW3_PATH}/lib
physcalc_LDADD += -lfftw3
endif

if LINK_STATICALLY
ray_trace_ellipse_LDFLAGS += -static
flatten_LDFLAGS += -static
physcalc_LDFLAGS += -static
makecluster_LDFLAGS += -static
endif