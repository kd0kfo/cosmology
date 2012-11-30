bin_PROGRAMS = ray_trace_ellipse flatten  physcalc makecluster plane2netcdf

ray_trace_ellipse_SOURCES = ray_trace_ellipse.cpp version.c
ray_trace_ellipse_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/ $(NETCDF_CPPFLAGS)
ray_trace_ellipse_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib $(NETCDF_LDFLAGS)
ray_trace_ellipse_LDADD = -lmygl -ldnstd $(NETCDF_LDADD)


flatten_SOURCES =  flattenmain.cpp flatten.cpp
flatten_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/ $(NETCDF_CPPFLAGS)
flatten_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib $(NETCDF_LDFLAGS)
flatten_LDADD = -lmygl -ldnstd $(NETCDF_LDADD)

AM_YFLAGS = -d
physcalc_SOURCES =  physcalc.yacc.ypp physcalc.cpp physcalc.lex.lpp functions.cpp version.c
physcalc_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/ $(NETCDF_CPPFLAGS)
physcalc_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib $(NETCDF_LDFLAGS)
physcalc_LDADD = -lmygl -ldnstd  $(NETCDF_LDADD)

makecluster_SOURCES =  makecluster.cpp create_cluster.cpp
makecluster_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/ $(NETCDF_CPPFLAGS)
makecluster_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib $(NETCDF_LDFLAGS)
makecluster_LDADD = -lmygl -ldnstd  $(NETCDF_LDADD)

plane2netcdf_SOURCES = plane2netcdf.cpp
plane2netcdf_CPPFLAGS = -I$(LIBMYGL_PATH)/include/  -I$(LIBDNSTD_PATH)/include/ $(NETCDF_CPPFLAGS)
plane2netcdf_LDFLAGS = -L${LIBMYGL_PATH}/lib -L${LIBDNSTD_PATH}/lib $(NETCDF_LDFLAGS)
plane2netcdf_LDADD = -lmygl -ldnstd $(NETCDF_LDADD)

if BUILD_WITH_MPI
ray_trace_ellipse_SOURCES +=  mpi_utils.cpp
ray_trace_ellipse_CPPFLAGS += -I${MPI_PATH}/include
ray_trace_ellipse_LDFLAGS += -L${MPI_PATH}/lib
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


version.c: 
	git rev-parse HEAD | awk ' BEGIN {print "#include \"version.h\""} {print "const char * build_git_sha = \"" $$0"\";"} END {}' > version.c
	date | awk 'BEGIN {} {print "const char * build_git_time = \""$$0"\";"} END {} ' >> version.c

.PHONY: version.c