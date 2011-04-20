#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(mycosmo, 1.1, david@davecoss.com)
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

#templates
AH_TEMPLATE([USE_MPI],[Define to compile with MPI.])

# Checks for programs.
AC_PROG_CXX([mpic++ icc g++])
AC_PROG_CC([mpicc icc gcc])
AC_PROG_YACC
AC_PROG_LEX

# Checks for libraries.
# FIXME: Replace `main' with a function in `-ldnstd':
AC_CHECK_LIB([dnstd], [main])
# FIXME: Replace `main' with a function in `-lfftw3':
AC_CHECK_LIB([fftw3], [main])
# FIXME: Replace `main' with a function in `-lmygl':
AC_CHECK_LIB([mygl], [main])

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h libintl.h malloc.h stddef.h stdlib.h string.h syslog.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AM_PROG_CC_C_O

#libdnstd
AC_ARG_WITH([libdnstd],AS_HELP_STRING([--with-libdnstd@<:@=dir@:>@],[Include libdnstd library @<:@Default=PREFIX/include@:>@]),with_libdnstd=${withval},AC_MSG_RESULT(Assuming libdnstd is installed in system paths))
if test x$with_libdnstd != x; then
   AC_SUBST(LIBDNSTD_PATH,$with_libdnstd)
fi

#libmygl
AC_ARG_WITH([libmygl],AS_HELP_STRING([--with-libmygl@<:@=dir@:>@],[Include libmygl library @<:@Default=PREFIX/include@:>@]),with_libmygl=${withval},AC_MSG_RESULT(Assuming libmygl is installed in system paths))
if test x$with_libmygl != x; then
   AC_SUBST(LIBMYGL_PATH,$with_libmygl)
fi

#fftw3
AC_ARG_WITH([fftw3],AS_HELP_STRING([--with-fftw3@<:@=dir@:>@],[Include fftw3 library @<:@Default=PREFIX/include@:>@]),with_fftw3=${withval},AC_MSG_RESULT(Assuming fftw3 is installed in system paths))
if test x$with_fftw3 != x; then
   AC_SUBST(FFTW3_PATH,$with_fftw3)
fi

#mpi
AC_ARG_WITH([mpi],AS_HELP_STRING([--with-mpi],[Includes the use of MPI.]),
[case "${withval}" in
      yes) mpi=. ;;
      no) mpi=false ;;
      *) mpi=${withval} ;;
esac]
AC_SUBST(MPI_PATH,$mpi)
AC_MSG_RESULT(Compiling with MPI)
AC_DEFINE([USE_MPI],[1])
,
[mpi=false]
)
AM_CONDITIONAL([BUILD_WITH_MPI],[test x$mpi != xfalse])



# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([alarm floor pow sqrt strerror])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
