#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(mycosmo, 1.0, david@davecoss.com)
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
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

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([alarm floor pow sqrt strerror])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
