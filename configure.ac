# Perlude.
AC_INIT([GTEnricher], [0.0.1], [yi.qiao@genetics.utah.edu])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

#Fixing file creation time
find ${srcdir} -name aclocal.m4 -exec touch {} \;
find ${srcdir} -name \*.in -exec touch {} \;
find ${srcdir} -name configure -exec touch {} \;

# Checks for programs
AC_PROG_CXX
AM_PROG_AR
LT_PREREQ([2.2])
LT_INIT()
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_14

AC_CONFIG_HEADERS([config.h])

## START: DEBUG / RELEASE ##
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging, default: no])],
              [case "${enableval}" in
                yes) debug=true ;;
                no)  debug=false ;;
                *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
               esac],
              [debug=false])
## END: DEBUG / RELEASE ##

AS_IF([test "x$debug" = xtrue], [AC_DEFINE([LOGLVL], [LOGLV_DEBUG], [Define log level to be debug])])


## Checks for libraries
AC_CHECK_LIB([hts], [bcf_read1], [], [AC_MSG_ERROR([Could not find htslib shared library])])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

echo \
"----------------------------------------------
${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'
Compiler: '${CXX}'
CFLAGS: '${CFLAGS}'
CPPFLAGS: '${CPPFLAGS}'
CXXFLAGS: '${CXXFLAGS}'

Extra variables:
LDFLAGS='${LDFLAGS}'
LDADD='${LDADD}'
LIBS='${LIBS}'

PTHREAD_LIBS='${PTHREAD_LIBS}'
PTHREAD_CFLAGS='${PTHREAD_CFLAGS}'
PTHREAD_CC='${PTHREAD_CC}'

DATADIR=${DATADIR}
debug=${debug}
"
