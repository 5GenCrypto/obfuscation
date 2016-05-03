m4_define([version_major], [0])
m4_define([version_minor], [2])
m4_define([version_micro], [0])

AC_INIT([libobf], [version_major.version_minor.version_micro], [amaloz@cs.umd.edu])

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([./])

AC_CONFIG_AUX_DIR([build/autoconf])
AC_CONFIG_MACRO_DIR([build/autoconf])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_DEFINE(VERSION_MAJOR, version_major, [libgarble major version])
AC_DEFINE(VERSION_MINOR, version_minor, [libgarble minor version])
AC_DEFINE(VERSION_MICRO, version_micro, [libgarble micro version])

AM_PROG_AR

LT_INIT

AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

AC_MSG_CHECKING([if debug option is enabled])
AC_ARG_ENABLE(debug,
  [AS_HELP_STRING([--enable-debug], [enable debugging, default: no])],
  [with_debug=y],
  [with_debug=n])

CFLAGS="-Wall -Wformat -Wformat-security -Wextra -Wunused \
-Wshadow -Wmissing-prototypes -Wfloat-equal -Wpointer-arith -Wcast-align \
-Wstrict-prototypes -Wredundant-decls -Wendif-labels -Wcast-qual \
-std=gnu11 -Wpedantic"
if test x"$with_debug" == x"y"; then
   CFLAGS="$CFLAGS -g3 -O0"
else
    CFLAGS="$CFLAGS -fomit-frame-pointer -O3"
fi

AC_FUNC_MALLOC

AC_SEARCH_LIBS(fmpz_init,flint)
if test "x$ac_cv_search_fmpz_init" = "xno"; then
  AC_MSG_ERROR([libflint not found])
fi
AC_SEARCH_LIBS(mmap_enc_mat_init,mmap)
if test "x$ac_cv_search_mmap_enc_mat_init" = "xno"; then
  AC_MSG_ERROR([libmmap not found])
fi

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
