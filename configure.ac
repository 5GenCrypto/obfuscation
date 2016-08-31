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

AC_ARG_ENABLE(debug,        [  --enable-debug          Enable assert() statements for debugging.], [enable_debug=yes])

CFLAGS=                         dnl get rid of default -g -O2
COMMON_CFLAGS="-Wall -Wformat -Wformat-security -Wextra -Wunused \
-Wshadow -Wmissing-prototypes -Wfloat-equal -Wpointer-arith -Wcast-align \
-Wstrict-prototypes -Wredundant-decls -Wendif-labels -Wcast-qual \
-std=gnu11 -Wpedantic"

if test "x$enable_debug" = x"yes"; then
  EXTRA_CFLAGS="-O0 -g"
else
  EXTRA_CFLAGS="-O3"
  AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi
AC_SUBST(COMMON_CFLAGS)
AC_SUBST(EXTRA_CFLAGS)

AC_FUNC_MALLOC

AC_CHECK_HEADERS([omp.h])

AC_SEARCH_LIBS(aes_randinit,aesrand)
if test "x$ac_cv_search_aes_randinit" = "xno"; then
  AC_MSG_ERROR([libaesrand not found])
fi
AC_SEARCH_LIBS(fmpz_init,flint)
if test "x$ac_cv_search_fmpz_init" = "xno"; then
  AC_MSG_ERROR([libflint not found])
fi
AC_SEARCH_LIBS(clt_state_new,clt13)
if test "x$ac_cv_search_clt_state_new" = "xno"; then
  AC_MSG_ERROR([libclt13 not found])
fi
AC_SEARCH_LIBS(mmap_enc_mat_init,mmap)
if test "x$ac_cv_search_mmap_enc_mat_init" = "xno"; then
  AC_MSG_ERROR([libmmap not found])
fi

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
