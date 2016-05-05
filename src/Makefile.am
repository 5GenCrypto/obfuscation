AUTOMAKE_OPTIONS = foreign -Wall
AM_CFLAGS  = $(COMMON_CFLAGS) $(EXTRA_CFLAGS) -I$(top_srcdir) -I$(top_srcdir)/dgs
AM_LDFLAGS = -lgomp

lib_LTLIBRARIES=libobf.la

libobf_la_SOURCES = circuit.c obfuscator.c thpool.c thpool_fns.c utils.c zobfuscator.c
libobf_la_LDFLAGS = -release 0.0.0 -no-undefined

pkgincludesubdir = $(includedir)/obf
pkgincludesub_HEADERS = obfuscator.h zobfuscator.h

