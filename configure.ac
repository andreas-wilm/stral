# CVS $Id: configure.in,v 1.6 2006-03-08 22:44:33 wilm Exp $

# Process this file with autoconf to produce a configure script.

# Autoconf-Manual, p. 15
AC_INIT(stral.c)
# Autoconf-Manual, p. 16
AC_CONFIG_SRCDIR(config.h.in)
AM_CONFIG_HEADER(config.h)
#

# Version changes must also be applied below !
AM_INIT_AUTOMAKE(stral, 0.5.4)
PACKAGE=stral
MAJOR_VERSION=0
MINOR_VERSION=5
PATCHLEVEL=4
#
VERSION=${MAJOR_VERSION}.${MINOR_VERSION}.${PATCHLEVEL}
NODOT_VERSION=${MAJOR_VERSION}_${MINOR_VERSION}_${PATCHLEVEL}
#
AC_SUBST(PACKAGE)
AC_SUBST(VERSION)
AC_SUBST(NODOT_VERSION)

#-----------------------------------------------------------------
#---------------------  Programs  --------------------------------
#-----------------------------------------------------------------

# Autoconf-Manual, p. 60
AC_PROG_CC
# Autoconf-Manual, p. 20, 60
CFLAGS="${CFLAGS} -ansi -Wall -O2"

# Autoconf-Manual, p. 35
AC_PROG_INSTALL

# Autoconf-Manual, p. 17
AC_PROG_MAKE_SET

# Autoconf-Manual, p. 36
AC_PROG_LN_S

# Autoconf-Manual, p. 3
AM_PROG_LIBTOOL 



#------------------------------------------------------------------------------
#------------------      LIBRARIES / EXT PROGRAMS  ----------------------------
#------------------------------------------------------------------------------
# Autoconf-Manual, p. 38
LIBS="-L."

# math-check
AC_CHECK_LIB(m, exp,, AC_MSG_ERROR("Could not find libm."))

# CAUTION: This checks for libs in local directory plus in $LIB_PATH
# lib from Vienna package
AC_CHECK_LIB(RNA, pf_fold, , AC_MSG_ERROR("Could not find libRNA (Vienna)."), -L. -lm)


# weighbor
AC_CHECK_LIB(weighbor, weighbor, , AC_MSG_ERROR("Could not find libweighbor."), -L. -lm)


#------------------------------------------------------------------------------
#------------------      Checks for header files   ----------------------------
#------------------------------------------------------------------------------

# Autoconf-Manual, p. 29
AH_VERBATIM([_GNU_SOURCE],
    [/* Enable GNU extensions on systems that have them. */
    #ifndef _GNU_SOURCE
        # define _GNU_SOURCE
    #endif])

# Autoconf-Manual, p. 49
AC_HEADER_DIRENT
# Autoconf-Manual, p. 51
AC_HEADER_STDC
# Autoconf-Manual, p. 53
AC_CHECK_HEADERS([fcntl.h float.h limits.h stddef.h stdlib.h string.h strings.h unistd.h squid.h])


###  squid
#
AC_CHECK_HEADERS([squid.h],,
   AC_MSG_WARN(" *** Couldn't find squid headers! Will use local header and library ***"); LIBS="${LIBS} -L. -Lsquid_precomp"; CFLAGS="${CFLAGS} -I./squid_precomp")
AC_CHECK_LIB(squid, MSAFileOpen,,
   AC_MSG_ERROR("Could not find libsquid."), -L. -lm)

### mhash
#
AC_CHECK_HEADERS([mhash.h],,
   AC_MSG_WARN(" *** Couldn't find mhash header! Will use local header and library ***"); LIBS="${LIBS} -L. -Lmhash_precomp"; CFLAGS="${CFLAGS} -I./mhash_precomp")
AC_CHECK_LIB(mhash, mhash_init, , AC_MSG_ERROR("Could not find libmhash."), -L.)

#------------------------------------------------------------------------------
#------- Checks for typedefs, structures, and compiler characteristics --------
#------------------------------------------------------------------------------

# Autoconf-Manual, p. 50
AC_HEADER_STDBOOL
# Autoconf-Manual, p. 61
AC_C_CONST
# Autoconf-Manual, p. 57
AC_TYPE_MODE_T
AC_TYPE_SIZE_T

#------------------------------------------------------------------------------
#--------------------- Checks for library functions ---------------------------
#------------------------------------------------------------------------------

# Autoconf-Manual, p. 44
AC_FUNC_MALLOC
# Autoconf-Manual, p. 45
AC_FUNC_REALLOC
# Autoconf-Manual, p. 57
AC_TYPE_SIGNAL

# Autoconf-Manual, p. 46
AC_CHECK_FUNCS(getcwd memset mkdir strchr strdup strrchr strstr)


AC_OUTPUT([
	Makefile doc/Makefile
])



