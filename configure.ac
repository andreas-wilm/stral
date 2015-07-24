# CVS $Id: configure.in,v 1.6 2006-03-08 22:44:33 wilm Exp $

# Process this file with autoconf to produce a configure script.

# Autoconf-Manual, p. 15
# http://www.gnu.org/software/automake/manual/automake.html#Modernize-AM_005fINIT_005fAUTOMAKE-invocation
AC_INIT([stral], [0.5.5])
PACKAGE=stral
MAJOR_VERSION=0
MINOR_VERSION=5
PATCHLEVEL=5

# Autoconf-Manual, p. 16
AC_CONFIG_SRCDIR(config.h.in)
AM_CONFIG_HEADER(config.h)
#

# Version changes must also be applied below !
AM_INIT_AUTOMAKE([foreign])
#
VERSION=${MAJOR_VERSION}.${MINOR_VERSION}.${PATCHLEVEL}
NODOT_VERSION=${MAJOR_VERSION}_${MINOR_VERSION}_${PATCHLEVEL}
AC_SUBST(PACKAGE)
AC_SUBST(VERSION)
AC_SUBST(NODOT_VERSION)

#-----------------------------------------------------------------
#---------------------  Programs  --------------------------------
#-----------------------------------------------------------------

# Autoconf-Manual, p. 60
AC_PROG_CC
# Autoconf-Manual, p. 20, 60
CFLAGS="${CFLAGS} -ansi -Wall -O3 -g"

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

# math-check
AC_CHECK_LIB(m, exp,, AC_MSG_ERROR("Could not find libm."))

AC_CHECK_HEADERS([ViennaRNA/part_func.h], ,
   AC_MSG_WARN(" *** Couldn't find ViennaRNA headers! 
   Please install from http://www.tbi.univie.ac.at/RNA/ (or see third-party folder) ***"))
AC_CHECK_LIB(RNA, pf_fold, , AC_MSG_ERROR("Could not find libRNA (ViennaRNA)."), -lm)


# weighbor
AC_CHECK_LIB(weighbor, weighbor, , AC_MSG_WARN("Could not find libweighbor."), -lm)


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
AC_CHECK_HEADERS([fcntl.h float.h limits.h stddef.h stdlib.h string.h strings.h unistd.h])


###  squid
#
AC_CHECK_HEADERS([squid.h], ,
   AC_MSG_WARN(" *** Couldn't find squid headers! 
   Please install from http://selab.janelia.org/software.html (or see third-party folder) ***"))
AC_CHECK_LIB(squid, MSAFileOpen, , AC_MSG_ERROR("Could not find libsquid."), -lm)

### mhash
#
AC_CHECK_HEADERS([mhash.h], ,
   AC_MSG_WARN(" *** Couldn't find mhash header! 
   Please install from http://mhash.sourceforge.net/ (or see third-party folder) ***"))
AC_CHECK_LIB(mhash, mhash_init, , AC_MSG_ERROR("Could not find libmhash."), )

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



