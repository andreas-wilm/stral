# CVS $Id: Makefile.am,v 1.6 2006-11-06 15:12:08 dalli Exp $

# configure.in! INCLUDES = -I$(includedir) -I.
# configure.in! CFLAGS   = -ansi -Wall -O2 ${CFLAGS}
# configure.in! LIBS     = -lRNA -lsquid -lweighbor -lm
# confgiure.in! LDFLAGS  = -L. ${LDFLAGS}
# confgiure.in! LDADD = ./libRNA.a ./libsquid.a ./libweighbor.a ./libmhash.a -lm



bin_PROGRAMS = stral
stral_SOURCES = align.c create_matrix.c energy_const.h iofuncs.h \
	neighbor.h part_func.c vectorops.h align.h create_matrix.h \
	energy_par.h iterate.c nrutils.c part_func.h stral.c  \
	weighted_nj.c BIONJ.c dist.c fold_vars.h iterate.h nrutils.h \
	phylip.c stral.h  weighted_nj.h bionj.h dist.h iofuncs.c neighbor.c \
	pair_mat.h phylip.h utils.h fft.h fft.c algebraic.h algebraic.c\
	alfuncs.h alfuncs.c aliterate.c aliterate.h wsp.h wsp.c \
	basic_funcs.c basic_funcs.h 


# These files should be added to the distribution
# CAUTION: These are COMPILED libraries!!!!
EXTRA_DIST = libRNA.a libweighbor.a mhash_precomp squid_precomp
SUBDIRS=doc

