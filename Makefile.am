# CVS $Id: Makefile.am,v 1.5 2006-03-08 22:44:30 wilm Exp $

# configure.in! INCLUDES = -I$(includedir) -I.
# configure.in! CFLAGS   = -ansi -Wall -O2 ${CFLAGS}
# configure.in! LIBS     = -lRNA -lsquid -lweighbor -lm
# confgiure.in! LDFLAGS  = -L. ${LDFLAGS}
# confgiure.in! LDADD = ./libRNA.a ./libsquid.a ./libweighbor.a ./libmhash.a -lm



bin_PROGRAMS = stral
stral_SOURCES = align.c create_matrix.c energy_const.h iofuncs.h \
	neighbor.h part_func.c vectorops.h align.h create_matrix.h \
	energy_par.h nrutils.c part_func.h stral.c  \
	weighted_nj.c BIONJ.c dist.c fold_vars.h  nrutils.h \
	phylip.c stral.h  weighted_nj.h bionj.h dist.h iofuncs.c neighbor.c \
	pair_mat.h phylip.h utils.h alfuncs.c alfuncs.h config.h

# These files should be added to the distribution
# CAUTION: These are COMPILED libraries!!!!
EXTRA_DIST = libRNA.a libweighbor.a mhash_precomp squid_precomp
SUBDIRS=doc
