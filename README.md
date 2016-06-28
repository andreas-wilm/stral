# StrAl: progressive alignment of non-coding RNA using base pairing probability vectors in quadratic time.

See [Dalli et al. (2006)](http://www.ncbi.nlm.nih.gov/pubmed/16613908) for more info

This is a near exact copy of StrAl 0.5.4 to make sure the code survives and is
compileable. The original is / used to be available at
http://www.biophys.uni-duesseldorf.de/stral/

Differences from the original 0.5.4:
- All precompiled libraries removed
- libweighbor only used if found at configure time


# Compilation

## Requirements:

StrAl requires ViennaRNA (version 1.8; not 2!), squid library 1.9g and
the mhash library to be installed beforehand. You can find copies of
versions which are known to work in the third-party folder.

If you install these libraries to non-standard directories let configure
know by modifying LDFLAGS and CFLAGS on the commandline.

Notes on ViennaRNA: On OsX, if you get `"ld: symbol(s) not found for
architecture x86_64"` while compiling/installing ViennaRNA try a newer
GCC (setting CC to gcc-mp-4.8 worked for me) or experiment with
`ARCHFLAGS="-arch x86_64"` and `--disable-openmp` as options to configure.


## First checkout

After a first checkout you will need to setup automake once:

    $ ./bootstrap

If you see correspondong warnings you might also need to run:

    $ automake --add-missing

Should ltmain.sh be mising please also run:

    $ [g]libtoolize && ./bootstrap


After that configure and Makefile.in's will have been generated and
you can do the GNU triple jump:

## Compile and install

    $ ./configure
    $ make
    $ make install

