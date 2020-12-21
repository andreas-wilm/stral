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

StrAl requires ViennaRNA (version 1.8; not 2!), Sean Eddy's Squid library
1.9g and the Mhash library to be installed beforehand. You can find copies of
versions which are known to work in the third-party folder of Stral's source code.
For each package, use the GNU triple jump to install it:
`./configure && make && sudo make install`. If you are not root, use
`./configure --prefix SOMEPATH` to install the packages to `SOMEPATH`
you have access to. In that case, let Stral's `configure`
know about the path by modifying LDFLAGS and CFLAGS on the commandline.
This is done by setting the CFLAGS and LDFLAGS environment variables.
For example, if you installed packages to `$HOME/local` run
`export CFLAGS="-I$HOME/local` and `export LDFLAGS="-L$HOME/local/`
before running `./configure`.

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

