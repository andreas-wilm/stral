# StrAl


See [Dalli et al. (2006)](http://www.ncbi.nlm.nih.gov/pubmed/16613908) for more info

This is a copy of StrAl taken from the SVN server  to make
sure the code survives. The original code used to
be available at
http://www.biophys.uni-duesseldorf.de/stral/


# Compilation

StrAl requires squid library 1.9g and the mhash library to be installed
beforehand.

After a first checkout you will need to setup automake once:

    $ ./reconf

and possibly

    $ [g]libtoolize

followed again by 

$ ./reconf

After that you can just use:

    $ ./configure
    $ make
    $ make install

# Welcome to StrAl.

StrAl is free software. Please see the file COPYING for details.
For documentation, please see the files in the doc subdirectory.
For building and installation instructions please see the INSTALL file.

For generic installation instructions se the file INSTALL.

The homepage of StrAl is
http://www.biophys.uni-duesseldorf.de/stral/

You will find the sources for download.
In the future we also plan to provide a webservice.

Contact, questions, bug reports:
    Gerhard Steger <steger@biophys.uni-duesseldorf.de>,
    Deniz Dalli <dalli@biophys.uni-duesseldorf.de>,
    or
    stral@biophys.uni-duesseldorf.de
