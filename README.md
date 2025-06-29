# DensToolKit

This is the version 2.0 of the suite DensToolKit.
If you like it, please consider it citing us in your work:

J.M. Solano-Altamirano, Julio M. Hernández-Pérez, J. Sandoval-Lira, and J. Barroso-Flores, *Jour. Chem. Phys.* (2024). [https://doi.org/10.1063/5.0239835](
https://doi.org/10.1063/5.0239835)

The first version was published in:

 J.M. Solano-Altamirano and Julio M. Hernández-Pérez, *Comput. Phys. Commun.* (2015), **196** (2015) 362-371.
 [http://dx.doi.org/10.1016/j.cpc.2015.07.005.](http://dx.doi.org/10.1016/j.cpc.2015.07.005.)
 
 Since June 2025, DensToolKit has been cloned the following number of times:
 
 <p align=center href='https://github.com/MShawon/github-clone-count-badge'><img alt='GitHub Clones' src='https://img.shields.io/badge/dynamic/json?color=success&label=Clone&query=count&url=https://gist.githubusercontent.com/jmsolano/ea34e49c1d80c50980410e1844b68535/raw/clone.json&logo=github'></p>

# Git instructions

You can get DensToolKit's source code as follows.
In your bash terminal type:

~~~~~~~~~~
$cd /local/path/to/dtk
$git clone https://github.com/jmsolano/denstoolkit.git
~~~~~~~~~~

After this, git will transfer the source files to ```/local/path/to/dtk/denstoolkit```


# Building the package

## Auxiliary programs

Optionally, you may want to first run

~~~~~~~~~~
$/local/path/to/dtk/src/checkdependencies
~~~~~~~~~~

This will check that auxiliary programs (to plot the data computed by DensToolKit) are installed
in your system. You can skip this step if you are only interested in obtaining the numerical
data. You may process the data later with your preferred plotting program.

## Compilation

For building the binaries, type:

~~~~~~~~~~
$cd denstoolkit/src
$make -j 4
$sudo make install
~~~~~~~~~~

This should compile and install the binaries into ```/usr/local/bin```

# Testing the suite

If you want to test the correct compilation of the program, type:

~~~~~~~~~~
$make runtest
~~~~~~~~~~

The above command will use the binaries present
in ```/local/path/to/dtk/denstoolkit/bin```, therefore, the
tests can be made without running "sudo make install".

# Updating DensToolKit

If you only need to get the latest version of DensToolKi, and you are not using
personal modifications to the source, in the local main directory (using the
above example should be ```/local/path/to/dtk/denstoolkit```):

~~~~~~~~~~
$git pull
$make distclean
$make
$sudo make install
~~~~~~~~~~

This should update your binaries in the local installation directory.


# Compiling with OpenMP

DensToolKit contains paralellized implementation of almost all basic functions. For using the parallel version (in this example with 4 cores), compile with the following commands:

~~~~~~~~~~
make SETDTKNPROC=4 -j4
~~~~~~~~~~

To use N processors, the above line, replace SETDTKNPROC=4 by SETDTKNPROC=N

~~~~~~~~~~
make SETDTKNPROC=N -j 4
~~~~~~~~~~

In the current version, it is not advised to use more than 4 Intel processors, nor more than 2 Apple M1 Silicon cores. There will be no speed improvement using more processors.

#DensToolKit manual
The manual, wherein some formal theory is developed, and also where more information about how to use the programs of the suite can be viewed, is ```/local/path/to/dtk/denstoolkit/tex/dtkmanual/dtk-manual.pdf```

#DensToolKitViewer
As of version 1.2.0, we provide an experimental graphical viewer. It is based on Qt, OpenGL, and GLUT. This program is under construction and is expected to change somewhat frequently in the near future. Unfortunately, this program requires a bit of extra-effort to compile it in Linux. For further instructions, please visit https://github.com/jmsolano/denstoolkit/tree/master/src/dtkview.
For MacOSX, we are temporarily distributing through JMSA's personal website: 

https://sites.google.com/site/jmsolanoalt/software/denstoolkit/downloads/DensToolKitViewer.dmg?attredirects=0&d=1

#Developer instructions

Please, notice that the described in the "Git instructions" section will provide read-only access to the repository.
If you would like to become a source code contributor of this project, then you would need a different approach. In essence you need to create your own github account, and then you fork this project (using the Fork button on the right upper corner of the github webpage). Your changes should be commited and pushed to your own repository.

After this, a pull request should be made. We will then
review your contribution, and consider it for inclusion in the official distribution. All contributions are welcome!

#Contributions and List of Contributors

We welcome all contributions to this project. If you would like to contribute to it, please tell us by creating a new issue (visit https://github.com/jmsolano/denstoolkit/issues and create a new issue).

For a complete list of contributors, please visit:
 
https://github.com/jmsolano/denstoolkit/blob/master/src/CONTRIBUTORS





