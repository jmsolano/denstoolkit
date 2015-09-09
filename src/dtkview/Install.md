#Compilation and installation instructions (DensToolKitViewer)
##Ubuntu 14
For compiling and installing DensToolKitViewer, please type in your terminal:

~~~~~~~~~~~~~
$sudo apt-get install aptitude build-essential make g++ git cmake freeglut3-dev qtdeclarative5-dev libxmu-dev libxi-dev
~~~~~~~~~~~~~

Then go to a directory wherein you plan to build dtk. In this example, replace ```/path/to/dtk/``` by the corresponding path.

~~~~~~~~~~~~~
$cd /path/to/dtk
$git clone git://github.com/jmsolano/denstoolkit.git builddtk
$cd builddtk/src/dtkview
$mkdir build
$cd build
$cmake ..
$make
$sudo cp bin/DensToolKitViewer /usr/local/bin
~~~~~~~~~~~~~

This should make available DensToolKitViewer from the command line.

