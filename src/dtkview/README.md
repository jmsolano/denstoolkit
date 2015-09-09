#Compilation and installation instructions (DensToolKitViewer)
##Ubuntu 14
For compiling and installing DensToolKitViewer, please type in your terminal:

~~~~~~~~~~~~~
$sudo apt-get install aptitude build-essential make g++ git cmake freeglut3-dev qtdeclarative5-dev libxmu-dev libxi-dev mesa-common-dev libgl1-mesa-dev
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

This should make available DensToolKitViewer from the command line. If this does not work, please check next section.

###Known probles/fixes
If there is a broken link to EGL and GL
libraries, execute:

~~~~~~~~~~~~~
$sudo rm /usr/lib/x86_64-linux-gnu/libEGL.so; sudo ln /usr/lib/x86_64-linux-gnu/libEGL.so.1 /usr/lib/x86_64-linux-gnu/libEGL.so

$sudo rm /usr/lib/x86_64-linux-gnu/libGL.so; sudo ln /usr/lib/x86_64-linux-gnu/libGL.so.1 /usr/lib/x86_64-linux-gnu/libGL.so
~~~~~~~~~~~~~

