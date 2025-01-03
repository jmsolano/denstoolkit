#Compilation and installation instructions (DensToolKitViewer)
##Install dependencies
Before building DensToolKitViewer, please install some dependecies as described below. After this, continue with the building of DensToolKitViewer as outlined in "Common procedure"
###Ubuntu 14 dependencies
For compiling and installing DensToolKitViewer in Ubuntu 14 (most likely any Debian-based distro), please type in your terminal:

~~~~~~~~~~~~~
$sudo apt-get install aptitude build-essential make g++ git cmake freeglut3-dev qtdeclarative5-dev libxmu-dev libxi-dev mesa-common-dev libgl1-mesa-dev
~~~~~~~~~~~~~

###Fedora 21

For compiling and installing DensToolKitViewer under Fedora (and most likely under any RPM-based distro), please type in your terminal:

~~~~~~~~~~~~~
sudo yum install gcc-g++ git cmake gnuplot epstool texlive-epstopdf povray GraphicsMagick freeglut-devel qt5-qtbase-devel libXmu-devel libXi-devel
~~~~~~~~~~~~~

##Common procedure

After installing the corresponding packages, go to a directory wherein you plan to build dtk. In this example, replace ```/path/to/dtk/``` by the corresponding path.

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
#MacOSX
The easiest way to obtain the viewer is for you to download a compiled version of it (https://sites.google.com/site/jmsolanoalt/software/denstoolkit/downloads/DensToolKitViewer.dmg?attredirects=0&d=1). You may also try to compile it by yourself using cmake. 
For compiling under MacOSX, Qt5, OpenGL (usually distributed with XCode), and freeglut packages/libraries are required. 

Also, perhaps you may need to adjust a few lines in the CMakeFileLists.txt in order to tell cmake where is the root Qt5 directory (_i.e._ you may need to redefine the ```CMAKE_PREFIX_PATH``` variable by hand).
The rest sould be the same as for Ubuntu 14 (see above), but instead of a binary, you will get a MacOSX application bundle (in /path/to/dtk/builddtk/src/dtkview/build/bin/DensToolKitViewer.app)

When compiling in mac, you may obtain the following message:

~~~~~~~~~~~~~
.../moc_dtkglwidget.cpp:...: error:
use of undeclared identifier 'QGLWidget';
did you mean 'QWidget'?
~~~~~~~~~~~~~

To correct this, correct the names of QGLWidget ---> QOpenGLWidget:

~~~~~~~~~~~~~
vim moc_dtkglwidget.cpp

~~~~~~~~~~~~~

and use:

~~~~~~~~~~~~~
%s/QGL/QOpenGL
~~~~~~~~~~~~~
Recompile:

~~~~~~~~~~~~~
$make
~~~~~~~~~~~~~

