#Compilation and installation instructions (DensToolKitViewer)

##Linux based systems

###Install dependencies
Before building DensToolKitViewer, please install some dependecies as described below. After this, continue with the building of DensToolKitViewer as outlined in "Common procedure"

####Ubuntu 14
For compiling and installing DensToolKitViewer in Ubuntu 14 (most likely any Debian-based distro), please type in your terminal:

~~~~~~~~~~~~~
sudo apt-get install aptitude build-essential make g++ git cmake freeglut3-dev qtdeclarative5-dev libxmu-dev libxi-dev mesa-common-dev libgl1-mesa-dev
~~~~~~~~~~~~~


####Fedora 21

For compiling and installing DensToolKitViewer under Fedora (and most likely under any RPM-based distro), please type in your terminal:

~~~~~~~~~~~~~
sudo yum install gcc-g++ git cmake gnuplot epstool texlive-epstopdf povray GraphicsMagick freeglut-devel qt5-qtbase-devel libXmu-devel libXi-devel
~~~~~~~~~~~~~

###Compilation in Linux

After installing the corresponding packages, go to a directory wherein you plan to build dtk. In this example, replace ```/path/to/dtk/``` by the corresponding path.

~~~~~~~~~~~~~
cd /path/to/dtk
git clone git://github.com/jmsolano/denstoolkit.git builddtk
cd builddtk/src/dtkview
mkdir build
cd build
cmake ..
make
sudo cp bin/DensToolKitViewer /usr/local/bin
~~~~~~~~~~~~~


This should make available DensToolKitViewer from the command line. If this does not work, please check next section.

###Known probles/fixes
If there is a broken link to EGL and GL
libraries, execute:

~~~~~~~~~~~~~
sudo rm /usr/lib/x86_64-linux-gnu/libEGL.so; sudo ln /usr/lib/x86_64-linux-gnu/libEGL.so.1 /usr/lib/x86_64-linux-gnu/libEGL.so

sudo rm /usr/lib/x86_64-linux-gnu/libGL.so; sudo ln /usr/lib/x86_64-linux-gnu/libGL.so.1 /usr/lib/x86_64-linux-gnu/libGL.so
~~~~~~~~~~~~~

##MacOSX

The easiest way to obtain the viewer is for you to download a compiled version of it (https://sites.google.com/site/jmsolanoalt/software/denstoolkit/downloads/DensToolKitViewer.dmg?attredirects=0&d=1). Unfortunately, this pre-compiled version only works for Macs with Intel processors.

For Macs with Apple Silicon processors (or if you want to compile your own version in your Mac), you may compile it using cmake, after installing Qt5 and FreeGlut.


### Qt5 libraries

For compiling under MacOSX, Qt5, OpenGL (which is usually distributed toghether with XCode), and freeglut packages/libraries are required. These libraries can be obtained *via* [HomeBrew](https://brew.sh) or [MacPorts](https://www.macports.org).

#### HomeBrew

~~~~~~~~~~~~~
brew install qt@5 freeglut
~~~~~~~~~~~~~

#### MacPorts

~~~~~~~~~~~~~
sudo port install qt5 freeglut
~~~~~~~~~~~~~

### Compilation in MacOSX

~~~~~~~~~~~~~
cd /path/to/dtk
git clone git://github.com/jmsolano/denstoolkit.git builddtk
cd builddtk/src/dtkview
mkdir build
cd build
cmake -DQT5_LIB_ROOT=/opt/homebrew/Cellar/qt\@5/5.15.16_1 ..
make
~~~~~~~~~~~~~

**Notice:** Replace the string ```/opt/homebrew/Cellar/qt\@5/5.15.16_1``` by the path wherein Qt5 libraries are installed. If you used Homebrew, this means to use the newest version of the library. Usually only the ```5.15.16_1``` needs to be replaced. If you installed Qt5 using macports, then the path is usually something like ```/opt/local/lib```. You can also install directly Qt5 thorugh Qt site; in this case, use the path wherein you installed Qt; e.g. ```/opt/local/Qt5/5.5/clang_64```.


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

(Alternately, open ```moc_dtkglwidget.cpp``` and replace these strings; 3 replacements are required.)

Recompile:

~~~~~~~~~~~~~
$make
~~~~~~~~~~~~~

### Installation in MacOSX

After a succesfull compilation, a MacOSX application bundle will be created as: ```/path/to/dtk/builddtk/src/dtkview/build/bin/DensToolKitViewer.app```

Finally, for installing the bundle in your system:

~~~~~~~~~~~~~
sudo cp -R bin/DensToolKitViewer.app /Applications/
~~~~~~~~~~~~~

