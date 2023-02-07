# suppose you have created a ui file called configwin.ui
cmake_minimum_required (VERSION 2.8)

project(DTKViewer)

set(the_usr_moc_hdrs
   ../common/bondnetwork.h 
   ../common/iofuncts-cpx.h 
   ../common/iofuncts-wfn.h 
   ../common/iofuncts-wfx.h 
   #../common/solpovtools.h 
   ../common/povraytools.h 
   #../common/solmemhand.h 
   ../common/mymemory.h 
   ../common/atomradiicust.h 
   #../common/solscrutils.h 
   ../common/screenutils.h 
   ../common/atomcolschjmol.h 
   #../common/solstringtools.h 
   ../common/stringtools.h 
   ../common/critptnetwork.h 
   ../common/gausswavefunction.h 
   #../common/eig2-4.h 
   ../common/eigendecompositionjama.h 
   ../common/figname.h
   #../common/solmath.h 
   ../common/mymath.h 
   #../common/solfileutils.h 
   ../common/fileutils.h 
   ../common/matrixvectoroperations3d.h
   ../common/atom.h
   ../common/molecule.h
   ../common/inputmolecule_xyz.h
   localdefs.h 
   dtkmainwindow.h
   dtkglwidget.h
   dtkglbondnetwork.h 
   dtkglutils.h
   dtkglcriticalpointnetwork.h
)

set(the_usr_bin_name "DensToolKitViewer")
set(usr_have_icons TRUE)

set(the_usr_rcs
   denstoolkitviewer.qrc
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
set(the_usr_include_dirs
   ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR}
)
set(the_usr_link_libraries
   ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES}
)

set(the_usr_uis
   dtkmainwindow.ui
)

set(the_usr_srcs
   main.cpp
   dtkmainwindow.cpp
   dtkglwidget.cpp
   ../common/bondnetwork.cpp 
   ../common/iofuncts-cpx.cpp 
   ../common/iofuncts-wfn.cpp 
   ../common/iofuncts-wfx.cpp 
   #../common/solpovtools.cpp 
   ../common/povraytools.cpp 
   #../common/solmemhand.cpp 
   ../common/mymemory.cpp 
   ../common/atomradiicust.cpp 
   #../common/solscrutils.cpp 
   ../common/screenutils.cpp 
   ../common/atomcolschjmol.cpp 
   #../common/solstringtools.cpp 
   ../common/stringtools.cpp 
   ../common/critptnetwork.cpp 
   ../common/gausswavefunction.cpp 
   #../common/eig2-4.cpp 
   ../common/eigendecompositionjama.cpp 
   #../common/solmath.cpp 
   ../common/mymath.cpp 
   #../common/solfileutils.cpp 
   ../common/fileutils.cpp 
   ../common/matrixvectoroperations3d.cpp
   ../common/atom.cpp
   ../common/molecule.cpp
   ../common/inputmolecule_xyz.cpp
   ../common/figname.cpp
   dtkglbondnetwork.cpp 
   dtkglutils.cpp
   dtkglcriticalpointnetwork.cpp
)

set(the_usr_qt5_modules Core Gui Widgets OpenGL)

set(CMAKE_AUTOMOC OFF)

if(UNIX)
   if(APPLE)
      set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}
         /Users/jmsolano/.binsol/Qt5/5.5/clang_64
      )
   else(APPLE)
      set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}
         /usr/lib/x86_64-linux-gnu/cmake
         /usr/lib/x86_64-linux-gnu/cmake/Qt5OpenGL
       )
      set(Qt5OpenGL_DIR 
         /usr/lib/x86_64-linux-gnu/cmake/Qt5OpenGL
      )
   endif(APPLE)
endif(UNIX)

find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt5Gui)
#if(UNIX AND NOT APPLE)
find_package(Qt5OpenGL)
#endif(UNIX AND NOT APPLE)



set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

QT5_ADD_RESOURCES(the_gen_rcs ${the_usr_rcs})

QT5_WRAP_UI(the_gen_uis ${the_usr_uis})

QT5_WRAP_CPP(the_gen_moc_hdrs ${the_usr_moc_hdrs})

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${the_usr_include_dirs})

# For Apple set the icns file containing icons
IF(APPLE)
   IF(${usr_have_icons})
# set how it shows up in the Info.plist file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      set(MACOSX_BUNDLE_ICON_FILE "${the_usr_bin_name}.icns") 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # set where in the bundle to put the icns file
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SET_SOURCE_FILES_PROPERTIES(
         ${CMAKE_CURRENT_SOURCE_DIR}/${the_usr_bin_name}.icns
         PROPERTIES MACOSX_PACKAGE_LOCATION Resources
      )
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # include the icns file in the target
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      set(the_usr_srcs 
         ${the_usr_srcs} 
         "${CMAKE_CURRENT_SOURCE_DIR}/${the_usr_bin_name}.icns"
      )
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ENDIF(${usr_have_icons})
ENDIF(APPLE)


# To put the application into a bin directory
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

#  CMake's way of creating an executable
add_executable( ${the_usr_bin_name} MACOSX_BUNDLE WIN32
  ${the_usr_srcs}
  ${the_gen_moc_hdrs} 
  ${the_gen_uis}
  ${the_gen_rcs}
)
#------------------------------------------------------------------------------


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Tell CMake which libraries we need to link our executable against.
TARGET_LINK_LIBRARIES (${the_usr_bin_name} 
   ${QT_LIBRARIES} ${the_usr_link_libraries}
)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Now the installation stuff below
#------------------------------------------------------------------------------
SET(plugin_dest_dir bin)
SET(qtconf_dest_dir bin)
SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/${the_usr_bin_name}")
IF(APPLE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SET(plugin_dest_dir ${the_usr_bin_name}.app/Contents/MacOS)
  SET(qtconf_dest_dir ${the_usr_bin_name}.app/Contents/Resources)
  SET(APPS "\${CMAKE_INSTALL_PREFIX}/${the_usr_bin_name}.app")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ENDIF(APPLE)
IF(WIN32)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/${the_usr_bin_name}.exe")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ENDIF(WIN32)

# Install the ${the_usr_bin_name} application, on Apple, the bundle is at the root of the
# install tree, and on other platforms it'll go into the bin directory.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
INSTALL(TARGETS ${the_usr_bin_name} 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
    )

# Install needed Qt plugins by copying directories from the qt installation
# One can cull what gets copied by using 'REGEX "..." EXCLUDE'
INSTALL(DIRECTORY "${QT_PLUGINS_DIR}/imageformats"
   DESTINATION ${plugin_dest_dir}/plugins COMPONENT Runtime
   )

# install a qt.conf file
# this inserts some cmake code into the install script to write the file
INSTALL(CODE "
    file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"\")
    " COMPONENT Runtime)


# Use BundleUtilities to get all other dependencies for the application to work.
# It takes a bundle or executable along with possible plugins and inspects it
# for dependencies.  If they are not system dependencies, they are copied.

# directories to look for dependencies
SET(DIRS ${QT_LIBRARY_DIRS})
qt5_use_modules(${the_usr_bin_name} ${the_usr_qt5_modules})
#execute_process(COMMAND "cat moc_dtkglwidget.cpp | sed 's/QGLWidget/QOpenGLWidget/' > moc_dtkglwidget.cpp"
#   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#   )

ADD_DEFINITIONS(-include localdefs.h -include ../common/globaldefs.h)

if (CMAKE_COMPILER_IS_GNUCXX)
   set(CMAKE_CXX_FLAGS "-Wall -O2 -m64 -funroll-loops -falign-loops=16")
   set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS})
endif()


# Now the work of copying dependencies into the bundle/package
# The quotes are escaped and variables to use at install time have their $ escaped
# An alternative is the do a configure_file() on a script and use install(SCRIPT  ...).
# Note that the image plugins depend on QtSvg and QtXml, and it got those copied
# over.
INSTALL(CODE 
   "file(GLOB_RECURSE QTPLUGINS
      \"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/plugins/*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
    include(BundleUtilities)
    fixup_bundle(\"${APPS}\" \"\${QTPLUGINS}\" \"${DIRS}\")"
    COMPONENT Runtime)


# To Create a package, one can run "cpack -G DragNDrop CPackConfig.cmake" on Mac OS X
# where CPackConfig.cmake is created by including CPack
# And then there's ways to customize this as well
set(CPACK_BINARY_DRAGNDROP ON)
include(CPack)


#add_executable(openglwindow ${the_usr_sources})
#qt5_use_modules(openglwindow Core Gui Widgets)
