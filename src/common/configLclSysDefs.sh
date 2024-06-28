#!/bin/bash
#
#*****************************************************************************************
#*****************************************************************************************
#The following part checks the existence of several programs needed for the proper running
# of the suite denstoolkit and save configuration macros into lclauxsoftwaredefs.h
#*****************************************************************************************
# 
if [ "$#" == "0" ]; then
   setverbose="F"
else
   setverbose="$1"
fi
have_gnuplot=1
if command -v gnuplot >/dev/null 2>&1 || command -v wgnuplot >/dev/null 2>&1; then
   if [ "$setverbose" == "T" ];then
      echo -e "We have gnuplot!\033[33m"
      echo "However, please ensure that the x11 terminal is available."
      echo -e "\033[m"
   fi
   have_gnuplot=0
   echo "#define _GNUPLOT_MAJ_VERSION_ $(gnuplot --version \
      | awk '{print $2}' | sed -e 's/\(.*\)\.\(.*\)/\1/g')" > lclauxsoftwaredefs.h
else
   if [ "$setverbose" == "T" ];then
      echo -e "\033[31m"
      echo -e "\nError: gnuplot must be installed in your system!\n"
      echo -e "\033[m"
   fi
fi
#
have_epstool=1
if command -v epstool >/dev/null 2>&1; then
   if [ "$setverbose" == "T" ];then
      echo "We have epstool!"
   fi
   have_epstool=0
else
   if [ "$setverbose" == "T" ];then
      echo -e "\033[31m"
      echo -e "\nError: epstool must be installed in your system!\n"
      echo -e "\033[m"
   fi
fi
#
have_epstopdf=1
if command -v epstopdf >/dev/null 2>&1; then
   if [ "$setverbose" == "T" ];then
      echo "We have epstopdf!"
   fi
   have_epstopdf=0
   else
   if [ "$setverbose" == "T" ];then
      echo -e "\033[31m"
      echo -e "\nError: epstopdf must be installed in your system!\n"
      echo -e "\033[m"
   fi
fi
#
have_povray=1
if command -v povray >/dev/null 2>&1 || command -v pvengine64.exe >/dev/null 2>&1 \
    || command -v pvengine.exe >/dev/null 2>&1; then
    if [ "$setverbose" == "T" ];then
       echo "We have povray!"
    fi
    have_povray=0
 else
    if [ "$setverbose" == "T" ];then
       echo -e "\033[31m"
       echo -e "\nError: povray must be installed in your system!\n"
       echo -e "\033[m"
    fi
fi
#
have_graphicsmagick=1
if command -v gm >/dev/null 2>&1; then
   if [ "$setverbose" == "T" ];then
      echo "We have gm!"
   fi
   have_graphicsmagick=0
else
   if [ "$setverbose" == "T" ];then
      echo -e "\033[31m"
      echo -e "\nError: gm must be installed in your system!\n"
      echo -e "\033[m"
   fi
fi
#
have_gzip=1
if command -v gzip >/dev/null 2>&1; then
   if [ "$setverbose" == "T" ];then
      echo "We have gzip!"
   fi
   have_gzip=0
else
   if [ "$setverbose" == "T" ];then
      echo -e "\033[31m"
      echo -e "\nError: gzip must be installed in your system!\n"
      echo -e "\033[m"
   fi
fi
have_all=$(($have_gnuplot + $have_epstool + $have_epstopdf + $have_povray\
   + $have_graphicsmagick + $have_gzip))
#echo -e "have_all=$have_all"
if [ "$setverbose" == "T" ];then
   if [[ $have_all -gt "0" ]];then
      echo -e "\033[31m"
      echo -e "\nError: Some of the packages are not installed in your system!\n"
      echo -e "       Please install them befor compiling DTK.\n"
      echo -e "\033[m"
   fi
   echo -e "\033[32m"
   echo -e "\nAll required packages are installed, you can proceed to"
   echo -e "build and install denstoolkit!\n"
   echo -e "\033[m"
fi
#*****************************************************************************************
#*****************************************************************************************
