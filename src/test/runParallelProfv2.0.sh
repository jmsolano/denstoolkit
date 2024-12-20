#!/bin/bash


#####################################################################
#set theNPts to the desired numver of evaluations.
#####################################################################
#
theNPts=4096
#
#####################################################################
#   READ THE FOLLOWING INSTRUCTIONS BEFOR RUNNING THIS SCRIPT.
#####################################################################
# Before runing this script, compile with the following steps:
#   make fullclean
#   make test_gausswavefunction.x
#   make test_gausswavefunction.x
#   mv test_gausswavefunction.x test_gausswavefunction01ncpu
#   
#Notice that for this test the first compilation must be executed twice.
# If you already had run the script ../checkdependencies, then
# the double initial compilation is not needed.
#
#Repeat with NCPU=2,4,6,8... (below, NN is NCPU padded with zeros):
#   make fullclean
#   make SETDTKNPROC=NCPU test_gausswavefunction.x
#   mv test_gausswavefunction.x test_gausswavefunctionNNncpu

#After compilation, run the script with:
# bash runParallelProfv2.0.sh
#
#####################################################################

echo "./test_gausswavefunction01ncpu -v -N $theNPts > profile01cpu.log"
./test_gausswavefunction01ncpu -v -N $theNPts > profile01cpu.log

echo "./test_gausswavefunction02ncpu -v -N $theNPts > profile02cpu.log"
./test_gausswavefunction02ncpu -v -N $theNPts > profile02cpu.log

echo "./test_gausswavefunction04ncpu -v -N $theNPts > profile04cpu.log"
./test_gausswavefunction04ncpu -v -N $theNPts > profile04cpu.log

echo "./test_gausswavefunction06ncpu -v -N $theNPts > profile06cpu.log"
./test_gausswavefunction06ncpu -v -N $theNPts > profile06cpu.log

echo "./test_gausswavefunction08ncpu -v -N $theNPts > profile08cpu.log"
./test_gausswavefunction08ncpu -v -N $theNPts > profile08cpu.log

echo "./test_gausswavefunction10ncpu -v -N $theNPts > profile10cpu.log"
./test_gausswavefunction10ncpu -v -N $theNPts > profile10cpu.log

#echo "Taring and gzipping..."
#tar cvf macm1ProfilesNcpu.tar profile*.log
#gzip -9 macm1ProfilesNcpu.tar

echo "Done!"

