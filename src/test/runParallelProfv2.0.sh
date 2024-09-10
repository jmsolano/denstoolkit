#!/bin/bash


#set theNPts to the desired numver of evaluations.
theNPts=4096

#Compile with the following steps:
#   make fullclean
#   make test_gausswavefunction.x
#   mv test_gausswavefunction.x test_gausswavefunction01ncpu
#Repeat with NCPU=2,4,6,8... (below, NN is NCPU padded with zeros):
#   make fullclean
#   make SETDTKNPROC=NCPU test_gausswavefunction.x
#   mv test_gausswavefunction.x test_gausswavefunctionNNncpu


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

