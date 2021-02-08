#!/bin/bash

for i in $(ls 37*.xyz);do
   ./dtkdrawer.x $i -A 5 4 11 -x 30 -Z 1.4
done

