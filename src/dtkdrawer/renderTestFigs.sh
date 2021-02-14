#!/bin/bash

for i in $(ls 37*.xyz);do
   ./dtkdrawer.x $i -A 5 4 11
done

