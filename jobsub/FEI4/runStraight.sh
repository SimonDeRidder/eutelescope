#!/bin/sh


first="3569"
last="3566"

RUNLIST="runlist-1fei4.csv"

for ((c=$first; c>=$last; c--));
do
 #jobsub.py -c config.cfg -csv $RUNLIST align  $c
 jobsub.py -c config.cfg -csv $RUNLIST trackdaf $c
done
 jobsub.py -c config.cfg -csv $RUNLIST trackdaf 3551
