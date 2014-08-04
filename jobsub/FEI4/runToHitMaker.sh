#!/bin/sh


first="3551"
last="3546"

RUNLIST="runlist-1fei4.csv"

for ((c=$first; c>=$last; c--));
do
 jobsub.py -c config.cfg -csv $RUNLIST converter  $c
 jobsub.py -c config.cfg -csv $RUNLIST clustering $c
 jobsub.py -c config.cfg -csv $RUNLIST hitmaker   $c
done
