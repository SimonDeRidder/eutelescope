#!/bin/sh

# parameters:
# 1 - runnumber         $RUN


while getopts a:r: option
do
        case "${option}"
        in
                r) RUN=${OPTARG};;
				a) ANGLE=${OPTARG};;
       esac
done

BeginX=0.
EndX=2.
BeginY=-1.
EndY=1.
Step=0.4

echo "Input recieved"
echo "Run: " $RUN
echo "Angle: " $ANGLE

gear="gear_ALIGNED/gear_150mm_1fei4_$ANGLE.xml"
echo $gear
outputfile="output/logs/trackdaf-00$RUN.zip"
echo $outputfile

tempgear="gear_temp_$ANGLE.xml"

tempgearhead="gear_temp_$ANGLE-head.xml"
tempgeartemptail="gear_temp_$ANGLE-temptail.xml"
tempgeartail="gear_temp_$ANGLE-tail.xml"
tempgearmid="gear_temp_$ANGLE-mid.xml"
tempgearline="gear_temp_$ANGLE-line.xml"

maxTracks="0"
bestX=0.
bestY=0.
nTracks="0"
loop="0"
converged="0"

head -n 71 $gear > $tempgearhead
tail -n 73 $gear > $tempgeartemptail
head -n 5 $tempgeartemptail > $tempgearmid
tail -n 67 $tempgeartemptail > $tempgeartail

while [ $converged -eq "0" ]
do
	echo "begin loop"
	if [[ $Step -le 0.001 ]]
	then
		converged="1"
	fi
	for (( x=$BeginX; x<=$EndX; x+=$Step ))
	do
		for (( y=$BeginY; y<=$EndY; y+=$Step ))
		do
			echo "							positionX=\"$x\"	positionY=\"$y\"positionZ=\"382.118\" " > $tempgearline
			cat $tempgearhead > $tempgear
			cat $tempgearline >> $tempgear
			cat $tempgearmid >> $tempgear
			cat $tempgearline >> $tempgear
			cat $tempgeartail >> $tempgear
			mv $tempgear $gear

			echo "gearfile edited: $x, $y"

			jobsub -c config.cfg -csv runlist-1fei4.csv -s hitmaker $RUN
			jobsub -c config.cfg -csv runlist-1fei4.csv -s trackdaf $RUN

			nTracks=`unzip  -p  $outputfile |grep "TTree with" |cut -d ' ' -f6`; 
			echo "$nTracks tracks created at $x, $y"

			if [ $nTracks -gt $maxTracks ]
			then
				maxTracks=$nTracks
				bestX=$x
				bestY=$y
			fi
		done
	done
	BeginX=$(( $bestX - $Step ))
	EndX=$(($bestX + $Step))
	BeginY=$(($bestY - $Step))
	EndY=$(($bestY + $Step))
	Step=$(($Step * 0.4))
	echo "step: $Step"
done
output="bestAligns/bestalign$RUN.txt"
echo $bestX > $output
echo $bestY >> $output
echo $bestX
echo $bestY

echo "							positionX=\"$bestX\"	positionY=\"$bestY\"	positionZ=\"382.118\" " > $tempgearline
cat $tempgearhead > $tempgear
cat $tempgearline >> $tempgear
cat $tempgearmid >> $tempgear
cat $tempgearline >> $tempgear
cat $tempgeartail >> $tempgear
mv $tempgear $gear

rm $tempgearhead
rm $tempgeartemptail
rm $tempgeartail
rm $tempgearmid
rm $tempgearline

