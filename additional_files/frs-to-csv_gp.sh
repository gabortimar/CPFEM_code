#! /bin/bash

for fname in *.frs
do

#echo $fname
#extension="${fname##*.}"
filename="${fname%.*}"

#echo $extension
#echo $filename

csvmkr <<STDIN -o other --options > /dev/null
$filename
2
$filename
STDIN


done
