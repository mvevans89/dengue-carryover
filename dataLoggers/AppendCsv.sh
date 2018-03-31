#!/bin/bash
#cd to folder storing all csv files and run this script by typing sh and then this script

OutFileName="All.csv"                       # Fix the output name
i=0                                       # Reset a counter
for filename in ./*.csv; do
awk 'NR=1{print $0 "",FILENAME""}NR>1{print $0 "",""FILENAME}' ""$filename"" > tempfile && mv tempfile ""$filename""  
 if [ "$filename"  != "$OutFileName" ] ;      # Avoid recursion 
 then 
   if [[ $i -eq 0 ]] ; then 
      head -1  $filename >   $OutFileName # Copy header if it is the first file
   fi
   tail -n +2  $filename >>  $OutFileName # Append from the 2nd line each file
   i=$(( $i + 1 ))                        # Increase the counter
 fi
done