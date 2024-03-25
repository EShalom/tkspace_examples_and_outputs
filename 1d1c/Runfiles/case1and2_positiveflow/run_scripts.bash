#!/bin/bash

# specify the list of scipt names to run
run_files="scriptrun_list.txt"

# read every line of the list to run
while IFS= read -r line
do
	# set each line to run with appropriately named log
	# tag & at end very important as moves process to background
	nohup python3 -u "$line" >& "${line%.py}"".log" & 

done < "$run_files"
