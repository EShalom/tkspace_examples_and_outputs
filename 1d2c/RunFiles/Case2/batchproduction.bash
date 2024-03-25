#!/bin/bash

# specify the list of run txt filenames
run="placeholder_value_files/batch_list.txt"

# specify the list of placeholders used for values
phs="placeholder_value_files/placeholders.txt"

# specify script with placeholders within
out="template_case2.py"

# extract placeholders in bash array
while IFS=\= read info
do
    places+=($info)
done < "$phs"

# read every line of the batch list
# this has the value associated with each placeholder
while IFS= read -r file
do
	# set the current line being scanned as output variable
	output="placeholder_value_files/$file"
	# set counter to loop thorugh placeholders
	COUNTER=0

	# copy scipt with placeholders and rename to specifc run
	cp -v "$out" "batch_python_files/${file%.txt}"".py"

	# read through each placeholder value in the run info .txt
	while IFS= read -r line
	do
		# run stream editor to replace each placeholder with corresponding value
		sed -i '' 's|'"${places[COUNTER]}"'|'"$line"'|g' "batch_python_files/${file%.txt}"".py"

		# increment counter
		((COUNTER++))
	done < "$output"
done < "$run"

