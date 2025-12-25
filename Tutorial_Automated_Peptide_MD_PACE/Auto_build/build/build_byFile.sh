#!/bin/bash

curdir=$(pwd)
pdir=$curdir/../

# Check if file is provided as argument
if [ $# -eq 0 ]; then
    echo "Please provide a text file as argument"
    exit 1
fi


cd $pdir
tgtdir=$pdir/setup


# Input file
input_file=$curdir/"$1"

# Check if file exists
if [ ! -f "$input_file" ]; then
    echo "File $input_file not found"
    exit 1
fi




# Read file line by line
while IFS= read -r word; do
    # Skip empty lines
    [ -z "$word" ] && continue
    
    # Do something with the word
    # Example: Print the word with some modification
    echo "Processing peptide: $word"
    cd $curdir    
    # Add your custom action here
    # For example:
    # - Convert to uppercase: echo "${word^^}"
    # - Check length: echo "${#word}"
    # - Search in another file: grep "$word" other_file.txt
    # - Make API call: curl "http://example.com/api/$word"

    bash prepare.sh $word 1 yes b
    mv $word-pace.top $tgtdir/"$word".top
    mv $word-pace-resid.pdb $tgtdir/"$word".pdb
    rm '#'*
    rm *"$word"*

    cd $tgtdir
    bash build_paceasm.sh $word
    
done < "$input_file"

