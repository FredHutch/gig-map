#!/bin/bash

# Make sure that there are no duplicate gene names in the input FASTA.

declare -A seen_names

while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        name=$(echo "$line" | cut -d' ' -f1 | cut -c2-)
        if [[ -v seen_names["$name"] ]]; then
            seen_names["$name"]=$((seen_names["$name"] + 1))
            name="${name}_${seen_names["$name"]}"
        else
            seen_names["$name"]=1
        fi
        echo ">$name"
    else
        echo "$line"
    fi
done < /dev/stdin
