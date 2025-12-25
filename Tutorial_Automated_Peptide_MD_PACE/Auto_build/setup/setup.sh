#!/bin/bash


# Store the input string
string="AAE AAK"

# Loop over each word in the string
for word in $string; do
    # Example action: Convert word to uppercase
# build and sim
#    bash sim_paceasm.sh $word simNoNeutr330
# build only
    bash build_paceasm.sh $word

    # Add your custom action here
done
