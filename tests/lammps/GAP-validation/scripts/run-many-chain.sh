#!/bin/bash

# NOTE pass in stride as argument 1, type as argument 2,  then the rest of the filenames after

stride="$1"
_type="$2"

module r QUIP_v3-2
module unload ASE

for ((i = 3; i <= $#; i++ )); do
    cd "${!i}/${_type}"
   
    # Delete old files
    /bin/rm *_chains_*txt

    printf 'Working on dir: %s\n' $(pwd)
    python ../../scripts/chain-analysis.py dump.3000K-equil.nc --stride ${stride}
    echo ""
    cd ../../
done

echo "Finished"
