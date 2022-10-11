#!/bin/bash

# NOTE pass in stride as argument 1, type as argument 2,  then the rest of the filenames after

stride="$1"
_type="$2"

if [ "${_type}" = "3000K" ] ; then
    filename="dump.3000K-equil.nc"
elif [ "${_type}" = "quench_010ps" ] ; then
    filename="dump.300K-equil.nc"
else
    echo "the type '${_type}' is not understood"
    exit 1
fi

module r QUIP_v3-2
module unload ASE

for ((i = 3; i <= $#; i++ )); do
    cd "${!i}/${_type}"
   
    printf 'Working on dir: %s\n' $(pwd)
    python ../../scripts/ring-analysis.py ${filename} --stride ${stride}
    echo ""
    cd ../../
done

echo "Finished"
