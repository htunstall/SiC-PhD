#!/bin/bash

touch check_mem.var

while [ -f check_mem.var ]; do
    let new_memory=$(free -m | awk 'NR==2{printf $3}')
    printf '%s,%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${new_memory}" >> time-memory.log
    sleep 5
done
