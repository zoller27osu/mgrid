#!/bin/bash

if [ -z "${1}" ]; then
    echo "Usage: $0 pbs_jobid"
    exit
fi

out_file=${1}.OU
while [ ! -f ${out_file} ]; do
    echo "Waiting for execution to begin."
    sleep 10s
done
echo
echo "Execution of job $1 has begun."

old_size=$(wc -c <"${out_file}")
while : ; do
    sleep 10s
    new_size=$(wc -c <"${out_file}")
    if [[ ${new_size}-${old_size} <= 0 ]]; then
        if [[ ${new_size}-${old_size} == 0 ]]; then
            echo "Execution stalled! Deleting job $1."
            qdel $(PBS_JOBID)
        fi # else, file was deleted - meaning execution finished naturally.
        break
    fi
    old_size=${new_size}
done
