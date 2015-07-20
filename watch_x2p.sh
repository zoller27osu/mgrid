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

#old_size=$(wc -c < ${out_file})
old_size=$(ls -n1 ${out_file} | awk '{print $5}')
while : ; do
    sleep 30s
    if [ ! -f ${out_file} ]; then
        # file was deleted, so...
        echo "Execution of job $1 finished normally!"
        exit
    fi
    new_size=$(ls -n1 ${out_file} | awk '{print $5}')
    if [ ${new_size} -eq ${old_size} ]; then
        echo "Execution stalled! Deleting job $1."
        qdel $1
        exit
    fi
    old_size=${new_size}
done
