#!/usr/bin/env bash

job_name=elastic_gcc

for n in 1 2 4 9 18 36 72 144 288 432 576 648 720; do
    if [ $n -le 72 ]; then
        QOS="--qos serial"
        node=1
        npn=$n
    else
        QOS="--qos parallel"
        node=$(( $n / 72 ))
        npn=72
    fi
    sbatch -J ${job_name} -e ${job_name}_%A.err -o ${job_name}_%A.out -N ${node} --ntasks-per-node ${npn} $QOS "$@" --exclusive ./run.sh
done
