#!/bin/bash
mkdir -p fasta
for f in fasta_compressed/*;
do
    tar xzf $f -C fasta/
done

mkdir -p get_objects
for f in get_objects_compressed/*;
do
    tar xzf $f -C get_objects
done

mkdir -p dRep_dir
for f in dRep_dir_compressed/*;
do
    tar xzf $f -C dRep_dir
done
