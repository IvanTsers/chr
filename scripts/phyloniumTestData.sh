#!/bin/bash

mkdir ../data/o

for case in identical tooDistant s1q1 s1q2a s1q2b s2q1a s2q1b s2q2 onEnds mutations
    do
	phylonium -p ../data/o/"$case".fasta -r ../data/s/"$case".fasta ../data/s/"$case".fasta ../data/q/"$case".fasta 2> /dev/null
    done
