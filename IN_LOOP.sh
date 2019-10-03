#!/bin/bash

for i in `cat list`
do
    p=`echo $i | cut -c 1-4`
    ch=`echo $i | cut -c 5`
#    gunzip $p.pdb.gz
    python ap_select_chains_from_pdb.py -f $p.pdb -o $p$ch.out -m '1' -c "$ch" > $p$ch.pdb
done