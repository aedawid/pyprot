#!/bin/bash


for i in `cat list`
do
    /Users/adawid/PROJECTS/repos/MSA/hmmer-3.2/src/phmmer -o $i-score_pair -A $i-MSA.aln --tblout $i-MSA_scores FASTA/$i.fasta /Users/adawid/Box/GROUP_data/LLPS-bioinformatics/databases/SEQ_DB/uniprot_swiss-prot.fasta
    /Users/adawid/PROJECTS/repos/MSA/hmmer-3.2/src/hmmbuild $i.profile $i-MSA.aln
    k=`cat $i.profile | grep "NSEQ" | awk '{print $2}'`
    touch $i-$k
    /Users/adawid/PROJECTS/repos/MSA/hmmer-3.2/src/hmmlogo $i.profile > $i.logo
    /Users/adawid/PROJECTS/repos/MSA/hmmer-3.2/src/hmmlogo --height relent abovebg $i.profile > $i.logo2
done