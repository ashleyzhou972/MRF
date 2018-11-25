#!/bin/bash
alias myRscript='/work/friedberg_lab/nzhou/configure/R_make/bin/Rscript'
OUT="FALSE"
round=2
param=3
while [ "$OUT" == "FALSE" ]; do
	if [[ $round -eq 3 ]] ; then
		break	
	fi
	echo "Round $round"
	outputfile=out_round"$round"_param"$param" 
	myRscript tuning_multiple_chains_1.R "$round"  "$param" >$outputfile
	#OUT='FALSE'
	#remove trailing whitespace
	OUT="$(tail -1 $outputfile | tr -d '[:space:]')"
	#echo $OUT
	round=$((round+1))	
	#echo $round
done
	

