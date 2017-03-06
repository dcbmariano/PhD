#!/bin/bash

############################################## GT ##############################################

#aCSM-ALL
i=1
j=1
while [ $i -lt 30 ]; do
	echo "Running cutoff $i"
	let i=i+1
	j=1
	while [ $j -lt 10 ]; do
		nohup perl aCSM.pl listas/mutantes.txt resultados/aCSM-ALL/mutantes_step_0_"$j"_cutoff_"$i".csv 0."$j" "$i".0 2 > logs/mutantes_step_0_"$j"_cutoff_"$i"_aCSM-ALL.log &
		let j=j+1
	done
done
