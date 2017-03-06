#!/bin/bash

############################################## GT ##############################################

#aCSM-ALL
nohup perl aCSM.pl listas/gt6.txt resultados/aCSM-ALL/gt6.csv 0.2 30.0 2 > logs/gt6_aCSM-ALL.log &
nohup perl aCSM.pl listas/gt6.5.txt resultados/aCSM-ALL/gt6.5.csv 0.2 30.0 2 > logs/gt6.5_aCSM-ALL.log &
nohup perl aCSM.pl listas/gt7.txt resultados/aCSM-ALL/gt7.csv 0.2 30.0 2 > logs/gt7_aCSM-ALL.log &
nohup perl aCSM.pl listas/gt7.5.txt resultados/aCSM-ALL/gt7.5.csv 0.2 30.0 2 > logs/gt7.5_aCSM-ALL.log &
nohup perl aCSM.pl listas/gt8.txt resultados/aCSM-ALL/gt8.csv 0.2 30.0 2 > logs/gt8_aCSM-ALL.log &

#aCSM-tradicional
nohup perl aCSM.pl listas/gt6.txt resultados/aCSM-tradicional/gt6.csv 0.2 30.0 0 > logs/gt6_aCSM-tradicional.log &
nohup perl aCSM.pl listas/gt6.5.txt resultados/aCSM-tradicional/gt6.5.csv 0.2 30.0 0 > logs/gt6.5_aCSM-tradicional.log &
nohup perl aCSM.pl listas/gt7.txt resultados/aCSM-tradicional/gt7.csv 0.2 30.0 0 > logs/gt7_aCSM-tradicional.log &
nohup perl aCSM.pl listas/gt7.5.txt resultados/aCSM-tradicional/gt7.5.csv 0.2 30.0 0 > logs/gt7.5_aCSM-tradicional.log &
nohup perl aCSM.pl listas/gt8.txt resultados/aCSM-tradicional/gt8.csv 0.2 30.0 0 > logs/gt8_aCSM-tradicional.log &

############################################## MOD ##############################################

#aCSM-ALL
nohup perl aCSM.pl listas/mod6.txt resultados/aCSM-ALL/mod6.csv 0.2 30.0 2 > logs/mod6_aCSM-ALL.log &
nohup perl aCSM.pl listas/mod6.5.txt resultados/aCSM-ALL/mod6.5.csv 0.2 30.0 2 > logs/mod6.5_aCSM-ALL.log &
nohup perl aCSM.pl listas/mod7.txt resultados/aCSM-ALL/mod7.csv 0.2 30.0 2 > logs/mod7_aCSM-ALL.log &
nohup perl aCSM.pl listas/mod7.5.txt resultados/aCSM-ALL/mod7.5.csv 0.2 30.0 2 > logs/mod7.5_aCSM-ALL.log &
nohup perl aCSM.pl listas/mod8.txt resultados/aCSM-ALL/mod8.csv 0.2 30.0 2 > logs/mod8_aCSM-ALL.log &

#aCSM-tradicional
nohup perl aCSM.pl listas/mod6.txt resultados/aCSM-tradicional/mod6.csv 0.2 30.0 0 > logs/mod6_aCSM-tradicional.log &
nohup perl aCSM.pl listas/mod6.5.txt resultados/aCSM-tradicional/mod6.5.csv 0.2 30.0 0 > logs/mod6.5_aCSM-tradicional.log &
nohup perl aCSM.pl listas/mod7.txt resultados/aCSM-tradicional/mod7.csv 0.2 30.0 0 > logs/mod7_aCSM-tradicional.log &
nohup perl aCSM.pl listas/mod7.5.txt resultados/aCSM-tradicional/mod7.5.csv 0.2 30.0 0 > logs/mod7.5_aCSM-tradicional.log &
nohup perl aCSM.pl listas/mod8.txt resultados/aCSM-tradicional/mod8.csv 0.2 30.0 0 > logs/mod8_aCSM-tradicional.log &
