#!/bin/bashi
#Este archivo crea las muestras de madgraph y las ntuplas de Delphes con los códigos creados para dichas funciones.
./home/jgodoy/PhenoThesisProject/Simulations/mg_files/create_sim.sh
./home/jgodoy/PhenoThesisProject/Simulations/mg_files/run_sim.sh
./home/jgodoy/PhenoThesisProject/Simulations/mg_files/clean_sim.sh
./home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/ntuple_creator.sh
./home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/organizer.sh
./home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/run_analysis.sh
./home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/get_files.sh
