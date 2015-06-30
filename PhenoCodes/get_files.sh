#!/bin/bash
#Se obtienen los archivos del análisis realizado. Escribir la carpeta donde está el archivo .cc
cd /home/jgodoy/PhenoAnalyzer/DelphesAnalyzer_staus
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
   echo copiando de $i
   mv ${i/.dat/}*_output.root /home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/${i/.dat/}
done
rm *.root
