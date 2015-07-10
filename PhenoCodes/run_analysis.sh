#!/bin/bash
#Este código realiza el análisis con el archivo en c++ creado para este propósito para todos los ParamCards en paralelo
CONTAR=1
#Ciclo en las ParamCards
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
  echo "Esta es la carpeta $CONTAR"
  COUNTER=1
  #Ciclo en las ntuplas existentes
  for j in $(ls ${i/.dat/}); do
    #Chequea si ya se realizó el análisis
    if test -e ${i/.dat/}/${i/.dat/}_run_${COUNTER}_output.root; then
    echo Ya existe esta simulación ${i/.dat/}_$COUNTER
    else 
      #Copia el archivo root a la carpeta donde está el archivo en c++ (PhenoAnalyzer.cc)
      cp ${i/.dat/}/*run_$COUNTER.root /home/jgodoy/PhenoAnalyzer/DelphesAnalyzer_staus
      #Se entra a la carpeta donde se encuetra el .cc
      cd /home/jgodoy/PhenoAnalyzer/DelphesAnalyzer_staus
      #Se ejecuta el análisis en paralelo
      ./PhenoAnalyzer ${i/.dat/}_run_$COUNTER &
        #Se vuelve a la carpeta de las ntuplas
      cd /home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes
    fi
    let COUNTER=$COUNTER+1
  done
  let CONTAR=$CONTAR+1
done
wait


#Se obtienen los archivos del análisis realizado. Escribir la carpeta donde está el archivo .cc
cd /home/jgodoy/PhenoAnalyzer/DelphesAnalyzer_staus
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
   echo copiando de $i
   mv ${i/.dat/}*_output.root /home/jgodoy/PhenoThesisProject/Simulations/ntuple_delphes/${i/.dat/}
done
rm *.root

echo "Análisis finalizado"
