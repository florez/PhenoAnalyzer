#!/bin/bash
#Este código crea las condiciones para la simulación en MadGraph. 
#Se crea un contador de prueba para verificar que se creen todas las simulaciones. Se ubica en la carpeta donde deben quedar las carpetas de las simulaciones.
COUNTER=0
echo El contador es $COUNTER
#Se realiza un ciclo que recorre todos las ParamCards. Acá se debe escribir la dirección donde se encuentran las ParamCards
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards); do
#Este comando ejecuta MadGraph. Se debe escribir la dirección de el ejecutable de MadGraph. El tomando toma como entrada un archivo "mgcode" que debe ser modificado de acuerdo a la colisión que se quiera emular. El archivo importa el modelo "import mssm", luego genera la simulación (ej: generate p p > ta1+ ta1-), crea un output ("output auto") y luego sale de Madgraph ("exit"). El símbolo "&" ejecuta madgraph en paralelo.
./home/jgodoy/PhenoThesisProject/Programs/MG5_aMC_v2_2_3/bin/mg5_aMC /home/jgodoy/PhenoThesisProject/Simulations/mg_files/mgcode &
let COUNTER=COUNTER+1
#Se imprime en pantalla las simulaciones creadas
echo El contador es $COUNTER
done 

