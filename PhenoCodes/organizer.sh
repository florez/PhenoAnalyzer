#!/bin/bash
#Mueve las ntuplas creadas dentro de su respectiva carpeta  elimina los archivos .hep. ste archivo debe estar en la carpeta donde están las carpetas de las ntuplas. Ejecutar después de ntuple_creator.sh
for i in $(ls /home/jgodoy/PhenoThesisProject/ParamCards);do
mv ${i/.dat/}*.root ${i/.dat/}
rm ${i/.dat/}/*.hep 
done

