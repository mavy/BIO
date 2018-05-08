import sys

import utils.analisis as analisis
import utils.file as futils

# Existe poblacionalmente relacion entre el fenotipo "asma" y "cancer de pulmon".
# Es decir, existe una proporcion significativa de individuos con los SNPs correspondientes simultaneamente.


if len(sys.argv) >= 2:
    inputFile = sys.argv[1]
else:
    print "Uso : <inputFile from filter_1000Genome>"

    sys.exit(-1)

# Leemos los identificadores de cada individuo
listInd = futils.llegir_order("individuosOrder.txt")
# Leemos el orden del cromosoma Y (hombres)
listIndY = futils.llegir_order("individuosOrder_Y.txt")
# Devuelve un diccionario [id] = lista  ( ssnp, value ) (pop, superpop y gender se encuentran en el mapping)
r = analisis.compute(open(inputFile))
# Devuelve el numero de individuos con todos los fenotipos activados,los que no, y una lista con sus identificadores
k = analisis.analisis_simultani(r)
# Imprimimos la salida
print k[0], k[1], [item for item in k[2]]
