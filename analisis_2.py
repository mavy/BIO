# coding=utf-8
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import utils.analisis as analisis
import utils.file as futils

# Existe relación entre unos determinados fenotipos y la poblacion (o superpoblacion).
# Es decir, hay una proporción
# significativa o una distribución no equitativa entre el fenotipo y el origen del individuo.

# Se muestra el snps en forma gráfica.

if len(sys.argv) >= 5:
    fenotypeData = sys.argv[1]
    dataDir = sys.argv[2]
    inputFile = sys.argv[3]
    outputFile = sys.argv[4]
else:
    print "Uso : <fenotypeData> <dataDir> <inputFile> <outputFile>"

    sys.exit(-1)


# Contamos los individuos con algún snp a 1, guardamos la entrada como [snp] = ( pob )
def analisis_sim(results):
    positives = dict()
    for result in results.keys():
        for (ssnp, value) in results[result]:
            if ssnp not in positives:
                positives[ssnp] = dict()
            if mapping[result][0] not in positives[ssnp]:
                positives[ssnp][mapping[result][0]] = 0
            if int(value) == 1:  # Es importante el cast a int
                positives[ssnp][mapping[result][0]] += 1
    return positives


# Cuenta cuantos individuos de cada fenotipo tenemos (genero, población y superpoblación)
def calcular_totales(mapping):
    total = dict()
    for (pop, super_pop, m_gender) in mapping.values():
        if pop not in total:
            total[pop] = 0
        if super_pop not in total:
            total[super_pop] = 0
        if m_gender not in total:
            total[m_gender] = 0
        total[pop] += 1
        total[super_pop] += 1
        total[m_gender] += 1
    return total


# Con los totales calculados por calcular_totales obtenemos la probabilidad de cada snp/población
def calcular_probabilidad(totales, positives):
    result = dict()
    for snp in positives:
        result[snp] = dict()
        for poblation in positives[snp]:
           result[snp][poblation] = (float(positives[snp][poblation]) / totales[poblation]) * 100.0
    return result


# Devolvemos la descripción (Mapped_trait) para el snp de entrada.
def encontrar_descripcion(fenotypes, rs):
    for feno in fenotypes:
        if feno[0] == rs:
            return feno[1]


# Generamos un pdf con las gráficas, este es el proceso más costoso ya que pueden haber muchos snps.
def generar_graficas(prob, fenotypeList, filename="analisis_2.pdf"):
    i = 0
    plots = list()
    pp = PdfPages(filename)
    for rs in prob.keys():
        x = prob[rs]
        colors = plt.cm.BuPu([float(val) / 100.0 for val in x.values()])
        # plots.append(plt.figure())
        p = plt.figure()
        plt.bar(x.keys(), x.values(), color=colors)
        plt.ylim((0, 100))
        plt.xticks(rotation=90)
        plt.title(rs + "\n" + encontrar_descripcion(fenotypeList, rs))
        i += 1
        pp.savefig(p)
        plt.close()

   # pp = PdfPages(filename)
   # for plot in plots:
   #     pp.savefig(plot)
    pp.close()


listInd = futils.llegir_order("individuosOrder.txt")
# Leemos el orden del cromosoma Y (hombres)
listIndY = futils.llegir_order("individuosOrder_Y.txt")

# Leemos el mapping individuo -> poblacion -> superpoblacion -> genero
mapping = futils.llegir_basic_fenotype(dataDir + "integrated_call_samples_v3.20130502.ALL.panel")

totales = calcular_totales(mapping)

fenotypeList = futils.llegir_selected_fenotype(open(fenotypeData))

r = analisis.compute(open(inputFile))

# individuos con algún snp activo
positives = analisis_sim(r)

# Convertimos numero de casos en probabilidad con el total de individuos.
prob = calcular_probabilidad(totales, positives)

# Generamos la salida
generar_graficas(prob, fenotypeList, outputFile)
