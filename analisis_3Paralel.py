# coding=utf-8
import cPickle as pickle
import operator
import sys
import tabix

from pyliftover import LiftOver

import utils.analisis as analisis
import utils.file as futils
import multiprocessing
from multiprocessing import Manager

FAST = True
FAST = False

lo = LiftOver('hg38', 'hg19')

# Dados una serie de fenotipos seleccionados, y utilizando los datos de 1000 genome encontrar otros fenotipos que suelen
#  aparecer juntos. Para realizarlo se filtrarán los datos de 1000 genomes con los fenotipos seleccionados,
# y se relacionarán el resto de fenotipos activos.
# Se mostrarán los fenotipos relacionados por probabilidad de aparición.


# Paso 1: filter_1000Genome nos da la lista de individuos con los fenotipos activados.
# Paso 2: Coger cada individuo, y ver que otros ssnps estan activados (GWAS + filter1000)
# Paso 2a: Leer todos los ssnps que faltan (All - fenotypeData)
# Paso 2b: De cada ssnp, mirar si el individuo lo tiene o no activado [Podemos optimizar]
# Paso 3: Sumar cada rs + 1 por activación | Separados por Pob ?
# Paso 4: Mostrar ranking.

if len(sys.argv) >= 8:
    fenotypeData = sys.argv[1]
    datadir = sys.argv[2]
    inputFile = sys.argv[3]
    GWASFile = sys.argv[4]
    number = sys.argv[5]
    procesos = int(sys.argv[6])
    sort = int(sys.argv[7])
else:
    print "Uso : <fenotypeData> <dataDir> <inputFile> <GWASFile> <elementos_a_tratar> <procesos>"

    sys.exit(-1)

# Paso 2b: De cada ssnp, mirar si el individuo lo tiene o no activado
# Paso 3: Sumar cada rs + 1 por activación | Separados por Pob
def contar_otros_snps(out_queue, fenotypeList, all_fenotypes, simultaneous, listInd, listIndY, number, index, total):
    SSNPsList = [fenotype[0] for fenotype in fenotypeList]
    humans = dict()
    humansY = dict()

    for human in simultaneous:
        humans[human] = listInd.index(human)
        try:
            humansY[human] = listIndY.index(human)
        except:
            None

    results = dict()
    # Recorremos por fenotipo, siempre habran mas que humanos
    treated = -1
    for fenotype in all_fenotypes:
        treated = treated + 1
        if treated == (number / total) * (index + 1):
            break
        if treated < (number / total) * (index):

            continue

        if fenotype['SNPS'] not in SSNPsList:
            # Debemos buscarlos

            if ':' in fenotype['SNPS']:
                try:
                    chromosome, position = fenotype['SNPS'].split(':')
                    chromosome = chromosome.lower().split('chr')[1]
                    position = position.split('_')[0]  # Casos especiales
                    position = position.split('-')[0]
                except:
                    continue
            else:
                position = fenotype['CHR_POS']
                chromosome = fenotype['CHR_ID']

            if 'x' in chromosome or chromosome == '' or ';' in chromosome:
                continue
            chromosome = chromosome.upper()

            if chromosome == 'X':
                url = datadir + 'ALL.chr' + chromosome + '.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
            elif chromosome == 'Y':
                url = datadir + 'ALL.chr' + chromosome + '.phase3_integrated_v2a.20130502.genotypes.vcf.gz'
            else:
                url = datadir + 'ALL.chr' + chromosome + '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

            coordinate = lo.convert_coordinate('chr' + chromosome, int(position))
            if coordinate is None or len(coordinate) == 0:
                continue

            tb = tabix.open(url)
            records = tb.query(chromosome, coordinate[0][1] - 1, coordinate[0][1])
            results[fenotype['SNPS']] = 0
            for record in records:
                if record[2] == fenotype['SNPS']:
                    list = ['1' in alelo for alelo in record[9:]]
                    if chromosome == 'Y':
                        for index in humansY.values():
                            if list[index]:
                                results[fenotype['SNPS']] += 1
                    else:
                        for index in humans.values():
                            if list[index]:
                                results[fenotype['SNPS']] += 1
            #print len(results), '/', len(all_fenotypes), treated

    out_queue.put(results)
    return results


def encontrar_descripcion(fenotypes, rs):
    for feno in fenotypes:
        if feno['SNPS'] == rs:
            return feno['MAPPED_TRAIT']


def merge_two_dicts(x, y):
    z = x.copy()  # start with x's keys and values
    for element in y:
        if element in z:
            z[element] = z[element]+y[element]
        else:
            z[element] = y[element]
    return z


listInd = futils.llegir_order("individuosOrder.txt")
# Leemos el orden del cromosoma Y (hombres)
listIndY = futils.llegir_order("individuosOrder_Y.txt")

# Leemos el mapping individuo -> poblacion -> superpoblacion -> genero
mapping = futils.llegir_basic_fenotype(datadir + "integrated_call_samples_v3.20130502.ALL.panel")

fenotypeList = futils.llegir_selected_fenotype(open(fenotypeData))

r = analisis.compute(open(inputFile))

simultani = analisis.analisis_simultani(r)

all_fenotypes = futils.llegir_GWAS(open(GWASFile))
# En simultani tenemos todos los individuos con los fenotipos activados simultaneamente
if sort==1:
    all_fenotypes = sorted(all_fenotypes, key=lambda tup: tup['CHR_ID'])

# Buscamos una lista de SNPs y ocurrencias
if FAST:
    results = pickle.load(open("analisis3_results.p", "rb"))
else:
    #results = contar_otros_snps(fenotypeList, all_fenotypes, simultani[2], listInd, listIndY)
    #pickle.dump( results, open( "analisis3_results.p", "wb" ))
    if __name__ == '__main__':
        jobs = []
        out_queue = multiprocessing.Queue()
        for i in range(procesos):
            p = multiprocessing.Process(target=contar_otros_snps, args=(
            out_queue, fenotypeList, all_fenotypes, simultani[2], listInd, listIndY, int(number), i, procesos))
            jobs.append(p)
            p.start()

        results = dict()
        for i in range(procesos):
            a = out_queue.get()
            results = merge_two_dicts(results, a)


results = sorted(results.items(), key=operator.itemgetter(1), reverse=True)

results2 = [(a, b) for (a, b) in results if b == len(simultani[2])]

desc = set()
for snp, value in results2:
    d = encontrar_descripcion(all_fenotypes, snp)
    if d is not None:
        desc.add((d,snp))

for d in desc:
    print d
