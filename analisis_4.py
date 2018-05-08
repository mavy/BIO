# coding=utf-8
import cPickle as pickle
import random
import sys
import tabix

from pyliftover import LiftOver
from sklearn.model_selection import KFold, cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.neural_network import MLPClassifier
import itertools
import utils.file as futils

lo = LiftOver('hg38', 'hg19')

# Podemos entrenar un modelo de ML para que dados un fenotipo nos dé la población más probable. Para ello
# entrenaremos con los datos de SNPs que tienen un fenotipo reconocido como entrada y utilizaremos la población como
# salida.

# Paso 1: Sacar una matriz Individuo - SNP reconocido
# Paso 1a: Se separará Training - Test aleatoriamente (k-fold)
# Paso 2: Usar como entrada en un ML
# Paso 3: Cogeremos una serie de SNPs y nos dará la población más probable...


# Nos permite saltar la fase de creación de la matriz leyendo un fichero pickle
FASTPhase1 = True

if len(sys.argv) >= 3:
    datadir = sys.argv[1]
    GWASFile = sys.argv[2]
else:
    print "Uso : <dataDir> <GWASFile>"

    sys.exit(-1)


# Devuelve listInd / fenotype (true, false)
def create_matrix(mapping, all_fenotypes, listInd):
    SSNPDict = dict()

    for fenotype in all_fenotypes:
        SSNPDict[fenotype['SNPS']] = (fenotype['CHR_ID'], fenotype['CHR_POS'])

    results = dict()

    for human in listInd:
        results[human] = []

    for SNP in SSNPDict.keys():

        # Debemos buscarlos

        if ':' in SNP:
            try:
                chromosome, position = SNP
                chromosome = chromosome.lower().split('chr')[1]
                position = position.split('_')[0]  # Casos especiales
                position = position.split('-')[0]
            except:
                continue
        else:
            chromosome, position = SSNPDict[SNP]

            if 'x' in chromosome or chromosome == '' or ';' in chromosome:
                continue
        chromosome = chromosome.upper()

        if chromosome == 'Y':  # Nos saltamos el chromosoma Y
            continue

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

        for record in records:
            if record[2] == SNP:
                list = ['1' in alelo for alelo in record[9:]]

                for idHuman in listInd:
                    position = listInd.index(idHuman)
                    results[idHuman].append(int(list[position]))
    return results

# Dividimos entre training y test, aunque la división la haremos con k-fold... este método nos preparara
# la lista de etiquetas.
def divide(matrix, mapping, listInd, trainingPer):
    random.seed()
    training = set()
    test = set()
    for i in range(0, len(listInd)):
        value = random.randint(0, 100)
        if value < trainingPer:
            training.add((listInd[i], mapping[listInd[i]]))
        else:
            test.add((listInd[i], mapping[listInd[i]]))
    print "Training ", len(training)
    print "Testing ", len(test)

    resultsTrainingX = []
    resultsTrainingY = []

    resultsTestX = []
    resultsTestY = []

    for (id, (pop, super_pop, genre)) in training:
        resultsTrainingY.append(pop)
        resultsTrainingX.append(matrix[id])

    for (id, (pop, super_pop, genre)) in test:
        resultsTestY.append(pop)
        resultsTestX.append(matrix[id])

    return resultsTrainingX, resultsTrainingY, resultsTestX, resultsTestY


listInd = futils.llegir_order("individuosOrder.txt")
# Leemos el orden del cromosoma Y (hombres)
listIndY = futils.llegir_order("individuosOrder_Y.txt")

# Leemos el mapping individuo -> poblacion -> superpoblacion -> genero
mapping = futils.llegir_basic_fenotype(datadir + "integrated_call_samples_v3.20130502.ALL.panel")

all_fenotypes = futils.llegir_GWAS(open(GWASFile))

if FASTPhase1:
    matrix = pickle.load(open("analisis4_allmatrix.p", "rb"))
else:
    matrix = create_matrix(mapping, all_fenotypes, listInd)
    # pickle.dump(matrix, open("analisis4_allmatrix.p", "wb"))

# No se necesita ya que podemos utilizar k-fold directamente.
(resultTrainingX, resultTrainingY, resultTestX, resultTestY) = divide(matrix, mapping, listInd, 75)

JoinX = resultTrainingX + resultTestX
JoinY = resultTrainingY + resultTestY

clf = MLPClassifier(solver="adam", alpha=1e-5, max_iter=200,
                    hidden_layer_sizes=(400, 300, 200), activation="relu", random_state=1, verbose=True)

results = cross_val_score(clf, JoinX, JoinY, cv=3, n_jobs=-1, verbose=1)


print results

