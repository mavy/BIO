import sys
import tabix
import utils.file as futils
from pyliftover import LiftOver

lo = LiftOver('hg38', 'hg19')

if len(sys.argv) >= 4:
    fenotypeData = sys.argv[1]
    dataDir = sys.argv[2]
    outputFile = sys.argv[3]
else:
    print "Uso : <fenotypeData> <dataDir> <outputFile>"

    sys.exit(-1)


def filter(fenotype, datadir, mapping):
    # 10	rs10995190	rs10995190-?	0.84	estrogen-receptor negative breast cancer	62518923	  --> 10      64278682        rs10995190      G       A

    results = dict()

    for feno in fenotype:
        if len(feno) < 1:
            continue
        if len(feno) < 4:  # Hay algun caso que se codifica el cromosoma y la posicion como chrx:posicion
            try:
                chromosome, position = feno[0].split(':')
                chromosome = chromosome.lower().split('chr')[1]
            except:
                continue  # There are some special cases, avoid them
        else:
            position = feno[4]
            chromosome = feno[3]

        if 'x' in chromosome or chromosome == '' or ';' in chromosome:
            continue
        chromosome = chromosome.upper()

        if chromosome == 'X':
            url = datadir + 'ALL.chr' + chromosome + '.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
        elif chromosome == 'Y':
            url = datadir + 'ALL.chr' + chromosome + '.phase3_integrated_v2a.20130502.genotypes.vcf.gz'
        else:
            url = datadir + 'ALL.chr' + chromosome + '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

        tb = tabix.open(url)
        try:
            coordinate = lo.convert_coordinate('chr' + chromosome, int(position))
            records = tb.query(chromosome, coordinate[0][1] - 1, coordinate[0][1])
        except:
            continue


        for record in records:
            if record[2] == feno[0]:  # Ignoramos los ssnp sin la misma id
                # Mostramos la mutacion
                lista_alelos = list()
                ones = 0
                zeros = 0
                i = 0

                for alelo in record[9:]:
                    if listInd[i] not in results:
                        results[listInd[i]] = list()
                    if '1' in alelo:
                        lista_alelos.append(1)
                        ones += 1
                        results[listInd[i]].append((record[2], 1))
                    else:
                        lista_alelos.append(0)
                        zeros += 1
                        results[listInd[i]].append((record[2], 0))
                    i += 1
    return results


# Leemos los identificadores de cada individuo

listInd = futils.llegir_order("individuosOrder.txt")
# Leemos el orden del cromosoma Y (hombres)
listIndY = futils.llegir_order("individuosOrder_Y.txt")

# Leemos el mapping individuo -> poblacion -> superpoblacion -> genero
mapping = futils.llegir_basic_fenotype(dataDir + "integrated_call_samples_v3.20130502.ALL.panel")

fenotypeList = futils.llegir_selected_fenotype(open(fenotypeData))

results = filter(fenotypeList, dataDir, mapping)

output = open(outputFile, "w")

for result in results:
    for (ssnp, value) in results[result]:
        # print result, ssnp, value, mapping[result][0], mapping[result][1], mapping[result][2]
        output.write(str(
            result + "\t" + ssnp + "\t" + str(value) + "\t" + mapping[result][0] + "\t" + mapping[result][1] + "\t" +
            mapping[result][2] + "\n"))

output.close()

# print summarizedTraits
