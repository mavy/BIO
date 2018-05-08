# coding=utf-8
import sys

import utils.file as futils

if len(sys.argv) >= 3:
    fileInput = sys.argv[1]
    command = sys.argv[2]
else:
    print "Uso : <ficheroGWAS.tsv> <mode>"
    print "mode : MostraKeywords"
    print "mode : Selecciona <keyword1> <keyword2> ..."
    sys.exit(-1)


# Busca el mapped trait dentro del diccionario de SNPs, devuelve el elemento
def search(mapped_trait, list_dic):
    return [element for element in list_dic if element['MAPPED_TRAIT'].lower() in mapped_trait]


# Guarda los fenotipos (sin repetidos) del diccionario de SNPs
def fenotypes(list_dic):
    return set([element['MAPPED_TRAIT'] for element in list_dic])


# Guarda todos los fenotipos (set)
def selectRegex(fenotypes, list_fenotypes):
    output = set()
    for feno in list_fenotypes:
        for regex in fenotypes:
            if regex.lower() in feno.lower():
                output.add(feno.lower())
    return list(output)


# Returns a unique, capitalized and sorted list of words found in the input (list of lines)
def keywords(fenotypes):
    output = set()
    for line in fenotypes:
        for word in line.split():
            if word[-1] == ',':
                word = word[0:-1]
            output.add(word.capitalize())
    return sorted(output)


output = futils.llegir_GWAS(open(fileInput))

traits = fenotypes(output)

if command == "MostraKeywords":
    summarizedTraits = keywords(traits)
    print summarizedTraits
elif command == "Selecciona":
    if len(sys.argv) > 3:
        listKeywords = sys.argv[3:]
        extendedTraits = selectRegex(listKeywords, traits)

        for lines in search(extendedTraits, output):
            cadena = ""
            for item in lines.itervalues():
                if item is None:
                    item = "NA"
                cadena += item + "\t"
            print cadena
    else:
        print "Falten paraules claus"

# print summarizedTraits
