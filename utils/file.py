# Llegim el fitxer de fenotypes de la fase 1
def llegir_selected_fenotype(file_object):
    output = []
    for linea in file_object:
        list_terms = []
        for value in linea.split('\t'):
            if len(value.strip()) != 0:
                list_terms.append(value.strip())
        output.append(list_terms)
    file_object.close()
    return output


def llegir_order(filename):
    with open(filename) as f:
        list = [word for line in f for word in line.split()]
        return list


def llegir_basic_fenotype(filename):
    # Leemos el mapping individuo -> poblacion -> superpoblacion -> genero
    mapping = dict()
    with open(filename) as f:
        for (id, pob, superpop, gender) in (line.split() for line in f):
            mapping[id] = (pob, superpop, gender)
    mapping.pop("sample")
    return mapping


# Llegim el fitxer que esta separat per tabuladors
# Molts snp estan repetits pero tenen un MAPPED_TRAIT diferent
# Farem un concat del text
def llegir_GWAS(file):
    header = []
    output = []
    keys = ("CHR_ID", "CHR_POS", "STRONGEST SNP-RISK ALLELE", "MAPPED_TRAIT", "SNPS")
    for linea in file:
        if len(header) == 0:
            for camp in linea.split('\t'):
                header.append(camp.strip())
        else:
            i = 0
            tdict = {"CHR_ID": 0, "CHR_POS": 0, "STRONGEST SNP-RISK ALLELE": 0, "MAPPED_TRAIT": 0, "SNPS": 0}
            for value in linea.split('\t'):
                if header[i] in keys:
                    tdict[header[i]] = value
                i += 1
            output.append(tdict.copy())
    snps = dict()
    for entrada in output:
        if entrada["SNPS"] not in snps:
            snps[entrada["SNPS"]] = {"CHR_ID": entrada["CHR_ID"], "CHR_POS": entrada["CHR_POS"], "STRONGEST SNP-RISK ALLELE": entrada["STRONGEST SNP-RISK ALLELE"], "MAPPED_TRAIT": entrada["MAPPED_TRAIT"], "SNPS" : entrada["SNPS"]}
        else:
            snps[entrada["SNPS"]]["MAPPED_TRAIT"] = snps[entrada["SNPS"]]["MAPPED_TRAIT"] + " " + entrada["MAPPED_TRAIT"]

    output = snps.values()
    return output