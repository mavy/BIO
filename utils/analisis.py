# Devuelve un diccionario [id] = lista  ( ssnp, value ) (pop, superpop y gender se encuentran en el mapping)
def compute(inp):
    result = dict()
    for (id, ssnp, value, pop, superpop, gender) in (line.split() for line in inp):
        if id not in result:
            result[id] = list()
        result[id].append((ssnp, value))
    return result


# Devuelve el numero de individuos con todos los fenotipos activados,los que no, y una lista con sus identificadores
# Entrada: Lista de [ID, {ssnp, valor}]
# Salida : positivos, negativos, set{ID}
def analisis_simultani(results):
    negative = 0
    positive = 0

    positives = set()
    for result in results.keys():
        found = True
        for (ssnp, value) in results[result]:
            if int(value) == 0:  # Es importante el cast a int
                negative = negative + 1
                found = False
                break
        if found:
            positive = positive + 1
            positives.add(result)
    return positive, negative, positives
