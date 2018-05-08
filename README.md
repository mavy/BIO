# BIO

Código del TFM : Conocimiento en 1000Genome y GWAS, análisis de rendimiento

Autor : Ramon Nou

### Descripción

* download.py : Descarga los datos de ambas bases de datos a un directorio
* read_GWAS.py : Permite seleccionar SNPs utilizando keywords del campo de Mapped_Treats
* filter_1000Genome.py : Filtra los datos de 1000Genome utilizando SNPs
* filter_1000GenomeParalel.py : Versión paralela del anterior con multiprocessing
* analisis_1.py : Muestra los individuos con todos los SNPs seleccionados activados
* analisis_2.py : Genera una gráfica de la distribución poblacional de los SNPs seleccionados
* analisis_3.py : Genera una lista de los SNPs activados simultaneamente con los seleccionados
* analisis_3Paralel.py : Versión paralela del anterior con multiprocessing
* analisis_4.py : Red neuronal para encontrar la población de un individuo dados sus SNPs
* analisis_4_dt.py : Árbol de decision para encontrar los SNPs más relevantes
* analisis_4_rf.py : Random Forest para encontrar la población de un individuo dados sus SNPs
* analisis_4_pca.py : Método para encontrar los SNPs más relevantes, comprobando su precisión con una red neuronal
* analisis_4_svm.py : SVM para encontrar la población de un individuo dados sus SNPs
