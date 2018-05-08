# Descarga de los ficheros del proyecto 1000 genes

# Accedemos a la url especificada y descargamos todos los ficheros al directorio especificado

import requests
import sys
import os
from BeautifulSoup import BeautifulSoup

urldownload = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
downloadPath = "./"

print sys.argv

if len(sys.argv)==3:
	urldownload = sys.argv[1]
	downloadPath = sys.argv[2]
else: 
	print "Uso : http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ ./"
	sys.exit(-1)

if (not os.path.exists(downloadPath)):
	os.mkdir(downloadPath)

r = requests.get(urldownload)

page = str(BeautifulSoup(r.content))

def getURLs(pagina):
	"""
	: pagina: pagina inicial
	: retorna: urls de la pagina , index processat
	"""

	link_inicial = pagina.find("a href")
	if link_inicial == -1:
		return None, 0 

	inici_link = page.find('"',link_inicial)
	stop_link = page.find('"', inici_link+1)
	url = page[inici_link+1: stop_link]
	return url, stop_link


while True:
	url, n = getURLs(page)
	page = page[n:]
	if url:
		if (url[-1]!='/'):
			if (url[0]!='?'):
				print "Starting Download [" + url
				down = requests.get(urldownload+url, stream=True)
				with open(downloadPath+url,'wb') as fd:
					for chunk in down.iter_content(chunk_size=1024*1024):
						fd.write(chunk)
				print "Finished Download [" + downloadPath+url+"]"
	else:
		break