import re
import urllib
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import json
from tqdm import tqdm
import datetime
from datetime import date
from collections import OrderedDict, defaultdict
import requests
import pandas as pd
import plotly.express as px
import io
from openpyxl import Workbook
import ipywidgets as widgets
from google.colab import files

import copy
#@title Run OmixPubMiner
#@markdown  Specify the title of the output file. If you have searched in the title only, it will be added to the filename
Filename = 'Suppl10_OPM_2024'#@param {type:"string"}

Protein_List = 'THBS2 CAV2 SCG2 SLC6A1 SAV1 SEZ6L2 ERO1A RAB3B OBSL1 CD109 PTPN14 MRPL35 LRPAP1' #@param {type:"string"}
#@markdown  Specify the type of the Input for the Protein List, set a tick either for the Accession or the Gene Name
Accession = False #@param {type:"boolean"}
Gene_name = True #@param {type:"boolean"}


if (Accession and Gene_name) or (not Accession and not Gene_name):
  raise RuntimeError(f'You need to enter either Accession or Gene_name')

#@markdown  Taxonomy ID (e.g. homo sapiens: 9606; mus musculus: 10090)
TaxID = "9606"#@param {type:"string"}

#@markdown  If you want to enter multiple keywords, please separate them by semicolons. If multiple keywords should be found together in the abstract or title, please tick the "Together" button. Else the tool will search whether one of the keywords was found.

Together = False #@param {type:"boolean"}
Keywords = 'migration' #@param {type:"string"}
Keywords = Keywords.replace(" ", "")


#@markdown  Checks if the keywords are mentioned  in the title, in the title and the abstract or in the full text
#@markdown The default is to search in the full text
Title = False #@param {type:"boolean"}
Title_Abstract = True #@param {type:"boolean"}
Full_Text = False #@param {type:"boolean"}

if (Title and Title_Abstract) or (Title and Full_Text) or (Full_Text and Title_Abstract):
  raise RuntimeError(f'You can only choose to serch in either the Title, the title and abstract or in the full text')
if (not Title and not Title_Abstract and not Full_Text):
  Full_Text = True

#@markdown  Specify the order for the publication output - standard is relevance
Publication_date = False #@param {type:"boolean"}
Relevance = True #@param {type:"boolean"}
if Publication_date and Relevance:
  raise RuntimeError(f'You can only sort by the publication date or the relevance')

protList = Protein_List.split()
keywords = Keywords

if Accession:
  IDType = "Accession"
else:
  IDType = "Gene"

def requestOnUniProtWebpage(uniprotID, idtype, taxid, proteinrequired):
  start = "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:"
  middle = '+AND+accession:'\
          if str(idtype).lower() == 'accession' \
          else "+AND+gene:"
  end = '&format=tsv'
  contents = urllib.request.urlopen('{}{}{}{}{}'.format(start,
                taxid, middle, uniprotID, end))\
                                                         .read().decode("utf-8")
  data = io.StringIO(contents)
  df = pd.read_csv(data, sep="\t", header=0)
  return df

def getUniProtSynonyms(uniprotID, idtype, taxid):
  contents = requestOnUniProtWebpage(uniprotID, idtype, taxid, False)
  if len(contents['Gene Names']) > 0:
    syns = [str(x) for x in contents['Gene Names']][0]
    syns = re.sub(r'\s', ',', syns)
    syns = syns.split(",")
  else:
    syns = ''
  if len(contents['Protein names']) > 0:
    prots = [str(x) for x in contents['Protein names']][0]
    contents = re.split(r'\)?\s[\(\[]', re.sub(r'\]|\)', '', prots))
    for idx in range(len(contents)):
      cur_str = re.sub(r'^\s', '', contents[idx])
      cur_str = re.sub(r'\s', '-', cur_str)
      contents[idx] = re.sub(r'\"|\;', '', cur_str)
      contents.append(uniprotID)
  else:
    contents = [uniprotID]
  return syns, contents

def uniProtQuery(uniprotID, idtype, taxid, keyword,  maxdate):

  synonyms, protnames = getUniProtSynonyms(uniprotID, idtype, taxid)


  searchsummary = OrderedDict()
  searchsummary['uniprotID'] = uniprotID
  searchsummary['idtype'] = idtype
  searchsummary['taxid'] = taxid
  searchsummary['synonyms'] = synonyms
  searchsummary['protein_names'] = protnames
  searchsummary['keywords'] = keyword
  searchsummary['totalResults'] = 0
  searchsummary['reviewResults'] = 0
  searchsummary['category'] = 0
  searchsummary['false'] = 0

  papercol = Papercollection(synonyms, protnames, keyword, maxdate)
  searchsummary['totalResults'] = papercol.resultcnt()
  searchsummary['reviewResults'] = papercol.revcnt()
  if searchsummary['reviewResults'] > 0:
    searchsummary['category'] = 1
  elif searchsummary['totalResults'] > 0:
    searchsummary['category'] = 2
  elif len(searchsummary['synonyms']) > 0:
    searchsummary['category'] = 3
  return searchsummary, papercol, papercol.resultcnt(), papercol.revcnt()


class Papercollection:
  def __init__(self, synonyms, protnames, keyword, maxdate):
    self._titlesPapercollection = ['Titles','Authors','Journal','Years',\
                             'PMID','DOI','Link']
    self._syns = synonyms
    self._prots = protnames
    self._papercollection = OrderedDict()
    self._bestcategory = 4

    self._keywordlist = keyword.split(';')
    self._keywordlist = list(self._keywordlist)

    if len(self._keywordlist) == 1:
      self._keyword = self._keywordlist
    #mehrere keywords durchgehen
    elif len(self._keywordlist) >1:
      if Together == True:
        placeholder = '%20AND%20'.join(self._keywordlist)
        self._keyword = '{}'.format(placeholder)
      else:
        self._keyword = self._keywordlist


    for title in self._titlesPapercollection:
      self._papercollection[title] = list()

    records, term = self.requestPubTator(self._syns, self._prots)
    self.parseResults(records)


  def requestPubTator(self, synonyms, proteinnames):
    start = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text="
    pubtator = "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text="
    combAnd = '%20AND%20'
    syns = []
    if len(synonyms) > 0 or len(proteinnames) > 0:
      synonyms = ["@GENE_" + gene for gene in synonyms]
      print
      syns =  synonyms + proteinnames
    synquery = []

    if len(self._keyword) >1:
      for keyword in self._keyword:
          for syn in syns:
            syn = '({}{}{})'.format(keyword, combAnd, syn)
            synquery.append(syn)
    else:
        for syn in syns:
          syn = '({}{}{})'.format(self._keyword, combAnd, syn)
          synquery.append(syn)
    protQuery = '%20OR%20'.join(synquery)
    time.sleep(0.5)

    url = '{}{}'.format(start,  protQuery)
    pubLink = '{}{}'.format(pubtator,  protQuery)
    if Title:
      title = "&sections=title"
      url = '{}{}'.format(url, title)
      pubLink = '{}{}'.format(pubLink, title)
    elif Title_Abstract:
      title_abstract = "&sections=title,abstract"
      url = '{}{}'.format(url, title_abstract)
      pubLink = '{}{}'.format(pubLink, title_abstract)
    if Publication_date:
      pub = "&sort=date%20desc"
      url = '{}{}'.format(url, pub)
      pubLink = '{}{}'.format(pubLink, pub)
    contents = requests.get(url)
    results = contents.json()
    self._url = pubLink
    return results, url

  def parseResults(self, records):
    self._resultcount = 0
    self._reviewcount = 0

    if "detail" in records.keys():
      if records["detail"] == "We are currently updating the Database. Please try again later":
        print("PubTator3 is updating the database. Not all genes could be searched. Please try again later")
        return
    results = records["results"]
    for result in results:
      self._papercollection['Titles'].append(result["title"])
      if "authors" in result:
        self._papercollection['Authors'].append(result["authors"])
      else:
        self._papercollection['Authors'].append("NaN")
      self._papercollection['PMID'].append(result["pmid"])
      self._papercollection['Journal'].append(result["journal"])
      if "doi" in result:
        self._papercollection['DOI'].append(result["doi"])
      else:
        self._papercollection['DOI'].append("NaN")
      date = result["meta_date_publication"]
      date_re = re.search(r'\b\d{4}\b', date)
      year = int(date_re.group())
      self._papercollection['Years'].append(year)
      self._papercollection['Link'].append(self._url)

    facets = records["facets"]
    if "facet_fields" in facets:
      facet_fields =facets["facet_fields"]
      facet_types = facet_fields["type"]
      for ftype in facet_types:
        if ftype["name"] == "Review":
          self._reviewcount =ftype["value"]
    else:
      self._reviewcount = 0

    self._resultcount = records["count"]

  def resultcnt(self):
    return(self._resultcount)

  def revcnt(self):
    return(self._reviewcount)

  def getPapercolletion(self):
    return self._papercollection



def runMain():
  #maxdate = datetime.date(2004,8,15)
  maxdate= date.today().strftime('%Y/%m/%d')
  papercol_titles = ['Titles','Authors','Journal','Years',\
                             'PMID','DOI','Link']
  resarray = list()
  paperSummary = Workbook()
  ps = paperSummary.active
  ps.append( ['UniprotID', 'Results', 'Reviews', 'Synonyms', 'Protein names', 'Category'])
  counter = 2
  for idx in tqdm(range(len(protList))):
    protquer, papercol, papercnt , revcnt= uniProtQuery(protList[idx], IDType, TaxID, \
                                                keywords, maxdate)
    resarray.append([protList[idx], papercnt, revcnt, protquer['synonyms'], \
                    protquer['protein_names'], protquer['category']])
    ps.append([protList[idx], papercnt,revcnt, str(protquer['synonyms']), \
                    str(protquer['protein_names']), protquer['category']])
    startcnt = counter+1
    counter+=1
    ps.append(papercol_titles)
    papers = papercol.getPapercolletion()
    for idx in range(len(papers['Titles'])):
      value = list()
      for key in papers.keys():
        value.append(str(papers[key][idx]))
      ps.append(value)
      counter +=1
    ps.row_dimensions.group(startcnt,counter, hidden=True)
    counter += 1

    allRes = pd.DataFrame(resarray, columns = ['UniprotID', 'Results', 'Reviews',  \
                                                'Synonyms', 'Protein names', \
                                                'Category'])

  if Title:
    paperSummary.save(f'{Filename}_inTitleOnly.xlsx')
  elif Title_Abstract:
    paperSummary.save(f'{Filename}_inTitleAbstract.xlsx')
  else:
    paperSummary.save(f'{Filename}.xlsx')

  return allRes

allRes = runMain()
display(allRes)


pivot_table = allRes.pivot_table(columns=['Category'], aggfunc='size')
color_map = {0: '#d4d4d4', 1: '#ff8000', 2: '#ffc080', 3: '#55a0fb'}
colors = [color_map[cat] for cat in pivot_table.index]
pivot_table.plot.pie(figsize=(6,6),
                     ylabel='',
                     colors=colors)
plt.show()


#Download button
def on_button_clicked(b):
  if Title:
    files.download(f'{Filename}_inTitleOnly.xlsx')
  elif Title_Abstract:
    files.download(f'{Filename}_inTitleAbstract.xlsx')
  else:
    files.download(f'{Filename}.xlsx')

button = widgets.Button(description='Download Excel File')
button.on_click(on_button_clicked)
display(button)


