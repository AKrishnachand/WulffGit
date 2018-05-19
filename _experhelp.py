from os import listdir
from os.path import isfile, join
import os as os

import matplotlib.pyplot as plt
import pickle as pl
import numpy as np

import ase.io as io
from ase.lattice.spacegroup import crystal
import ase.utils.geometry as geometry
from ase.visualize import view

import pymatgen as mg
from _exper import WulffShape
from pymatgen import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import Slab, SlabGenerator, get_symmetrically_distinct_miller_indices

import numpy 
import math
import os
import copy

import ast
#####################################

def sorter(sma):
	# list in order of ['($001$) : 0.0051', '($100$) : 0.0281',......
	# output
	#   1-list sorted in order of highest to lowest area
	#   2-swirly brackets not ()
	#   3-an index used to sort list and not use miller with area = 0
	
	#swig=['{$001$} : 0.0051', '($100$) : 0.0281', '($101$) : 0.0', '($102$) : 0.003', '($110$) : 0.0', '($111$) : 0.035', '($2\\overline{1}2$) : 0.0454', '($201$) : 0.0507', '($210$) : 0.0', '($211$) : 0.3123', '($212$) : 0.0028', '($221$) : 0.0']
	
	milist=[] #list of miller index part
	arlist=[] #list of areas
	
	for i in sma:
		sp=i.split(' ')
		mi=list(sp[0])
		ar=float(sp[2])
		l=len(mi)
		mi[0]='{'
		mi[l-1]='}'
		mi=''.join(mi)
		milist.append(mi)
		arlist.append(ar)
		
	arsum=sum(arlist)
	#print (arlist,'arlist sorter')
	pelist=[] #percentage of area
	for ar in arlist:
		pe=ar/arsum*100.00000000000000
		pelist.append(pe)
	#print (pelist,'pelist sorter')
	pelist=np.array(pelist)	
	#print (pelist,'pelist in sorter')
	sorti1=numpy.argsort(pelist)[::-1]
	pelist=pelist[sorti1]
	sorti2=np.where(pelist > 0.0000000000000000000000000001)
	pelist=pelist[sorti2]
	milist=np.array(milist,dtype=object)
	milist=milist[sorti1]
	milist=milist[sorti2]
	swig=[]
	ii=len(milist)
	i=0
	while(i<ii):
		mi=milist[i]
		pe=int(round(pelist[i]))
		#print (pelist[i])
		st=mi+' : '+str(pe)
		#print (st)
		swig.append(st)
		i=i+1
	
	
	return swig,sorti1,sorti2
	
		
	
	

def extractor(arlist,mi,en):
				
	arsum=sum(arlist)
	#print (arlist,'arlist')
	pelist=[] #percentage of area
	for ar in arlist:
		pe=ar/arsum*100.0000000000
		pelist.append(pe)
	pelist=np.array(pelist)	
	extrac=np.where(pelist > 0.5)
	
	pelist=pelist[extrac]
	#print (pelist,extrac,'pelist,extrac')
	
	mi=np.array(mi,dtype=object)
	mi=mi[extrac]
	
	en=np.array(en,dtype=object)
	en=en[extrac]
	#print (en, 'extractor')
	return mi,en
	
	
	