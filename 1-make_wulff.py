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

import _experhelp
import ast
#####################################
dataf='data/'
izaif='IZA-input/'
outpf=''
#####################################
pwd=''
mypath=pwd+dataf
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
######################################

def extr(f,ext):
	## take file name
	## take out surfaceE,millerIND
	## regardless if reaxff or estimate
	co=open(f,'r')
	line=co.readline()
	print(line)
	line=line.rstrip('\n')
	co.close()
	
	spli=line.split(' ')
	spli=spli[1:]
	l=len(spli)
	print(spli)
	p=int(l/2)
	mi=spli[0:p]
	en=spli[p:]
	print('youuu')
	print(en)
	if ext == 'e':
		nen=[]
		for e in en:
			nen.append(float(e)*100.0)
		en=nen
	
			
	return mi,en

	
	
def WulffAssem(mi,en,struct):
	## assemble info to make wulff
	print(mi)
	print(en)
	millerIND=[ast.literal_eval(m) for m in mi]
	surfaceE=[float(e) for e in en]
	dictionary = dict(zip(millerIND, surfaceE))
	miller_list = dictionary.keys()
	e_surf_list = dictionary.values()
	wulffshape = WulffShape(struct.lattice, miller_list, e_surf_list)
	return wulffshape

	
#['CAS.e','CFI.e','EEI.e','EON.r','ETL.e','IRR.r','ITE.e','ITH.e','ITR.e','JOZ.e']
#'JRY.e','JSW.e'
longt=['CON.e','CON.r','OSO.r']	
invam=['CAS.e'] #cas.e invalid when minimum stuff when n.n
invad=['CFI.e','ZON.e','ZON.r']  #invalid double scalars 360 %360
nowul=['EEI.e','EON.r','ETL.e','IRR.r','ITE.e','ITE.r','ITR.r','JOZ.e','JRY.e','JSW.e','MSE.e','MVY.r','MWW.e','MWW.r','NES.e','NES.r','PON.e','PUN.r','SEW.r','SYY.r','STI.e'] # no wulff no error
#pun.r has negative surf energy

ii=len(onlyfiles)
i=0
while(i<10):
	f=onlyfiles[i]
	
	
	## get pymatgen stucture object from iza framework
	frm=f.split('_')[1]
	frmc=frm+'.cif'
	## to print name out asap
	ext=f.split('_')[2]
	picknom=frm+'.'+ext
	## get miller indexes and energies from data file
	fd=dataf+f
	mi,en=extr(fd,ext)
	
	
	print (picknom,' ',i)
	if picknom in longt :
		print ('skipped ',picknom)
	#elif picknom == 'JOZ.e':
	else:
		#----------
		frmcd=izaif+frmc
		zeo=io.read(frmcd)
		struct=AseAtomsAdaptor.get_structure(zeo)
		fwulffshape=WulffAssem(mi,en,struct)
		
		miad=fwulffshape.miller_area_dict
		#print (miad.keys(),miad.values())
		mi,en=_experhelp.extractor(miad.values(),mi,en)
		#print (mi,en,'mien')
		wulffshape=WulffAssem(mi,en,struct)		
		
		##make a pickle file to store for future
		fl = open(picknom,'wb')
		pl.dump(wulffshape,fl)
		print('lesseee')
		## show wulff
		wulffshape.show(show_area=True,bar_on=True,legend_on=True)
	
	i=i+1
	
'''
zeo=io.read('LTA.cif')
struct=AseAtomsAdaptor.get_structure(zeo)


millerIND=[(1,0,0),(0,1,0),(0,0,1),(1,0,1)]
surfaceE=[10,10,10,9.5]

dictionary = dict(zip(millerIND, surfaceE))

miller_list = dictionary.keys()
e_surf_list = dictionary.values()
wulffshape = WulffShape(struct.lattice, miller_list, e_surf_list)
wulffshape.show()

name='LTA'
picklename=name+'.pickle'
fl = open(picklename,'wb')
pl.dump(wulffshape,fl)
'''