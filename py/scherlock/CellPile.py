import subprocess
import atexit
import sys
import time
import pathlib
import os

import os
jardir = os.path.join(os.path.dirname(__file__),'data')
jarfile = os.path.join(jardir,'pycellpile.jar')
#jarfile = '/home/mahogny/javaproj/isoform/py/scherlock/pycellpile.jar'

########### Need to set these things before importing jnius
import jnius_config
jnius_config.add_options('-Xrs', '-Xmx1024M')
jnius_config.set_classpath('.', jarfile)
import jnius
from jnius import autoclass


from IPython.core.display import SVG




###########################################
## One CellPile file
class CellPile:

	def __init__(self):
		#Prepare metapile
		self.cp = autoclass('isoform.cellpile.CellPileManager')()

	def addPile(self,fname,pileName=""):
		"""Load one pile and add to the list"""
		jFile = autoclass("java.io.File")
		fname = os.path.abspath(fname)
		if pileName=="" and self.cp.getNumCellPiles()>1:
			print("WARNING: No pile name given. You most likely want to use this if you work with multiple piles")
		f = jFile(fname)
		onepile = autoclass('isoform.cellpile.CellPileFile').open(f)
		self.cp.addCellPile(onepile,pileName)

	def addGTF(self, fname, trackName="gtf"):
		"""Load a GTF file to display features"""
		jFile = autoclass("java.io.File")
		fname = os.path.abspath(fname)
		fname = jFile(fname)
		self.cp.addTrackGTF(trackName, fname)


	def getView(self, gene):
		"""Return the view that will cover the span of a gene"""
		grange = self.cp.getRangeForGene(gene)
		if grange is None:
			raise Exception('Gene does not exist in GTF')
		return (grange.getSource(), grange.getFrom(), grange.getTo())




	def getBarcodes(self, index):
		"""Get the barcodes in one pileup file"""
		return self.cp.getBarcodes(index)

	def _toStringArray(self,pa):
		"""Turn python array into java array"""
		jString = autoclass("java.lang.String")
		jArray = autoclass("java.lang.reflect.Array")
		str_array = jArray.newInstance(jString, len(pa))
		for i in range(0,len(pa)):
			str_array[i] = jString(pa[i])
		return str_array

	def pileup(self, viewrange, cellBC=None, cellCluster=None, cellFile=None, numdiv=1000):
		"""Build a pileup"""

		#Split up the viewrange. [sequence, from, to]  where sequence is a string
		seq, sfrom, sto = viewrange

		#If no cells given then everything into one group
		if cellBC is None:
			pileup=self.cp.buildPileup(seq, sfrom, sto, numdiv)
			renderer=self.cp.render(seq, sfrom, sto, pileup)
			return Pileup(renderer)

		#If no cell cluster given, put everything in one group
		if cellCluster is None:
			cellCluster = [""] * len(cellBC)

		#If no cell file given, assume all ""
		if cellFile is None:
			cellFile = [""] * len(cellBC)

		#Convenience conversions, especially for scanpy
		if hasattr(cellFile, 'tolist'):
			#if isinstance(cellFile, pandas.core.series.Series):
			cellFile = cellFile.tolist()
		if hasattr(cellBC, 'tolist'):
			#if isinstance(cellBC, pandas.core.series.Series):
			cellBC = cellBC.tolist()
		if hasattr(cellCluster, 'tolist'):
			#if isinstance(cellFile, pandas.core.series.Series):
			cellCluster = cellCluster.tolist()

		pileup=self.cp.buildPileup(seq, sfrom, sto, numdiv, 
			self._toStringArray(cellBC), 
			self._toStringArray(cellFile), 
			self._toStringArray(cellCluster))

		renderer=self.cp.render(seq, sfrom, sto, pileup)

		return Pileup(renderer)


############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################




########################################
## Wrapper around a pileup
class Pileup:
	def __init__(self,renderer):
		self.renderer=renderer
	
	def plot(self, save=None, trackWidth=800, labelsWidth=200, showLog=False):
		"""Produce the plot, optionally save to file"""
		self.renderer.trackWidth=trackWidth
		self.renderer.labelsWidth=labelsWidth  #does not work?
	
		self.renderer.setShowLog(showLog)
		svg=self.renderer.toSVG()
		if save is not None:
			with open(save, "w") as text_file:
				text_file.write(svg)
		
		return SVG(data=svg)

