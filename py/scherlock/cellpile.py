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
jnius_config.add_options('-Xrs', '-Xmx100000M')
jnius_config.set_classpath('.', jarfile)
import jnius
from jnius import autoclass


from IPython.core.display import SVG




###########################################
## One CellPile file
class cellpile:

	def __init__(self):
		#Prepare metapile
		self.cp = autoclass('isoform.cellpile.CellPileManager')()



	def add_cellpile(self,fname,cellpile_name=""):
		"""Load one pile and add to the list"""
		jFile = autoclass("java.io.File")
		fname = os.path.abspath(fname)
		if cellpile_name=="" and self.cp.getNumCellPiles()>1:
			print("WARNING: No pile name given. You most likely want to use this if you work with multiple piles")
		f = jFile(fname)
		onepile = autoclass('isoform.cellpile.CellPileFile').open(f)
		self.cp.addCellPile(onepile,cellpile_name)



	def add_gtf(self, fname, trackName="gtf"):
		"""Load a GTF file to display features"""
		jFile = autoclass("java.io.File")
		fname = os.path.abspath(fname)
		fname = jFile(fname)
		self.cp.addTrackGTF(trackName, fname)

	def dump_gtf(self, fname, trackName="gtf"):
		"""Dump a parsed gtf file on user"""
		jFile = autoclass("java.io.File")
		fname = os.path.abspath(fname)
		fname = jFile(fname)
		return self.cp.dumpTrackGTF(trackName, fname)



	def add_bed(self, fname, trackName="bed"):
		"""Load a Bed file to display features"""
		jFile = autoclass("java.io.File")
		# fname = os.path.abspath(fname)
		# fname = jFile(fname)
		self.cp.addTrackBed(trackName, fname)



	def get_view(self, gene):
		"""Return the view that will cover the span of a gene"""
		grange = self.cp.getRangeForGene(gene)
		if grange is None:
			raise Exception('Gene does not exist in GTF')
		return (grange.getSource(), grange.getFrom(), grange.getTo())

	def get_barcodes(self, index):
		"""Get the barcodes in one pileup file"""
		return self.cp.getBarcodes(index)

	def get_number_of_cellpiles(self):
		"""Get the number of pileups"""
		return self.cp.getNumCellPiles()



	def _toStringArray(self,pa):
		"""Turn python array into java array"""
		jString = autoclass("java.lang.String")
		jArray = autoclass("java.lang.reflect.Array")
		str_array = jArray.newInstance(jString, len(pa))
		for i in range(0,len(pa)):
			str_array[i] = jString(pa[i])
		return str_array



	def pileup(self, viewrange, 
				barcodes=None, track_labels=None, cellpile_names=None, 
				numdiv=1000, 
				show_inbetweens=False, individual_track_scaling=False):
		"""Build a pileup"""

		#Split up the viewrange. [sequence, from, to]  where sequence is a string
		seq, sfrom, sto = viewrange

		#If no cells given then everything into one group
		if barcodes is None:
			pileup=self.cp.buildPileup(seq, sfrom, sto, numdiv)
			renderer=self.cp.render(seq, sfrom, sto, pileup)
			return Pileup(renderer)

		#If no cell cluster given, put everything in one group
		if track_labels is None:
			track_labels = [""] * len(barcodes)

		#If no cell file given, assume all ""
		if cellpile_names is None:
			cellpile_names = [""] * len(barcodes)

		#Convenience conversions, especially for scanpy
		if hasattr(cellpile_names, 'tolist'):
			#if isinstance(cellpile_names, pandas.core.series.Series):
			cellpile_names = cellpile_names.tolist()
		if hasattr(barcodes, 'tolist'):
			#if isinstance(barcodes, pandas.core.series.Series):
			barcodes = barcodes.tolist()
		if hasattr(track_labels, 'tolist'):
			#if isinstance(cellpile_names, pandas.core.series.Series):
			track_labels = track_labels.tolist()

		pileup=self.cp.buildPileup(seq, sfrom, sto, numdiv, 
			self._toStringArray(barcodes), 
			self._toStringArray(cellpile_names), 
			self._toStringArray(track_labels))

		renderer=self.cp.render(seq, sfrom, sto, pileup, 
			show_inbetweens, individual_track_scaling)

		return Pileup(renderer)






	def pileup_raw(self, viewrange, barcodes=None, track_labels=None, cellpile_names=None, numdiv=1000):
		"""Build a pileup"""

		#Split up the viewrange. [sequence, from, to]  where sequence is a string
		seq, sfrom, sto = viewrange

		#If no cells given then everything into one group
		if barcodes is None:
			pileup=self.cp.buildPileup(seq, sfrom, sto, numdiv)
			renderer=self.cp.render(seq, sfrom, sto, pileup)
			return Pileup(renderer)

		#If no cell cluster given, put everything in one group
		if track_labels is None:
			track_labels = [""] * len(barcodes)

		#If no cell file given, assume all ""
		if cellpile_names is None:
			cellpile_names = [""] * len(barcodes)

		#Convenience conversions, especially for scanpy
		if hasattr(cellpile_names, 'tolist'):
			#if isinstance(cellpile_names, pandas.core.series.Series):
			cellpile_names = cellpile_names.tolist()
		if hasattr(barcodes, 'tolist'):
			#if isinstance(barcodes, pandas.core.series.Series):
			barcodes = barcodes.tolist()
		if hasattr(track_labels, 'tolist'):
			#if isinstance(cellpile_names, pandas.core.series.Series):
			track_labels = track_labels.tolist()

		pileup=self.cp.buildPileup(seq, sfrom, sto, numdiv, 
			self._toStringArray(barcodes), 
			self._toStringArray(cellpile_names), 
			self._toStringArray(track_labels))


		return(pileup)


		# renderer=self.cp.render(seq, sfrom, sto, pileup)

		# return Pileup(renderer)	






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

