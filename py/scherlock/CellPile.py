import subprocess
import atexit
import sys
import time
import pathlib
import os
from py4j.java_gateway import JavaGateway
from py4j.java_gateway import GatewayParameters
from py4j.java_collections import SetConverter, MapConverter, ListConverter

from IPython.core.display import SVG


### This could be an alternative pipeline. does not rely on sockets which can be a pain
#http://jpype.sourceforge.net/doc/user-guide/userguide.html



###########################################
## One CellPile file
class CellPile:

	#static variables - the JVM gateway
	_gateway=None
	_pid=None

	
	def __init__(self):
		"""Open a file and return a handle"""
		#self.fname=fname
		
		#Figure out where the jar-file is
		##jardir=str(pathlib.Path(__file__).parent.absolute())
		jardir = os.path.join(os.path.dirname(__file__),'data')
		jarfile = os.path.join(jardir,'pycellpile.jar')
		print("---------jar-----------")
		print(jarfile)

		##If there is no JVM running already, ensure there is
		if CellPile._gateway is None:
			self.port = 25336
			atexit.register(cellpile_cleanup)
		
			#Start the JVM - could pass a port here!
			#CellPile._pid = subprocess.Popen(["java","-jar",jardir+"/pycellpile.jar",str(self.port)])
			CellPile._pid = subprocess.Popen(["java","-jar",jarfile,str(self.port)])
			time.sleep(2)
			
			#Connect to the JVM
			CellPile._gateway = JavaGateway(gateway_parameters=GatewayParameters(port=self.port))
		
		#Prepare metapile
		self.cp = CellPile._gateway.jvm.isoform.cellpile.CellPileManager()



	def addPile(self,fname,pileName=""):
		"""Load one pile and add to the list"""
		fname = os.path.abspath(fname)
		if pileName=="" and self.cp.getNumFiles():
			print("WARNING: No pile name given. You most likely want to use this if you work with multiple piles")
		f = CellPile._gateway.jvm.java.io.File(fname)
		onepile = CellPile._gateway.jvm.isoform.cellpile.CellPileFile.open(f)
		self.cp.addCellPile(onepile,pileName)
		

	def addGTF(self, fname, trackName="gtf"):
		"""Load a GTF file to display features"""
		fname = os.path.abspath(fname)
		fname = CellPile._gateway.jvm.java.io.File(fname)
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
		str_array = CellPile._gateway.new_array(CellPile._gateway.jvm.java.lang.String,len(pa))
		for i in range(0,len(pa)):
			str_array[i] = pa[i]   #Can be 50k elements easily. Any faster method?
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
		
		
		

###########################################
## post-process cleanup. working so-so
def cellpile_cleanup_blahnotworking():
	print("--------- cleanup")
	if not CellPile._gateway is None:
		#First stop the java gateway as it seems to want to do some final communication otherwise
		print("--------- cleanup1")
		CellPile._gateway.jvm.java.lang.System.exit(0)
		CellPile._gateway.close()


###########################################
## post-process cleanup. working so-so
def cellpile_cleanup():
	if not CellPile._gateway is None:
		#First stop the java gateway as it seems to want to do some final communication otherwise
		CellPile._gateway.close()
		
		#Now kill the process. not sure why so complicated? but a timeout might be needed
		timeout_sec = 1
		p=CellPile._pid
		p_sec = 0
		for second in range(timeout_sec):
			if p.poll() == None:
				time.sleep(1)
				p_sec += 1
			if p_sec >= timeout_sec:
				p.kill() # supported from python 2.6
				################################################### The problem appears to be copies of CellPile._gateway ... maybe
	


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







