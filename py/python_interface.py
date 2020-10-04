import subprocess
import atexit
import sys
import time
from py4j.java_gateway import JavaGateway
from py4j.java_gateway import GatewayParameters
from py4j.java_collections import SetConverter, MapConverter, ListConverter


### This could be an alternative pipeline. does not rely on sockets which can be a pain
#http://jpype.sourceforge.net/doc/user-guide/userguide.html


###########################################
## One CellPile file
class CellPile:

	#static variables - the JVM gateway
	_gateway=None
	_pid=None

	
	def __init__(self,fname):
		"""Open a file and return a handle"""
		self.fname=fname

		##If there is no JVM running already, ensure there is
		if CellPile._gateway is None:
			self.port = 25336
			atexit.register(cellpile_cleanup)
		
			#Start the JVM - could pass a port here!
			CellPile._pid = subprocess.Popen(["java","-jar","pycellpile.jar",str(self.port)])
			time.sleep(2)
			
			#Connect to the JVM
			CellPile._gateway = JavaGateway(gateway_parameters=GatewayParameters(port=self.port))


		#Actually open the file	
		f = CellPile._gateway.jvm.java.io.File(fname)
		self.cp = CellPile._gateway.jvm.isoform.cellpile.CellPileFile.open(f)
		
		#Create an empty GTF
		self.gtf = CellPile._gateway.jvm.isoform.util.GtfParser()
		

	def getBarcodes(self):
		return self.cp.getListBarcodes()

	def getSequences(self):
		return self.cp.getListSequences()

	def getView(self, gene):
		"""Return the view that will cover the span of a gene"""
		grange = self.gtf.getRangeForGene(gene)
		if grange is None:
			raise Exception('Gene does not exist in GTF')
		return (grange.source, grange.getFrom(), grange.getTo())

	def loadGTF(self, fname):
		"""Load a GTF file to display features"""
		fname = CellPile._gateway.jvm.io.File(fname)
		self.gtf = CellPile._gateway.jvm.isoform.util.GtfParser(fname)



	def pileup(self, viewrange, cellgroups, numdiv=1000):
		"""Build a pileup"""
	
		#Split up the viewrange. [sequence, from, to]  where sequence is a string
		seq, sfrom, sto = viewrange
	
		#Format barcodes by cluster the way it is expected
		
		#{"group1":[list of cells]}
	
#		public int[][] cellGroups;
	
		return Pileup(self.cp.buildPileup(seq, sfrom, sto, numdiv, cellgroups))


###########################################
## post-process cleanup. working so-so
def cellpile_cleanup():
	if not CellPile._gateway is None:
		#First stop the java gateway as it seems to want to do some final communication otherwise
		CellPile._gateway.close()
		
		#Now kill the process. not sure why so complicated? but a timeout might be needed
		timeout_sec = 2
		p=CellPile._pid
		p_sec = 0
		for second in range(timeout_sec):
			if p.poll() == None:
				time.sleep(1)
				p_sec += 1
			if p_sec >= timeout_sec:
				p.kill() # supported from python 2.6
	

########################################
## Wrapper around a pileup
class Pileup:
	def __init__(self,pileup,gtf):
		selp.pileup=pileup
		self.gtf=gtf
	
	def plot(self, save=None):
		"""Produce the plot, optionally save to file"""
		svg=self.pileup.tosvg(cp.gtf)
		if save is not None:
			print("should save image TODO")
		
		return svg







