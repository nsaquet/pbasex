import numpy as np
import scipy as sp
import ctypes as ct

from scipy.ndimage.measurements import center_of_mass as com
import pyfits

def Generate_Basis(lmax,odd):
    _libpBasis=np.ctypeslib.load_library('libCoreBasis','.')
    _libpBasis.write_forward.argtypes=[ct.c_int,ct.c_int]
    _libpBasis.write_forward.restype=ct.c_void_p
    _libpBasis.write_forward(ct.c_int(lmax),ct.c_int(odd))

def CenterImage2D(image):
    return com(image)
    

class Datas():
    def __init__(self):
        self.lmax=2
        self.odd=0
        self.raw=2.*np.random.normal(0.5,size=(256,256))
        self.datas=0.5*self.raw
        self.center=(0.,0.)
        self.r=0.
        
    def OpenFile(self,filepath):
        if filepath[-4:]=='.fit':
           hdulist = pyfits.open(filepath)
           #hdulist.info()
           scidata = hdulist[0].data
           hdulist.close()
           #print scidata.shape, scidata.dtype.name
           self.raw=scidata
           self.datas=scidata
        elif filepath[-4:]=='.dat' or filepath[-4:]=='.txt':
            scidata = np.loadtxt(filepath,int)
            #print scidata.shape
            self.raw=scidata
            self.datas=scidata
        elif filepath[-4:]=='.jpg' or filepath[-4:]=='.bmp' or filepath[-5:]=='.tiff' or filepath[-4:]=='.png':
            scidata = sp.ndimage.imread(filepath)
            self.raw=scidata
            self.datas=scidata
        else: print "Incorrect data file format"
        
    def get_com(self):
        self.center=CenterImage2D(self.datas)
        
    def Symmetrize(self):
        """
    		Symmetrise a 2_D circular selection vertically (ie about a horizontal axis). 
    		Assume that the centre is mid-pixel (x0,y0) rather than at lower left corner
    		of pixel x0,y0. Assume destination array bb[][] is pre-zeroed. Note that no 
    		symmetrisation is needed horizontally since the Legendre polynomials are 
    		already symmetric along the vertical axis. (Vertically being the polarisation 
    		axis of the ligth, or the direction of propagation in the case of cpl).
    	"""
    	#Need to build up the selected indexes within self.r
        xindex=np.arange(self.center[0]-self.r,self.center[0]+self.r,dtype=int)
        yindex=np.arange(self.center[1]-self.r,self.center[1]+self.r,dtype=int)
        for k in xindex: self.datas[yindex,k]=0.5*(self.datas[yindex,k]+self.datas[yindex[::-1],k])
        