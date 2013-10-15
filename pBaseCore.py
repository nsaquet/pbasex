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
        self.raw=0.*np.random.normal(0.5,size=(256,256))
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
           
            
        
        