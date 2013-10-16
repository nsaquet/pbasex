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

class Datas():
    def __init__(self):
        self.lmax=2
        self.odd=0
        self.raw=2.*np.random.normal(0.5,size=(256,256))
        x=np.arange(0,10,0.01)
        y=np.arange(0,10,0.01)
        X,Y=np.meshgrid(x,y)
        self.datas=np.exp(-((X-5)**2+(Y-4)**2))
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
        datmax=self.datas.max()
        mask=(self.datas>datmax*0.25)
        self.center=com(self.datas.T,mask.T)
        
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
        
    def AutoCenter(self):
        """
            Using the Bordas TT* criterion walk across the image from trial 
            centre in a (crude) search for that which maximises TT* function.
            NB there ought to be checking that the region (x0,y0) and radius dr 
            being searched always falls totally within image area, else array 
            bound problems are possible

        """
        print self.center,self.r
        Cmax=0
        center,Cn=self.Newcenter(10)
        for i in np.arange(20):
        	if Cn>Cmax:
        		self.center=center
        		Cmax=Cn
        		print Cn, center
        		center,Cn=self.Newcenter(10)
        	else: break
        
    def Newcenter(self,dr):
    	xl=np.arange(self.center[0]-dr,self.center[0]+dr,dtype=int)
    	yl=np.arange(self.center[1]-dr,self.center[1]+dr,dtype=int)
    	crit=np.array([self.CenteringCriterion(i,j) for i in xl for j in yl]).reshape((2*dr,2*dr))
    	critmax=np.where(crit==crit.max())
    	print critmax
    	return (critmax[0][0]+self.center[0]-dr,critmax[1][0]+self.center[1]-dr), crit.max()
    	
    def CenteringCriterion(self,x0,y0):
        """
            Evaluates Bordas centring criterion: Sum T_ij bt T*_ij 
            (where * indicate reflection through the trial centre (x0,y0).
            Implicitly assumes centre is mid-pixel
        """        
        #Indexes within  the circle 
        listx=np.arange(x0-self.r,x0+self.r,dtype=int)  
        listy=np.arange(y0-self.r,y0+self.r,dtype=int)
        xx,yy=np.meshgrid(listx,listy)   
        listsquare=((xx-x0)**2 + (yy-y0)**2)<self.r**2
        xind=xx[listsquare]
        yind=yy[listsquare]
        r_xind=int(2*x0)-xind
        r_yind=int(2*y0)-yind
        
        return sum(self.datas[xind,yind]*self.datas[r_xind,r_yind])
                
"""
	Need the auto center methode
	Need Invert
	Need Save
	Need status bar
	Need stability
"""