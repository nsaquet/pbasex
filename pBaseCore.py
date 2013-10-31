import numpy as np
import scipy as sp
import ctypes as ct
from os.path import isfile
from scipy.ndimage.interpolation import map_coordinates as mpc
from scipy.ndimage.interpolation import spline_filter
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
        x=np.arange(0,700)
        y=np.arange(0,500)
        X,Y=np.meshgrid(x,y)
        self.datas=np.exp(-(((X-255)/np.sqrt(2)/20.)**2+((Y-255)/np.sqrt(2)/20.)**2))
        self.center=(0.,0.)
        self.r=0.
        
    def get_NumberPoly(self):
    	if not self.odd: return np.arange(0,self.lmax+1,2).shape[0]
    	else: return np.arange(0,self.lmax+1,1).shape[0]
        
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
        
    def Symmetrize(self,data):
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
        for k in xindex: data[yindex,k]=0.5*(data[yindex,k]+data[yindex[::-1],k])
        
    def AutoCenter(self):
        """
            Using the Bordas TT* criterion walk across the image from trial 
            centre in a (crude) search for that which maximises TT* function.
            NB there ought to be checking that the region (x0,y0) and radius dr 
            being searched always falls totally within image area, else array 
            bound problems are possible

        """
        #print self.center,self.r
        Cmax=0
        center,Cn=self.Newcenter(10)
        for i in np.arange(20):
        	if Cn>Cmax:
        		self.center=center
        		Cmax=Cn
        		#print Cn, center
        		center,Cn=self.Newcenter(10)
        	else: break
        
    def Newcenter(self,dr):
    	xl=np.arange(self.center[0]-dr,self.center[0]+dr,dtype=int)
    	yl=np.arange(self.center[1]-dr,self.center[1]+dr,dtype=int)
    	crit=np.array([self.CenteringCriterion(j,i) for i in xl for j in yl]).reshape((2*dr,2*dr))
    	critmax=np.where(crit==crit.max())
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
    
    def LoadBasis(self):
    	"""
    		Look for Basis file and load it. Creation if necessary.
    	"""
    	#Build basis file name
    	path='./BasisFiles/'
    	if not self.odd: filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1,2)])
    	else:filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1)])
    	fileextension='.dat'
    	Sfilename=path+'S'+filenamescore+fileextension
    	Ufilename=path+'U'+filenamescore+fileextension
    	Vfilename=path+'V'+filenamescore+fileextension
    	#Generate Basis file if not existing
    	if not (isfile(Sfilename) and isfile(Ufilename) and isfile(Vfilename)): 
    	    print 'files missing'
    	    return -1
    	
    	#Loading file
    	s=np.diag(np.fromfile(Sfilename))
    	temp=np.fromfile(Ufilename)
    	u=temp.reshape((temp.shape[0]/s.shape[0],s.shape[0]))
    	v=np.fromfile(Vfilename).reshape((s.shape[0],s.shape[0]))
    	
    	return s,u,v
    
    def Invert(self,s,u,v):
    	"""
    		Load Basis and do the inversion.
    		Generate PES.
    	"""
    	rdata=self.raw
    	#Symmetrize the problem
    	self.Symmetrize(rdata)
    	data=rdata[self.center[1]-self.r:self.center[1]+self.r,self.center[0]-self.r:self.center[0]+self.r]
    	#Convert into polar coordinates
        rad=np.arange(256)
        theta=np.pi*np.arange(256)/255.
        fact=np.min(rdata.shape)/2./rad.shape[0]	#Scaling parameter
        #x_step=rad.max()/(0.5*data.shape[0])
        #y_step=rad.max()/(0.5*data.shape[1])
        #print x_step,y_step
        polar_dat=cart_to_polar(data,1./fact,1./fact,rad,theta,(0.5*data.shape[0],0.5*data.shape[0]),3)
        polar_dat=np.nan_to_num(polar_dat)	#Otherwise there is no solution
        """
        	Resolve Basis*coeff=polar_dat
        	s is diagonal of size Legendre polynioms x number of function
        	u is orthogonal of size (Radial binning x Angular bining).(Legendre polynioms x number of function)
        """
        #A=np.dot(u,np.dot(s,v))
        coefficients,res,siz,singular=np.linalg.lstsq(u,polar_dat.T.ravel())
        #Need to build up angular distribution and PES
        #Angular Matrix is of size NLxNR where NL is the number of Legendre polynoms and Nr the radial binning
        angular=np.zeros((s.shape[0]/128,rad.shape[0]))
        polradius=int(self.r/fact)
        kvector=[np.arange(np.maximum(0,r/2-5),np.minimum(r/2+6,127)+1) for r in rad[:polradius]]
        radfunc=[np.exp(-0.5*(i-2*j)**2) for i,j in zip(rad,kvector)]
        NL=self.get_NumberPoly()
        angfunc=[coefficients[i*NL+j] for i in kvector for j in np.arange(NL)]
        angular[:,:polradius]=np.array([sum(radfunc[i]*angfunc[j::NL][i]) for i in np.arange(len(radfunc)) for j in np.arange(NL)]).reshape(NL,len(radfunc))
        ang0=angular[0,:]
        ind=np.where(ang0<0)[0]
        if len(ind)>0: 
    		ang0[ind]=0.0
    		angular[:,ind]=0.0
    	ind=np.where(ang0>0)[0]
    	angular[1:,ind]/=ang0[ind]
        pmax=(ang0*rad**2).max()
        normed_PES=ang0*rad**2 /pmax
        
        return data,angular,normed_PES,fact   

# Auxiliary function to map polar data to a cartesian plane
def polar_to_cart(polar_data, theta_step, range_step, x, y, order=2):

    # "x" and "y" are numpy arrays with the desired cartesian coordinates
    # we make a meshgrid with them
    X, Y = np.meshgrid(x, y)

    # Now that we have the X and Y coordinates of each point in the output plane
    # we can calculate their corresponding theta and range
    Tc = np.degrees(np.arctan2(Y, X)).ravel()
    Rc = (np.sqrt(X**2 + Y**2)).ravel()

    # Negative angles are corrected
    Tc[Tc < 0] = 360 + Tc[Tc < 0]

    # Using the known theta and range steps, the coordinates are mapped to
    # those of the data grid
    Tc = Tc / theta_step
    Rc = Rc / range_step

    # An array of polar coordinates is created stacking the previous arrays
    coords = np.vstack((Rc, Tc))

    # To avoid holes in the 360 - 0 boundary, the last column of the data
    # copied in the begining
    polar_data = np.vstack((polar_data, polar_data[-1,:]))

    # The data is mapped to the new coordinates
    # Values outside range are substituted with zeros
    filt_polar_data=spline_filter(polar_data,3)
    cart_data = mpc(filt_polar_data, coords, order=order, mode='constant', cval=0.0)

    # The data is reshaped and returned
    return(cart_data.reshape(len(y), len(x)).T)


# Auxiliary function to map cartesian data to a polar plane
def cart_to_polar(cart_data, x_step, y_step, r, theta, origin, order=3):

    # "r" and "t" are numpy arrays with the desired cartesian coordinates
    # we make a meshgrid with them
    R, T = np.meshgrid(r, theta)

    # Now that we have the R and T coordinates of each point in the output plane
    # we can calculate their corresponding theta and range
    xind = (R*np.sin(T)).ravel()
    yind = (R*np.cos(T)).ravel()

    # Using the known x and y steps, the coordinates are mapped to
    # those of the data grid
    xind = origin[0]- (xind / x_step)
    yind = origin[1]+ (yind / y_step)

    # An array of polar coordinates is created stacking the previous arrays
    coords = np.vstack((xind, yind))

    # The data is mapped to the new coordinates
    # Values outside range are substituted with nans
    polar_data = mpc(cart_data, coords,order=order, mode='constant', cval=0.0)
    
    # The data is reshaped and returned
    return(polar_data.reshape(len(theta), len(r)).T)


"""
	To polar coordinates
	Invert
	---> PES sort of OK but need smoothing
	
	Need to build display of the inversion result
	Need Save
	Need status bar
	Need stability
"""
