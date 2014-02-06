import numpy as np
import scipy as sp
from os.path import isfile
from scipy.ndimage.interpolation import map_coordinates as mpc
from scipy.ndimage.interpolation import spline_filter
from scipy.ndimage.measurements import center_of_mass as com
import pyfits

"""
	Definition of constants that are linked to the calculated basis set
"""
Rbin=440
Angbin=440
Funcnumber=220
Bspace=2.	#Basis space is currently set to 2
Bwidth=1.	#Basis width is currently set to 1

"""
	Structure containing the datas informations (X,Y)
	It generate as well scaling parameters for (X,Y)->(R,T)
"""
class ArrayInfos():
	def __init__(self,data):
		self.nX=data.shape[0]
		if len(data.shape)>1: self.nY=data.shape[1]
		else: self.nY=0
		self.Xfact=1
		self.Yfact=1
		self.Rfact=max(self.Xfact,self.Yfact)
		self.nR=1
		self.update_nR()
	
	def update_nR(self):
		nr=min(0.5*self.Xfact*(self.nX-1)/self.Rfact,0.5*self.Yfact*(self.nY-1)/self.Rfact)
		if nr>Rbin:
			nr=Rbin
			self.Rfact*=nr/float(Rbin)
		self.nR=int(nr)

"""
	Core class of the program. 
	Contains all parameters: Datas (raw and displayed),Max Legendre polynoms, Odd boolean, 
	Center of the image, size of the area of interest in radius and scaling parameters
	
	Contains all function to:
		Get the center
		Get the number of Legendre polynoms for the fit
		Dealing with files
		Symmetrize the picture
		Loading the Basis set
		Do the inversion
		
"""
class Datas():
    def __init__(self):
        self.lmax=2
        self.odd=0
        self.raw=2.*np.random.normal(0.5,size=(256,256))
        x=np.arange(0,700)
        y=np.arange(0,500)
        X,Y=np.meshgrid(x,y)
        self.datas=np.exp(-(((X-255)/np.sqrt(2)/20.)**2+((Y-205)/np.sqrt(2)/20.)**2))
        self.center=(0.,0.)
        self.get_com()
        self.r=10.
        self.scale=ArrayInfos(self.datas)
        
    def get_NumberPoly(self):
    	if not self.odd: return np.arange(0,self.lmax+1,2).shape[0]
    	else: return np.arange(0,self.lmax+1,1).shape[0]
        
    def OpenFile(self,filepath):
        if filepath[-4:]=='.fit':
           hdulist = pyfits.open(filepath)
           #hdulist.info()
           scidata = hdulist[0].data
           hdulist.close()
           ind=np.where(scidata.ravel()<0)[0]
           scidata.ravel()[ind]=0
           #print scidata.shape, scidata.dtype.name
           self.raw=scidata
           self.datas=scidata
           self.scale=ArrayInfos(self.datas)
        elif filepath[-4:]=='.dat' or filepath[-4:]=='.txt':
            scidata = np.loadtxt(filepath,int)
            ind=np.where(scidata.ravel()<0)[0]
            scidata.ravel()[ind]=0
            #print scidata.shape
            self.raw=scidata
            self.datas=scidata
            self.get_com()
            self.scale=ArrayInfos(self.datas)
        elif filepath[-4:]=='.jpg' or filepath[-4:]=='.bmp' or filepath[-5:]=='.tiff' or filepath[-4:]=='.png':
            scidata = sp.ndimage.imread(filepath)
            ind=np.where(scidata.ravel()<0)[0]
            scidata.ravel()[ind]=0
            self.raw=scidata
            self.datas=scidata
            self.get_com()
            self.scale=ArrayInfos(self.datas)
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
        yindex=np.arange(self.center[1]-self.r,self.center[1]+self.r,dtype=int)
        xindex=np.arange(self.center[0]-self.r,self.center[0]+self.r,dtype=int)
        for k,l in zip(xindex[round(len(xindex)/2.):],xindex[len(xindex)/2 -1::-1]): 
        	yind=np.where((k-self.center[0])**2+(yindex-self.center[1])**2<self.r**2)[0]
        	data.T[k,yindex[yind]]=0.5*(data.T[k,yindex[yind]]+data.T[l,yindex[yind]])
        	data.T[l,yindex[yind]]=data.T[k,yindex[yind]]
        #if len(xindex)%2: data.T[xindex[len(xindex)/2],yindex]+=data.T[xindex[len(xindex)/2],yindex]
        
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
    		Look for Basis file and load it. We need only the initial basis for leastsq
    	"""
    	#Build basis file name
    	path='./BasisFiles/'
    	if not self.odd: filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1,2)])
    	else:filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1)])
    	fileextension='.dat'
    	Ufilename=path+'U'+filenamescore+fileextension
    	#Generate Basis file if not existing
    	if not (isfile(Ufilename)): 
    	    print 'files missing'
    	    return -1
    	
    	#Loading file as a vector as in PBaseX.
    	M=440*440   #correspond to the polar array dimension in the basis file NR=440 NTH=440
    	N=self.get_NumberPoly()*220 #Number of polynoms by the number of function set to NR/2.
    	u=np.fromfile(Ufilename,dtype=float64).reshape((M,N))
    	
    	return u
    
    def Invert(self,u):
    	"""
    		Load Basis and do the inversion.
    		Generate PES.
    	"""
    	rdata=self.raw
    	#Symmetrize the problem
    	self.Symmetrize(rdata)
    	
    	#Convert into polar coordinates
        polar,rad_dist =cart2pol(rdata,self.scale,self.center)
        polar_dat=np.concatenate((polar,np.zeros(u.shape[0]-len(polar),dtype=polar.dtype)))
        
        
        #Resolve Basis*coeff=polar_dat
        coefficients,residuts,siz,singulars=np.linalg.lstsq(u,polar_dat)       #Need to build up angular distribution and PES
        
        width=2*Bwidth**2
        irmax=int((self.r/self.scale.Rfact)**2 /self.scale.nR )
        if irmax>Rbin: irmax=Rbin
        rad=np.sqrt(np.arange(irmax)*self.scale.nR*1.)
        kvector=np.arange(Funcnumber)
        
        #Angular Matrix is of size NLxNR where NL is the number of Legendre polynoms and Nr the radial binning
        NL=self.get_NumberPoly()
        angular=np.zeros((NL,Angbin),dtype='float')
        coefs=coefficients.reshape((Funcnumber,NL))
        
        
        for r in rad:
            fradial=np.exp(-(r-Bspace*kvector)**2/width)
            for l in np.arange(NL):
                angular[l][r]=sum(fradial*coefficients[l::NL])
        
        ang0=angular[0,:]
        ind=np.where(ang0<0)[0]
        if len(ind)>0: 
    		ang0[ind]=0.0
    		angular[:,ind]=0.0
    	ind=np.where(ang0>0)[0]
    	#angular[1:,ind]/=ang0[ind]
        #pmax=(ang0*rad**2).max()
        normed_PES=ang0#*rad**2 /pmax
        
        return normed_PES 
        """
        #polradius=int(self.r/fact)
        #kvector=[np.arange(np.maximum(0,r/2-5),np.minimum(r/2+6,127)) for r in rad[:polradius]]
        radfunc=[np.exp(-0.5*(i-2*j)**2) for i,j in zip(rad,kvector)]
        
        
        
        for i,r in enumerate(radfunc):
            angular[:,i]=np.sum(r*coefs[:,kvector[i]],axis=1)
        #angfunc=[coefficients[i*NL+j] for i in kvector for j in np.arange(NL)]
        #angular[:,:polradius]=np.array([sum(radfunc[i]*angfunc[j::NL][i]) for i in np.arange(len(radfunc)) for j in np.arange(NL)]).reshape(NL,len(radfunc))
        ang0=angular[0,:]
        ind=np.where(ang0<0)[0]
        if len(ind)>0: 
    		ang0[ind]=0.0
    		angular[:,ind]=0.0
    	ind=np.where(ang0>0)[0]
    	angular[1:,ind]/=ang0[ind]
        pmax=(ang0*rad**2).max()
        normed_PES=ang0*rad**2 /pmax
        
        return data,polar_dat,normed_PES,fact,ang0 
        """
# Auxiliary function to map polar data to a cartesian plane
def polar_to_cart(polar_data, theta_step, range_step, x, y, order=3):

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
    filt_polar_data=spline_filter(polar_data,order)
    cart_data = mpc(polar_data, coords, order=order, mode='constant', cval=0.0)

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
    filt_cart_data=spline_filter(cart_data,2)
    polar_data = mpc(filt_cart_data, coords,order=order, mode='constant', cval=0.0)
    
    # The data is reshaped and returned
    return(polar_data.reshape(len(theta), len(r)).T)

#Change to G. Garcia function from IGOR
def cart2pol(data,scale,center):
	Rfact=scale.Rfact
	nR=scale.nR
	
	if nR>Rbin: nr=Rbin #Security
	polar=np.array([])
	rfunc=[]
	for l in np.arange(nR):
		nth=2*l+1	#define a r dependant angular binning 
		
		rad=l*Rfact	#Radius in image unit
		
		theta=np.pi*np.arange(nth)/nth
		x=rad*np.cos(theta) + center[1]
		y=-rad*np.sin(theta) + center[0]
		
		#Cubic interpolation
		xpix=x/scale.Xfact
		ypix=y/scale.Yfact
		ix=xpix.astype('int')
		iy=ypix.astype('int')
		dx=xpix-ix
		dy=ypix-iy
		
		Pol=np.zeros_like(theta)
		
		"""
			Fastest implementation of the cubic interpolation so far
			Faster than max(0,(x+1)**3) by 20% or any of the factorization by bool (x+2>0) by 60-65%
		"""
		def cubic(x):
			p0=lambda y: (y+2)**3 if (y+2>0) else 0
			p1=lambda y: (y+1)**3 if (y+1>0) else 0
			p2=lambda y: (y)**3 if (y>0) else 0
			p3=lambda y: (y-1)**3 if (y-1>0) else 0
			return (p0(x)-4*p1(x)+6*p2(x)-4*p3(x))/6.
		for index in np.arange(nth):
			for i in np.arange(-1,3):
				for j in np.arange(-1,3):
					condx=np.logical_and((ix[index]+i)<scale.nX,(ix[index]+i)>=0)
					condy=np.logical_and((iy[index]+j)<scale.nY,(iy[index]+j)>=0)
					cub=cubic(i-dx[index])*cubic(dy[index]-j)*np.logical_and(condx,condy)
					#print i-dx[index],cubic(i-dx[index]) ,cub
					Pol[index]+=data.ravel()[(ix[index]+i)*scale.nY+(iy[index]+j)]*cub	
		polar=np.concatenate((polar,Pol))
		rfunc.append(Pol.sum()/Pol.shape[0])
	return polar,rfunc
    
"""
	To polar coordinates
	Invert
	---> PES sort of OK but need smoothing
	
	Need to build display of the inversion result
	Need Save
	Need status bar
	Need stability
"""
def testinvert(coef):
	r=250
	deltar=1.
	nR=383.
	NL=2
	pmax=0.
	width=2.*(Bwidth**2)
	irmax=round((r/deltar)**2/nR)
	angular=np.zeros((NL,irmax),dtype='float')
	pes=np.zeros((irmax),dtype='float')
	
	print angular.shape
	for ir in np.arange(irmax):
		rad=np.sqrt(ir*nR)
		kmin=0#max(0,int(rad/Bspace)-5)
		kmax=Funcnumber#min(Funcnumber,int(rad/Bspace)+6)
		for ik in np.arange(kmin,kmax):
			func=np.exp(-(rad-ik*Bspace)**2 /width)
			for il in np.arange(NL):
				angular[il][ir]+=coef[ik*NL+il]*func
		a0=angular[0][ir]
		
		if a0<=0.: 
			a0=0.
			for il in np.arange(NL):
				angular[il][ir]=0.
		else:
			for il in np.arange(1,NL):
				angular[il][ir]/=a0
		
		if a0*rad>pmax: pmax=a0*rad
		pes[ir]=a0*rad
	
	return angular,pes/pmax
		
	