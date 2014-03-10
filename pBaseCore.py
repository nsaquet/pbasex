import numpy as np
import scipy as sp
from os.path import isfile
from scipy.ndimage.interpolation import map_coordinates as mpc
from scipy.ndimage.interpolation import spline_filter
from scipy.ndimage.measurements import center_of_mass as com
from numpy.polynomial.legendre import legval

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
			self.Rfact=nr/float(Rbin)
			nr=Rbin
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
        #Simulate a ring image
        x=np.arange(0,500)
        y=np.arange(0,500)
        X,Y=np.meshgrid(x,y) #Create a 2d map
        X-=250 #Center the ring
        Y-=250 #Center the ring
        self.theta=theta_f(X,Y)
        r=np.sqrt(X**2+Y**2)
        self.datas=np.exp(-(r-80)**2/50)*legval(np.cos(self.theta),[1,0,1])
        self.datas/=self.datas.max()
        self.raw=self.datas
    
        self.center=(0.,0.)
        self.get_com()
        self.r=100.
        self.scale=ArrayInfos(self.datas)
        
        #Output
        self.normed_pes=np.zeros(Rbin)
        self.ang=np.zeros((Angbin,self.get_NumberPoly()))
        
    
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
        
    def SaveFileFits(self,filepath):
    	hdu_pes=pyfits.PrimaryHDU(zip(np.arange(self.normed_pes.shape[0]),self.normed_pes))
    	hdu_ang=pyfits.PrimaryHDU(self.normed_ang)
    	col1=pyfits.Column(name='Max Legendre',format='E',array=np.array([self.l]))
    	col2=pyfits.Column(name='Image Center',format='E',array=np.array([self.center]))
    	col3=pyfits.Column(name='Selection radius',format='E',array=np.array([self.r]))
        col4=pyfits.Column(name='odd',format='E',array=np.array([self.odd]))
        cols=pyfits.ColDefs([col1,col4,col2,col3])
        tbhdu=pyfits.new_table(cols)
        hdulist=pyfits.HDUList([hdu_pes,hdu_ang,tbhdu])
        hdulist.writeto(filepath)
        
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
    	M=Rbin*Angbin  #correspond to the polar array dimension in the basis file NR=440 NTH=440
    	N=self.get_NumberPoly()*Funcnumber #Number of polynoms by the number of function set to NR/2.
    	u=np.fromfile(Ufilename,dtype=np.float64).reshape((M,N))
    	
    	return u
    
    def Invert(self,u):
    	"""
    		Load Basis and do the inversion.
    		Generate PES.
    	"""
    	rdata=self.raw
    	#Symmetrize the problem
    	#self.Symmetrize(rdata)
    	
    	#Convert into polar coordinates
        polar =cart2pol(rdata,self.scale,self.center,self.r)
        #polar_dat=np.concatenate((polar,np.zeros(u.shape[0]-len(polar),dtype=polar.dtype)))
        
        
        #Resolve Basis*coeff=polar_dat
        coefficients,residuts,siz,singulars=np.linalg.lstsq(u,polar)       #Need to build up angular distribution and PES
                
        width=2*Bwidth**2
        #irmax=min(int((self.r/self.scale.Rfact)**2 /self.scale.nR ),Rbin)
        #if irmax>Rbin: irmax=Rbin
        rad=np.arange(Rbin)#np.sqrt(np.arange(irmax)*self.scale.nR*1.)
        kvector=np.arange(Funcnumber)
        
        #Angular Matrix is of size NLxNR where NL is the number of Legendre polynoms and Nr the radial binning
        NL=self.get_NumberPoly()
        angular=np.zeros((NL,Angbin),dtype='float')
        coefs=coefficients.reshape((Funcnumber,NL))
        
        
        for r in rad:
            fradial=np.exp(-(r*Rbin/float(self.r)-Bspace*kvector)**2/width)
            for l in np.arange(NL):
                angular[l][r]=sum(fradial*coefs[:,l])
        
        ang0=angular[0,:]
        ind=np.where(ang0<0)[0]
        if len(ind)>0: 
    		ang0[ind]=0.0
    		angular[:,ind]=0.0
    	ind=np.where(ang0>0)[0]
    	angular[1:,ind]/=ang0[ind]
        pmax=(ang0*rad**2).max()
        self.normed_pes=ang0*rad**2 /pmax
        self.ang=angular
        
        


# Auxiliary function to map cartesian data to a polar plane
#Change to G. Garcia function from IGOR
def cart2pol(data,scale,center,r):

	"""
		Cubic interpolation
		Fastest implementation of the cubic interpolation so far for an array
	"""
	def cubic(y):
		p0=(y+2)**3
    		p0[p0<0]=0
    		p1 = (y+1)**3
    		p1[p1<0]=0
    		p2 = y**3
    		p2[p2<0]=0
    		p3 =(y-1)**3
    		p3[p3<0]=0
    		return (p0-4.*p1+6.*p2-4.*p3)/6.
	
	"""
		Adapt the  selected area size to polar basis size
	"""
	nR=min(r,scale.nR)
	if nR>Rbin: nR=Rbin #Security
	
	Rfact=nR/float(Rbin)
	scale.Rfact=Rfact
	
	polar=np.array([])
	
	rad=Rfact*np.concatenate([r*np.ones(2*r+1) for r in np.arange(Rbin)])
	theta=np.concatenate([np.pi*np.arange(t)/t for t in 2*np.arange(Rbin)+1])
	x=rad*np.cos(theta) + center[1]
	y=-rad*np.sin(theta) + center[0]
	#Cubic interpolation
	xpix=x/scale.Xfact
	ypix=y/scale.Yfact
	ix=xpix.astype('int')
	iy=ypix.astype('int')
	dx=xpix-ix
	dy=ypix-iy
	polar=np.zeros_like(rad)
	
	i,j=np.meshgrid(np.arange(-1,3),np.arange(-1,3))
	i=i.ravel()
	j=j.ravel()
	I,IX=np.meshgrid(i,ix)
	J,IY=np.meshgrid(i,iy)
	IXI=IX+I
	IYJ=IY+J
	rule1=(IXI<scale.nX) & (IXI>=0)
	rule2=(IYJ<scale.nY) & (IYJ>=0)
	cub=cubic(I-dx.reshape((dx.shape[0],1)))*cubic(dy.reshape((dy.shape[0],1))-J)
	dat=data.ravel()[(IXI*scale.nY+IYJ).ravel()].reshape(IYJ.shape)
	polar+=(dat*cub).sum(axis=1)
	return polar


def theta_f(x,y):
	ang=np.arctan(np.fabs(x)/np.fabs(y))
	#Check for singularity and remap angles over [0 2pi]
	indx=np.where(np.ravel(x)==0)[0]
	indy=np.where(np.ravel(y)==0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]=2.*np.pi
	
	indx=np.where(np.ravel(x)==0)[0]
	indy=np.where(np.ravel(y)>0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]=0.
	
	indx=np.where(np.ravel(x)==0)[0]
	indy=np.where(np.ravel(y)<0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]=np.pi
	
	indx=np.where(np.ravel(x)>0)[0]
	indy=np.where(np.ravel(y)==0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]=np.pi/2.
	
	indx=np.where(np.ravel(x)<0)[0]
	indy=np.where(np.ravel(y)==0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]=3.*np.pi/2.
	
	indx=np.where(np.ravel(x)<0)[0]
	indy=np.where(np.ravel(y)<0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]+=np.pi
	
	indx=np.where(np.ravel(x)<0)[0]
	indy=np.where(np.ravel(y)>0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]*=-1
	np.ravel(ang)[np.intersect1d(indx,indy)]+=2.*np.pi
	
	indx=np.where(np.ravel(x)>0)[0]
	indy=np.where(np.ravel(y)<0)[0]
	np.ravel(ang)[np.intersect1d(indx,indy)]*=-1
	np.ravel(ang)[np.intersect1d(indx,indy)]+=np.pi
	return ang
    
"""
	Need to build display of the inversion result
	Need Save
	Need status bar
	Need stability
"""




#Old function slow but working
"""
def cart2pol_old(data,scale,center,r):
		#Cubic interpolation
		#Fastest implementation of the cubic interpolation so far
		#Faster than max(0,(x+1)**3) by 20% or any of the factorization by bool (x+2>0) by 60-65%
	def cubic(x):
		p0=lambda y: (y+2)**3 if (y+2>0) else 0
		p1=lambda y: (y+1)**3 if (y+1>0) else 0
		p2=lambda y: (y)**3 if (y>0) else 0
		p3=lambda y: (y-1)**3 if (y-1>0) else 0
		return (p0(x)-4*p1(x)+6*p2(x)-4*p3(x))/6.
	
	
	#Adapt the  selected area size to polar basis size
	
	nR=min(r,scale.nR)
	if nR>Rbin: nR=Rbin #Security
	
	Rfact=nR/float(Rbin)
	scale.Rfact=Rfact
	print Rfact
	polar=np.array([])
	rfunc=[]
	
	#Calculate the polar map
	for l in np.arange(Rbin):
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
		
		for index in np.arange(nth):
			for i in np.arange(-1,3):
				for j in np.arange(-1,3):
					condx=np.logical_and((ix[index]+i)<scale.nX,(ix[index]+i)>=0)
					condy=np.logical_and((iy[index]+j)<scale.nY,(iy[index]+j)>=0)
					if np.logical_and(condx,condy):
						cub=cubic(i-dx[index])*cubic(dy[index]-j)
						#print i-dx[index],cubic(i-dx[index]) ,cub
						Pol[index]+=data.ravel()[(ix[index]+i)*scale.nY+(iy[index]+j)]*cub	
					
		polar=np.concatenate((polar,Pol))
		rfunc.append(Pol.sum()/Pol.shape[0])
	return polar,rfunc
"""