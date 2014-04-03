import numpy as np
import scipy as sp
import pyfits
from os.path import isfile,join
from scipy.ndimage.measurements import center_of_mass as com
from scipy.special import eval_legendre

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
        r=np.sqrt(X**2+Y**2)
        X-=250 #Center the ring
        Y-=250 #Center the ring
        theta=theta_f(X,Y)
        r=np.sqrt(X**2+Y**2)
        self.datas=np.exp(-(r-80)**2/50)*eval_legendre(2,np.cos(theta))
        self.datas[self.datas<0]=0.
        self.datas/=self.datas.max()
        #self.datas[self.datas<0.]=0.
        self.raw=self.datas
    
        self.center=(128.,128.)
        self.get_com()
        self.r=100.
        self.scale=ArrayInfos(self.datas)
        
        #Output
        self.normed_pes=np.zeros(Rbin)
        self.ang=np.zeros((self.get_NumberPoly(),Angbin))
        self.output=np.zeros_like(self.datas)
        self.pes_error=np.zeros(Rbin)
        self.ang_var=np.zeros((self.get_NumberPoly(),Angbin))
    
    def reset(self):
    	self.raw=np.zeros((256,256))
    	self.datas=self.raw
    	self.center=(128.,128.)
    	self.r=100.
    	self.normed_pes=np.zeros(Rbin)
        self.ang=np.zeros((self.get_NumberPoly(),Angbin))
        self.output=np.zeros_like(self.datas)
        self.pes_error=np.zeros(Rbin)
        self.ang_var=np.zeros((self.get_NumberPoly(),Angbin))
    	
    
    def get_NumberPoly(self):
    	if not self.odd: return np.arange(0,self.lmax+1,2).shape[0]
    	else: return np.arange(0,self.lmax+1,1).shape[0]
        
    def OpenFile(self,filepath):
    	self.reset()
        if filepath[-4:]=='.fit':
           hdulist = pyfits.open(filepath)
           #hdulist.info()
           scidata = hdulist[0].data
           hdulist.close()
           ind=np.where(scidata.ravel()<0)[0]
           scidata.ravel()[ind]=0
           self.raw=scidata
           self.datas=scidata
           self.get_com()
           self.scale=ArrayInfos(self.datas)
        elif filepath[-4:]=='.dat' or filepath[-4:]=='.txt':
            scidata = np.loadtxt(filepath)
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
    	hdu_pes=pyfits.ImageHDU(np.array([np.arange(self.normed_pes.shape[0]),self.normed_pes,self.pes_error]))
    	hdu_ang=pyfits.ImageHDU(self.ang)
    	hdu_ang_var=pyfits.ImageHDU(self.ang_var)
    	hdu_out=pyfits.PrimaryHDU(self.output)
    	col1=pyfits.Column(name='Max Legendre',format='I',array=np.array([self.lmax]))
    	col2=pyfits.Column(name='Image Center',format='2I',array=np.array(self.center))
    	col3=pyfits.Column(name='Selection radius',format='I',array=np.array([self.r]))
        col4=pyfits.Column(name='odd',format='I',array=np.array([self.odd]))
        tbhdu=pyfits.new_table([col1,col4,col2,col3])#, tbtype='TableHDU')
        hdulist=pyfits.HDUList([hdu_out,hdu_pes,hdu_ang,hdu_ang_var,tbhdu])
        hdulist.writeto(filepath,clobber=True)
    
    def SaveFileDat(self,filepath):
    	root=filepath[:-4]
    	np.savetxt(root+'_pes.dat',np.hstack((np.arange(Rbin),self.normed_pes,self.pes_error)).reshape((3,Rbin)),fmt='%f')
    	np.savetxt(root+'_img_inv.dat',self.output,fmt='%f')
    	for beta in np.arange(1,self.get_NumberPoly()):
    		if self.odd: i=beta
    		else: i=beta+1
    		np.savetxt(root+'_ang_b'+str(i)+'.dat',np.hstack((np.arange(Rbin),self.ang[beta,:],self.ang_var[beta,:])).reshape((3,Rbin)))
        
        
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
    
    def LoadBasis(self,path):
    	"""
    		Look for Basis file and load it. We need only the initial basis for leastsq
    	"""
    	#Build basis file name
    	if not self.odd: filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1,2)])
    	else:filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1)])
    	fileextension='.dat'
    	fname='U'+filenamescore+fileextension
    	Ufilename=join(path,fname)
    	#Generate Basis file if not existing
    	if not (isfile(Ufilename)): 
    	    print 'files missing'
    	    return -1
    	
    	#Loading file as a vector as in PBaseX.
    	M=Rbin*Angbin  #correspond to the polar array dimension in the basis file NR=440 NTH=440
    	N=self.get_NumberPoly()*Funcnumber #Number of polynoms by the number of function set to NR/2.
    	u=np.fromfile(Ufilename,dtype=np.float64).reshape((M,N))
    	
    	return u
    
    def to_polar(self,data):
    	return (cart2pol(data,self.scale,self.center,self.r),cart2pol_var(data,self.scale,self.center,self.r))
    	
    def Invert(self,polar,u):
        #Resolve Basis*coeff=polar_dat
        self.coefficients,residuts,siz,singulars=np.linalg.lstsq(u,polar[0])       #Need to build up angular distribution and PES
        self.coefficients_var,residuts,siz,singulars=np.linalg.lstsq(u*u,polar[1])
        del residuts,siz,singulars,polar,u
                         
        width=2*Bwidth**2
        rad=np.arange(Rbin)
        kvector=np.arange(Funcnumber)
        
        #Angular Matrix is of size NLxNR where NL is the number of Legendre polynoms and Nr the radial binning
        NL=self.get_NumberPoly()
        angular=np.zeros((NL,Angbin),dtype='float')
        angular_var=np.zeros((NL,Angbin),dtype='float')
        coefs=self.coefficients.reshape((Funcnumber,NL))
        coefs_var=self.coefficients_var.reshape((Funcnumber,NL))

        for r in rad:
            fradial=np.exp(-(r*Rbin/float(self.r)-Bspace*kvector)**2/width)
            for l in np.arange(NL):
                angular[l][r]=sum(fradial*coefs[:,l])
                angular_var[l][r]=sum(fradial*fradial*coefs_var[:,l])
        
        ang0=angular[0,:]
        ang0_var=angular_var[0,:]
        angular[:,ang0<=1e-8]=0.0
        angular_var[:,ang0<=1e-8]=0.0
        ang0[ang0<=1e-8]=0.0
        
        a0=ang0[ang0>1e-6]
    	angular[1:,ang0>1e-6]/=a0
        angular_var[1:,ang0>1e-6]/=(a0**2)
        angular_var[1:,ang0>1e-6]+=(angular[1:,ang0>1e-6]**2)*ang0_var[ang0>1e-6]/(a0**4)
        
        pmax=(ang0*rad**2).max()
        self.normed_pes=ang0*rad**2 /pmax
        angular_var*=(rad**2 /pmax)**2
        self.ang=angular
        self.ang_var=np.sqrt(np.abs(angular_var))
    	self.pes_error=self.ang_var[0,:]
        	
    def image_for_display(self):
    	dim=int(self.r+1)
    	#Calculate new image in cartesian coordinates and return it for display
    	X,Y=np.meshgrid(np.arange(-dim,dim+1),np.arange(-dim,dim+1))
    	new_r=np.sqrt(X**2+Y**2).ravel()
    	new_t=theta_f(X,Y).ravel()
        del X,Y
        #List of legendre polynoms
        leg=ldist(new_t,self.lmax,self.odd,self.coefficients)
        temp=np.concatenate([k*np.ones(self.get_NumberPoly()) for k in np.arange(Funcnumber)])
        K,Rad=np.meshgrid(temp,new_r)
        del temp,new_t
    	
    	Rad*=Rbin/float(self.r)
    	func=np.exp(-(Rad-K*Bspace)**2/(2*Bwidth**2))
    	del Rad,K
    	outdata=(func*leg).sum(axis=1)*new_r
    	outdata[outdata<0]=0
    	self.output=np.zeros_like(self.raw)
    	self.output[self.center[1]-dim:self.center[1]+dim+1,self.center[0]-dim:self.center[0]+dim+1]=outdata.reshape((2*dim+1,2*dim+1))
    	
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

def cart2pol_var(data,scale,center,r):

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
	polar+=(dat*cub*cub).sum(axis=1)
	return polar

def theta_f(x,y):
	ang=np.arctan2(np.fabs(x),np.fabs(y))
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

def ldist(X,lmax,odd,coeff):
	if(odd): NL=lmax+1
	else: NL=lmax/2+1
	x=np.cos(X)
	
	K,T=np.meshgrid(np.arange(Funcnumber),x)
	del x,K
	pl=np.zeros((T.shape[0],T.shape[1],lmax+1))
	
	pl[:,:,0]=np.ones_like(T)
	pl[:,:,1]=T
	
	if lmax >= 2:
		twox=2.*T
		f2=T
		d=1.
		f1=d
		for j in np.arange(2,lmax+1):
			d+=1.
			f2 += twox
			pl[:,:,j]=(f2*pl[:,:,j-1]-f1*pl[:,:,j-2])/d
			f1+=d
		if odd: return pl.reshape((X.shape[0],NL*Funcnumber))*coeff
		else: return pl[:,:,::2].reshape((X.shape[0],NL*Funcnumber))*coeff
	else: 
		if odd: return pl.reshape((X.shape[0],NL*Funcnumber))*coeff
		else: return pl[:,:,0].reshape((X.shape[0],NL*Funcnumber))*coeff