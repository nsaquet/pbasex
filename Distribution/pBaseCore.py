import numpy as np
import scipy as sp
import pyfits
from os.path import isfile,join
from scipy.ndimage.measurements import center_of_mass as com
from scipy.special import eval_legendre

from scipy.weave import converters
from scipy import weave

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
		self.ellipticity=1.0
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
        self.half=0
        self.raw=2.*np.random.normal(0.5,size=(256,256))
        #Simulate a ring image
        x=np.arange(0,600)
        y=np.arange(0,500)
        X,Y=np.meshgrid(x,y) #Create a 2d map
        r=np.sqrt(X**2+Y**2)
        X-=280 #Center the ring
        Y-=236 #Center the ring
        theta=theta_f(X,Y)
        r=np.sqrt(X**2+Y**2)
        self.datas=np.zeros_like(r)
        #self.datas[np.all([(r>75), (r<95)], axis=0)]=1-1.*eval_legendre(2,np.cos(theta[np.all([(r>75), (r<95)], axis=0)])) #Test distribution to convert to polar
        self.datas=np.exp(-(r-140)**2/30)*(1.+1*eval_legendre(2,np.cos(theta))) +np.exp(-(r-60)**2/5)*(1.-1*eval_legendre(2,np.cos(theta))) +np.exp(-(r-80)**2/20)*(1.+2*eval_legendre(2,np.cos(theta)))
        self.datas[self.datas<0]=0.
        self.datas/=self.datas.max()
        #self.datas[self.datas<0.]=0.
        self.raw=self.datas
    
        self.center=(128.,128.)
        self.get_com()
        self.r=100.
        self.scale=ArrayInfos(self.datas)
        
        
        #Output
        self.radial=np.zeros(Rbin)
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
        if filepath[-4:]=='.fit' or filepath[-5:]=='.fits' :
        
           hdulist = pyfits.open(filepath)
           #hdulist.info()
           scidata = hdulist[0].data
           try:
           	len(scidata)
           except:
           	print('Entry 0 not an image')
           	scidata = hdulist[1].data
           	len(scidata)
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
    	hdu_pes=pyfits.ImageHDU(np.array([self.radial,self.normed_pes,self.pes_error]))
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
    	np.savetxt(root+'_pes.dat',np.hstack((self.radial,self.normed_pes,self.pes_error)).reshape((3,Rbin)),fmt='%f')
    	np.savetxt(root+'_img_inv.dat',self.output,fmt='%f')
    	for beta in np.arange(1,self.get_NumberPoly()):
    		if self.odd: i=beta
    		else: i=2*beta
    		np.savetxt(root+'_ang_b'+str(i)+'.dat',np.hstack((self.radial,self.ang[beta,:],self.ang_var[beta,:])).reshape((3,Rbin)))
              
    def get_com(self):
        datmax=self.datas.max()
        mask=(self.datas>datmax*0.25)
        self.center=com(self.datas.T,mask.T)
        
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
    	
    def LoadBasis_svd(self,path):
    	"""
    		Look for Basis file and load it. We need only the initial basis for svd
    	"""
    	#Build basis file name
    	if not self.odd: filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1,2)])
    	else:filenamescore=''.join(['P%i' %i for i in np.arange(0,self.lmax+1)])
    	fileextension='.dat'
    	fnames=[i+filenamescore+fileextension for i in ['U','V','S']]
    	filenames=[join(path,i) for i in fnames]
    	#Generate Basis file if not existing
    	if not (isfile(filenames[1])): 
    	    print 'files missing'
    	    return -1
    	
    	#Loading file as a vector as in PBaseX.
    	M=Rbin*Angbin  #correspond to the polar array dimension in the basis file NR=440 NTH=440
    	N=self.get_NumberPoly()*Funcnumber #Number of polynoms by the number of function set to NR/2.
    	u=np.fromfile(filenames[0],dtype=np.float64).reshape((M,N))
    	v=np.fromfile(filenames[1],dtype=np.float64).reshape((N,N))
    	s=np.fromfile(filenames[2],dtype=np.float64)
    	
    	return u,v,s
    
    def to_polar(self,data):
    	return cart2pol(self,data)
    	
    def Invert(self,polar,u):
        #Resolve Basis*coeff=polar_dat
        self.coefficients,residuts,siz,singulars=sp.linalg.lstsq(u,polar[0],lapack_driver='gelss')       #Need to build up angular distribution and PES
        self.coefficients_var,residuts,siz,singulars=sp.linalg.lstsq(u**2,polar[1],lapack_driver='gelss')
        del residuts,siz,singulars,polar,u
                         
        width=2*Bwidth**2
        rad=np.arange(Rbin,dtype='f16')
        rad=np.sqrt(rad*self.r**2/float(Rbin))
        self.radial=rad
        kvector=np.arange(Funcnumber)
        
        #Angular Matrix is of size NLxNR where NL is the number of Legendre polynoms and Nr the radial binning
        NL=self.get_NumberPoly()
        angular=np.zeros((NL,Angbin))
        angular_var=np.zeros((NL,Angbin),dtype='f16')
        coefs=self.coefficients.reshape((Funcnumber,NL))
        coefs_var=self.coefficients_var.reshape((Funcnumber,NL))

        for ir in np.arange(Rbin,dtype='f16'):
            fradial=np.exp(-(rad[ir]*Rbin/float(self.r)-Bspace*kvector)**2/width)
            for l in np.arange(NL):
                angular[l][ir]=sum(fradial*coefs[:,l])
                angular_var[l][ir]=sum(fradial*fradial*coefs_var[:,l])

        ang0=angular[0,:]
        angular[:,ang0<=1e-3]=0.0
        angular_var[:,ang0<=1e-3]=0.0
        a0=angular[0,:]
        ind=np.where(a0>0)[0]

        a0_var=angular_var[0,ind]
    	angular[1:,ind]/=a0[ind]
        angular_var[1:,ind]= angular_var[1:,ind]/a0[ind]/a0[ind] + (angular[1:,ind]**2)*a0_var/(a0[ind]**4) 
        
        pmax=(a0*rad).max()
        angular[0,:]*=rad/pmax
        angular_var=np.sqrt(np.abs(angular_var))
        angular_var[0,:]*=rad/pmax
        
        self.normed_pes=angular[0,:]
        self.ang=angular
        self.ang_var=angular_var
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
def cart2pol(struct,data):

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
	nR=min(struct.r,struct.scale.nR)
	if nR>Rbin: nR=Rbin #Security
	
	#Norm datas before cubic interpolation:
	data/=data.max()
	
	Rfact=nR/float(Rbin)
	struct.scale.Rfact=Rfact
	
	rad=Rfact*np.concatenate([r*np.ones(2*r+1) for r in np.arange(Rbin)])
	theta=np.concatenate([np.pi*np.arange(t)/t for t in 2*np.arange(Rbin)+1])
	y=rad*np.cos(theta) + struct.center[1]
	x=-rad*np.sin(theta)*struct.scale.ellipticity + struct.center[0]
	#Cubic interpolation
	xpix=x/struct.scale.Xfact
	ypix=y/struct.scale.Yfact
	ix=xpix.astype('int')
	iy=ypix.astype('int')
	dx=xpix-ix
	dy=ypix-iy
	
	polar=np.zeros_like(rad)
	polar_var=np.zeros_like(rad,dtype='f8')
	
	i,j=np.meshgrid(np.arange(-1,3),np.arange(-1,3))
	i=i.ravel()
	j=j.ravel()
	J,IX=np.meshgrid(j,ix)
	I,IY=np.meshgrid(i,iy)
	IXJ=IX+J
	IYI=IY+I
	rule1=(IXJ<struct.scale.nX) & (IXJ>=0)
	rule2=(IYI<struct.scale.nY) & (IYI>=0)
	mask=np.all(rule1 & rule2,axis=1)
	cub=cubic(J-dx.reshape((dx.shape[0],1)))*cubic(dy.reshape((dy.shape[0],1))-I)
	dat=data.ravel()[(IYI*struct.scale.nY+IXJ).ravel()].reshape(IXJ.shape)
	polar+=(dat*cub).sum(axis=1)
	polar_var+=(dat*cub*cub).sum(axis=1)
	return (polar,polar_var)

def weave_cart2pol(data,scale,center,r):
	nR=min(r,scale.nR)
	if nR>Rbin: nR=Rbin #Security
	Rfact=(nR/float(Rbin))
	scale.Rfact=Rfact
	
	#Norm datas before cubic interpolation:
	data/=data.max()
	
	yc=center[1]
	xc=center[0]
	polar=np.zeros(Rbin**2)
	polar_var=np.zeros(Rbin**2,dtype='f8')
	ell=scale.ellipticity
	Xfact=scale.Xfact
	Yfact=scale.Yfact
	nX=scale.nX
	nY=scale.nY
	dat=np.copy(data.ravel(order='C').astype('float'))
	support_code = 	"""
						inline double cubic(double x){
    						double p0,p1,p2,p3;
    						if((x+2)>0) p0=(x+2)*(x+2)*(x+2);
   							else p0=0.0;
    						if((x+1)>0) p1=(x+1)*(x+1)*(x+1);
    						else p1=0.0;
    						if((x)>0) p2=(x)*(x)*(x);
    						else p2=0.0;
    						if((x-1)>0) p3=(x-1)*(x-1)*(x-1);
    						else p3=0.0;
    						return 1.0/6.0*(p0-4.0*p1+6.0*p2-4.0*p3);
						};
					"""
	
	code =	"""
				using namespace std;
				int l,nth,k,i,j,lowx,lowy,vi;
				double p,pp,x,y,rad,thta,xpix,ypix,dx,dy,cub;

				vi=0;
				for(l=0;l<Rbin;l++){
					nth=2*l+1;
        			rad=(double)l*Rfact; // Radius in image units
        			
        			for(k=0;k<nth;k++){
            			thta=(double)k/(double)nth*M_PI; // Theta in radians
            			y=(double)yc+rad*cos(thta)*ell;
            			x=(double)xc-rad*sin(thta);
            			// Cubic interpolation
            			xpix=(double)x;
            			ypix=(double)y;
            			lowx=(int)(xpix);
            			lowy=(int)(ypix);
            			dx=xpix-lowx;
            			dy=ypix-lowy;
            			
            			
            			p=0.0;
            			pp=0.0;
            			for (i=-1;i<=2;i++){
                			for (j=-1;j<=2;j++){
                    			if(((lowy+i)<nY)&&((lowx+j)<nX)&&((lowy+i)>=0)&&((lowx+j)>=0)){
                        			cub=cubic(j-dx)*cubic(dy-i);		
                        			p+=dat((lowy+i)*nY+(lowx+j))*cub;
                        			pp+=dat((lowy+i)*nY+(lowx+j))*cub*cub;
                    			}
                			}
            			}  
            			    
            			polar(vi)=p;
            			polar_var(vi)=pp;
            			vi++;
					}
				}
			"""
	weave.inline(code,arg_names=['polar','polar_var','dat','Rbin','Rfact','ell','Xfact','Yfact','nX','nY','yc','xc'], type_converters=converters.blitz, extra_compile_args =['-stdlib=libc++ -std=gnu++11'], compiler='gcc',support_code=support_code,verbose=1)
	return (polar,polar_var)

def weave_invert_svd(data,path,polar):
	u,v,s = data.LoadBasis_svd(path)
	wi=(1/s)*(np.dot(u.T,polar[0]))
	data.coefficients = np.dot(v,wi)
	wi_var=(1/(s*s))*(np.dot(u.T*u.T,polar[1]))
	data.coefficients_var = np.dot(v*v,wi_var) 
	
	width=2*Bwidth**2
	Rfact=data.r**2/float(Rbin)**2
	print data.r, Rbin,Rfact
	#Angular Matrix is of size NLxNR where NL is the number of Legendre polynoms and Nr the radial binning
	NL=data.get_NumberPoly()
	angular=np.zeros((NL,Rbin))
	angular_var=np.zeros((NL,Rbin))
	coeff=data.coefficients
	coeff_var=data.coefficients_var
	code =	"""
				using namespace std;
				int kmin,kmax,ir,ik,il;
				double pmax,rad,func,a0,a0_var,ai,vi;
				
				pmax=0.0;
				for(ir=0;ir<Rbin;ir++){
					rad=sqrt((double)(ir*Rfact*Rbin));
					kmin=(int)(rad/Bspace/sqrt(Rfact))-5;
					if(kmin<0) kmin=0;
					kmax=(int)(rad/Bspace/sqrt(Rfact))+6;
					if(kmax>Funcnumber) kmax=Funcnumber;
					
					for(ik=kmin;ik<kmax;ik++){
						func=exp(-(rad/sqrt(Rfact)-ik*Bspace)*(rad/sqrt(Rfact)-ik*Bspace)/width);						
						for(il=0;il<NL;il++){
							angular(il,ir)+=coeff(ik*NL+il)*func;
							angular_var(il,ir)+=coeff_var(ik*NL+il)*func*func;
						}
					}
					a0=angular(0,ir);
					a0_var=angular_var(0,ir);
					if(a0<=1e-3){
						a0=0.0;
						a0_var=0.0;
						for(il=0;il<NL;il++){
							angular(il,ir)=0.0; 
							angular_var(il,ir)=0.0;							
						}
					}
					else{
						for(il=1;il<NL;il++){
							ai=angular(il,ir);
							angular(il,ir)=ai/a0;
							vi=angular_var(il,ir);
							angular_var(il,ir)=vi/a0/a0+ai*ai/a0/a0/a0/a0*a0_var;
						}
					}
		
					if((a0*rad)>pmax) pmax=a0*rad;
				}		
				for(ir=0;ir<Rbin;ir++){
					rad=sqrt((double)(ir*Rfact*Rbin));
					angular(0,ir)*=(rad/pmax);
					for(il=0;il<NL;il++){
						angular_var(il,ir)=sqrt(angular_var(il,ir));
						}
					angular_var(0,ir)*=(rad/pmax);
				}
			"""
	weave.inline(code,arg_names=['Rbin','Rfact','Bspace','Funcnumber','width','NL','angular','angular_var','coeff','coeff_var'], type_converters=converters.blitz, extra_compile_args =['-stdlib=libc++ -std=gnu++11'], compiler='gcc',verbose=1)

	data.normed_pes=angular[0,:]
	data.ang=angular
	data.ang_var=angular_var
	data.pes_error=data.ang_var[0,:]
	return data

def weave_image_for_display(data):
	data.output=np.zeros_like(data.raw)
	nX=data.scale.nX
	nY=data.scale.nY
	yc=data.center[1]
	xc=data.center[0]
	dr=data.r
	width=2*Bwidth**2
	NL=data.get_NumberPoly()
	odd=data.odd
	coeff=data.coefficients
	Imrz=np.zeros_like(data.raw).ravel()
	pl=np.zeros(NL+1)
	Rfact=(dr/float(Rbin))
	
	support_code = 	"""
						inline double theta_f(double x,double y){
							double thta;
							thta=atan(fabs(x)/fabs(y));
							if (x<0&&y>0) thta=2.0*M_PI-thta;
							if (x>0&&y<0) thta=M_PI-thta;
							if (x<0&&y<0) thta=M_PI+thta;
							if(x==0&&y>0) thta=0.0;
							if(x==0&&y<0) thta=M_PI;
							if(y==0&&x>0) thta=M_PI/2.0;
							if(y==0&&x<0) thta=3.0*M_PI/2.0;
							if(x==0&&y==0) thta=2.0*M_PI;
							return thta;
						};
						
						inline double ldist(double X, blitz::Array<double, 1> coeff, double pl[], int npl, int odd, int ik){
							double v=0.0;
							int i,j,index,NL;
							double twox,f2,f1,d,x;
							
							if(odd) NL=npl+1;
							else NL=npl/2+1;
							
							x=(float)cos(X);
							pl[0]=1.0;
							pl[1]=x;
							if (npl >= 2) {
								twox=2.0*x;
								f2=x;
								d=1.0;
								for (j=2;j <= npl ;j++) {
									f1=d++;
									f2 += twox;
									pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
								}
							}
							for (i=0; i<=npl; i++){
								if((i%2)){
									if(odd)	v += coeff(ik*NL+i)*pl[i];
								}
								else{
									if(odd)	v += coeff(ik*NL+i)*pl[i];
									else{
										index=i/2;
										v += coeff(ik*NL+index)*pl[i];
									}
								}
							}
							return v;
						};
					"""
	code =	"""
				using namespace std;
				int kmin,kmax,ir,ik,il;
				double rad,func,dy,dx,thta,leg,sum,*pl;
				pl=(double*)calloc(NL+1,sizeof(double)); // Allocate Legendre polynomials
				
				for(int iz=0;iz<nY;iz++){
					dy=iz-(double)yc;
					for(int ix=xc;ix<nX;ix++){
						dx=ix-(double)xc;
						rad=sqrt(dx*dx+dy*dy);
						sum=0.0;
						if (rad<(int)dr){
							thta=theta_f(dx,dy);
							kmin=(int)(rad/(double)Rfact/Bspace)-5;
							if(kmin<0) kmin=0;
							kmax=(int)(rad/(double)Rfact/Bspace)+6;
							if(kmax>Funcnumber) kmax=Funcnumber;
							for(ik=kmin;ik<kmax;ik++){
								func=exp(-(rad/(double)Rfact-ik*Bspace)*(rad/(double)Rfact-ik*Bspace)/width);
								leg=ldist(thta,coeff,pl,NL,odd,ik);
								leg*=func;
								sum+=leg;
							}
							if(sum>=0) Imrz(iz*nY+ix)=sum*rad;
							else Imrz(iz*nY+ix)=0.0;
							Imrz(iz*nY+2*(double)xc-ix)=Imrz(iz*nY+ix);
						}
						else Imrz(iz*nY+ix)=0.0;
					}
				}
				free(pl);
			"""
	weave.inline(code,arg_names=['nY','nX','yc','xc','dr','Rfact','Bspace','Funcnumber','width','coeff','NL','odd','Imrz'], type_converters=converters.blitz, extra_compile_args =['-stdlib=libc++ -std=gnu++11'], compiler='gcc',support_code=support_code,verbose=1)
	
	data.output=Imrz.reshape(data.raw.shape)
	return data

	
	


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
		
def symmetrize(data,center,r):
        """
    		Symmetrise a 2_D circular selection vertically (ie about a horizontal axis). 
    		Assume that the centre is mid-pixel (x0,y0) rather than at lower left corner
    		of pixel x0,y0. Assume destination array bb[][] is pre-zeroed. Note that no 
    		symmetrisation is needed horizontally since the Legendre polynomials are 
    		already symmetric along the vertical axis. (Vertically being the polarisation 
    		axis of the ligth, or the direction of propagation in the case of cpl).
    	"""
    	#Need to build up the selected indexes within self.r
        yindex=np.arange(center[1]-r,center[1]+r,dtype=int)
        xindex=np.arange(center[0]-r,center[0]+r,dtype=int)
        for k,l in zip(xindex[round(len(xindex)/2.):],xindex[len(xindex)/2 -1::-1]): 
        	yind=np.where((k-center[0])**2+(yindex-center[1])**2<r**2)[0]
        	data.T[k,yindex[yind]]=0.5*(data.T[k,yindex[yind]]+data.T[l,yindex[yind]])
        	data.T[l,yindex[yind]]=data.T[k,yindex[yind]]
        return data
        #if len(xindex)%2: data.T[xindex[len(xindex)/2],yindex]+=data.T[xindex[len(xindex)/2],yindex]