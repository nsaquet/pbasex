import numpy as np
from scipy.special import eval_legendre
from scipy.weave import converters
from PySide import QtCore, QtGui

def ldist(X,lmax, odd):
	if(odd): NL=lmax+1
	else: NL=lmax/2+1
	
	pl=np.zeros((lmax+1,X.shape[0]))
	x=np.cos(X)
	pl[0]=np.ones_like(X)
	pl[1]=x
	
	if lmax >= 2:
		twox=2.*x
		f2=x
		d=1.
		f1=d
		for j in np.arange(2,lmax+1):
			d+=1.
			f2 += twox
			pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d
			f1+=d
	if odd: return pl
	else: return pl[::2,:]

def ldist_coeff(X,lmax,odd,Funcnumber):
	if(odd): NL=lmax+1
	else: NL=lmax/2+1
	x=np.cos(X)
	
	#K,T=np.meshgrid(np.arange(Funcnumber),x)
	T,K=np.meshgrid(x,np.arange(Funcnumber))
	
	pl=np.zeros((lmax+1,T.shape[0],T.shape[1]))
	
	pl[0]=np.ones_like(T)
	pl[1]=T
	
	if lmax >= 2:
		twox=2.*T
		f2=T
		d=1.
		f1=d
		for j in np.arange(2,lmax+1):
			d+=1.
			f2 += twox
			pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d
			f1+=d
		
		if odd: return pl.reshape((NL*Funcnumber,X.shape[0]))
		else: return pl[::2].reshape((NL*Funcnumber,X.shape[0]))
	else: 
		if odd: return pl.reshape((NL*Funcnumber,X.shape[0]))
		else: return pl[0].reshape((NL*Funcnumber,X.shape[0]))
		
def inlineTimeStep(self, dt=0.0):
        """Takes a time step using inlined C code -- this version uses
        blitz arrays."""
        g = self.grid
        nx, ny = g.u.shape
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u

        code = """
               #line 120 "laplace.py" (This is only useful for debugging)
               double tmp, err, diff;
               err = 0.0;
               for (int i=1; i<nx-1; ++i) {
                   for (int j=1; j<ny-1; ++j) {
                       tmp = u(i,j);
                       u(i,j) = ((u(i-1,j) + u(i+1,j))*dy2 +
                                 (u(i,j-1) + u(i,j+1))*dx2)*dnr_inv;
                       diff = u(i,j) - tmp;
                       err += diff*diff;
                   }
               }
               return_val = sqrt(err);
               """
        # compiler keyword only needed on windows with MSVC installed
        err = weave.inline(code,
                           ['u', 'dx2', 'dy2', 'dnr_inv', 'nx', 'ny'],
                           type_converters=converters.blitz,
                           compiler = 'gcc')
        return err