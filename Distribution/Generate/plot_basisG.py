# plot basis
import numpy as np
import pylab

N=440
rad=np.concatenate([r*np.ones(2*r+1) for r in np.arange(N)])
theta=np.concatenate([np.pi*np.arange(t)/t for t in 2*np.arange(N)+1])
x=rad*np.cos(theta)
y=rad*np.sin(theta)

def plot_basis(file,nl,k,coord):
	basis=np.fromfile(file,np.float64).reshape((N**2,N*nl/2))
	
	if coord=='pol':
		pylab.figure()
		pylab.scatter(theta,rad,marker='.',c=basis[:,k],alpha=0.5,hold=False,edgecolor='none')
		pylab.xlim(0,np.pi)
		pylab.ylim(0,440)
		pylab.show()
	
	else:
		pylab.figure()
		pylab.scatter(x,y,marker='.',c=basis[:,k],alpha=0.5,edgecolor='none')
		pylab.scatter(x,-y,marker='.',c=basis[:,k],alpha=0.5,edgecolor='none')
		pylab.xlim(-rad.max(),rad.max())
		pylab.ylim(-rad.max(),rad.max())
		pylab.axis('equal')
		pylab.show()