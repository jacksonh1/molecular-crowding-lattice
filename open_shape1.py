'''
created by: Jackson H
started 11/09/2016
'''


import numpy as np
import matplotlib.pyplot as plt
np.random.seed(123)
#np.set_printoptions(threshold='nan',precision=0)

#initialization should be up here
#the export situation is ridiculous

def setup(n,k,c):#generates nxn sized lattice with randomly distributed cxc sized crowders 
	lat=np.zeros(n*n).reshape(n,n)
	while np.sum(lat) != k:
		x=np.random.randint(0,n-c+1)#draws randomly between int. not including high int in distri
		y=np.random.randint(0,n-c+1)
		if np.all(lat[x:x+c,y:y+c]==np.zeros([c,c])):
			lat[x:x+c,y:y+c]=np.ones([c,c])
		#else:
			#print 'blocked'
	return lat

def latprint(lat):# prints the crowded lattice nicely
	for i in lat:
		for t in i:
			print int(t),
		print ''

def fitO(lat,O1,O2,H90o,n):# counts how many times open state fits in crowded lattice
	Wo=0.0
	Ri=n-O1+1#Ri and Ci are the row/column positions that are available to the protein without hitting lattice edge
	Ci=n-O2+1
	for i in range(0,Ri):
		for t in range(0,Ci):
			check=lat[i:i+O1,t:t+O2]+H90o#sum the section of lattice and H90o shape. a value of 2 is a block
			if np.all(check!=2):
				Wo=Wo+1.0
	return Wo

def fitC(lat,C1,C2,H90c,n):# counts how many times closed state fits in crowded lattice
	Wc=0.0
	Ri=n-C1+1
	Ci=n-C2+1
	for i in range(0,Ri):
		for t in range(0,Ci):
			check=lat[i:i+C1,t:t+C2]+H90c#sum the section of lattice and H90o shape. a value of 2 is a block
			if np.all(check!=2):
				Wc=Wc+1.0
	return Wc

def calc(Wc,Wo,n,k):
	R=float(Wc)/float(Wo)
	dS=np.log(Wc/Wo)# The 'entropy' difference from the # of configurations of the two states
	Fvc=float(k)/(n*n)# fraction of 'volume' that is excluded by the crowders
	return R,dS,Fvc

n=200#size of the matrix: nxn matrix
c=np.arange(1,8)# size of the crowder: cxc crowder


#open Shape 1
O1=6#vertical dimension of open state
O2=10#horizontal dimension of open state
H90o=np.zeros([O1,O2])#open state shape
H90o[4:6,4:6]=np.ones([2,2])
m=0
for i in range(0,4):
    H90o[m:m+2,m]=1
    m=m+1
m=0
m2=3
for i in range(0,4):
    H90o[m2:m2+2,m+6]=1
    m2=m2-1
    m=m+1
latprint(H90o)

C1=6#vertical dimension of closed state
C2=4#horizontal dimension of closed state
H90c=np.zeros([C1,C2])
one=np.ones(4)
H90c[:4]=H90c[:4]+one
H90c[4:6,1:3]=np.ones([2,2])
latprint(H90c)

loops=25#number of iterations of crowded lattice creation, open/closed counting, and calculations for each crowder size

filename='open_shape1'

fracV=np.array([0,0.025,0.05,0.075,0.1,0.125,0.15,0.2])#fraction excluded volumes

expo=np.zeros([len(c),len(fracV)*4])#size series for each concentration
expo2=np.zeros([len(fracV),len(c)*4])#concentration series for each size
l2=0
for j in range(0,len(fracV)):
	
	expo3=np.zeros([loops,len(c)*6])

	k=np.zeros([len(c)])
	k=[(int(fracV[j]*(n*n))/(i*i))*(i*i) for i in c]#k is the # of lattice sites that are excluded.
	# the integer calculation of k makes sure that k is a number divisible by the crowder size. otherwise
	#setup will loop infinitely
	
	l=0
	l3=0
	for t in range(0,len(c)):
		dS=np.zeros(loops)
		R=np.zeros(loops)
		Wo=np.zeros(loops)
		Wc=np.zeros(loops)
		print c[t]
		for i in range(0,loops):#this should be a function
			lat=setup(n,k[t],c[t])
			Wo[i]=fitO(lat,O1,O2,H90o,n)
			Wc[i]=fitC(lat,C1,C2,H90c,n)
			R[i],dS[i],Fvc=calc(Wc[i],Wo[i],n,k[t])
			aveR=np.sum(R)/float(loops)
		print R,aveR,Fvc
		expo3[:,l3]=c[t]
		expo3[:,l3+1]=np.arange(loops)
		expo3[:,l3+2]=Wo
		expo3[:,l3+3]=Wc
		expo3[:,l3+4]=R
		expo3[:,l3+5]=dS
		expo[t,l2]=Fvc#size series for each concentration
		expo[t,l2+1]=c[t]
		expo[t,l2+2]=aveR
		expo[t,l2+3]=np.std(R)
		expo2[j,l]=c[t]#concentration series for each size
		expo2[j,l+1]=Fvc
		expo2[j,l+2]=aveR
		expo2[j,l+3]=np.std(R)
		l3=l3+6
		l=l+4
	plt.plot(expo[:,l2+1],expo[:,l2+2],'o')
	plt.savefig(filename+str(fracV[j])+'.png', dpi=100,format='png')
	
	l2=l2+4

	np.savetxt(filename+'_'+str(fracV[j])+'trials'+'.CSV',expo3,delimiter=',',header='c,trial #,Wo,Wc,R,dS,'*len(c))
	#plt.show()
np.savetxt(filename+'size_series'+'.CSV',expo,delimiter=',',header='Fvc=%r,c,AveR,SD_R,'*len(fracV)%(0,0.025,0.05,0.075,0.1,0.125,0.15,0.2))
np.savetxt(filename+'ex_vol_series'+'.CSV',expo2,delimiter=',',header='c=%r,Fvc,AveR,SD_R,'*len(c)%(1,2,3,4,5,6,7))
