#!/usr/bin/python

import sys
import functools
import numpy as np
import matplotlib.pyplot as plt

class appinputs:
	minmaxJobsPerSec=[]
	numOfApps=0
	numOfAvailImplPerApp=[]
	availImplPerApp={}
	minmaxInsPerApp={}
	minmaxVMPerApp={}
	minmaxProcPerVM={}
	minmaxMemPerVM={}
	minmaxStoPerVM={}
	minmaxNetPerVM={}
	minmaxActP={}
	minmaxActM={}
	minmaxActN={}
	accelerator={}
	rhoAcc={}
	typeOfActP=[]
	typeOfActM=[]
	typeOfActN=[]
	
	

class brinputs:
	pollInterval=0

class powinputs:
	accelerator=0
	typeCpu=0
	cpuPmin=0.0
	cpuPmax=0.0
	numOfPoints=0
	cpubins=[]
	cpuP=[]
	cpuC=0.0
	typeAcc=0
	accPmin=0.0
	accPmax=0.0
	accC=0.0

class resinputs:
	rtype=-1
	numOfProcUnits=0.0
	totMem=0.0
	totSto=0.0
	overCommitProc=0.0
	overCommitMem=0.0
	compCap=0.0
	accelerator=0
	totAcc=0
	accCompCap=0.0
	
class netinputs:
	netBW=0
	overCommitNet=0	

class cellinputs:
	numOfTypes=0
	types=[]
	numOfResourcesPerType=[]
	ni=[]
	ri=[]
	pi=[]
	bi=[]

class siminputs:
	maxTime=0.0
	upInterval=0.0
	numOfCells=0
	ci=[]
	ai=[]

print "Loading Simulation Inputs and Parameters..."

si=siminputs();

f = open('../input/main','r')
f2 = open('../input/broker','r')
f3 = open('../input/app','r')
si.maxTime=float(f.readline())
si.upInterval=float(f.readline())
si.numOfCells=int(f.readline())

for i in range(0,si.numOfCells):
	si.ci.append(cellinputs())
	si.ci[i].numOfTypes=int(f.readline())
	si.ci[i].ni.append(netinputs())
	si.ci[i].ni[0].netBW=float(f.readline())
	si.ci[i].ni[0].overCommitNet=float(f.readline())
	si.ci[i].bi.append(brinputs())
	si.ci[i].bi[0]=float(f2.readline())
	for j in range(0,si.ci[i].numOfTypes):
		si.ci[i].types.append(int(f.readline()))
		si.ci[i].numOfResourcesPerType.append(int(f.readline()))
		si.ci[i].ri.append(resinputs())
		si.ci[i].ri[j].rtype=si.ci[i].types[j]
		si.ci[i].ri[j].numOfProcUnits=float(f.readline())
		si.ci[i].ri[j].totMem=float(f.readline())
		si.ci[i].ri[j].totSto=float(f.readline())
		si.ci[i].ri[j].overCommitProc=float(f.readline())
		si.ci[i].ri[j].overCommitMem=float(f.readline())
		si.ci[i].ri[j].compCap=float(f.readline())
		si.ci[i].ri[j].accelerator=int(f.readline())
		si.ci[i].ri[j].totAcc=int(f.readline())
		si.ci[i].ri[j].accCompCap=float(f.readline())
		si.ci[i].pi.append(powinputs())
		si.ci[i].pi[j].accelerator=si.ci[i].ri[j].accelerator
		si.ci[i].pi[j].typeCpu=int(f.readline())
		si.ci[i].pi[j].cpuPmin=float(f.readline())
		si.ci[i].pi[j].cpuPmax=float(f.readline())
		si.ci[i].pi[j].numOfPoints=int(f.readline())
		if (si.ci[i].pi[j].numOfPoints>0):
			si.ci[i].pi[j].cpubins=(map((lambda x: float(x)),(f.readline()).split()))
			si.ci[i].pi[j].cpuP=(map((lambda x: float(x)),(f.readline()).split()))
		else:
			dummy=float(f.readline())
			dummy=float(f.readline())
		si.ci[i].pi[j].cpuC=float(f.readline())
		si.ci[i].pi[j].typeAcc=int(f.readline())
		si.ci[i].pi[j].accPmin=float(f.readline())
		si.ci[i].pi[j].accPmax=float(f.readline())
		si.ci[i].pi[j].accC=float(f.readline())

f.close()
f2.close()

si.ai.append(appinputs);
si.ai[0].minmaxJobsPerSec=(map((lambda x: float(x)),(f3.readline()).split()))
si.ai[0].numOfApps=int(f3.readline())

for i in range(0,si.ai[0].numOfApps):
	si.ai[0].numOfAvailImplPerApp.append(int(f3.readline()))
	si.ai[0].availImplPerApp[i,0]=(map((lambda x: int(x)),(f3.readline()).split()))
	si.ai[0].minmaxInsPerApp[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	si.ai[0].minmaxVMPerApp[i,0]=(map((lambda x: int(x)),(f3.readline()).split()))
	si.ai[0].minmaxProcPerVM[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	si.ai[0].minmaxMemPerVM[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	si.ai[0].minmaxStoPerVM[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	si.ai[0].minmaxNetPerVM[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	t=(map((lambda x: int(x)),(f3.readline()).split()))
	si.ai[0].typeOfActP.append(t[0])
	si.ai[0].typeOfActM.append(t[1])
	si.ai[0].typeOfActN.append(t[2])
	si.ai[0].minmaxActP[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	si.ai[0].minmaxActM[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	si.ai[0].minmaxActN[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))
	si.ai[0].accelerator[i,0]=(map((lambda x: int(x)),(f3.readline()).split()))
	si.ai[0].rhoAcc[i,0]=(map((lambda x: float(x)),(f3.readline()).split()))

f3.close()

print "Loading Simulation Results..."

f = open('../output/output11','r')
data = [tuple(float(n) for n in line.split()) for line in f]
data = np.array(data)


#[:,0]: Time Steps
#[:,1]: activeSrvs
#[:,2]: Active VMs
#[:,3]: accTasks
#[:,4]: rejTasks
#[:,5]: availProc
#[:,6]: utilProc
#[:,7]: autilProc
#[:,8]: totProc
#[:,9]: phyProc
#[:,10]: availMem
#[:,11]: utilMem
#[:,12]: autilMem
#[:,13]: totMem
#[:,14]: phyMem
#[:,15]: availSto
#[:,16]: utilSto
#[:,17]: availSto
#[:,18]: phySto
#[:,19]: availNetw
#[:,20]: utilNetw
#[:,21]: autilNetw
#[:,22]: totNetw
#[:,23]: phyNetw
#[:,24]: availAcc
#[:,25]: utilAcc
#[:,26]: totAcc
#[:,27]: totPcons




plt.plot(data[:,0],np.divide(data[:,6],data[:,1])/(44*1.5))
#plt.plot(data[:,0],data[:,8],'r')
plt.show()

