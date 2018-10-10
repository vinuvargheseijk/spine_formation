import matplotlib.pyplot as plt
import numpy as np
IB=np.arange(0.001,40,0.005)
IBAR=[0.001,0.002,0.003,0.004,0.005,0.006,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.08,1,1.1,1.5,1.8,2,3,4,5,6,7,8,9,10,11,13,15,18,20,22,24,26,28,30,32,34,35,40]
curvature=[0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.1]
curv_IR=[0,1,2,3,4,5,6,7,8,9]
clust=[]
clust_dose=[1.,2.,3.,4.,5.,6.,7.,8.,9.]
curv_IBAR=[]
curv_clust=[]
curv=[]
Rad=[]
respto_curv=[]
respto_clust=[]
respto_clust.append(0.0)
curv.append(1.0/np.sqrt(85/(2*(0.05))))
Rad.append(np.sqrt(85/(2*(0.05))))
for i in IB:
   print i
   clust.append(((0.027*3*i**2)/(1+0.027*3*i**2))*10)
##   curv_IR.append(((0.027*3*i**2)/(1+0.027*3*i**2))*np.exp(0.055*c-0.5*c*c))
for j in range(1,len(clust)):
   curv.append(1.0/np.sqrt(85/(2*(0.05-0.08*np.log(10.0*clust[j]/100.0)))))
   Rad.append(np.sqrt(85/(2*(0.05-0.08*np.log(clust[j]/100.0)))))
for k in np.sort(curv):
   curv_IBAR.append(np.exp(40*54*(0.055*0.055-0.5*0.055**2))-np.exp(40*54*(0.055*k-0.5*k**2)))
   curv_clust.append(np.exp(40*54*(0.055*k-0.5*k**2)))
for l in curvature:
   respto_curv.append(np.exp(40*54*(0.055*l-0.5*l**2))-1.0)
   ##respto_curv.append(np.exp(40*54*(0.055*l-0.5*l**2)))
   print l, np.exp(40*54*(0.055*l-0.5*l**2))
dos=np.sort(respto_curv)
for m in range(1,len(dos)):
   respto_clust.append(1.0/np.sqrt(((40+0.01*np.exp(max(dos)-dos[m]))*4)/(2*(0.05-0.08*np.log((dos[m])*0.0068)))))
   print 'rigidity, max of dos, dos',(15+np.exp(max(dos)-dos[m])), max(dos), dos[m]
plt.xlabel('number of curv_IRSp53')
plt.ylabel('Curvature')
plt.plot(dos,respto_clust,'*-',label='curvature in response to curv_IRSp53')
plt.plot(respto_curv,curvature,label='curv_IRSp53 recruited in response to curvature')
plt.legend()
plt.show()

##n=max(curv)
##m=max(clust)
##for k in range(len(clust)):
##   clust[k]=clust[k]/m
##   curv[k]=curv[k]/n


