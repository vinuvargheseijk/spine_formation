#################################################################
## This program is part of 'MOOSE', the
## Messaging Object Oriented Simulation Environment.
##           Copyright (C) 2015 Upinder S. Bhalla. and NCBS
## It is made available under the terms of the
## GNU Lesser General Public License version 2.1
## See the file COPYING.LIB for the full notice.
##
## rxdSpineSize.py: Builds a cell with spines and a propagating reaction
## wave. Products diffuse into the spine and cause it to get bigger.
##################################################################
## Cdc42 activation profile (spatial and time) has to be accounted, 
## Tension modification
import math
import pylab
import numpy as np
import matplotlib.pyplot as plt
import moose
import sys
sys.path.append( '../util' )
import rdesigneur as rd
import matplotlib
import moose.SBML

PI = 3.141592653
ScalingForTesting = 10
RM = 1.0 / ScalingForTesting
RA = 1.0 * ScalingForTesting
CM = 0.01 * ScalingForTesting
runtime = 250
frameruntime = 1.0
diffConst = 5e-10
dendLen = 10e-6
diffLen = 2.0e-7
dendDia = 1e-6
somaDia = 1e-6
concInit = 0.1    # 10 micromolar
spineSpacing = 7e-6
spineSpacingDistrib = 1e-6
spineSize = 1.0
spineSizeDistrib = 0.5
spineAngle= np.pi / 2.0
spineAngleDistrib = 0.0

def makeCellProto( name ):
    elec = moose.Neuron( '/library/' + name )
    ecompt = []
    ###soma = rd.buildCompt( elec, 'soma', dx=somaDia, dia=somaDia, x=-somaDia, RM=RM, RA=RA, CM=CM )
    dend = rd.buildCompt( elec, 'dend', dx=dendLen, dia=dendDia, x=0, RM=RM, RA=RA, CM=CM )
    ###moose.connect( soma, 'axial', dend, 'raxial' )
    elec.buildSegmentTree()

def makeChemProto(name='hydra'):
        moose.Neutral('/library/')
        meshName='/library/'+name
        moose.SBML.mooseReadSBML('./rec_mech2_5_4.xml',meshName)
        Rad=moose.element('/library/hydra/kinetics/Rad')
        I_C=moose.element(meshName+'/kinetics/IRSp53')
        Cdc42=moose.element(meshName+'/kinetics/Cdc42')
        Eps8=moose.element(meshName+'/kinetics/Eps8')
        I_Cact=moose.element(meshName+'/kinetics/IRSp53_a')
        I_M=moose.element(meshName+'/kinetics/IRSp53_m')	
        Idimer=moose.element(meshName+'/kinetics/IRSp53_dimer')
        PIP2=moose.element(meshName+'/kinetics/PIP2')
        cPIP2=moose.element(meshName+'/kinetics/cPIP2')
        cip2=moose.element(meshName+'/kinetics/cip2')
        cip3=moose.element(meshName+'/kinetics/cip3')
        height=moose.element(meshName+'/kinetics/height')
        recBAR=moose.Function(meshName+'/kinetics/recBAR')
        gradBAR=moose.Function(meshName+'/kinetics/gradBAR')
        ###Radfun=moose.element(meshName+'/kinetics/Rad/func')
        Radfun=moose.Function(meshName+'/kinetics/Radfun')
        Ipool=moose.element(meshName+'/kinetics/Ipool')
        sort=moose.element(meshName+'/kinetics/sort')
        curv_IRSp53=moose.element(meshName+'/kinetics/curv_IRSp53')
        curv=moose.element(meshName+'/kinetics/curv')
        curv_pip=moose.element(meshName+'/kinetics/curv_pip')
        simpleFun=moose.Function(meshName+'/kinetics/simpleFun')
        recFun=moose.Function(meshName+'/kinetics/recFun')
        enzn=moose.element(meshName+'/kinetics/curv/Enz')
        enzn1=moose.element(meshName+'/kinetics/curv_pip/Enz')
        enzncplx=moose.element(meshName+'/kinetics/curv/Enz/Enz2_cplx')
        ##tube_reac=moose.element(meshName+'/kinetics/tube_reac')
        enzn2=moose.element(meshName+'/kinetics/sort/Enz')
        enzn3=moose.element(meshName+'/kinetics/detach/Enz')
        gradBAR.x.num=4
        Radfun.x.num=4
        recBAR.x.num=2
        simpleFun.x.num=1
        recFun.x.num=1
        #curv.x.num=1
        #curv.expr="0*x0+x0^-1/(0.2*1e-6)"
        #moose.connect(Rad,'nOut',curv.x[0],'input')
        gradBAR.expr="0*x0+2.1*10^-28*2*((x1+x2+x3)-x0)*(1/(20.0*20.0*1e-18))*(x1+x2+x3)"##+(1/(x2))*2*10^-14*(x2>0)"
        ###simpleFun.expr="(1.0/10^10)*sqrt(x0^2)"
        simpleFun.expr="sqrt((55*4)*(1/(x0+1)))"

### I assume that deformation happens only after a threshold in number of IBAR is reached. Therefore a condition is given in Rad calculation

##        Radfun.expr="sqrt((55*4*1e-18)/(2*(0.1+0.08*(ln(1/(1*(((3.0*x0+1.0*x1)/(2*3.14*(0.05*0.45)))*1e-6*54*100)))*(x0>1.0)))))" ##+  sqrt((55*4*1e-18)/(2*0.03))*(x0<0.1)"
        recFun.expr="0.0*x0+(exp(70*54*((0.055/(x0*1e9))-(1/(2.0*(x0*1e9)*(x0*1e9))))))"
        Radfun.expr="sqrt((75*4*1e-18)/(2*(0.04+0.08*(ln(1/(1*(((1.0*x0+1.0*x1+1.0*x3)/(4*3.14*(0.05*0.05)))*1e-6*54)))*(((1.0*x0+1.0*x1+1.0*x2+x3)*54*1e-6*100)/(4.0*3.14*(0.05*0.05))>2.0)))))" ##+  sqrt((55*4*1e-18)/(2*0.03))*(x0<0.1)"
        ###recBAR.expr="0.0004+0.0*x0+30.0*(exp((1/(x0))*0.02*1e-6)-1)"
        ##recBAR.expr="0.0763249212343+exp(-(20*1e-9)/x0)"
        ###recBAR.expr="0.0004+0.0*x0+30.0*((1/(x0+1))*2-1)*(x0>0)"
        #gradBAR_expr=gradBAR_expr.lstrip("0 +")
        #gradBAR_expr.replace(" ","")
        #gradBAR.expr = gradBAR_expr
        #gradBAR.expr="x0+x1+0.0*x2"
        #Radfun.expr="sqrt((55*4*10^-18)/(2*(0.002+0.08*ln(1/0.15))))"
        #recBAR.expr="x0/2"
        #recBAR.expr="0*x0+190*x1*x2-4*x3"
        moose.connect(I_M,'nOut',gradBAR.x[0],'input')
        moose.connect(curv_IRSp53,'nOut',gradBAR.x[1],'input')
        moose.connect(cip2,'nOut',gradBAR.x[2],'input')
        moose.connect(cip3,'nOut',gradBAR.x[3],'input')
        moose.connect(gradBAR,'valueOut',height,'increment')

        moose.connect(enzn1,'sub',PIP2,'reac') 
        moose.connect(enzn1,'sub',PIP2,'reac') 
        moose.connect(enzn1,'sub',PIP2,'reac') 
        moose.connect(enzn1,'sub',PIP2,'reac') 
        moose.connect(enzn1,'sub',PIP2,'reac') 

        ##moose.connect(enzn,'sub',I_M,'reac')
        ##moose.connect(enzn,'sub',I_M,'reac')
        

        ##moose.connect(enzn3,'prd',I_M,'reac')
        ##moose.connect(enzn3,'prd',I_M,'reac')

        moose.connect(curv_IRSp53,'nOut',Radfun.x[0],'input')
        moose.connect(I_M,'nOut',Radfun.x[1],'input')
        moose.connect(cip2,'nOut',Radfun.x[2],'input')
        moose.connect(enzncplx,'nOut',Radfun.x[3],'input')
        moose.connect(Radfun,'valueOut',Rad,'setN')
        print gradBAR.expr
        print Radfun.expr
        
     
        moose.connect(Rad,'nOut',recFun.x[0],'input')
        moose.connect(recFun,'valueOut',Ipool,'setN')

### Ipool is further copied to sort to participate in the reaction tube_reac or enz (model version 2.5.3)

        I_C.concInit=0.0012
        

        ## Diffconsts are set here

        I_C.diffConst=10e-12
        I_Cact.diffConst=10e-12
        Eps8.diffConst=10e-12
        Cdc42.diffConst=10e-12
        I_M.diffConst=1e-12
        ##I_M.diffConst=1e-12
        Idimer.diffConst=10e-12
        ##I_Mdim.diffConst=1e-12
        PIP2.diffConst=0.8e-12
        cPIP2.diffConst=1e-14
        ##IP_clx.diffConst=1e-12
        height.diffConst=2e-14
        curv_IRSp53.diffConst=1e-14
        cip2.diffConst=0.0
        #moose.connect(curv,'valueOut',height,'increment')
        #moose.element(meshName+'/Compartment/Rad').nInit=1e-9
        #moose.element(meshName+'/Compartment/height').concInit=0.001
        ###print moose.le(meshName+'/Compartment')
        return meshName

def makeModel():
    moose.Neutral( '/library' )
    makeCellProto( 'cellProto' )
    makeChemProto()
    #moose.SBML.readSBML('/home/vinu/Documents/Biophysics/new_mod/detailed_models/oct01_volScale/mv1.xml','/model')
    rdes = rd.rdesigneur( 
            useGssa = False,
            turnOffElec = True,
            chemDt = 0.01,
            diffDt = 0.01,
            chemPlotDt = 0.1,
            combineSegments = False,
            stealCellFromLibrary = True,
            diffusionLength = 50e-9,
            cellProto = [['cellProto', 'elec' ]] ,
            chemProto = [['hydra','chem']] ,
            chemDistrib = [[ "chem", "#", "install", "1" ]]
            #plotList = [
                #['psd/z', '1', '.', 'n', 'z n' ]
                #['psd/z', '1', '.', 'n', 'z n' ]
            #],
        )
    moose.seed(100)
    rdes.buildModel( '/model' )
    moose.element('/model/chem/dend/Rad').nInit=30e-9
   ## moose.element('/model/chem/dend/stoich').allowNegative=1
    print 'built model'

def makePlot( name, srcVec, field ):
    tab = moose.Table2('/graphs/' + name + 'Tab', len( srcVec ) ).vec
    print 'source vec is', srcVec
    for i in zip(srcVec, tab):
        moose.connect(i[1], 'requestOut', i[0], field)
    return tab


def displayPlots():
    fig = plt.figure( figsize = ( 20, 20 ) )
    #plt.suptitle( 'Demo: Spine diameter changes as per z conc, y just diffuses. Pool concs are reported unless otherwise indicated.', fontsize = 18)
    wild = moose.wildcardFind( '/graphs/#[0]' )
    j = 1
    for x in wild:
        plt.subplot( len(wild)/3, 3, j )
        j += 1
        tab = moose.vec( x )
        for i in range( len( tab ) ):
            v = np.array( tab[i].vector )
            plt.plot( v )
            #plt.plot( v, label=x.name + " " + str( i ) )
        #plt.title( x.name )
        plt.ylabel( x.name[:-3], fontsize = 16 )
        #pylab.legend()
        #pylab.figure()

def radCalc():
    model=moose.Neutral('/model')
    hello=moose.PyRun('/model/hello')
    hello.initString="""print 'initialising'
count=0
print 'initial count',count"""
    hello.runString="""print 'running'
count+=1
print 'count',count"""
    moose.useClock(0,hello.path,'process')
    moose.reinit()
    moose.start(0.001)

def main():
    """
    This illustrates the use of rdesigneur to build a simple dendrite with
    spines, and then to resize them using spine fields. These are the
    fields that would be changed dynamically in a simulation with reactions
    that affect spine geometry.
    In this simulation there is a propagating reaction wave using a
    highly abstracted equation, whose product diffuses into the spines and
    makes them bigger.
    """
    radCalc()
    moose.delete('/model')
    makeModel()
    for i in range( 11, 18 ): 
        moose.setClock( i, 0.01 )
    moose.setClock( 18, 0.01 )
    elec = moose.element( '/model/elec' )
    graphs = moose.Neutral( '/graphs' )
    Epsinit=moose.vec('/model/chem/dend/Eps8')
    Cdcinit=moose.vec('/model/chem/dend/Cdc42')
    IRSpinit=moose.vec('/model/chem/dend/IRSp53')
    PIPinit=moose.vec('/model/chem/dend/PIP2')
    Cdcinit.concInit=0.0
##    Cdcinit[100].concInit=0.01
    ##Epsinit.concInit=0.0
    ##Epsinit[50].concInit=0.010

###    IRSpinit.nInit=0.0
    for ci in range(90,110):
         Cdcinit[ci].concInit=np.exp(-((ci-100)**2)/10.**2)*0.010
##        Cdcinit[ci].concInit=0.010
    PIPinit.nInit=565

    print moose.vec('/model/chem/dend/Cdc42').nInit
    print moose.vec('/model/chem/dend/IRSp53').nInit
    makePlot( 'IRSp53', moose.vec( '/model/chem/dend/IRSp53' ), 'getN' )
    makePlot( 'height', moose.vec( '/model/chem/dend/height' ), 'getN' )
    makePlot( 'IM', moose.vec( '/model/chem/dend/IRSp53_m' ), 'getN' )
    makePlot( 'Rad', moose.vec( '/model/chem/dend/Rad' ), 'getN' )
    makePlot( 'curv_IRSp53', moose.vec( '/model/chem/dend/curv_IRSp53' ), 'getN' )
    makePlot( 'Ipool', moose.vec( '/model/chem/dend/Ipool' ), 'getN' )
    makePlot( 'sort', moose.vec( '/model/chem/dend/sort' ), 'getN' )
    rad=moose.element('/model/chem/dend/Rad')
   ##curv[gri]=(radv[gri+1]-2*radv[gri]+radv[gri-1])/(200*200*1e-9)
##   curv[gric]=1e-9*(heightvec[gric+1]-2*heightvec[gric]+heightvec[gric-1])/(0.2*0.2*1e-12)
    mypyrun=moose.PyRun('/model/mypyrun')
##    curv[gri]=np.exp(6*curv_IRSp53vec[gri]+6*cip2vec[gri]+6*cip3vec[gri]+2*IRSp53_mvec[gri])
##cPIP2.diffConst=abs(0.33e-12-(1.0/np.exp(abs(max(radv)*1e9-20.)))*0.33e-12)
##curv_IRSp53.diffConst=cPIP2.diffConst
   ##curv[gric]=1e-9*(-heightvec[gric+2]+16*heightvec[gric+1]-30*heightvec[gric]+16*heightvec[gric-1]-heightvec[gric-2])/(12*20.0*20.0*1e-18)
   ###curv[gric]=1e-9*(-radv[gric+2]+16*radv[gric+1]-30*radv[gric]+16*radv[gric-1]-radv[gric-2])/(12*50.0*50.0*1e-18)
   ##curv[gric]=1e-9*(radv[gric+1]-2*radv[gric]+radv[gric-1])/(0.05*0.05*1e-12)
##for gric in range(0,len(curv)):
##   curv[gric]=1.0/(100.0-2*1e9*radv[gric])
##   curv[gric]=curv[gric]

##   curv[gri].conc=(3*27000*IRSp53_mvec[gri].conc**2)/(1+3*27000*IRSp53_mvec[gri].conc**2)
##   curv[gri]=1030*((3*0.027*(3*curv_IRSp53vec[gri]+3*cip2vec[gri]+3*cip3vec[gri]+IRSp53_mvec[gri])**2)/(1.0+3.0*0.027*(3.0*curv_IRSp53vec[gri]+3.0*cip2vec[gri]+3.0*cip3vec[gri]+IRSp53_mvec[gri])**2))
    print moose.showfield('/model/chem/dend/curv/Enz'), 'number of substrates'
    mypyrun.initString="""count=0"""
    mypyrun.runString="""count=count+1
radv=moose.vec('/model/chem/dend/Rad').n
radm=radv[10]
output=radm
IRvec=moose.vec( '/model/chem/dend/IRSp53' ).n
curv_IRSp53=moose.element( '/model/chem/dend/curv_IRSp53')
Cdcvec=moose.vec( '/model/chem/dend/Cdc42' ).n
heightvec=moose.vec( '/model/chem/dend/height' ).n
curv_IRSp53vec=moose.vec('/model/chem/dend/curv_IRSp53').n
IRSp53_mvec=moose.vec('/model/chem/dend/IRSp53_m').n
cip2vec=moose.vec('/model/chem/dend/cip2').n
cip3vec=moose.vec('/model/chem/dend/cip3').n
Ipoolvec=moose.vec('/model/chem/dend/Ipool').n
curv_ir_vec=moose.vec('/model/chem/dend/curv/Enz/Enz2_cplx').n
curv=moose.vec('/model/chem/dend/curv').n
cPIP2=moose.element('/model/chem/dend/cPIP2')
##moose.vec('/model/chem/dend/sort').n=(-Ipoolvec+max(Ipoolvec))
moose.vec('/model/chem/dend/sort').n=(Ipoolvec-min(Ipoolvec))
moose.vec('/model/chem/dend/detach').n=-moose.vec('/model/chem/dend/sort').n+max(moose.vec('/model/chem/dend/sort').n)
##if max(Ipoolvec)>min(Ipoolvec):
##    plt.plot(Ipoolvec)
##    plt.plot(moose.vec('/model/chem/dend/sort').n)
##    plt.figure(2)
##    plt.plot(moose.vec('/model/chem/dend/Rad').n)

for gri in range(0,len(Cdcvec)):
   curv[gri]=1.0*((3*0.027*(curv_IRSp53vec[gri]+IRSp53_mvec[gri])**2)/(1.0+3.0*0.027*(curv_IRSp53vec[gri]+curv_ir_vec[gri]+IRSp53_mvec[gri])**2))*0.055*np.exp(-(moose.vec('/model/chem/dend/mesh')[gri].Coordinates[0]-5e-6)**2/(2*0.1e-6**2))
if count==1:
   maxCIR=[]
   maxIR=[]
   max_curv=[]
moose.vec('/model/chem/dend/curv').n=curv
moose.vec('/model/chem/dend/curv_pip').n=curv
maxCIR.append(3*max(moose.vec('/model/chem/dend/curv_IRSp53').n))
maxIR.append(max(moose.vec('/model/chem/dend/IRSp53_m').n))
max_curv.append(max(moose.vec('/model/chem/dend/curv').n))
pref_curv=0.055
PIP2vec=moose.vec('/model/chem/dend/PIP2').n
plt.xlabel('Number of IRSp53 dimers')
plt.ylabel('Deformation strength')
if count%500==0:
   print 'Dimer n', moose.vec('/model/chem/dend/IRSp53_dimer').n
   print 'PIP2', moose.vec('/model/chem/dend/PIP2').n
   print 'Radius', radv
   print 'Curvature', moose.vec('/model/chem/dend/curv').n
   print 'CurvIRSp53', curv_IRSp53vec
   plt.plot(maxCIR,max_curv)
radmax=max(radv)
Rn=moose.element( '/model/chem/dend/Rad' ).n
"""
    '''
    debug = moose.PyRun( '/pyrun' )
    debug.tick = 10
    debug.runString = """print "RUNNING: ", moose.element( '/model/chem/psd/z' ).n, moose.element( '/model/elec/head0' ).diameter"""
    '''
    ###moose.connect(rad,'nOut',mypyrun,'trigger')
    outputTab=moose.Table('/model/output')
    moose.connect(mypyrun,'output',outputTab,'input')
    moose.setClock(30, 0.01)
    moose.setClock(10, 0.01)
    moose.setClock(15, 0.01)
    moose.setClock(16, 0.01)
    
    mypyrun.mode=1
    ###radCalc()
    moose.reinit()
    ###moose.start( runtime )
    test = 'False'
    step=0


    if test=='True':
       while step<int(runtime/0.001):
          I_M = moose.element('/model/chem/dend/IRSp53_m').n
          I_C = moose.element('/model/chem/dend/IRSp53').n
          Rad = moose.element('/model/chem/dend/Rad').n
          height = moose.element('/model/chem/dend/height').n
          gradBAR = moose.element('/model/chem/dend/gradBAR').expr
          Radfun = moose.element('/model/chem/dend/Radfun').expr
          recBAR = moose.element('/model/chem/dend/recBAR').expr
          gradBAReval = 2.1*1e-28*(I_M-I_C)*(1/(0.2*0.2*1e-12))*I_C+0*Rad
          print 'gradBAR', gradBAReval
          print I_C,I_M,I_P,Rad
          print moose.element('/model/chem/dend/IPreac').numKf
          if I_P > 0:
            Radfuneval=np.sqrt((55*4*1e-18)/(2*(0.02-0.08*(np.log(1/(I_P*50*100+1))))))
            print np.log(1/(I_P*50*100+1))
            print (2*(0.02-0.08*(np.log(1/(I_P*50*100+1)))))
            print (55*4*1e-18)/(2*(0.02-0.08*(np.log(1/(I_P*50*100+1)))))
            print 'Radfun', Radfuneval
          raw_input()
          moose.start(0.001)
          step=step+1
          print step
    else:
       moose.start(0.001) 
    print 'IMTab', moose.element('/graphs/IMTab').vector
    print 'IRSp53Tab', moose.element('/graphs/IRSp53Tab').vector
    print 'Radius', moose.element('/graphs/RadTab').vector
    print 'height', moose.element('/graphs/heightTab').vector
    print 'curv_IRSp53', moose.element('/graphs/curv_IRSp53Tab').vector
    print 'Ipool', moose.element('/graphs/IpoolTab').vector
    print 'sort', moose.element('/graphs/sortTab').vector
    ##print 'recruited IP', moose.vec('/model/chem/dend/IP_clx').n
    ##plt.plot(moose.vec('/model/chem/dend/IP_clx').n,label='IP at time 150')
    ##moose.start(250)
    ##if height=='True':
    ##    plt.plot(moose.vec('/model/chem/dend/height').n,label='At time 300')
    ##else:
    ##    plt.plot(moose.vec('/model/chem/dend/Rad').n,label='At time 300')
    plt.legend()
    plt.show()
       
    #displayPlots()
    #plt.legend()

# Run the 'main' if this script is executed standalone.

if __name__ == '__main__':
    main()
