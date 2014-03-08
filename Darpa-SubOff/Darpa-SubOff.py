# -*- coding: iso-8859-1 -*-

###
### This file is generated automatically by SALOME v7.3.0 with dump python functionality
###

import math
import salome
import collections

## INIT THE SALOME PART
salome.salome_init()
theStudy = salome.myStudy
from salome.geom import geomBuilder
geompy = geomBuilder.New(theStudy)

## STD GLOBAL VECTORS, CENTER AND PLANES
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Plane_X = geompy.MakePlane(O, OX, 1)
Plane_Y = geompy.MakePlane(O, OY, 1)
Plane_Z = geompy.MakePlane(O, OZ, 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Plane_X, 'Plane_X' )
geompy.addToStudy( Plane_Y, 'Plane_Y' )
geompy.addToStudy( Plane_Z, 'Plane_Z' )


## Create counters and dictionaries to hold the objects
vertCounter=0
lineCounter=0
curveCounter=0
revolutionCounter=0
extrusionCounter=0
compoundCounter=0
mirrorCounter=0
shellCounter=0
supressHolesCounter=0
multiRotationCounter=0
v={}
li={}
c={}
r={}
e={}
co={}
fi={}
mi={}
sh={}
sp={}
mr={}

## These two lists are used to hold the coordinates for the sail
## This makes it possible to make the Sail from 1 Face and not 4
zList=[]
xCapList=[]

## Global DARPA parameters
RMAX=0.8333333
XB= 3.333333
XM= 10.645833
XA= 13.979167
XC= 14.291667


def linspace(a, b, n=50):
    '''Function used to make a linspace between two floats\
    as numpy cannot be imported easily into Salome because of\
    some hooks.py import errors'''
    if n < 2:
        return b
    diff = (float(b) - a)/(n - 1)
    return [diff * i + a  for i in range(n)]

def createPoint(x,y,z,showInStudy=False):
    global vertCounter
    v["Vertex_{0}".format(vertCounter)]=geompy.MakeVertex(x, y, z)
    if showInStudy==True:
        geompy.addToStudy(v["Vertex_{0}".format(vertCounter)], 'Vertex_'+str(vertCounter))
    vertCounter+=1

def createLine(obj1,obj2,showInStudy=False):
    global lineCounter
    li["Line_{0}".format(lineCounter)]=geompy.MakeLineTwoPnt(v["Vertex_{0}".format(obj1)], v["Vertex_{0}".format(obj2)])
    if showInStudy==True:
        geompy.addToStudy(v["Line_{0}".format(lineCounter)], 'Line_'+str(lineCounter))
    lineCounter+=1


def createCurve(start,end,showInStudy=False):
    global curveCounter
    vertexList=[]
    od = collections.OrderedDict(sorted(v.items()))
    for x in xrange(start,end):
        vertexList.append(od["Vertex_{0}".format(x)])
    c["Curve_{0}".format(curveCounter)] = geompy.MakeInterpol(vertexList, False, False)
    if showInStudy==True:
        geompy.addToStudy(c["Curve_{0}".format(curveCounter)], 'Curve_'+str(curveCounter))
    curveCounter+=1


def createRevolution(curve,showInStudy=False):
    global revolutionCounter
    r["Revolution_{0}".format(revolutionCounter)] = geompy.MakeRevolution(curve, OX, 360*math.pi/180.0)
    if showInStudy==True:
        geompy.addToStudy(r["Revolution_{0}".format(revolutionCounter)], 'Revolution_'+str(revolutionCounter))
    revolutionCounter+=1


def createExtrusion(curve,showInStudy=False):
    global extrusionCounter
    e["Extrusion_{0}".format(extrusionCounter)] = geompy.MakePrismVecH(curve, OY, -1)
    if showInStudy==True:
        geompy.addToStudy(e["Extrusion_{0}".format(extrusionCounter)], 'Extrusion_'+str(extrusionCounter))
    extrusionCounter+=1


def createFilling(start,end,showInStudy=False):
    global compoundCounter
    curveList=[]
    od = collections.OrderedDict(sorted(c.items()))
    for x in xrange(start,end):
        curveList.append(od["Curve_{0}".format(x)])
    co["Compound_{0}".format(compoundCounter)] = geompy.MakeCompound(curveList)
    fi["Filling_{0}".format(compoundCounter)] = geompy.MakeFilling(co["Compound_{0}".format(compoundCounter)], theMinDeg=3, theMaxDeg=6, theNbIter=1)
    if showInStudy==True:
        geompy.addToStudy(co["Compound_{0}".format(compoundCounter)], 'Compound_'+str(compoundCounter))
        geompy.addToStudy(fi["Filling_{0}".format(compoundCounter)], 'Filling_'+str(compoundCounter))
    compoundCounter+=1

def createMirror(face,showInStudy=False):
    global mirrorCounter
    mi["Mirror_{0}".format(mirrorCounter)] = geompy.MakeMirrorByPlane(face, Plane_Z)
    if showInStudy==True:
        geompy.addToStudy(mi["Mirror_{0}".format(mirrorCounter)], 'Mirror_'+str(mirrorCounter))
    mirrorCounter+=1


def createShell(shellList,name=None,showInStudy=False):
    global shellCounter
    sh["Shell_{0}".format(shellCounter)] = geompy.MakeShell(shellList)
    if showInStudy==True:
        if name!=None:
            geompy.addToStudy(sh["Shell_{0}".format(shellCounter)], name)
        else:
            geompy.addToStudy(sh["Shell_{0}".format(shellCounter)], 'Shell_'+str(shellCounter))
    shellCounter+=1

def createSupressHoles(obj,name=None,showInStudy=False):
    global supressHolesCounter
    sp["SupressHoles_{0}".format(supressHolesCounter)] = geompy.SuppressHoles(obj, [])
    if showInStudy==True:
        if name!=None:
            geompy.addToStudy(sp["SupressHoles_{0}".format(supressHolesCounter)], name)
        else:
            geompy.addToStudy(sp["SupressHoles_{0}".format(supressHolesCounter)], 'SupressHoles_'+str(supressHolesCounter))
    supressHolesCounter+=1

def createMultiRotation(obj,name=None,showInStudy=False):
    global multiRotationCounter
    mr["Multi_Rotation_{0}".format(multiRotationCounter)] = geompy.MultiRotate1DNbTimes(obj, OX, 4)
    if showInStudy==True:
        if name!=None:
            geompy.addToStudy(mr["Multi_Rotation_{0}".format(multiRotationCounter)], name)
        else:
            geompy.addToStudy(mr["Multi_Rotation_{0}".format(multiRotationCounter)], 'Multi_Rotation_'+str(multiRotationCounter))
    multiRotationCounter+=1


def makeBow():
    listStart=vertCounter
    rang = linspace(0,XB)
    for x in rang:
        y = RMAX*math.pow((1.126395101*x*math.pow((0.3*x-1.),4)\
        +0.442874707*math.pow(x,2)*math.pow((0.3*x-1.),3)\
        +1-math.pow((0.3*x-1.),4)*(1.2*x+1.)),(1./2.1))
        createPoint(x,y,0)
    listEnd=vertCounter
    createCurve(listStart,listEnd)
    createRevolution(c["Curve_{0}".format(curveCounter-1)])


def makeMiddle():
    listStart=vertCounter
    rang = linspace(XB,XM)
    for x in rang:
        y = RMAX
        createPoint(x,y,0)
    listEnd=vertCounter
    createCurve(listStart,listEnd)
    createRevolution(c["Curve_{0}".format(curveCounter-1)])
    
def makeAfterBody():
    rh=0.1175
    k0=10
    k1=44.6244
    
    listStart=vertCounter
    rang = linspace(XM,XA)
    for x in rang:
        e=(13.979167-x)/3.333333
        y = RMAX*math.pow((\
        math.pow(rh,2)+rh*k0*math.pow(e,2)\
        +(20-20*math.pow(rh,2)-4*rh*k0-1.0/3.0*k1)*math.pow(e,3)\
        +(-45+45*math.pow(rh,2)+6*rh*k0+k1)*math.pow(e,4)\
        +(36-36*math.pow(rh,2)-4*rh*k0-k1)*math.pow(e,5)\
        +(-10+10*math.pow(rh,2)+rh*k0+(1./3.)*k1)*math.pow(e,6)\
        ),(1./2.))
        createPoint(x,y,0)
    listEnd=vertCounter
    createCurve(listStart,listEnd)
    createRevolution(c["Curve_{0}".format(curveCounter-1)])
    

def makeAfterBodyCap():
    listStart=vertCounter
    rang = linspace(XA,XC)
    for x in rang:
        try:
            yint = 1.-math.pow((3.2*x-44.733333),2)
            y=0.1175*RMAX*math.pow(yint,0.5)
        except:
            y=0
        createPoint(x,y,0)
    listEnd=vertCounter
    createCurve(listStart,listEnd)
    createRevolution(c["Curve_{0}".format(curveCounter-1)])


def makeSailForemostbody():
#    listStart=vertCounter
    Zmax=0.109375
    rang = linspace(3.032986,3.069155,10)
    y = 1.507813
    for x in rang:
        D=3.0720*(x-3.032986)
        C=1-math.pow((D-1),4)*(4*D+1)
        B=1./3.*(D)**2*(D-1)**3
        A=2*D*(D-1)**4
        z = Zmax * math.pow((2.094759*A+0.2071781*B+C),0.5)
        zList.append(z)
        xCapList.append(x)
        createPoint(x,y,z)
#    listEnd=vertCounter
#    createCurve(listStart,listEnd)
    
def makeSailForebody():
#    listStart=vertCounter
    Zmax=0.109375
    rang = linspace(3.069155,3.358507,10)
    y = 1.507813
    for x in rang:
        D=3.0720*(x-3.032986)
        C=1-math.pow((D-1),4)*(4*D+1)
        B=1./3.*(D)**2*(D-1)**3
        A=2*D*(D-1)**4
        z = Zmax * math.pow((2.094759*A+0.2071781*B+C),0.5)
        
        if x==3.069155:
            continue
        createPoint(x,y,z)
        xCapList.append(x)
        zList.append(z)
#    listEnd=vertCounter
#    createCurve(listStart,listEnd)

def makeSailMiddlebody():
#    listStart=vertCounter
    Zmax=0.109375
    rang = linspace(3.358507,3.559028,10)
    y = 1.507813
    for x in rang:
        z = Zmax
        
        
        if x==3.358507:
            continue
        createPoint(x,y,z)
        xCapList.append(x)
        zList.append(z)
#    listEnd=vertCounter
#    createCurve(listStart,listEnd)

def makeSailAfterbody():
#    listStart=vertCounter
    Zmax=0.109375
    rang = linspace(3.559028,4.241319,10)
    y = 1.507813
    for x in rang:
        E=(4.241319-x)/0.6822917
        F=E-1
        G=2.238361*E*F**4
        H=3.106529*(E**2)*(F**3)
        P=1-(F**4)*(4*E+1)
        z = Zmax * (G+H+P)
        if x==3.559028:
            continue
        createPoint(x,y,z)
        xCapList.append(x)
        zList.append(z)

def makeSailCap():
    xlist=[]
    ylist=[]
    zlist=[]
    curveListStart = curveCounter
    
    for i,x in enumerate(xCapList):
        z1=zList[i]
        xlist.append(x)
        yrang = linspace(1.507813,(z1/2)+1.507813,10)
        for y in yrang:
            value = z1**2-(2*(y-1.507813))**2
            try:
                z=math.sqrt(value)
            except:
                z=0.
            
            ylist.append(y)
            zlist.append(z)
    
    for i in xrange(0,10):
        listStart=vertCounter
        counter=i

        for x in xlist:
            y=ylist[counter]
            z=zlist[counter]
            counter+=10
            createPoint(x,y,z)
            
        listEnd=vertCounter
        createCurve(listStart,listEnd)
    curveListEnd = curveCounter
    createFilling(curveListStart,curveListEnd)
    createMirror(fi["Filling_{0}".format(compoundCounter-1)])


def makeHull():
    makeBow()
    makeMiddle()
    makeAfterBody()
    makeAfterBodyCap()
    
    createShell([r["Revolution_{0}".format(revolutionCounter-1)],\
    r["Revolution_{0}".format(revolutionCounter-2)],\
    r["Revolution_{0}".format(revolutionCounter-3)],\
    r["Revolution_{0}".format(revolutionCounter-4)]],showInStudy=True,name="Hull")

def makeSail():
    listStart=vertCounter
    makeSailForemostbody()
    makeSailForebody()
    makeSailMiddlebody()
    makeSailAfterbody()
    listEnd=vertCounter
    createCurve(listStart,listEnd)
    createExtrusion(c["Curve_{0}".format(curveCounter-1)])
    createMirror(e["Extrusion_{0}".format(extrusionCounter-1)])
    makeSailCap()
    
    createShell([e["Extrusion_{0}".format(extrusionCounter-1)],\
    fi["Filling_{0}".format(compoundCounter-1)],\
    mi["Mirror_{0}".format(mirrorCounter-1)],\
    mi["Mirror_{0}".format(mirrorCounter-2)]])
    createSupressHoles(sh["Shell_{0}".format(shellCounter-1)],showInStudy=True,name="Sail")


def makeSternAppendage():
    H=[12.729618,13.146284,13.552950]
    
    DELR=0.025
    RR=0.075
    
    XXI=[0.0,0.005,0.0125,0.025,0.05\
    ,0.075,0.1,0.15,0.2,0.25,0.3,0.4\
    ,0.5,0.6,0.7,0.8,0.9,0.95,1.]
    
    for nu,h in enumerate(H):
        curveListStart = curveCounter
        for i in xrange(1,31):
            CY=-0.466308*RR+0.88859
            listStart=vertCounter
            for j in XXI:
                XXX = (j-1)*CY+h
                a=0.2969*math.sqrt(j)
                b=0.126*j
                c=0.3516*j*j
                d=0.2852*j**3
                e=0.1045*j**4
                value = a - b - c + d - e 
                Z=value*CY
                createPoint(XXX,RR,Z)
            listEnd=vertCounter
            createCurve(listStart,listEnd)
            RR=RR+DELR
        curveListEnd = curveCounter
        createFilling(curveListStart,curveListEnd)
        createMirror(fi["Filling_{0}".format(compoundCounter-1)])
        
        createShell([fi["Filling_{0}".format(compoundCounter-1)],\
        mi["Mirror_{0}".format(mirrorCounter-1)]])
        
        createSupressHoles(sh["Shell_{0}".format(shellCounter-1)])
        
        createMultiRotation(sp["SupressHoles_{0}".format(supressHolesCounter-1)]\
        ,showInStudy=True,name="Appendage_"+str(nu+1))
        RR=0.075
    
def makeRingWing():
    
    YC=[0]
    YCP=[0]
    YT=[0]
    XU=[0,0]
    YU=[0,0]
    XL=[0,0]
    YL=[0,0]
    XDLE=[13.46990,13.46990]
    YDLE=[0.43004,0.47681]
    XDTE=[14.21661,14.2074]
    YDTE=[0.35659,0.33856]
    
    XC=[0.,0.0,0.005,0.0075,0.0125,0.025\
    ,0.05,0.075,0.1,0.15,0.2,0.25,0.3\
    ,0.35,0.4,0.45,0.5,0.55,0.6,0.65\
    ,0.7,0.75,0.8,0.85,0.9,0.95,1.]
    
    B=[0,0.43756,-0.08136,-0.06496,-0.01926,-0.00185\
    ,0.00348,0.00156,-0.00113,-0.00058,0.00027\
    ,0.00080,0.00006,-0.00027,-0.00033,0.00005\
    ,0.00014,0.00008]
    
    for i in xrange(1,27):
        X=XC[i]
        D=0.4-X
        E=1.0-X
        if abs(X-0.0)<1.0e-20:
            X=1.0e-30
        if abs(D)<1.0e-20:
            D=1.0e-30
        if abs(E)<1.0e-20:
            E=1.0e-30
        ycVal1=-0.049921*(0.5*D*D*math.log(abs(D))-0.5*E*E*math.log(E)\
        +0.25*E*E-0.25*D*D)
        ycVal2=ycVal1+0.029953*(X*math.log(X)+0.227828-0.531076*X)
        YC.append(ycVal2)
        
        ycpVal1=-0.049921*(E*math.log(E)-D*math.log(abs(D)))\
        +0.02995253*(math.log(X)+0.4689244)
        YCP.append(ycpVal1)
    
    for i in xrange(1,27):
        X=XC[i]
        
        if i >= 16:
            XC1=1.0-X
            yVal=0.03333+1.696969*XC1-1.441945*XC1**2\
            -0.366363*XC1**3+0.333049*XC1**4
#            print yVal
#            print yVal*0.1
            YT.append(yVal*0.1)
        else:
            OM=math.acos(2.0*X-1.)
            YY=0.
            for j in xrange(1,18):
#                print j
                YY=YY+(B[j]*math.sin(j*OM))
#            print YY
#            print YY*0.1
            YT.append(YY*0.1)
    
    for i in xrange(2,27):
        TH=math.atan(YCP[i])
#        print "TH =",math.degrees(TH)
        SINTH=math.sin(TH)
        COSTH=math.cos(TH)
        
        XU.append(XC[i]-YT[i]*SINTH)
        YU.append(YC[i]+YT[i]*COSTH)
        
        XL.append(XC[i]+YT[i]*SINTH)
        YL.append(YC[i]-YT[i]*COSTH)
    
    
    for k in xrange(0,2):
        PHI=math.atan2((YDTE[k]-YDLE[k]),(XDTE[k]-XDLE[k]))
        CS=math.cos(PHI)
        SN=math.sin(PHI)
        CHORD=math.sqrt((YDTE[k]-YDLE[k])**2+(XDTE[k]-XDLE[k])**2)
        listStart=vertCounter
        for i in xrange(1,27):
            x=XDLE[k]+(CHORD*(XU[i]*CS-YU[i]*SN))
            y=YDLE[k]+(CHORD*(XU[i]*SN+YU[i]*CS))
            createPoint(x,y,0)
        lp1=vertCounter-1
        listEnd=vertCounter
        createCurve(listStart,listEnd)
        createRevolution(c["Curve_{0}".format(curveCounter-1)])
        
        listStart=vertCounter
        for i in xrange(1,27):
#            XLL=XL[i]
            x=XDLE[k]+(CHORD*(XL[i]*CS-YL[i]*SN))
            y=YDLE[k]+(CHORD*(XL[i]*SN+YL[i]*CS))
            createPoint(x,y,0)
        listEnd=vertCounter
        lp2=vertCounter-1
        createCurve(listStart,listEnd)
        createRevolution(c["Curve_{0}".format(curveCounter-1)])
        
        createLine(lp1,lp2)
        createRevolution(li["Line_{0}".format(lineCounter-1)])
        
        createShell([r["Revolution_{0}".format(revolutionCounter-1)],\
        r["Revolution_{0}".format(revolutionCounter-2)],\
        r["Revolution_{0}".format(revolutionCounter-3)]]\
        ,showInStudy=True,name="Ringwing_"+str(k+1))


def makeRingWingStrut():
    
    XHLE=13.589
    XHTE=13.83582
    YHLE=0.14726
    YHTE=0.10547
    
    XDLE=[13.63845,13.64487]
    YDLE=[0.36886,0.39755]
    XDTE=[13.88818,13.89023]
    YDTE=[0.34002,0.34932]
    
    
    XC=[0.,0.0,0.005,0.0125,0.025,0.05\
    ,0.075,0.1,0.15,0.2,0.25,0.3\
    ,0.4,0.5,0.6,0.7,0.8,0.9,1.]

    for k in xrange(0,2):
        R1=YHLE
        R2=YDLE[k]
        DELR=(R2-R1)/9.
        R=R2+DELR*2.
        curveListStart = curveCounter
        for l in xrange(1,14):
            X0=0.223221*R+13.556128
            listStart=vertCounter
            for i in xrange(1,len(XC)):
                XI = XC[i]
                x=(X0+0.243995*XI)
                y=(R-0.054465*XI)
                z=0.15*(0.29690*math.sqrt(XI)-0.126*XI-0.3516*XI**2\
                + 0.28520*XI**3 - 0.10450*XI**4)
                createPoint(x,y,z)
            listEnd=vertCounter
            createCurve(listStart,listEnd)
            R=R-DELR
        curveListEnd = curveCounter
        createFilling(curveListStart,curveListEnd)
        createMirror(fi["Filling_{0}".format(compoundCounter-1)])
        
        createShell([fi["Filling_{0}".format(compoundCounter-1)],\
        mi["Mirror_{0}".format(mirrorCounter-1)]])
        
        createSupressHoles(sh["Shell_{0}".format(shellCounter-1)])
        
        createMultiRotation(sp["SupressHoles_{0}".format(supressHolesCounter-1)]\
        ,showInStudy=True,name="Ringwing_Strut_"+str(k+1))
        
        geompy.Rotate(mr["Multi_Rotation_{0}".format(multiRotationCounter-1)], OX, 45*math.pi/180.0)
        

def main():
    ## START OF PROGRAM

    makeHull() ## Make the HULL
    makeSail() ## Make the SAIL
    makeSternAppendage() ## Make the 3 Stern Appendages
    makeRingWing() ## Make the 2 Ring Wings
    makeRingWingStrut() ## Make the 2 Struts for the Ring Wings

if __name__ == '__main__':
    main()