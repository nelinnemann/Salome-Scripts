# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v7.5.1 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
geompy = geomBuilder.New(theStudy)
import math
import SALOMEDS

points = 50
radius = 0.018*1000
pipeRadius = 5
turnLength = (0.877/10.)*1000
numberOfturns = 10

def Helix(a, t):
	return [a*math.cos(t), a*math.sin(t), turnLength*t/(2*math.pi)]

def dHelix(a, t):
	return [a*math.sin(t), a*math.cos(t), turnLength/(2*math.pi)]

# create vertices
pointList=[]
for i in range(numberOfturns*points+1):
	Helixx = Helix(radius, 2*math.pi/points * i)
	pointList.append(geompy.MakeVertex(Helixx.pop(0), Helixx.pop(0), Helixx.pop(0)))

HelixCurve = geompy.MakePolyline(pointList,False)

Helixx = tuple(Helix(radius, 0))
dHelixx = dHelix(radius, 0)

P1forVector = geompy.MakeVertex(Helixx[0], Helixx[1], Helixx[2])
P2forVector = geompy.MakeVertex(Helixx[0]+dHelixx.pop(0), Helixx[1]+dHelixx.pop(0), Helixx[2]+dHelixx.pop(0))

VectorForBaseCircle = geompy.MakeVector(P1forVector, P2forVector)

Divided_Disk_1 = geompy.MakeDividedDiskPntVecR(P1forVector, VectorForBaseCircle, pipeRadius, GEOM.HEXAGON)

# add objects in the study
geompy.addToStudy(HelixCurve, "HelixCurve")
geompy.addToStudy(Divided_Disk_1, 'Divided_Disk_1')

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Divided_Disk_1)
Regular_1D = Mesh_1.Segment()
Nb_Segments_1 = Regular_1D.NumberOfSegments(5)
Nb_Segments_1.SetDistrType( 0 )
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_1.Compute()
Nb_Segments_2 = smesh.CreateHypothesis('NumberOfSegments')
Nb_Segments_2.SetNumberOfSegments( 2 )
Nb_Segments_2.SetDistrType( 0 )
Mesh_2 = smesh.Mesh(HelixCurve)
status = Mesh_2.AddHypothesis(Nb_Segments_2)
status = Mesh_2.AddHypothesis(Regular_1D)
isDone = Mesh_2.Compute()
error = Mesh_1.ExtrusionAlongPathX( Mesh_1, Mesh_2, 1, 0, [  ], 0, 0, [ 0, 0, 0 ], 0, SMESH.FACE )


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Nb_Segments_2, 'Nb. Segments_2')
smesh.SetName(Nb_Segments_1, 'Nb. Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
