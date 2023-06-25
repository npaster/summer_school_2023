#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.4.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

pt_r = 200
gd_r = 1000

# pt_r = 100
# gd_r = 200

nb_seg_ep = 4
nb_seg_arc = 4

geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Sphere_1 = geompy.MakeSphereR(gd_r)
Sphere_2 = geompy.MakeSphereR(pt_r)
Sph_tot = geompy.MakeCutList(Sphere_1, [Sphere_2], True)
Face_1 = geompy.MakeFaceObjHW(OY, gd_r*2, gd_r*2)
Face_2 = geompy.MakeFaceObjHW(OX, gd_r*2, gd_r*2)
Face_3 = geompy.MakeFaceObjHW(OZ, gd_r*2, gd_r*2)
Partition_1 = geompy.MakePartition([Sph_tot], [Face_1, Face_2, Face_3], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
[Solid_1,Solid_2,Solid_3,Solid_4,Solid_5,Solid_6,Solid_7,Solid_8] = geompy.ExtractShapes(Partition_1, geompy.ShapeType["SOLID"], True)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

Box_1 = geompy.MakeBoxDXDYDZ(pt_r/2, pt_r/2, pt_r/2)
Box_1_vertex_19 = geompy.GetSubShape(Box_1, [19])
Vertex_1 = geompy.MakeVertex(gd_r, 0, gd_r)
Vertex_2 = geompy.MakeVertex(gd_r, 0, 0)
Box_1_vertex_17 = geompy.GetSubShape(Box_1, [17])
Line_1 = geompy.MakeLineTwoPnt(Vertex_2, Box_1_vertex_17)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_1)
Box_1_vertex_16 = geompy.GetSubShape(Box_1, [16])
Line_3 = geompy.MakeLineTwoPnt(Vertex_1, Box_1_vertex_16)
Vertex_3 = geompy.MakeVertex(gd_r, gd_r, 0)
Vertex_4 = geompy.MakeVertex(gd_r, gd_r, gd_r)
Vertex_5 = geompy.MakeVertex(0, gd_r, gd_r)
Vertex_6 = geompy.MakeVertex(0, gd_r, 0)
Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Box_1_vertex_19)
Line_4_vertex_2 = geompy.GetSubShape(Line_4, [2])
Line_5 = geompy.MakeLineTwoPnt(Vertex_3, Line_4_vertex_2)
Box_1_vertex_21 = geompy.GetSubShape(Box_1, [21])
Line_5_vertex_2 = geompy.GetSubShape(Line_5, [2])
Line_6 = geompy.MakeLineTwoPnt(Box_1_vertex_21, Line_5_vertex_2)
Line_7 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_6)
Box_1_vertex_11 = geompy.GetSubShape(Box_1, [11])
Line_8 = geompy.MakeLineTwoPnt(Vertex_6, Box_1_vertex_11)
Line_7_vertex_2 = geompy.GetSubShape(Line_7, [2])
Box_1_vertex_9 = geompy.GetSubShape(Box_1, [9])
Line_9 = geompy.MakeLineTwoPnt(Line_7_vertex_2, Box_1_vertex_9)
Line_3_vertex_2 = geompy.GetSubShape(Line_3, [2])
Line_10 = geompy.MakeLineTwoPnt(Line_4_vertex_2, Line_3_vertex_2)
Line_11 = geompy.MakeLineTwoPnt(Line_7_vertex_2, Line_4_vertex_2)
Line_1_vertex_2 = geompy.GetSubShape(Line_1, [2])
Line_12 = geompy.MakeLineTwoPnt(Box_1_vertex_17, Line_1_vertex_2)
Box_1_edge_15 = geompy.GetSubShape(Box_1, [15])
Wire_1 = geompy.MakeWire([Box_1_edge_15, Line_1, Line_2, Line_3], 1e-07)
Box_1_edge_18 = geompy.GetSubShape(Box_1, [18])
Wire_1_edge_5 = geompy.GetSubShape(Wire_1, [5])
Wire_2 = geompy.MakeWire([Box_1_edge_18, Line_4, Line_10, Wire_1_edge_5], 1e-07)
Box_1_edge_20 = geompy.GetSubShape(Box_1, [20])
Wire_2_edge_2 = geompy.GetSubShape(Wire_2, [2])
Wire_3 = geompy.MakeWire([Box_1_edge_20, Line_5, Line_6, Wire_2_edge_2], 1e-07)
Box_1_edge_30 = geompy.GetSubShape(Box_1, [30])
Wire_3_edge_5 = geompy.GetSubShape(Wire_3, [5])
Wire_4 = geompy.MakeWire([Box_1_edge_30, Line_9, Line_11, Wire_3_edge_5], 1e-07)
Face_2 = geompy.MakeFaceWires([Wire_2], 1)
Face_3 = geompy.MakeFaceWires([Wire_3], 1)
Face_4 = geompy.MakeFaceWires([Wire_4], 1)
Face_5 = geompy.MakeFaceWires([Wire_4], 1)

Partition_2 = geompy.MakePartition([Solid_8], [Face_2, Face_3, Face_4, Face_5], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
Group_1 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Group_1, [47])
Sym_y = geompy.CreateGroup(Partition_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Sym_y, [68, 50])
Sym_x = geompy.CreateGroup(Partition_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Sym_x, [73, 31])
Sym_z = geompy.CreateGroup(Partition_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Sym_z, [14, 45])
Press = geompy.CreateGroup(Partition_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Press, [76, 58, 34])
Test_dz = geompy.CreateGroup(Partition_2, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(Test_dz, [71])
Test_dx = geompy.CreateGroup(Partition_2, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(Test_dx, [48])
Test_dy = geompy.CreateGroup(Partition_2, geompy.ShapeType["VERTEX"])
geompy.UnionIDs(Test_dy, [19])

geompy.addToStudy( Partition_2, 'Partition_2' )
geompy.addToStudyInFather( Partition_2, Group_1, 'Group_1' )
geompy.addToStudyInFather( Partition_2, Sym_y, 'Sym_y' )
geompy.addToStudyInFather( Partition_2, Sym_x, 'Sym_x' )
geompy.addToStudyInFather( Partition_2, Sym_z, 'Sym_z' )
geompy.addToStudyInFather( Partition_2, Press, 'Press' )
geompy.addToStudyInFather( Partition_2, Test_dz, 'Test_dz' )
geompy.addToStudyInFather( Partition_2, Test_dx, 'Test_dx' )
geompy.addToStudyInFather( Partition_2, Test_dy, 'Test_dy' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)


Mesh_2 = smesh.Mesh(Partition_2)
Regular_1D_1 = Mesh_2.Segment()

Number_of_Segments_1 = Regular_1D_1.NumberOfSegments(nb_seg_arc)
status = Mesh_2.AddHypothesis(Number_of_Segments_1)

Quadrangle_2D_1 = Mesh_2.Quadrangle(algo=smeshBuilder.QUADRANGLE)

Hexa_3D_1 = Mesh_2.Hexahedron(algo=smeshBuilder.Hexa)

Quadratic_Mesh_1 = Regular_1D_1.QuadraticMesh()

Regular_1D_2 = Mesh_2.Segment(geom=Group_1)
Number_of_Segments_2 = Regular_1D_2.NumberOfSegments(nb_seg_ep)
status = Mesh_2.AddHypothesis(Quadratic_Mesh_1,Group_1)
Propagation_of_1D_Hyp = Regular_1D_2.Propagation()
isDone = Mesh_2.Compute()
Sub_mesh_1 = Regular_1D_2.GetSubMesh()
coincident_nodes_on_part = Mesh_2.FindCoincidentNodesOnPart( [ Mesh_2 ], 1e-05, [], 0 )
Mesh_2.MergeNodes(coincident_nodes_on_part, [], 0)

Sym_y_1 = Mesh_2.GroupOnGeom(Sym_y,'Sym_y',SMESH.FACE)
Sym_x_1 = Mesh_2.GroupOnGeom(Sym_x,'Sym_x',SMESH.FACE)
Sym_z_1 = Mesh_2.GroupOnGeom(Sym_z,'Sym_z',SMESH.FACE)
Press_1 = Mesh_2.GroupOnGeom(Press,'Press',SMESH.FACE)
Test_dx_1 = Mesh_2.GroupOnGeom(Test_dx,'Test_dx',SMESH.NODE)
Test_dy_1 = Mesh_2.GroupOnGeom(Test_dy,'Test_dy',SMESH.NODE)
Test_dz_1 = Mesh_2.GroupOnGeom(Test_dz,'Test_dz',SMESH.NODE)
## Set names of Mesh objects

smesh.SetName(Mesh_2.GetMesh(), 'Sphere')
smesh.SetName(Sym_y_1, 'Sym_y')
smesh.SetName(Sym_x_1, 'Sym_x')
smesh.SetName(Sym_z_1, 'Sym_z')
smesh.SetName(Test_dx_1, 'Test_dx')
smesh.SetName(Test_dy_1, 'Test_dy')
smesh.SetName(Test_dz_1, 'Test_dz')
smesh.SetName(Press_1, 'Press')

# try:
#   Mesh_2.ExportMED(r'/home/A21173/Travail/2021/3M/LP/sph.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=1)
#   pass
# except:
#   print('ExportMED() failed. Invalid file name?')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()


