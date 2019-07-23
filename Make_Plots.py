#!/usr/bin/env python

import sys
import os.path
from os import path
import os
import ROOT
from optparse import OptionParser
from array import array
from math import sqrt
from math import acos
from math import pi
import subprocess
import random


class NeutronCandidate:
    def __init__(self, entry, index, source):
        self.pos = ROOT.TVector3(entry.nPosX[index], entry.nPosY[index], entry.nPosZ[index])
        self.time = entry.nPosT[index]
        self.E = entry.nE[index]
        self.Emax = entry.nEmax[index]
        self.PDG = entry.nTruePDG[index]
        self.isPrimary = entry.isPrimary[index]
        self.TrueKE = entry.nTrueKE[index]
        self.index = index
        self.source = source
    def getPos(self):
        return self.pos
    def getTime(self):
        return self.time
    def getE(self):
        return self.E
    def getEmax(self):
        return self.Emax
    def getPDG(self):
        return self.PDG
    def getisPrimary(self):
        return self.isPrimary
    def getTrueKE(self):
        return self.TrueKE
    def getIndex(self):
        return self.index
    def getSource(self):
        return self.source



def Initialize_Plot_Dicts(str1 = '', str2 = '', includevtx = True): #String 1 can be whatever but string 2 should be "True Neutron: " or "True Photon: "
    hDict = {}
    xmin = -1000; ymin = -1000; zmin = -1000; Tmin = 0
    xmax = 1000; ymax = 1000; zmax = 1000; Tmax = 1000

    if includevtx:
        hDict['vtxPosX'] = ROOT.TH1D(str1 + "Nu vtx X", "Nu Vertex X Distribution; X-Position (cm);", 200, xmin, xmax  )
        hDict['vtxPosY'] = ROOT.TH1D(str1 + "Nu vtx Y", "Nu Vertex Y Distribution; Y-Position (cm);", 200, ymin, ymax  )
        hDict['vtxPosZ'] = ROOT.TH1D(str1 + "Nu vtx Z", "Nu Vertex Z Distribution; Z-Position (cm);", 200, zmin, zmax  )
        hDict['vtxPosXY']= ROOT.TH2D(str1 + "Nu vtx XvY", "Nu Vertex XY Distribution; Y-Position (cm); X-Position (cm)", 200, xmin, xmax, 200, ymin, ymax )
        hDict['vtxPosXZ']= ROOT.TH2D(str1 +  "Nu vtx XvZ", "Nu Vertex XZ Distribution; Z-Position (cm); X-Position (cm)", 200, xmin, xmax, 200, zmin, zmax )
        hDict['vtxPosYZ']= ROOT.TH2D(str1 + "Nu vtx YvZ", "Nu Vertex YZ Distribution; Z-Position (cm); Y-Position (cm)", 200, ymin, ymax, 200, zmin, zmax )
        #hDict['vtxTime'] = ROOT.TH1D(str1 + "Nu vtx Time", "Nu Vertex Time Distribution; Time (ns);" 100, tmin, tmax)

    xmin = -1000; ymin = -1000; zmin = -1000; tmin = 0
    xmax = 1000; ymax = 1000; zmax = 1000; tmax = 1000
    rmin = 0; rmax = 1000
    hDict['CX'] = ROOT.TH1D(str1 +  "Candidate X", str2 + "Neutron Candidate X Distribution; X-Position (cm);", 200, xmin, xmax  )
    hDict['CY'] = ROOT.TH1D(str1 +  "Candidate Y", str2 + "Neutron Candidate Y Distribution; Y-Position (cm);", 200, ymin, ymax  )
    hDict['CZ'] = ROOT.TH1D(str1 +  "Candidate Z", str2 + "Neutron Candidate Z Distribution; Z-Position (cm);", 200, zmin, zmax  )
    hDict['CXY']= ROOT.TH2D(str1 +  "Candidate XvY", str2 + "Neutron Candidate XY Distribution; Y-Position (cm); X-Position (cm)", 200, xmin, xmax, 200, ymin, ymax )
    hDict['CXZ']= ROOT.TH2D(str1 +  "Candidate XvZ", str2 + "Neutron Candidate XZ Distribution; Z-Position (cm); X-Position (cm)", 200, xmin, xmax, 200, zmin, zmax )
    hDict['CYZ']= ROOT.TH2D(str1 +  "Candidate YvZ", str2 + "Neutron Candidate YZ Distribution; Z-Position (cm); Y-Position (cm)", 200, ymin, ymax, 200, zmin, zmax )
    hDict['CT'] = ROOT.TH1D(str1 +  "Candidate RTime", str2 + "Neutron Candidate Time Distribution (Relative to Nu Vertex); Time (ns);", 100, tmin, tmax)
    hDict['CR'] = ROOT.TH1D(str1 +  "Candidate Distance from Nu VtX",str2 + "Candidate Distance from Nu Vertex Distribution; Distance (cm);", 100, rmin, rmax  )

    hDict['CVel'] = ROOT.TH1D(str1 + 'Candidate Velocity', 'Candidate Beta Distribution; Beta;', 1000, 0, 1)
    Nmin = 0; Nmax = 0
    TNmin = 0; TNmax = 0
#    hDict['TrueNCvNC'] = ROOT.TH2D(str1 + 'NoTrueNCvNC', 'Number of True Neutron Candidates vs Neutron Candidates; Number of Neutron Candidates; Number of True Neutrons', 100, Tmin, Tmax, 100, TNmin, TNmax)

    return hDict

def Fill_Vertex_Info(Plot_Dict, vtx):
    PosX = vertex.X(); PosY = vtx.Y(); PosZ = vertex.Z()
    Plot_Dict['vtxPosX'].Fill(PosX)
    Plot_Dict['vtxPosY'].Fill(PosY)
    Plot_Dict['vtxPosZ'].Fill(PosZ)
    Plot_Dict['vtxPosXY'].Fill(PosX, PosY)
    Plot_Dict['vtxPosXZ'].Fill(PosX, PosZ)
    Plot_Dict['vtxPosYZ'].Fill(PosY, PosZ)

def Fill_Candidate_Info(Plot_Dict, Candidates, vertex):
    c = 29.9792# in cm/ns
    vtx = ROOT.TVector3(vertex.Position.X()/10., vertex.Position.Y()/10., vertex.Position.Z()/10.)
    for cluster in Candidates:
        cposx = cluster.getPos().X(); cposy = cluster.getPos().Y(); cposz = cluster.getPos().Z()
        cpost = cluster.getTime(); diff = cluster.getPos() - vtx; cposr = sqrt(diff.Dot(diff))
        Plot_Dict['CX'].Fill(cposx); Plot_Dict['CY'].Fill(cposy); Plot_Dict['CZ'].Fill(cposz)
        Plot_Dict['CXY'].Fill(cposx,cposy); Plot_Dict['CXZ'].Fill(cposx, cposz); Plot_Dict['CYZ'].Fill(cposy, cposz)
        Plot_Dict['CT'].Fill(cpost); Plot_Dict['CR'].Fill(cposr)
        Plot_Dict['CVel'].Fill(abs(cposr/cpost)/c)

def loop(tree):
    NPlot_Dict = Initialize_Plot_Dicts(str1 = '1', str2 = 'True Neutron: ')
    GPlot_Dict = Initialize_Plot_Dicts(str1 = '2', str2 = 'True Photon: ', includevtx = False )
    for entry in tree:
        vtx = ROOT.TVector3(entry.vtxX, entry.vtxY, entry.vtxZ)
        Fill_Vertex_Info(NPlot_Dict)
        Candidates = [NeutronCandidate(entry, index) for index in range(entry.nCandidates)]
        TrueNeutronC = [cluster for cluster in Candidates if cluster.PDG == 2112]
        TruePhotonC = [cluster for cluster in Candidates if cluster.PDG == 22]
        Fill_Candidate_Info(NPlot_Dict, TrueNeutronC)
        Fill_Candidate_Info(GPlot_Dict, TruePhotonC)
    NOut = ROOT.TFile('TrueNeutron.root', "RECREATE")
    c = ROOT.TCanvas()
    for key in NPlot_Dict:
        if 'vtx' not in key:
            print(key)
            NPlot_Dict[key].Write()
            NPlot_Dict[key].Draw()
            c.Print('TrueNeutron_%s.eps'%(key))
    del NOut

    NuOut = ROOT.TFile('Neutrino.root', "RECREATE")
    for key in NPlot_Dict:
        if 'vtx' in key:
            NPlot_Dict[key].Write()
            NPlot_Dict[key].Draw()
            c.Print('%s.eps'%key)
    del NuOut

    GOut = ROOT.TFile('TruePhoton.root', "RECREATE")
    for key in GPlot_Dict:
        GPlot_Dict[key].Write()
        GPlot_Dict[key].Draw()
        c.Print('TruePhoton_%s.eps'%key)
    del GOut



if __name__ == "__main__":
    #Edit Defaults
    parser = OptionParser()
    parser.add_option('--outdir', help="Output Directory", default='/pnfs/dune/persistent/users/rsahay/Test/')
    parser.add_option('--infile', help="Input Directory", default='/dune/app/users/rsahay/ND_ECAL_neutrons/out.root')
    (args, dummy) = parser.parse_args()
    File_Exists = path.exists(args.infile)
    if File_Exists:
        InFile = ROOT.TFile('%s.root'%File_Exists)
        tree = InFile.Get('tree')
        loop(tree)
