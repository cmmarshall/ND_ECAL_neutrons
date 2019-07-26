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

def dist(TVec1, TVec2):
    diff = TVec1 - TVec2
    return sqrt(diff.Dot(diff))





def Cone_Reject(Candidates, vtx): #Cand
    Candidates = sorted(Candidates, key = lambda cluster: dist(vtx, cluster.getPos()))
    output = []; rejec_clusters = []
    for i, clusteri in enumerate(Candidates):
        if i not in rejec_clusters:
            for j, clusterj in enumerate(Candidates):
                vtx_2_i = clusteri.getPos() - vtx; magi = sqrt(vtx_2_i.Dot(vtx_2_i))
                vtx_2_j = clusterj.getPos() - vtx; magj = sqrt(vtx_2_j.Dot(vtx_2_j))
                angle = acos(vtx_2_i.Dot(vtx_2_j)/(magi*magj))
                if i != j and angle <= pi/6:
                    rejec_clusters.append(j)
            output.append(clusteri)
    return output


def GetRecoE(Candidate,vtx):
    c = 1 #Figure out what the speed of light is...
    mn = 939.565
    beta = dist(Candidate.getPos(), vtx)/(c*Candidate.getTime())
    gamma = 1/sqrt(1 - beta**2);
    return mn*(gamma - 1)


def GetRockEvts(t0):
    rando = ROOT.TRandom3(random.randint(0,int(1E5)))
    NoRockin10mus = 1
    Mean_Rock_Evts = ((t0 + 100)/1000)*NoRockin10mus
    Poisson_Rock = rando.Poisson(Mean_Rock_Evts)
    return Poisson_Rock

def GetHallEvts(t0):
    rando = ROOT.TRandom3(random.randint(0,int(1E5)))
    NoHallin10mus = 1
    Mean_Hall_Evts = ((t0 + 100)/1000)*NoHallin10mus
    Poisson_Hall = rando.Poisson(Mean_Hall_Evts)
    return Poisson_Hall


def Initialize_Plot_Dict():
    h_dict = {}
    h_dict['Reco'] = ROOT.TH1D('Reco', 'Reconstructed Neutron Energy Distribution; Reconstructed Neutron Energy (MeV);', 5000,0, 700)
    h_dict['FReco'] = ROOT.TH1D('Fractional_Reco', 'Fractional Reconstructed Energy Residual; Fractional Residual;', 500, 0, 5)
    hEffic['Eff'] =  ROOT.TH1D('Efficiency', 'Reconstruction Efficiency; Efficiency;', 100, 0, 1)


def loop(GTree, RTree, HTree, Emin, dt, Plot_dict):
    for entry in GTree: #Each Entry corresponds to a Neutrino
        vtx = ROOT.TVector3(entry.vtxX, entry.vtxY, entry.vtxZ)
        vtxA = entry.vtxA
        #pick a time from a flat distribution from 0, 10 microseconds
        t0 = random.uniform(0,1000)
        GCandidates = [NeutronCandidate(entry, index) for index in range(entry.nCandidates)]
        for cand in GCandidates:
            cand.time+=t0
            cand.time = random.normalvariate(cand.time, dt)

        Candidates = GCandidates; NoMax = GetRockEvts(t0)
        for No, Rentry in enumerate(RTree):
            if No >= NoMax:
                break;
            RCandidates = [NeutronCandidate(Rentry, index) for index in range(Rentry.nCandidates)]
            t1 = random.uniform(0,t0+100)
            for cand in RCandidates:
                cand.time+=t1
                cand.time = random.normalvariate(cand.time, dt)
            Candidates += RCandidates

        NoMax = GetHallEvts(t0)
        for No, Hentry in enumerate(HTree):
            if No >= NoMax:
                break;
            HCandidates = [NeutronCandidate(Hentry, index) for index in range(Hentry.nCandidates)]]
            t1 = random.uniform(0,t0+100)
            for cand in HCandidates:
                cand.time+=t1
                cand.time = random.normalvariate(cand.time, dt)
            Candidates += HCandidates

        Candidates = sorted(Cone_Reject(Candidates, vtx), key = lambda cluster: random.normalvariate(cluster.getTime(), dt))
        RecoE = GetRecoE(Candidates[0], vtx)











if __name__ == "__main__":

    #Edit Defaults
    parser = OptionParser()
    parser.add_option('--outfile', help='Output File Name', default='FSout.root')

    parser.add_option('--topdir', help='Directory containing Input', default='/pnfs/dune/persistent/users/rsahay/neutronCandidates')
    parser.add_option('--Gfile_str', help='Input Files for Gas Argon File String', default='GArOut.root')
    parser.add_option('--Rfile_str', help='Input Files for Rock File String', default='WorldOut.root')
    parser.add_option('--Hfile_str', help='Input Files for Hall File String', default='DetOut.root')
    parser.add_option('--Emin', help='Energy Threshhold for Neutron Candidates', default= 5)
    parser.add_option('--dt', help='Time Resolution of our Detector', default=0.1)

    (args, dummy) = parser.parse_args()
    Gfile_str = args.Gfile_str;
    Rfile_str = args.Rfile_str;
    Hfile_str = args.Hfile_str;
    thresh = args.Emin; dt = args.dt

    GExists = path.exists('%s/%s'%(topdir, Gfile_str))
    RExists = path.exists('%s/%s'%(topdir, Rfile_str))
    HExists = path.exists('%s/%s'%(topdir, Hfile_str))

    if GExists and RExists and HExists:
        Gfile = ROOT.TFile('%s/%s'%(topdir, Gfile_str))
        Rfile = ROOT.TFile('%s/%s'%(topdir, Rfile_str))
        Hfile = ROOT.TFile('%s/%s'%(topdir, Hfile_str))
    else:
        print('Does GAr File Exist:%r'%(GExists))
        print('Does Rock File Exist:%r'%(RExists))
        print('Does Hall File Exist:%r'%(HExists))
        sys.exit()

    GTree =  Gfile.Get('argon')
    RTree = Rfile.Get('argon')
    HTree = Hfile.Get('argon')

    Plot_Dict = Initialize_Plot_Dict()
