#!/usr/bin/env python

import ROOT
from math import sqrt

Geo_Pos = {}
Geo_Pos['TPC'] = ROOT.TVector3(0, -215.106 , 585.300)
Geo_Pos['BarrelECal_stave01'] = ROOT.TVector3(0, -447.736 , 394.210)
Geo_Pos['BarrelECal_stave02'] = ROOT.TVector3(0, -514.169 , 616.074)
Geo_Pos['BarrelECal_stave03'] = ROOT.TVector3(0, -404.492 , 818.049)
Geo_Pos['BarrelECal_stave04'] = ROOT.TVector3(0, -184.466 , 883.584)
Geo_Pos['BarrelECal_stave05'] = ROOT.TVector3(0, 18.503 , 773.929)
Geo_Pos['BarrelECal_stave06'] = ROOT.TVector3(0, 84.193 , 552.650)
Geo_Pos['BarrelECal_stave07'] = ROOT.TVector3(0, -25.400 , 351.026)
Geo_Pos['BarrelECal_stave08'] = ROOT.TVector3(0, -246.454 , 285.059)

# build coordinate transformations for each stave, so that z axis points directly into the stave
# this simplifies the process of voxylizing within each stave
matrices = [ROOT.TMatrixD(3,3) for i in range(9)]

for i in range(3):
    for j in range(3):
        matrices[0][i][j] = (i == (j+2)%3) # rotate z to point along x
for m in range(1,9):
    newcoord = [None, None, None]
    newcoord[2] = (Geo_Pos["BarrelECal_stave%02d" % m] - Geo_Pos["TPC"]).Unit()
    newcoord[0] = ROOT.TVector3(1,0,0)
    newcoord[1] = newcoord[2].Cross(newcoord[0])

    for i in range(3):
        matrices[m][i][0] = newcoord[i].x()
        matrices[m][i][1] = newcoord[i].y()
        matrices[m][i][2] = newcoord[i].z()

class Hit:

    def __init__(self, pos, volName, energy, time, pdg, neutral_tid, primary_tid):
        self.pos = pos
        self.volName = volName
        self.energy = energy
        self.time = time
        self.pdg = pdg
        self.neutral_tid = neutral_tid # TID of parent photon/neutron
        self.primary_tid = primary_tid # TID of primary parent    

    def getPos(self):
        return self.pos
    def getVolName(self):
        return self.volName
    def getEnergy(self):
        return self.energy
    def getTime(self):
        return self.time
    def getPDG(self):
        return self.pdg
    def getParentTID(self):
        return self.neutral_tid
    def getPrimaryTID(self):
        return self.primary_tid

class Cluster:

    def __init__(self, dist_gap = 5.):
        # Cluster consists of Hit objects
        # Every quantity about the cluster can be determined from its constituent hits
        self.hits = []
        self.dist_gap = dist_gap # allowed gap between hits to make different cluster

        # These are cluster-specific variables that will be None until CalcStuff() is called
        self.voxyls = [{} for i in range(10)] # sorted by stave; 0 is left endcap 10 is right endcap
        self.energy = None
        self.time = None
        self.centroid = None
        self.true_pdg = None
        self.parent_tid = None
        self.primary_tid = None
        self.energyFracPar = None
        self.energyFracPrim = None
        self.sigmaX = None
        self.sigmaY = None
        self.nvoxyls = 0
        self.layers = []
        self.maxvoxyl = 0.

    def Voxylize(self):

        for hit in self.hits:

            volName = hit.getVolName()
            isEndcap = ("Endcap" in volName)

            # Get stave number and layer number from string
            staveidx = volName.find("stave")
            stave_number = int(volName[staveidx+5:staveidx+7])

            layeridx = volName.find("layer_")
            layer_number = int(volName[layeridx+6:layeridx+8])

            matrix = matrices[0] if isEndcap else matrices[stave_number]
            vect = ROOT.TVectorD(3)
            vect[0] = hit.getPos().x()
            vect[1] = hit.getPos().y()
            vect[2] = hit.getPos().z()
            vect *= matrix # this actually does v = M*v
            stave_pos = ROOT.TVector3( vect[0], vect[1], vect[2] )

            voxX = int(stave_pos.x() / 2.2)
            voxY = int(stave_pos.y() / 2.2)

            if isEndcap:
                if hit.getPos().x() < 0.: stave_number = 0
                else: stave_number = 9

            # see if this hit is in an existing voxyl or a new one
            if layer_number in self.voxyls[stave_number]:
                if (voxX,voxY) in self.voxyls[stave_number][layer_number]: # existing voxyl
                    self.voxyls[stave_number][layer_number][(voxX,voxY)] += hit.energy
                else: # new voxyl in existing layer
                    self.voxyls[stave_number][layer_number][(voxX,voxY)] = hit.energy
            else: # new layer
                self.voxyls[stave_number][layer_number] = { (voxX,voxY):hit.energy }

    def calcSigma(self):
        self.sigmaX = 0.
        self.sigmaY = 0.
        self.layers = []
        self.nvoxyls = 0
        self.maxvoxyl = 0.
        # loop over voxyls
        for s,stave in enumerate(self.voxyls):
            for layer_number in stave:
                self.layers.append(layer_number)
                this_layer = self.voxyls[s][layer_number]
                sigX = 0.
                sigY = 0.
                meanX = 0.
                meanY = 0.
                layere = 0.
                for xy in this_layer:
                    e = this_layer[xy]
                    layere += e
                    if e > 0.1: self.nvoxyls += 1
                    if e > self.maxvoxyl: self.maxvoxyl = e
                    meanX += xy[0]*e
                    meanY += xy[1]*e
                meanX /= layere
                meanY /= layere
                for xy in this_layer:
                    sigX += ((xy[0]-meanX)*this_layer[xy])**2/layere
                    sigY += ((xy[1]-meanY)*this_layer[xy])**2/layere
                self.sigmaX += sigX
                self.sigmaY += sigY
        self.sigmaX *= 2.2/len(self.layers) # convert to cm and normalize per layer
        self.sigmaY *= 2.2/len(self.layers)

    def addHit(self, hit):
        # if htis are added, cluster-level variables will be wrong
        # ensure that wrong info is not used by setting everything to garbage
        if self.energy is not None: 
            self.setToBS()

        # add the hit to the cluster
        self.hits.append(hit)

    def getHits(self):
        return self.hits

    # find if a hit should be part of this cluster
    def isInCluster(self, newHit):
        for i,hit in enumerate(self.hits):
            if (hit.getPos()-newHit.getPos()).Mag() < self.dist_gap: # this hit should be in this cluster
                return True
        return False

    def mergeCluster(self, new_cluster):
        for hit in new_cluster.getHits():
            self.addHit(hit)

    # set the cluster-level variables to BS
    def setToBS(self):
        self.cells = {}
        self.energy = None
        self.time = None
        self.centroid = None
        self.true_pdg = None
        self.parent_tid = None
        self.primary_tid = None

    # Calculate stuff about this cluster
    # This is intended to be called before filling ntuple, and isn't fast
    # Don't call this function until you need to know the centroid, cluster energy, etc.
    def CalcStuff(self):
        self.centroid = ROOT.TVector3(0.,0.,0.) # initialize so we can += it later
        self.energy = 0.
        self.time = 999999999999.9

        self.nCylinder = 1

        pdg_energy = {} # map from PDG making energy deposits (i.e. electron, proton) --> energy
        par_energy = {} # map of neutral parent TID (i.e. neutron, photon) --> energy
        prim_energy = {} # map of primary parent TID --> energy
        self.Voxylize()
        for hit in self.hits:
            self.centroid += (hit.getEnergy() * hit.getPos()) # scale position 3-vector by energy
            self.energy += hit.getEnergy()

            if hit.getTime() < self.time:
                self.time = hit.getTime()

            pdg = abs(hit.getPDG()) # don't distinguish between electrons and positrons; there aren't antiphotons or antiprotons
            par = hit.getParentTID()
            prim = hit.getPrimaryTID()
            if pdg in pdg_energy:
                pdg_energy[pdg] += hit.getEnergy()
            else:
                pdg_energy[pdg] = hit.getEnergy()

            if par in par_energy:
                par_energy[par] += hit.getEnergy()
            else:
                par_energy[par] = hit.getEnergy()

            if prim in prim_energy:
                prim_energy[prim] += hit.getEnergy()
            else:
                prim_energy[prim] = hit.getEnergy()                

        self.centroid *= (1. / self.energy) # scale to make meaningful energy-weighted position

        max_pdg = 0.
        max_par = 0.
        max_prim = 0.
        for pdg in pdg_energy:
            if pdg_energy[pdg] > max_pdg:
                max_pdg = pdg_energy[pdg]
                self.true_pdg = pdg
        for par in par_energy:
            if par_energy[par] > max_par:
                max_par = par_energy[par]
                self.parent_tid = par
                self.energyFracPar = max_par / self.energy
        for prim in prim_energy:
            if prim_energy[prim] > max_prim:
                max_prim = prim_energy[prim]
                self.primary_tid = prim
                self.energyFracPrim = max_prim / self.energy
        self.calcSigma()

    # simple getters for the calculated things
    def getEnergy(self):
        return self.energy

    def getTime(self):
        return self.time

    def getCells(self):
        return self.cells

    def getCentroid(self):
        return self.centroid

    def getTruePDG(self):
        return self.true_pdg

    def getTrueParent(self):
        return self.parent_tid

    def getPrimary(self):
        return self.primary_tid

    def getNcell(self):
        return self.nvoxyls

    def getMaxCell(self):
        return self.maxvoxyl

    def getSigma(self):
        return self.sigmaX, self.sigmaY

    def getNlayers(self):
        return len(self.layers)

    def getLayerGaps(self):
        ordered = sorted(self.layers)
        gap = 0
        last = ordered[0]
        for layer in ordered:
            if layer - last > 1:
                gap += layer - last - 1
            last = layer
        return gap

    def getFracMaxPar(self):
        return self.energyFracPar

    def getFracMaxPrim(self):
        return self.energyFracPrim



