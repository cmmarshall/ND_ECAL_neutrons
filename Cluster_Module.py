#!/usr/bin/env python

import ROOT
from math import sqrt

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
        self.cells = {}
        self.energy = None
        self.time = None
        self.centroid = None
        self.true_pdg = None
        self.parent_tid = None
        self.primary_tid = None
        self.sigmas = None

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
        self.sigmas = None

    # Calculate stuff about this cluster
    # This is intended to be called before filling ntuple, and isn't fast
    # Don't call this function until you need to know the centroid, cluster energy, etc.
    def CalcStuff(self):
        self.centroid = ROOT.TVector3(0.,0.,0.) # initialize so we can += it later
        self.energy = 0.
        self.time = 999999999999.9

        pdg_energy = {} # map from PDG making energy deposits (i.e. electron, proton) --> energy
        par_energy = {} # map of neutral parent TID (i.e. neutron, photon) --> energy
        prim_energy = {} # map of primary parent TID --> energy

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

            # see if this hit is in a cell already included in the cluster
            if hit.getVolName() not in self.cells:
                self.cells[hit.getVolName()] = hit.getEnergy()
            else:
                self.cells[hit.getVolName()] += hit.getEnergy()

        self.centroid *= (1. / self.energy) # scale to make meaningful energy-weighted position

        # sigmas from the centroid in each dimension
        self.sigmas = [ 0., 0., 0. ]
        for hit in self.hits:
            self.sigmas[0] += (hit.getPos().x() - self.centroid.x())**2 * hit.getEnergy()
            self.sigmas[1] += (hit.getPos().y() - self.centroid.y())**2 * hit.getEnergy()
            self.sigmas[2] += (hit.getPos().z() - self.centroid.z())**2 * hit.getEnergy()
        for i in range(3): self.sigmas[i] = sqrt(self.sigmas[i] / self.energy)

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
        for prim in prim_energy:
            if prim_energy[prim] > max_prim:
                max_prim = prim_energy[prim]
                self.primary_tid = prim

    # simple getters for the calculated things
    def getEnergy(self):
        if self.energy is None:
            CalcStuff()
        return self.energy

    def getTime(self):
        if self.time is None:
            CalcStuff()
        return self.time

    def getCells(self):
        if not len(self.cells):
            CalcStuff()
        return self.cells

    def getCentroid(self):
        if self.centroid is None:
            CalcStuff()
        return self.centroid

    def getTruePDG(self):
        if self.true_pdg is None:
            CalcStuff()
        return self.true_pdg

    def getTrueParent(self):
        if self.parent_tid is None:
            CalcStuff()
        return self.parent_tid

    def getPrimary(self):
        if self.primary_tid is None:
            CalcStuff()
        return self.primary_tid

    def getNcell(self, cut):
        n = 0
        for cell in self.cells:
            if self.cells[cell] > cut:
                n += 1
        return n

    def getMaxCell(self):
        m = 0.
        for cell in self.cells:
            if self.cells[cell] > m:
                m = self.cells[cell]
        return m

    def getSigmas(self):
        if self.sigmas is None:
            CalcStuff()
        return self.sigmas



