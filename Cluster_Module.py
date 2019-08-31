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

def Matrix_Multiply(M, vec): #Matrix Multiplication Method because I can't figure out what is wrong with TMatrixD
    vect = [vec.X(), vec.Y(), vec.Z()]
    output = [0,0,0]
    for i in range(len(M)):
        for j in range(3):
            output[i] += M[i][j]*vect[j]
    return ROOT.TVector3(output[0],output[1],output[2])


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
        self.energyFracPar = None
        self.energyFracPrim = None
        self.voxhits = None
        self.SigmaVX = None
        self.SigmaVY = None
    def Voxelize(self):
        self.voxhits = {} #Thing keyed by voxel
        hits = self.sortHits()
        Cart = {}; Vec = {}
        Cart[0] = ROOT.TVector3(1,0,0)
        Cart[1] = ROOT.TVector3(0,1,0)
        Cart[2] = ROOT.TVector3(0,0,1)
        for key in hits:
            if 'Barrel' in key:
                Gkey = 'BarrelECal_stave0'
                for i in range(1, 9):
                    if 'stave0%d'%i in key:
                        Gkey += str(i) 
                TempMag = (Geo_Pos[Gkey] - Geo_Pos['TPC']).Mag()            
                Vec[2] = (Geo_Pos[Gkey] - Geo_Pos['TPC'])*(1/TempMag)
                Vec[0] = ROOT.TVector3(1,0,0)
                Vec[1] = Vec[2].Cross(Vec[0])
            elif 'Endcap' in key:
                Vec[2] = ROOT.TVector3(1,0,0)
                Vec[0] = ROOT.TVector3(0,1,0)
                Vec[1] = ROOT.TVector3(0,0,1)

            M = [[0,0,0],[0,0,0],[0,0,0]]
            for i in range(3):
                for j in range(3):
                    M[i][j] = Vec[i].Dot(Cart[j])
            
            for hit in hits[key]: 
                Pos = hit.getPos()
                Pos = Matrix_Multiply(M,Pos)
                VoxX = int(Pos.X()/2.2); VoxY = int(Pos.Y()/2.2)
                Vox = '%d,%d,%s'%(VoxX, VoxY, key)
                if Vox not in self.voxhits:
                    self.voxhits[Vox] = []
                self.voxhits[Vox].append(hit)
           
    def getSigmaV(self):
        self.SigmaVX = 0
        self.SigmaVY = 0
        if self.voxhits is None:
            self.Voxelize()
        if 1 ==1:
            VPos = [[int(key.split(',')[0]), int(key.split(',')[1]), key.split(',')[2]] for key in self.voxhits]
            if len(VPos) == 1:
                self.SigmaVX = 0
                self.SigmaVY = 0
            else:
                aux_arr = []
                for i in range(0, len(VPos)-1):
                    for j in range(i+1, len(VPos)):
                        SigmaX = abs(VPos[j][0] - VPos[i][0]); SigmaY = abs(VPos[j][1] - VPos[i][1])
                        aux_arr.append([SigmaX, SigmaY])
                self.SigmaVX = max(aux_arr, key = lambda x: x[0])[0]
                self.SigmaVY = max(aux_arr, key = lambda x: x[1])[1]


    def sortHits(self):
        output_dict = {}
        for hit in self.hits:
	    volName = hit.getVolName()
            key = ''
            if 'Barrel' in volName:
                key += 'Barrel'
            else:
                key += 'Endcap'

            if key == '':
                #print('What the actual hell. The key is %s. The volName is %s.'%(key, volName))
                continue

            key += 'ECal'

            for i in range(1, 9):
                if 'stave0%d'%i in volName:
                     key += '_stave0%d'%i
            
            
            for i in range(1, 61):
                temp_str = 'layer_%d'%i
                if i < 10:
                    temp_str = 'layer_0%d'%i
                #print(temp_str, volName)
                if temp_str in volName:
                    key += '_'; key+= temp_str
            
            if 'layer' not in key or 'stave' not in key:
                #print('What the actual hell. The key is %s. The volName is %s' %(key, volName))
                continue   
            
            if key not in output_dict:
                output_dict[key] = []
            
            output_dict[key].append(hit)
    	return output_dict
       	
    

#    def sortHits(self):
#        output_dict = {}
#        for hit in self.hits:
#            key1 = None ; key2 = None
#            for i in range(1, 9):
#                if 'stave0%d'%i in hit.getVolName():
#                    key1 = 'stave0%d'%i
#            for i in range(1,10):
#                if 'layer_0%d'%i in hit.getVolName():
#                    key2 = 'layer%d'%i
#            for i in range(10, 61):
#                if 'layer_%d'%i in hit.getVolName():
#                    key2 = 'layer%d'%i
#
#            if key1 is None or key2 is None:
#                print('ALERT: Hit is in %s')
#
#            key = key1 + key2
           # if key1 not in output_dict:
           #     output_dict[key1] = {}
           # if key2 not in output_dict:
           #     output_dict[key1][key2] = []
           # if key1 is None or key2 is None:
           #     continue
           # else:
#            output_dict[key1][key2].append(hit)
#        return output_dict


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
        #print('Hello') 
        self.Voxelize()
        self.getSigmaV()
        #print(self.SigmaVX, self.SigmaVY) 
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
                self.energyFracPar = max_par / self.energy
        for prim in prim_energy:
            if prim_energy[prim] > max_prim:
                max_prim = prim_energy[prim]
                self.primary_tid = prim
                self.energyFracPrim = max_prim / self.energy

    # simple getters for the calculated things

    def getVoxHits(self):
        if self.voxhits is None:
            CalcStuff()
        return self.voxhits

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

    def getSigmaVX(self):
        if self.SigmaVX is None:
            CalcStuff()
        return self.SigmaVX

    def getSigmaVY(self):
        if self.SigmaVY is None:
            self.CalcStuff()
        return self.SigmaVY

    def getFracMaxPar(self):
        return self.energyFracPar

    def getFracMaxPrim(self):
        return self.energyFracPrim



