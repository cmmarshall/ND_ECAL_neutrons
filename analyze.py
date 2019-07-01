#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array
from math import sqrt
import numpy as np

MAXCANDIDATES = 1000

# Output TTree array variables
t_run = array( 'i', [0] )
t_event = array( 'i', [0] )
t_vtxX = array( 'd', [0.] )
t_vtxY = array( 'd', [0.] )
t_vtxZ = array( 'd', [0.] )
t_vtxA = array( 'i', [0] )
t_nCandidates = array( 'i', [0] )
t_nPosX = array( 'd', MAXCANDIDATES*[0.] )
t_nPosY = array( 'd', MAXCANDIDATES*[0.] )
t_nPosZ = array( 'd', MAXCANDIDATES*[0.] )
t_nPosT = array( 'd', MAXCANDIDATES*[0.] )
t_nE = array( 'd', MAXCANDIDATES*[0.] )
t_nEmax = array( 'd', MAXCANDIDATES*[0.] )
t_nNcell = array( 'i', MAXCANDIDATES*[0] )
t_nNcell1MeV = array( 'i', MAXCANDIDATES*[0] )
t_nNcell3MeV = array( 'i', MAXCANDIDATES*[0] )
t_nTruePDG = array( 'i', MAXCANDIDATES*[0] )
t_nIsPrimary = array( 'i', MAXCANDIDATES*[0] )
t_nTrueKE = array( 'd', MAXCANDIDATES*[0.] )

class NeutronCandidate:

    def __init__(self, tid, pdg, ke):
        self.intpt = None
        self.inttime = None
        self.hits = []
        self.energy = 0.
        self.tid = tid
        self.cells = {}
        self.nhits = 0
        self.truePDG = pdg
        self.trueKE = ke

    def addHit(self, pos, volName, hEnergy, hTime, parent, neutral_tid):
        self.nhits += 1
        self.hits.append([pos, parent, neutral_tid])

        if volName not in self.cells:
            self.cells[volName] = hEnergy
        else:
            self.cells[volName] += hEnergy

        self.energy += hEnergy

        if self.inttime is None:
            self.inttime = hTime
        elif hTime < self.inttime:
            self.inttime = hTime

        if self.intpt is None:
            self.intpt = pos*hEnergy
        else:
            self.intpt += pos*hEnergy # energy-weighted average

    def getHits(self):
        return self.hits

    def getPos(self):
        return (1./self.energy) * self.intpt

    def getTime(self):
        return self.inttime

    def getEnergy(self):
        return self.energy

    def getNcell(self):
        return len(self.cells)

    def getNcellCut(self, cut):
        count = 0
        for c in self.cells:
            if self.cells[c] > cut:
                count += 1
        return count

    def getMaxCell(self):
        m = 0.
        for c in self.cells:
            if self.cells[c] > m:
                m = self.cells[c]
        return m

    def getTruePDG(self):
        return self.truePDG

    def getTrueKE(self):
        return self.trueKE

    def printNC(self):
        print "  Neutron candidate from track ID %d at time %1.1f" % (self.tid, self.inttime)
        print "     %d hits on %d cells, total energy %1.1f MeV, centroid (%1.1f, %1.1f, %1.1f)" % (self.nhits, len(self.cells), self.energy, self.intpt.x()/self.energy, self.intpt.y()/self.energy, self.intpt.z()/self.energy)



def isinCluster(candidate, cluster, thresh):
    #candidate should be a T3Vector
    #cluster is a NeutronCandidate object
    hits = cluster.getHits()
    for hit in hits:
        diff = hit[0] - candidate
        distance = sqrt(diff.dot(diff))
        if distance <= thresh:
            return True
    return False



def setToBogus():
    t_run[0] = 0
    t_event[0] = 0
    t_vtxX[0] = 0.
    t_vtxY[0] = 0.
    t_vtxZ[0] = 0.
    t_vtxA[0] = 0
    t_nCandidates[0] = 0

# return track ID of most immediate neutron parent, if there is one (otherwise return 0)
def neutronParent( event, tid ):
    traj = event.Trajectories[tid]
    if traj.PDGCode == 2112: # this hit is from a neutron
        return tid
    parentid = traj.ParentId
    while parentid != -1: # primary particles have parentid -1
        parent = event.Trajectories[parentid]
        if parent.PDGCode == 2112:
            return parentid
        parentid = parent.ParentId # recursively go through ancestry to look for neutrons
    return -1 # no neutron parent

# return track ID of earliest photon parent; this is different from the neutron case because EM showers
# will have zillions of photons
def photonParent( event, tid ):

    traj = event.Trajectories[tid]
    if traj.PDGCode not in [11, -11, 22]:
        return -1 # not an EM particle, so it won't have a photon parent

    parentid = traj.ParentId
    parent = event.Trajectories[parentid]
    # recursively go back through the ancestry of EM shower to find the track ID of the first EM particle
    while parentid != -1 and parent.PDGCode in [11, -11, 22]: # primary particles have parentid -1
        tmpid = parent.ParentId
        if event.Trajectories[tmpid].PDGCode not in [11, -11, 22]: # non-EM parent, we've found the first one
            break
        parentid = tmpid # parent is EM, so keep searching back
        parent = event.Trajectories[parentid]

    if parent.PDGCode == 22:
        return parentid
    else: # if earliest EM particle is electron/positron, then you have a track entering the detector and we don't need to worry about it
        return -1

def loop( events, tgeo, tree, Cluster_Threshold = 1 ): # ** CHRIS: WHAT SHOULD I SET THE DEFAULT THRESHOLD TO **

    offset = [ 0., 305., 5. ]

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()
    for ient in range(N):
        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)

        events.GetEntry(ient)
        for ivtx,vertex in enumerate(event.Primaries):

            setToBogus()

            reaction = vertex.Reaction #vertex.GetReaction()
            t_vtxX[0] = vertex.Position.X()/10. #vertex.GetPosition().X()/10. # cm
            t_vtxY[0] = vertex.Position.Y()/10. #vertex.GetPosition().Y()/10.
            t_vtxZ[0] = vertex.Position.Z()/10. #vertex.GetPosition().Z()/10.
            # atomic mass of target
            t_vtxA[0] = int((reaction.split(";")[1].split(":")[1])[6:9])

            #node = tgeo.FindNode( vertex.Position[0], vertex.Position[1], vertex.Position[2] )
            #print "Interaction vertex in %s at (%1.1f, %1.1f, %1.1f)" % (node.GetName(), vertex.Position[0]/10., vertex.Position[1]/10., vertex.Position[2]/10.)

            # Loop over primary particles only, i.e. direct products of neutrino interaction
            #neutron_tids = []
            #for ipart,particle in enumerate(vertex.Particles):
            #    e = particle.GetMomentum().E()
            #    p = (particle.GetMomentum().X()**2 + particle.GetMomentum().Y()**2 + particle.GetMomentum().Z()**2)**0.5
            #    m = (e**2 - p**2)**0.5
            #    pdg = particle.GetPDGCode()

            # Loop over all true trajectories
            #for traj in event.Trajectories:
            #    mom = traj.GetParentId()
            #    tid = traj.GetTrackId()
            #    pid = traj.GetPDGCode()

            # Get ECal hits
            ecal_hits = []
            for det in event.SegmentDetectors:
                if "ECal" in det.first: # there are several sub-detectors that make up the ECal
                    ecal_hits += det.second

            candidates = []
            for hit in ecal_hits:
                if hit.EnergyDeposit < 0.01: # skip tiny deposits, this cut needs to be tuned
                    continue

                hStart = ROOT.TVector3( hit.Start.X()/10., hit.Start.Y()/10., hit.Start.Z()/10. )
                hStop = ROOT.TVector3( hit.Stop.X()/10., hit.Stop.Y()/10., hit.Stop.Z()/10. )

                #Truth-match the hits
                #neutral_tid = -1
                parent = None; neutral_tid = -1
                for i in range(len(hit.Contrib)):
                    tid = hit.Contrib[i]
                    photon_tid = photonParent(event, tid); neutron_tid = neutronParent( event, tid )
                    if neutron_tid > 0:
                        parent = 'n'; neutral_tid = neutron_tid
                    else:
                        if photon_tid > 1:
                            parent = 'g'; neutral_tid = photon_tid
                        else:
                            parent = None; neutral_tid = photon_tid

                    if parent is not None:
                        inCluster = False
                        for cluster in candidates:
                            if isinCluster(hStart, cluster, Cluster_Threshold) and not inCluster:
                                cluster.addHit(hStart, node.GetName(), hit.EnergyDeposit, hit.Start[3], parent)
                                inCluster = True
                        if not inCluster:
                            truePDG = event.Trajectories[neutral_tid].PDGCode
                            mom = event.Trajectories[neutral_tid].InitialMomentum
                            c = NeutronCandidate(neutral_tid, truePDG, mom.E()-mom.M())
                            c.addHit(hStart, node.GetName(), hit.EnergyDeposit, hit.Start[3], parent, int(neutral_tid))
                            candidates.append(c)


            """
            candidates = {}
            for hit in ecal_hits:

                if hit.EnergyDeposit < 0.01: # skip tiny deposits,if this cut needs to be tuned
                    continue

                #hStart = ROOT.TVector3( hit.GetStart().X()/10., hit.GetStart().Y()/10., hit.GetStart().Z()/10. )
                #hStop = ROOT.TVector3( hit.GetStop().X()/10., hit.GetStop().Y()/10., hit.GetStop().Z()/10. )
                hStart = ROOT.TVector3( hit.Start.X()/10., hit.Start.Y()/10., hit.Start.Z()/10. )
                hStop = ROOT.TVector3( hit.Stop.X()/10., hit.Stop.Y()/10., hit.Stop.Z()/10. )

                # Truth-match the hits
                neutral_tid = -1
                for i in range(len(hit.Contrib)): # most hits have only one contributor
                    tid = hit.Contrib[i]
                    # see if this hit is from a neutron or photon, including one that was produced in secondary interactions
                    # this basically requires that all contributors be neutral descendents, which is fine
                    photon_tid = photonParent( event, tid )
                    neutron_tid = neutronParent( event, tid )
                    #print "ECal Hit due to %d has photon %d neutron %d" % (tid, photon_tid, neutron_tid)
                    if neutron_tid > 0:
                        neutral_tid = neutron_tid
                    else:
                        neutral_tid = photon_tid

                if neutral_tid > 1: # neutron- or photon-iduced thing, make a candidate or add to an existing one
                    if neutral_tid not in candidates: # first hit from this particular neutron/photon
                        truePDG = event.Trajectories[neutral_tid].PDGCode
                        mom = event.Trajectories[neutral_tid].InitialMomentum
                        c = NeutronCandidate( neutral_tid, truePDG, mom.E()-mom.M() )
                        candidates[neutral_tid] = c

                    # get detector element where this hit occured
                    #node = tgeo.FindNode( hit.GetStart.X(), hit.GetStart.Y(), hit.GetStart.Z() )
                    node = tgeo.FindNode( hit.Start.X(), hit.Start.Y(), hit.Start.Z() )
                    #print "  adding hit to tid %d: (%1.1f, %1.1f, %1.1f), %s, %1.2f MeV, time = %1.2f" % (neutral_tid, hStart.x(), hStart.y(), hStart.z(), node.GetName(), hit.EnergyDeposit, hit.Start[3])
                    candidates[neutral_tid].addHit( hStart, node.GetName(), hit.EnergyDeposit, hit.Start[3] )
            """;
            t_nCandidates[0] = 0
            for key in range(len(candidates)):
                #isPrimary = (event.Trajectories[key].ParentId == -1)
                #for hit in candidate[key].getHits():
                neutral_tids = list(np.array(candidate[key].getHits())[:,2])
                largestContrib = max(set(neutral_tids), key=neutral_tids.count)
                isPrimary = (event.Trajectories[largestContrib].ParentId == -1)
                t_nPosX[t_nCandidates[0]] = candidates[key].getPos().x()
                t_nPosY[t_nCandidates[0]] = candidates[key].getPos().y()
                t_nPosZ[t_nCandidates[0]] = candidates[key].getPos().z()
                t_nPosT[t_nCandidates[0]] = candidates[key].getTime()
                t_nE[t_nCandidates[0]] = candidates[key].getEnergy()
                t_nEmax[t_nCandidates[0]] = candidates[key].getMaxCell()
                t_nNcell[t_nCandidates[0]] = candidates[key].getNcell()
                t_nNcell1MeV[t_nCandidates[0]] = candidates[key].getNcellCut(1.)
                t_nNcell3MeV[t_nCandidates[0]] = candidates[key].getNcellCut(3.)
                t_nTruePDG[t_nCandidates[0]] = candidates[key].getTruePDG()
                t_nIsPrimary[t_nCandidates[0]] = isPrimary
                t_nTrueKE[t_nCandidates[0]] = candidates[key].getTrueKE()
                t_nCandidates[0] += 1

                if t_nCandidates[0] == MAXCANDIDATES:
                    print "Event has more than maximum %d neutron candidates" % MAXCANDIDATES
                    break

            tree.Fill()

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="/pnfs/dune/persistent/users/marshalc/neutronSim/EDep")
    parser.add_option('--first_run', type=int, help='First run number', default=0)
    parser.add_option('--last_run', type=int, help='Last run number', default=0)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--geom',help='top volume of interactions', default="DetEnclosure")

    # python analyze --topdir /pnfs/dune/persistent/users/marshalc/neutronSim/EDep --first_run 0 --last_run 0 --geom DetEnclosure --outfile out.root

    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tree = ROOT.TTree( "tree","tree" )
    tree.Branch( "run", t_run, "run/I" )
    tree.Branch( "event", t_event, "event/I" )
    tree.Branch( "vtxX", t_vtxX, "vtxX/D" )
    tree.Branch( "vtxY", t_vtxY, "vtxY/D" )
    tree.Branch( "vtxZ", t_vtxZ, "vtxZ/D" )
    tree.Branch( "vtxA", t_vtxA, "vtxA/I" )
    tree.Branch( "nCandidates", t_nCandidates, "nCandidates/I" )
    tree.Branch( "nPosX", t_nPosX, "nPosX[nCandidates]/D" )
    tree.Branch( "nPosY", t_nPosY, "nPosY[nCandidates]/D" )
    tree.Branch( "nPosZ", t_nPosZ, "nPosZ[nCandidates]/D" )
    tree.Branch( "nPosT", t_nPosT, "nPosT[nCandidates]/D" )
    tree.Branch( "nE", t_nE, "nE[nCandidates]/D" )
    tree.Branch( "nEmax", t_nEmax, "nEmax[nCandidates]/D" )
    tree.Branch( "nNcell", t_nNcell, "nNcell[nCandidates]/I" )
    tree.Branch( "nNcell1MeV", t_nNcell1MeV, "nNcell1MeV[nCandidates]/I" )
    tree.Branch( "nNcell3MeV", t_nNcell3MeV, "nNcell3MeV[nCandidates]/I" )
    tree.Branch( "nTruePDG", t_nTruePDG, "nTruePDG[nCandidates]/I" )
    tree.Branch( "nIsPrimary", t_nIsPrimary, "nIsPrimary[nCandidates]/I" )
    tree.Branch( "nTrueKE", t_nTrueKE, "nTrueKE[nCandidates]/D" )

    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino" if not args.rhc else "antineutrino"
    horn = "FHC" if not args.rhc else "RHC"

    for run in range( args.first_run, args.last_run+1 ):
        fname = "%s/%s/%s/%s.%d.edep.root" % (args.topdir, horn, args.geom, neutrino, run)
        if not os.access( fname, os.R_OK ):
            print "Can't access file: %s" % fname
            continue
        tf = ROOT.TFile( fname )
        if tf.TestBit(ROOT.TFile.kRecovered):
            print "File is crap: %s" % fname
            continue

        if tgeo is None:
            tf.MakeProject("EDepSimEvents","*","RECREATE++")
            tgeo = tf.Get( "EDepSimGeometry" )

        print "Adding: %s" % fname
        events.Add( fname )

    loop( events, tgeo, tree )
    fout.cd()
    tree.Write()
