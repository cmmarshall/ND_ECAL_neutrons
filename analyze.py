#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array
from math import sqrt
import subprocess
import random
from Cluster_Module import *

MAXCANDIDATES = 1000
MAXNEUTRONS = 100

# Output TTree array variables
t_run = array( 'i', [0] )
t_event = array( 'i', [0] )
t_vtxX = array( 'd', [0.] )
t_vtxY = array( 'd', [0.] )
t_vtxZ = array( 'd', [0.] )
t_vtxA = array( 'i', [0] )
t_ECAL_visE = array( 'd', [0.] )
t_ECAL_vetoT = array( 'd', [0.] )
t_ECAL_vetoX = array( 'd', [0.] )
t_ECAL_vetoY = array( 'd', [0.] )
t_ECAL_vetoZ = array( 'd', [0.] )
t_nCandidates = array( 'i', [0] )
t_nPosX = array( 'd', MAXCANDIDATES*[0.] )
t_nPosY = array( 'd', MAXCANDIDATES*[0.] )
t_nPosZ = array( 'd', MAXCANDIDATES*[0.] )
t_nPosT = array( 'd', MAXCANDIDATES*[0.] )
t_nSigmaX = array( 'd', MAXCANDIDATES*[0.] )
t_nSigmaY = array( 'd', MAXCANDIDATES*[0.] )
t_nSigmaZ = array( 'd', MAXCANDIDATES*[0.] )
t_nE = array( 'd', MAXCANDIDATES*[0.] )
t_nIso = array( 'd', MAXCANDIDATES*[0.] )
t_nNcell = array( 'i', MAXCANDIDATES*[0] )
t_nMaxCell = array( 'd', MAXCANDIDATES*[0.] )
t_nTruePDG = array( 'i', MAXCANDIDATES*[0] )
t_nTrueKE = array( 'd', MAXCANDIDATES*[0.] )
t_nParTID = array('i', MAXCANDIDATES*[0])
t_nPrimTID = array('i', MAXCANDIDATES*[0])
t_pParticles = array( 'i', [0] )
t_pTID = array( 'i', MAXNEUTRONS*[0] )
t_pPDG = array( 'i', MAXNEUTRONS*[0] )
t_pKE = array( 'd', MAXNEUTRONS*[0.] )
t_pExitX = array( 'd', MAXNEUTRONS*[0.] )
t_pExitY = array( 'd', MAXNEUTRONS*[0.] )
t_pExitZ = array( 'd', MAXNEUTRONS*[0.] )
t_pExitdX = array( 'd', MAXNEUTRONS*[0.] )
t_pExitdY = array( 'd', MAXNEUTRONS*[0.] )
t_pExitdZ = array( 'd', MAXNEUTRONS*[0.] )

def setToBogus():
    t_run[0] = 0
    t_event[0] = 0
    t_vtxX[0] = 0.
    t_vtxY[0] = 0.
    t_vtxZ[0] = 0.
    t_vtxA[0] = 0
    t_ECAL_visE[0] = 0.
    t_ECAL_vetoT[0] = 0.
    t_ECAL_vetoX[0] = 0.
    t_ECAL_vetoY[0] = 0.
    t_ECAL_vetoZ[0] = 0.
    t_nCandidates[0] = 0
    t_pParticles[0] = 0

def isCharged(pdg):
    if abs(pdg) in [11, 13, 211, 321, 2212]:
        return True
    return False

def primaryParent( event, tid ):
    otid = tid
    parentid = event.Trajectories[otid].ParentId
    while parentid != -1: # primaries have parent ID -1
        otid = parentid
        parentid = event.Trajectories[parentid].ParentId
    return otid

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

def loop( events, tgeo, tree, cluster_gap = 10 ):

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()
    for ient in range(0, N):
        events.GetEntry(ient)

        info = ROOT.ProcInfo_t()
        ROOT.gSystem.GetProcInfo(info)
        vm = info.fMemVirtual
        rm = info.fMemResident

        if ient % 100 == 0:
            print "Event %d of %d, VM = %1.2f MB, RM = %1.2f MB" % (ient,N,vm/1000.,rm/1000.)

        t_event[0] = ient
        clusters = []
        tid_idx = {}

        setToBogus()

        for ivtx,vertex in enumerate(event.Primaries):

            reaction = vertex.Reaction #vertex.GetReaction()
            t_vtxX[0] = vertex.Position.X()/10. #vertex.GetPosition().X()/10. # cm
            t_vtxY[0] = vertex.Position.Y()/10. #vertex.GetPosition().Y()/10.
            t_vtxZ[0] = vertex.Position.Z()/10. #vertex.GetPosition().Z()/10.
            # atomic mass of target
            t_vtxA[0] = int((reaction.split(";")[1].split(":")[1])[6:9])

            # Save true primary neutrons for efficiency determination
            for ipart,particle in enumerate(vertex.Particles):
                mom = particle.Momentum
                pdg = particle.PDGCode
                tid = particle.TrackId
                t_pTID[t_pParticles[0]] = tid
                t_pPDG[t_pParticles[0]] = pdg
                t_pKE[t_pParticles[0]] = mom.E() - mom.M()
                t_pExitX[t_pParticles[0]] = -9999.9
                t_pExitY[t_pParticles[0]] = -9999.9
                t_pExitZ[t_pParticles[0]] = -9999.9
                t_pExitdX[t_pParticles[0]] = -9999.9
                t_pExitdY[t_pParticles[0]] = -9999.9
                t_pExitdZ[t_pParticles[0]] = -9999.9

                if isCharged(pdg): tid_idx[tid] = t_pParticles[0]

                t_pParticles[0] += 1

            for tid in tid_idx:
                pts = event.Trajectories[tid].Points
                lastPt = None
                for pt in pts:
                    node = tgeo.FindNode( pt.Position.X(), pt.Position.Y(), pt.Position.Z() )
                    volName = node.GetName()
                    if "TPC" in volName: # still in Gas TPC
                        lastPt = pt.Position.Vect()
                    else: # we just poked out of the Gas TPC
                        if lastPt is None: # vertex isn't in the gas TPC, this variable is irrelevant for hall/rock background
                            break
                        idx = tid_idx[tid]
                        thisPt = pt.Position.Vect()
                        t_pExitX[idx] = lastPt.x()
                        t_pExitY[idx] = lastPt.y()
                        t_pExitZ[idx] = lastPt.z()
                        t_pExitdX[idx] = (thisPt-lastPt).Unit().x()
                        t_pExitdY[idx] = (thisPt-lastPt).Unit().y()
                        t_pExitdZ[idx] = (thisPt-lastPt).Unit().z()

            ecal_hits = []
            for det in event.SegmentDetectors:
                if "ECal" in det.first: # there are several sub-detectors that make up the ECal
                    ecal_hits += det.second

            t_ECAL_visE[0] = 0.
            t_ECAL_vetoT[0] = 9999999.

            charged_hits = []

            # Loop over edep-sim "hits", call it edep to avoid confusin with Hit class
            for kk, edep in enumerate(ecal_hits):
                if edep.EnergyDeposit < 0.01: # skip tiny deposits, this cut needs to be tuned
                    continue
                hStart = ROOT.TVector3( edep.Start.X()/10., edep.Start.Y()/10., edep.Start.Z()/10. )

                #--------------DETERMINE PARENT PARTICLE--------------------#
                neutral_tid = -1
                primary_tid = -1
                edep_tid = 0
                for ii in range(len(edep.Contrib)):
                    tid = edep.Contrib[ii]
                    photon_tid = photonParent(event, tid)
                    neutron_tid = neutronParent(event, tid)
                    if neutron_tid > 0:
                        neutral_tid = neutron_tid
                        edep_tid = tid
                        primary_tid = primaryParent(event, tid)
                    elif neutral_tid <= 0:
                        neutral_tid = photon_tid
                        edep_tid = tid
                        primary_tid = primaryParent(event, tid)

                #----------------------------------------------------------#
                #-------------------\\\\\Clustering\\\\\-------------------#
                #----------------------------------------------------------#
                if neutral_tid > 0: # ensure neutral parent

                    node = tgeo.FindNode( edep.Start.X(), edep.Start.Y(), edep.Start.Z())
                    this_hit = Hit(hStart, node.GetName(), edep.EnergyDeposit, edep.Start.T(), event.Trajectories[edep_tid].PDGCode, neutral_tid, primary_tid)

                    # Should this hit be added to any existing clusters?
                    # It is possible that this hit can be in multiple clusters, in which case they should be merged
                    inWhichClusters = [] # cluster indices of clusters that this hit matches
                    for i,cluster in enumerate(clusters):
                        if cluster.isInCluster(this_hit):
                            inWhichClusters.append(i)

                    # If this hit belongs in multiple clusters, they all need to be merged
                    # order of hits in clusters is totally irrelevant, so just merge them into the first one
                    if len(inWhichClusters) > 1:
                        idx = inWhichClusters[0]
                        for i in range(len(inWhichClusters)-1, 0, -1):
                            clusters[idx].mergeCluster( clusters[inWhichClusters[i]] )
                            clusters.pop(inWhichClusters[i]) # reverse order loop assures indices don't get messed up
                        clusters[idx].addHit(this_hit) # finally, add the hit that caused the merger
                    elif len(inWhichClusters) == 1:
                        # hit is in only one cluster, so just add the hit
                        idx = inWhichClusters[0]
                        clusters[idx].addHit(this_hit)
                    elif len(inWhichClusters) == 0:
                        # make a new cluster with this hit
                        this_cluster = Cluster(cluster_gap)
                        this_cluster.addHit(this_hit)
                        clusters.append(this_cluster)
                else: # charged parent
                    t_ECAL_visE[0] += edep.EnergyDeposit
                    charged_hits.append( ROOT.TVector3(edep.Start.X()/10., edep.Start.Y()/10., edep.Start.Z()/10.) )
                    reco_t = random.normalvariate(edep.Start.T(), 0.7)
                    if reco_t < t_ECAL_vetoT[0]:
                        t_ECAL_vetoT[0] = reco_t
                        t_ECAL_vetoX[0] = edep.Start.X()/10.
                        t_ECAL_vetoY[0] = edep.Start.Y()/10.
                        t_ECAL_vetoZ[0] = edep.Start.Z()/10.

        # Fill the output ntuple
        t_nCandidates[0] = 0
        for cluster in clusters:

            cluster.CalcStuff()

            cpos = cluster.getCentroid()
            t_nIso[t_nCandidates[0]] = 999999.9
            for hit in charged_hits:
                if (cpos-hit).Mag() < t_nIso[t_nCandidates[0]]:
                    t_nIso[t_nCandidates[0]] = (cpos-hit).Mag()

            t_nPosX[t_nCandidates[0]] = cluster.getCentroid().x()
            t_nPosY[t_nCandidates[0]] = cluster.getCentroid().y()
            t_nPosZ[t_nCandidates[0]] = cluster.getCentroid().z()
            t_nSigmaX[t_nCandidates[0]] = cluster.getSigmas()[0]
            t_nSigmaY[t_nCandidates[0]] = cluster.getSigmas()[1]
            t_nSigmaZ[t_nCandidates[0]] = cluster.getSigmas()[2]
            t_nPosT[t_nCandidates[0]] = cluster.getTime()
            t_nE[t_nCandidates[0]] = cluster.getEnergy()
            t_nTruePDG[t_nCandidates[0]] = cluster.getTruePDG()
            t_nParTID[t_nCandidates[0]] = cluster.getTrueParent()
            t_nPrimTID[t_nCandidates[0]] = cluster.getPrimary()
            t_nMaxCell[t_nCandidates[0]] = cluster.getMaxCell()
            t_nNcell[t_nCandidates[0]] = cluster.getNcell(0.5)
            parent = event.Trajectories[cluster.getTrueParent()]
            t_nTrueKE[t_nCandidates[0]] = parent.InitialMomentum.E() - parent.InitialMomentum.M()

            #vtx = ROOT.TVector3( t_vtxX[0], t_vtxY[0], t_vtxZ[0] )
            #dt = cluster.getTime()
            #dx = (cluster.getCentroid() - vtx).Mag()
            #beta = dx/(29.9792*dt)
            #recoe = -1.
            #if beta > 0. and beta < 1.:
            #    recoe = 939.565*( 1./sqrt(1.-beta**2) - 1. )

            #print "Neutron cluster due to %d, dt %1.1f dx %1.1f, reco E %1.1f, primary %d of %1.1f MeV" % (cluster.getTruePDG(), 
            #                                                                                               dt, 
            #                                                                                               dx, 
            #                                                                                               recoe, 
            #                                                                                               tid_pdg_ke[cluster.getPrimary()][0],
            #                                                                                               tid_pdg_ke[cluster.getPrimary()][1])

            t_nCandidates[0] += 1

            if t_nCandidates[0] == MAXCANDIDATES:
                print "Event has more than maximum %d neutron candidates" % MAXCANDIDATES
                break

        tree.Fill()




if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="/pnfs/dune/persistent/users/marshalc/neutronSim")
    parser.add_option('--first_run', type=int, help='First run number', default=1001)
    parser.add_option('--last_run', type=int, help='Last run number', default=1001)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--geom',help='top volume of interactions', default="GArTPC")
    parser.add_option('--cgap',help='Set Cluster Gap', default=5.)
    parser.add_option('--grid',action='store_true', help='Grid mode')
    # python analyze --topdir /pnfs/dune/persistent/users/marshalc/neutronSim/EDep --first_run 0 --last_run 0 --geom DetEnclosure --outfile out.root

    (args, dummy) = parser.parse_args()

    rhcarg = "--rhc" if args.rhc else ""
    gridarg = "--grid" if args.grid else ""
    cppopts = ['./getPOT', '--topdir', args.topdir, '--first', str(args.first_run), '--last', str(args.last_run), '--geom', args.geom, rhcarg, gridarg]
    sp = subprocess.Popen(cppopts, stdout=subprocess.PIPE, stderr=None)
    the_POT = float(sp.communicate()[0])
   

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tree = ROOT.TTree( "tree","tree" )
    tree.Branch( "run", t_run, "run/I" )
    tree.Branch( "event", t_event, "event/I" )
    tree.Branch( "vtxX", t_vtxX, "vtxX/D" )
    tree.Branch( "vtxY", t_vtxY, "vtxY/D" )
    tree.Branch( "vtxZ", t_vtxZ, "vtxZ/D" )
    tree.Branch( "vtxA", t_vtxA, "vtxA/I" )
    tree.Branch( "ECAL_visE", t_ECAL_visE, "ECAL_visE/D" )
    tree.Branch( "ECAL_vetoT", t_ECAL_vetoT, "ECAL_vetoT/D" )
    tree.Branch( "ECAL_vetoX", t_ECAL_vetoX, "ECAL_vetoX/D" )
    tree.Branch( "ECAL_vetoY", t_ECAL_vetoY, "ECAL_vetoY/D" )
    tree.Branch( "ECAL_vetoZ", t_ECAL_vetoZ, "ECAL_vetoZ/D" )
    tree.Branch( "nCandidates", t_nCandidates, "nCandidates/I" )
    tree.Branch( "nPosX", t_nPosX, "nPosX[nCandidates]/D" )
    tree.Branch( "nPosY", t_nPosY, "nPosY[nCandidates]/D" )
    tree.Branch( "nPosZ", t_nPosZ, "nPosZ[nCandidates]/D" )
    tree.Branch( "nPosT", t_nPosT, "nPosT[nCandidates]/D" )
    tree.Branch( "nSigmaX", t_nSigmaX, "nSigmaX[nCandidates]/D" )
    tree.Branch( "nSigmaY", t_nSigmaY, "nSigmaY[nCandidates]/D" )
    tree.Branch( "nSigmaZ", t_nSigmaZ, "nSigmaZ[nCandidates]/D" )
    tree.Branch( "nE", t_nE, "nE[nCandidates]/D" )
    tree.Branch( "nIso", t_nIso, "nIso[nCandidates]/D" )
    tree.Branch( "nNcell", t_nNcell, "nNcell[nCandidates]/I" )
    tree.Branch( "nMaxCell", t_nMaxCell, "nMaxCell[nCandidates]/D" )
    tree.Branch( "nTruePDG", t_nTruePDG, "nTruePDG[nCandidates]/I" )
    tree.Branch( "nTrueKE", t_nTrueKE, "nTrueKE[nCandidates]/D" )
    tree.Branch( "nParTID", t_nParTID, "nParTID[nCandidates]/I")
    tree.Branch( "nPrimTID", t_nPrimTID, "nPrimTID[nCandidates]/I")
    tree.Branch( "pParticles", t_pParticles, "pParticles/I" )
    tree.Branch( "pTID", t_pTID, "pTID[pParticles]/I")
    tree.Branch( "pPDG", t_pPDG, "pPDG[pParticles]/I")
    tree.Branch( "pKE", t_pKE, "pKE[pParticles]/D")
    tree.Branch( "pExitX", t_pExitX, "pExitX[pParticles]/D")
    tree.Branch( "pExitY", t_pExitY, "pExitY[pParticles]/D")
    tree.Branch( "pExitZ", t_pExitZ, "pExitZ[pParticles]/D")
    tree.Branch( "pExitdX", t_pExitdX, "pExitdX[pParticles]/D")
    tree.Branch( "pExitdY", t_pExitdY, "pExitdY[pParticles]/D")
    tree.Branch( "pExitdZ", t_pExitdZ, "pExitdZ[pParticles]/D")

    meta = ROOT.TTree( "potTree", "potTree" )
    t_pot = array( 'd', [0.] )
    meta.Branch( "pot", t_pot, "pot/D" )

    tgeo = None

    neutrino = "neutrino" if not args.rhc else "antineutrino"
    horn = "FHC" if not args.rhc else "RHC"

    t_pot[0] = 0.
    for run in range( args.first_run, args.last_run+1 ):
        fname = None
        if args.grid:
            fname = "%s.%d.edepsim.root" % (neutrino, run)
        else:
            fname = "%s/EDep/%s/%s/%s.%d.edepsim.root" % (args.topdir, horn, args.geom, neutrino, run)
        if not os.access( fname, os.R_OK ):
            print "Can't access file: %s" % fname
            continue
        tf = ROOT.TFile( fname )
        events = tf.Get( "EDepSimEvents" )

        if tgeo is None:
            tf.MakeProject("EDepSimEvents","*","RECREATE++")
            tgeo = tf.Get( "EDepSimGeometry" )

        print "Looping over: %s" % fname
        fout.cd()
        loop( events, tgeo, tree, cluster_gap=float(args.cgap))
        tf.Close()

    fout.cd()
    t_pot[0] = the_POT
    meta.Fill()

    tree.Write()
    meta.Write()
