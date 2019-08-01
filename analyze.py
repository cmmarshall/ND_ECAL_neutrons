#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array
from math import sqrt
import subprocess
from Cluster_Module import *
#import psutil
#import numpy as np

MAXCANDIDATES = 1000
MAXNEUTRONS = 100

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
t_nTruePDG = array( 'i', MAXCANDIDATES*[0] )
t_nTrueKE = array( 'd', MAXCANDIDATES*[0.] )
t_nParTID = array('i', MAXCANDIDATES*[0])
t_pNeutrons = array( 'i', [0] )
t_pTID = array( 'i', MAXNEUTRONS*[0] )
t_pKE = array( 'd', MAXNEUTRONS*[0.] )


def intersection(list1, list2):
    return list(set(list1)&set(list2))


#def GetLayer(hit):
#    node = tgeo.FindNode( hit.Start.X(), hit.Start.Y(), hit.Start.Z())
#    Cell_Name = node.GetName()
#    decomp = [int(s) for s in Cell_Name.split('_') if s.isdigit()]
#    return decomp[0]



def setToBogus():
    t_run[0] = 0
    t_event[0] = 0
    t_vtxX[0] = 0.
    t_vtxY[0] = 0.
    t_vtxZ[0] = 0.
    t_vtxA[0] = 0
    t_nCandidates[0] = 0
    t_pNeutrons[0] = 0

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

        if ient % 10 == 0:
            print "Event %d of %d, VM = %1.2f MB, RM = %1.2f MB" % (ient,N,vm/1000.,rm/1000.)

        for ivtx,vertex in enumerate(event.Primaries):
            clusters = []

            setToBogus()

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
                if pdg == 2112:
                    t_pTID[t_pNeutrons[0]] = tid
                    t_pKE[t_pNeutrons[0]] = mom.E() - mom.M()
                    t_pNeutrons[0] += 1


            ecal_hits = []
            for det in event.SegmentDetectors:
                if "ECal" in det.first: # there are several sub-detectors that make up the ECal
                    ecal_hits += det.second

            # Sort deposits by ECAL layer so that cluster merge logic occurs in predictable order?
#            ecal_hits = sorted(ecal_hits, key = lambda edep: GetLayer(edep))

            # Loop over edep-sim "hits", call it edep to avoid confusin with Hit class
            for kk, edep in enumerate(ecal_hits):
                if edep.EnergyDeposit < 0.01: # skip tiny deposits, this cut needs to be tuned
                    continue
                hStart = ROOT.TVector3( edep.Start.X()/10., edep.Start.Y()/10., edep.Start.Z()/10. )
                #hStop = ROOT.TVector3( edep.Stop.X()/10., edep.Stop.Y()/10., edep.Stop.Z()/10. )

                #--------------DETERMINE PARENT PARTICLE--------------------#
                neutral_tid = -1
                edep_tid = 0
                for ii in range(len(edep.Contrib)):
                    tid = edep.Contrib[ii]
                    photon_tid = photonParent(event, tid)
                    neutron_tid = neutronParent(event, tid)
                    if neutron_tid > 0:
                        neutral_tid = neutron_tid
                        edep_tid = tid
                    elif neutral_tid <= 0:
                        neutral_tid = photon_tid
                        edep_tid = tid

                #----------------------------------------------------------#
                #-------------------\\\\\Clustering\\\\\-------------------#
                #----------------------------------------------------------#
                if neutral_tid > 0: # ensure neutral parent

                    node = tgeo.FindNode( edep.Start.X(), edep.Start.Y(), edep.Start.Z())
                    this_hit = Hit(hStart, node.GetName(), edep.EnergyDeposit, edep.Start[3], event.Trajectories[edep_tid].PDGCode, neutral_tid)

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

        # Fill the output ntuple
        t_nCandidates[0] = 0
        for cluster in clusters:

            cluster.CalcStuff()

            t_nPosX[t_nCandidates[0]] = cluster.getCentroid().x()
            t_nPosY[t_nCandidates[0]] = cluster.getCentroid().y()
            t_nPosZ[t_nCandidates[0]] = cluster.getCentroid().z()
            t_nPosT[t_nCandidates[0]] = cluster.getTime()
            t_nE[t_nCandidates[0]] = cluster.getEnergy()
            t_nTruePDG[t_nCandidates[0]] = cluster.getTruePDG()
            t_nParTID[t_nCandidates[0]] = cluster.getTrueParent()
            parent = event.Trajectories[cluster.getTrueParent()]
            t_nTrueKE[t_nCandidates[0]] = parent.InitialMomentum.E() - parent.InitialMomentum.M()

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
    parser.add_option('--cgap',help='Set Cluster Gap', default=10.)
    # python analyze --topdir /pnfs/dune/persistent/users/marshalc/neutronSim/EDep --first_run 0 --last_run 0 --geom DetEnclosure --outfile out.root

    (args, dummy) = parser.parse_args()

    rhcarg = "--rhc" if args.rhc else ""
    cppopts = ['./getPOT', '--topdir', args.topdir, '--first', str(args.first_run), '--last', str(args.last_run), '--geom', args.geom, rhcarg]
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
    tree.Branch( "nCandidates", t_nCandidates, "nCandidates/I" )
    tree.Branch( "nPosX", t_nPosX, "nPosX[nCandidates]/D" )
    tree.Branch( "nPosY", t_nPosY, "nPosY[nCandidates]/D" )
    tree.Branch( "nPosZ", t_nPosZ, "nPosZ[nCandidates]/D" )
    tree.Branch( "nPosT", t_nPosT, "nPosT[nCandidates]/D" )
    tree.Branch( "nE", t_nE, "nE[nCandidates]/D" )
    tree.Branch( "nTruePDG", t_nTruePDG, "nTruePDG[nCandidates]/I" )
    tree.Branch( "nTrueKE", t_nTrueKE, "nTrueKE[nCandidates]/D" )
    tree.Branch( "nParTID", t_nParTID, "nParTID[nCandidates]/I")
    tree.Branch( "pNeutrons", t_pNeutrons, "pNeutrons/I" )
    tree.Branch( "pTID", t_pTID, "pTID[pNeutrons]/I")
    tree.Branch( "pKE", t_pKE, "pKE[pNeutrons]/D")

    meta = ROOT.TTree( "potTree", "potTree" )
    t_pot = array( 'd', [0.] )
    meta.Branch( "pot", t_pot, "pot/D" )

    tgeo = None

    neutrino = "neutrino" if not args.rhc else "antineutrino"
    horn = "FHC" if not args.rhc else "RHC"

    t_pot[0] = 0.
    for run in range( args.first_run, args.last_run+1 ):
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
