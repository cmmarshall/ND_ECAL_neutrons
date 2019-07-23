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
#import numpy as np

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


def intersection(list1, list2):
    return list(set(list1)&set(list2))


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

def Fill_Vertex_Info(Plot_Dict, vertex):
    PosX = vertex.Position.X()/10.; PosY = vertex.Position.Y()/10.; PosZ = vertex.Position.Z()/10.
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

def GetLayer(hit):
    node = tgeo.FindNode( hit.Start.X(), hit.Start.Y(), hit.Start.Z())
    Cell_Name = node.GetName()
    decomp = [int(s) for s in Cell_Name.split('_') if s.isdigit()]
    return decomp[0]

def MergeClusters(candidates, Cluster_Threshold):
    merge_dict = {}
    for i in range(len(candidates)):
        #print('Candidate %d of %d'%(i, len(candidates)))
        clusteri = candidates[i]#; posi = clusteri.getPos()
        hitsi = clusteri.getHits()
        merge_dict[i] = []
        for j in range(len(candidates)):
            clusterj = candidates[j];# posj = clusterj.getPos()
            hitsj = clusterj.getHits()
            for hiti in hitsi:
                for hitj in hitsj:
                    posi = hiti.getPos(); posj = hitj.getPos()
                    diff = posi - posj; distance = sqrt(diff.Dot(diff))
                    if distance < Cluster_Threshold:
                        merge_dict[i].append(j)
                        break; break;

    #STEP 2
    reduced_merges = []; reduced_keys = []
    for key1 in merge_dict:
        if key1 not in reduced_keys:
            reduced_merges.append(merge_dict[key1])
            reduced_keys.append(key1)
            for key2 in merge_dict:
                if len(intersection(reduced_merges[len(reduced_merges)-1], merge_dict[key2])) > 0:
                    reduced_keys.append(key2)
                    reduced_merges[len(reduced_merges)-1] = list( set(reduced_merges[len(reduced_merges)-1]) | set(merge_dict[key2]) )

    #STEP 3
    new_candidates = []
    for thing in reduced_merges:
        if len(thing) == 1:
            new_candidates.append(candidates[thing[0]])
        else:
            output_cluster = candidates[thing[0]]
            for i in range(1,len(thing)):
#               print('Length of Cluster:%d'%(len(candidates[thing[i]].getHits())))
                output_cluster+=candidates[thing[i]]
            new_candidates.append(output_cluster)
    return new_candidates




def loop( events, tgeo, tree, Cluster_Threshold = 10 ): # ** CHRIS: WHAT SHOULD I SET THE DEFAULT THRESHOLD TO **
    offset = [ 0., 305., 5. ]

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

#    NPlot_Dict = Initialize_Plot_Dicts(str1 = '1', str2 = 'True Neutron: ')
#    GPlot_Dict = Initialize_Plot_Dicts(str1 = '2', str2 = 'True Photon: ', includevtx = False )

    N = events.GetEntries()
    for ient in range(0, N):
        if ient % 1 == 0:
            print "Event %d of %d..." % (ient,N)
	    if ient > 1000:
	        break;
        events.GetEntry(ient)
        candidates = []
        for ivtx,vertex in enumerate(event.Primaries):
            setToBogus()
            reaction = vertex.Reaction #vertex.GetReaction()
            t_vtxX[0] = vertex.Position.X()/10. #vertex.GetPosition().X()/10. # cm
            t_vtxY[0] = vertex.Position.Y()/10. #vertex.GetPosition().Y()/10.
            t_vtxZ[0] = vertex.Position.Z()/10. #vertex.GetPosition().Z()/10.
            # atomic mass of target
            t_vtxA[0] = int((reaction.split(";")[1].split(":")[1])[6:9])

            ecal_hits = []
            for det in event.SegmentDetectors:
                if "ECal" in det.first: # there are several sub-detectors that make up the ECal
                    ecal_hits += det.second

            ecal_hits = sorted(ecal_hits, key = lambda hit: GetLayer(hit))

            for kk, hit in enumerate(ecal_hits):
                #print('hit %d of %d'%(kk, len(ecal_hits)))
                #print(tgeo.FindNode( hit.Start.X(), hit.Start.Y(), hit.Start.Z()))
                if hit.EnergyDeposit < 0.01: # skip tiny deposits, this cut needs to be tuned
                    continue
                hStart = ROOT.TVector3( hit.Start.X()/10., hit.Start.Y()/10., hit.Start.Z()/10. )
                hStop = ROOT.TVector3( hit.Stop.X()/10., hit.Stop.Y()/10., hit.Stop.Z()/10. )

                #--------------DETERMINE PARENT PARTICLE--------------------#
                parent = None; neutral_tid = -1
                for ii in range(len(hit.Contrib)):
                    tid = hit.Contrib[ii]
                    photon_tid = photonParent(event, tid); neutron_tid = neutronParent( event, tid )
                    if neutron_tid > 0:
                        parent = 'n'; neutral_tid = neutron_tid
                    else:
                        if photon_tid > 1:
                            parent = 'g'; neutral_tid = photon_tid
                        else:
                            parent = None; neutral_tid = photon_tid

                    #----------------------------------------------------------#
                    #-------------------\\\\\Clustering\\\\\-------------------#
                    #----------------------------------------------------------#
                    if parent is not None:
                        #Define Proposed Hit
                        node = tgeo.FindNode( hit.Start.X(), hit.Start.Y(), hit.Start.Z())
                        mom = event.Trajectories[neutral_tid].InitialMomentum
                        New_Hit = Hit(hStart, node.GetName(), hit.EnergyDeposit, hit.Start[3], parent, int(neutral_tid), mom.E() - mom.M())


                        #Figure out what clusters our hit is in
                        inCluster = False; inwhichClusters= []
                        for i, cluster in enumerate(candidates):
                            if New_Hit in cluster:
                                inwhichClusters.append(i)


                        #If hit belongs to multiple clusters merge the clusters
                        if len(inwhichClusters) > 1:
                            output_cluster = candidates[inwhichClusters[0]]
                            for i in range(1, len(inwhichClusters)):
                                index = inwhichClusters[i]
                                this_cluster = candidates[index]
                                output_cluster += this_cluster
                            output_cluster.addHit(New_Hit)

                            #Make sure to remove old clusters
                            inwhichClusters = sorted(inwhichClusters); pindex = 0
                            for i in range(len(inwhichClusters)):
                                cindex = inwhichClusters[i] - pindex
                                candidates.pop(cindex)
                                pindex+=1
                            #Append New Cluster
                            candidates.append(output_cluster)


                        #If hit only belongs to one cluster
                        elif len(inwhichClusters) == 1:
                            index = inwhichClusters[0]
                            candidates[index].addHit(New_Hit)


                        #If hit doesn't belong to any clusters
                        else:
                            c = NeutronCandidate()
                            c.addHit(New_Hit)
                            candidates.append(c)

                #---------------------------------------#
                #-------------MERGE CLUSTERS------------#
                #---------------------------------------#

            candidates = MergeClusters(candidates, Cluster_Threshold)
            


        GetPurityData(candidates, ient )
        Closest_Cluster_Distribution(candidates, ient)



        t_nCandidates[0] = 0
        CTrueN = []; CTrueG =[]
        for key in range(len(candidates)):
            #isPrimary = (event.Trajectories[key].ParentId == -1)
            #for hit in candidate[key].getHits():

            candidates[key].UpdateTruth()
            if candidates[key].getTruePDG() == 2112:
                CTrueN.append(candidates[key])
            else:
                CTrueG.append(candidates[key])

            neutral_tids = list([thing.getNeutralTID() for thing in candidates[key].getHits()])
            #neutral_tids = list(np.array(candidate[key].getHits())[:,2])
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
        #Fill_Candidate_Info(NPlot_Dict, CTrueN, vertex)
        #Fill_Candidate_Info(GPlot_Dict, CTrueG, vertex)
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

    meta = ROOT.TTree( "potTree", "potTree" )
    t_pot = array( 'd', [0.] )
    meta.Branch( "pot", t_pot, "pot/D" )

    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino" if not args.rhc else "antineutrino"
    horn = "FHC" if not args.rhc else "RHC"

    t_pot[0] = 0.
    for run in range( args.first_run, args.last_run+1 ):
        fname = "%s/EDep/%s/%s/%s.%d.edepsim.root" % (args.topdir, horn, args.geom, neutrino, run)
        print('Does the file %s exist: %r'%(fname, os.path.exists(fname)))
        print('Do you have access to %s:%r'%(fname, os.access( fname, os.R_OK )))
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

    print('I am here')
    rhcarg = "--rhc" if args.rhc else ""
    cppopts = ['./getPOT', '--topdir', args.topdir, '--first', str(args.first_run), '--last', str(args.last_run), '--geom', args.geom, rhcarg]
    sp = subprocess.Popen(cppopts, stdout=subprocess.PIPE, stderr=None)
    t_pot[0] = float(sp.communicate()[0])
    meta.Fill()

    loop( events, tgeo, tree )
    fout.cd()
    tree.Write()
    meta.Write()
