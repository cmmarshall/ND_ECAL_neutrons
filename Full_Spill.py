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
from NTuple import *

def dist(TVec1, TVec2):
    return (TVec1-TVec2).Mag()

# distance from point p to line specified by point a and unit vector n
def distFromLine( a, n, p ):
    return ((a-p) - n*((a-p).Dot(n))).Mag()

def neutron_cut(cans, vtx):

    # Apply 2 cuts. 1) neutron recoil cut, 2) cone cut
    candidates = sorted(cans, key = lambda cluster: dist(vtx, cluster.getPos()))

    i = 0
    while i < len(candidates):
        can = candidates[i]
        nCut = (can.nMaxCell > 10. or (can.nMaxCell>3. and can.nNcell <= 4))
        if not nCut:
            candidates.pop(i)
        else:
            i += 1

    # cylinder rejection, start with nearest candidate to vertex
    i = 0
    while i < len(candidates):
        a = candidates[i].getPos()
        n = (vtx-a).Unit()
        j = i+1
        while j < len(candidates):
            distance = distFromLine( a, n, candidates[j].getPos() )
            if distance < 20.:
                candidates.pop(j)
            else:
                j += 1
        i += 1
    return candidates

    # old cone
    #i = 0
    #while i < len(candidates):
    #    vecti = candidates[i].getPos() - vtx
    #    j = i+1
    #    while j < len(candidates):
    #        vectj = candidates[j].getPos() - vtx
    #        angle = vecti.Angle(vectj)
    #        if angle <= pi/6:
    #            candidates.pop(j)
    #        else:
    #            j += 1
    #    i += 1
    #return candidates

def GetRecoE(np,vtx,tof):
    c = 29.9792 # cm/ns
    mn = 939.565 # MeV/c2
    beta = dist(np, vtx)/(c*tof)
    if beta < 0. or beta >= 1.:
        return None
    gamma = 1/sqrt(1 - beta**2);
    return mn*(gamma - 1)

def inChargedCylinder( cylinders, cand ):
    p = cand.getPos()
    for cylinder in cylinders:
        a = cylinder[1]
        n = cylinder[2].Unit()
        d = distFromLine(a, n, p)

        #print "Candidate (%1.0f, %1.0f, %1.0f), exit point (%1.0f, %1.0f, %1.0f) dir (%3.3f, %3.3f, %3.3f) d = %1.1fcm" % (p.x(), p.y(), p.z(), a.x(), a.y(), a.z(), n.x(), n.y(), n.z(), d)

        if d < 20.: # neutron candidate is in line with charged track from gas TPC
            #print "Rejecting neutron %1.1f cm from axis of %d" % (d, cylinder[0])
            return True
    return False

def isVeto( veto_time_vtx, cand ):
    p = cand.getPos()
    nt = cand.nPosT
    for vtxt,vtx in veto_time_vtx:
        dt = nt - vtxt
        d = (p-vtx).Mag()

    return False
        

def Initialize_Plot_Dict(h_dict):
    h_dict["reco"] = {}
    h_dict["eres"] = ROOT.TH2D( "eres", ";True neutron KE (MeV);Fractional energy resolution", 70, 0., 700., 100, -1., 1. )
    h_dict["eff"] = ROOT.TH1D( "eff", ";True neutron KE (MeV)", 70, 0., 700. )
    h_dict["denom"] = ROOT.TH1D( "denom", ";True neutron KE (MeV)", 70, 0., 700. )

    cats = ["signal", "gasn", "gasg", "halln", "hallg", "rockn", "rockg"]
    cuts = ["none", "gamma", "iso", "veto"]
    for cat in cats:
        h_dict["reco"][cat] = {}
        for cut in cuts:
            h_dict["reco"][cat][cut] = ROOT.TH1D( "reco_%s_%s" % (cat,cut), ";Reco neutron KE (MeV)", 70, 0., 700. )

def loop(GTree, RTree, HTree, Emin, dt, muRock, muHall, Plot_dict, Ngas=None):

    rando = ROOT.TRandom3(12345)

    # save indices for rock, hall
    irock = 0
    ihall = 0

    # Loop over GAr interactions
    if Ngas is None:
        Ngas = GTree.GetEntries()
    Nrock = RTree.GetEntries()
    Nhall = HTree.GetEntries()
    for igas in range(Ngas):
        if igas % 1 == 0:
            print "Gas event %d of %d..." % (igas, Ngas)
        GasEvent = Event(GTree, igas)

        vtx = ROOT.TVector3(GasEvent.vtxX, GasEvent.vtxY, GasEvent.vtxZ)
        # impose a fiducial volume cut; for simplicity
        R = sqrt((vtx.y()+217.)**2 + (vtx.z()-585.)**2) # YZ circle centered at (-217, 585)
        if R > 200. or abs(vtx.x()) > 200.: 
            continue # x is centered -250, 250

        # Signal time
        t0 = random.uniform(0.,10000.) # ns
        reco_t0 = random.normalvariate(t0, dt)

        n_rock = rando.Poisson(muRock)
        n_hall = rando.Poisson(muHall)

        # generate potential veto vertex times
        veto_time_vtx = []
        rock_times = [ random.uniform(0., 10000.) for i in range(n_rock) ]
        hall_times = [ random.uniform(0., 10000.) for i in range(n_hall) ]

        # Generate vertex time stamps and add to the veto list, if applicable
        for i,t in enumerate(rock_times):
            # only events just prior to signal neutron candidates matter, i.e. 200 ns around is sufficient
            if abs(t-t0) < 300.:
                idx = i + irock
                if idx > Nrock: idx -= Nrock # reset the rock counter
                RockEvent = Event(Rtree, idx)
                if RockEvent.ECAL_visE > 10.:
                    veto_time_vtx.append( (RockEvent.ECAL_vetoT+t, RockEvent.ECAL_vetoP) )

        for i,t in enumerate(hall_times):
            # only events just prior to signal neutron candidates matter, i.e. 200 ns around is sufficient
            if abs(t-t0) < 300.:
                idx = i + irock
                if idx > Nhall: idx -= Nhall # reset the hall counter
                HallEvent = Event(Htree, idx)
                if HallEvent.ECAL_visE > 10.:
                    veto_time_vtx.append( (HallEvent.ECAL_vetoT+t, HallEvent.ECAL_vetoP) )

        # cone-filter the signal neutron candidates
        gas_candidates = neutron_cut(GasEvent.candidates, vtx)

        # Get true neutrons, and other primaries that might produce backgrounds
        tid_KE = {}
        cylinders = []
        for pn in GasEvent.primaries:
            if pn.pPDG == 2112:
                #print "True neutron KE %1.1f" % (pn.pKE)
                tid_KE[pn.pTID] = pn.pKE
                h_dict["denom"].Fill( pn.pKE )

            if pn.pPos is not None:
                cylinders.append( (pn.pPDG, pn.pPos, pn.pDir) )

        reco_tids = [] # save track IDs of reco neutrons to check for duplicates, i.e. 1 neutron --> 2 or more candidates

        for cand in gas_candidates:
            cand.nPosT += t0
            cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

            np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
            reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
            true_neutron = (cand.nTruePDG == 2212 or cand.nTruePDG == 2112)

            # series of cuts to remove background and clean sample

            # remove super-luminal neutrons, or neutrons so delayed w.r.t. neutrino that they have very low reco energy
            if reco_KE is None or reco_KE < 5.: continue

            # check the ECAL cylinders of entering charged particles, which are likely due to those particles
            passCyl = (not inChargedCylinder(cylinders, cand))

            # Veto the event if it could have come from background
            passVeto = (not isVeto(veto_time_vtx, cand))

            tid = cand.nPrimTID
            if tid in tid_KE: # signal neutron!!
                trueKE = tid_KE[tid]
                print "Signal neutron true %1.1f reco %1.1f" % (trueKE, reco_KE)
                #if trueKE < 1.E-6: continue # is this possible? it shouldn't be
                if tid in reco_tids: # duplicate neutron; count the second one as background
                    #print "   Whoops! This neutron was already reconstructed"
                    h_dict["reco"]["gasnbkg"].Fill( reco_KE )
                    if passCyl: h_dict["reco"]["passgas"].Fill( reco_KE )
                else: # Good signal neutron
                    h_dict["reco"]["signal"].Fill( reco_KE )
                    if passCyl:
                        h_dict["reco"]["passsignal"].Fill( reco_KE )
                        h_dict["eff"].Fill( trueKE ) # eff numerator fill with true KE
                        h_dict["eres"].Fill( trueKE, (reco_KE-trueKE)/trueKE )
                    #else: print "   Whoops! rejected signal by cylinder cut"
                    reco_tids.append(tid)
            else: # not true primary neutron, i.e. gas background
                if true_neutron: h_dict["reco"]["gasnbkg"].Fill( reco_KE )
                else:            h_dict["reco"]["gasgbkg"].Fill( reco_KE )
                if passCyl:
                    h_dict["reco"]["passgas"].Fill( reco_KE )
                else:
                    print "    rejected by cylinder!"


        # add rock background
        for trock in rock_times:
            RockEvent = Event(RTree, irock)

            rock_candidates = neutron_cut(RockEvent.candidates, vtx)

            for cand in rock_candidates:
                true_neutron = (cand.nTruePDG == 2212 or cand.nTruePDG == 2112)
                cand.nPosT += trock
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

                np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
                reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
                if reco_KE is not None: # physical neutron
                    if true_neutron: h_dict["reco"]["rockn"].Fill( reco_KE )
                    else:            h_dict["reco"]["rockg"].Fill( reco_KE )

            irock += 1
            if irock == Nrock: # we've run out of rock events; loop again
                print "We've run through %d rock events; resetting" % Nrock
                n_rock -= Nrock
                irock = 0

        # add hall background
        while ihall < n_hall:
            HallEvent = Event(GTree, ihall)

            hall_candidates = neutron_cut(HallEvent.candidates, vtx)

            # pick a time
            thall = random.uniform(0.,bkg_end) # ns

            # Veto around "connected" ECAL energy from charged ancestors
            veto_vtx = None
            if HallEvent.ECAL_visE > 10.:
                veto_vtx = HallEvent.ECAL_centroid

            for cand in hall_candidates:
                true_neutron = (cand.nTruePDG == 2212 or cand.nTruePDG == 2112)
                cand.nPosT += thall
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

                np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
                reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
                if reco_KE is not None: # physical neutron
                    if true_neutron: h_dict["reco"]["halln"].Fill( reco_KE )
                    else:            h_dict["reco"]["hallg"].Fill( reco_KE )


            ihall += 1
            if ihall == Nhall: # we've run out of hall events; loop again
                print "We've run through %d hall events; resetting" % Nhall
                n_hall -= Nhall
                ihall = 0


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--outfile', help='Output File Name', default='FSout.root')
    parser.add_option('--topdir', help='Directory containing Input', default='/pnfs/dune/persistent/users/marshalc/neutronSim/ntuple')
    parser.add_option('--Gfile_str', help='Input Files for Gas Argon File String', default='outGArTPCv3.root')
    parser.add_option('--Rfile_str', help='Input Files for Rock File String', default='outRock.root')
    parser.add_option('--Hfile_str', help='Input Files for Hall File String', default='outDetEnclosure.root')
    parser.add_option('--Emin', type=float, help='Energy Threshhold for Neutron Candidates', default= 5)
    parser.add_option('--dt', type=float, help='Time Resolution of our Detector', default=0.7)
    parser.add_option('--spill_pot', type=float, help='POT per spill', default=7.5E13)

    (args, dummy) = parser.parse_args()
    topdir = args.topdir
    Gfile_str = args.Gfile_str
    Rfile_str = args.Rfile_str
    Hfile_str = args.Hfile_str
    thresh = float(args.Emin )
    dt = float(args.dt)
    spill_pot = float(args.spill_pot)

    #Gfile = ROOT.TFile('%s/%s' % (topdir, Gfile_str))
    Gfile = ROOT.TFile("out.root")
    Rfile = ROOT.TFile('%s/%s' % (topdir, Rfile_str))
    Hfile = ROOT.TFile('%s/%s' % (topdir, Hfile_str))

    GTree =  Gfile.Get('tree')
    RTree = Rfile.Get('tree')
    HTree = Hfile.Get('tree')

    setBranches(GTree)
    setBranches(RTree)
    setBranches(HTree)

    rock_meta = Rfile.Get('potTree')
    rock_pot = 0.
    for entry in rock_meta: rock_pot += entry.pot
    muRock = RTree.GetEntries() * spill_pot / rock_pot

    enc_meta = Hfile.Get('potTree')
    enc_pot = 0.
    for entry in enc_meta: enc_pot += entry.pot
    muHall = HTree.GetEntries() * spill_pot / enc_pot

    h_dict = {}
    Initialize_Plot_Dict(h_dict)

    print "Got %1.3g POT rock and %1.3g POT hall" % (rock_pot, enc_pot)
    print "Mean events per spill (%1.3g POT): rock %1.1f hall %1.1f" % (spill_pot, muRock, muHall)

    loop( GTree, RTree, HTree, thresh, dt, muRock, muHall, h_dict )

    tfout = ROOT.TFile( args.outfile, "RECREATE" )
    for key in h_dict:
        if key != "reco":
            h_dict[key].Write()
        else:
            for key2 in h_dict[key]:
                h_dict[key][key2].Write()


