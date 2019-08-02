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

    # cone rejection
    i = 0
    while i < len(candidates):
        vecti = candidates[i].getPos() - vtx
        j = i+1
        while j < len(candidates):
            vectj = candidates[j].getPos() - vtx
            angle = vecti.Angle(vectj)
            if angle <= pi/6:
                candidates.pop(j)
            else:
                j += 1
        i += 1
    return candidates

def GetRecoE(np,vtx,tof):
    c = 29.9792 # cm/ns
    mn = 939.565 # MeV/c2
    beta = dist(np, vtx)/(c*tof)
    if beta < 0. or beta >= 1.:
        return None
    gamma = 1/sqrt(1 - beta**2);
    return mn*(gamma - 1)

def Initialize_Plot_Dict(h_dict):
    h_dict["reco"] = {}
    h_dict["eres"] = ROOT.TH2D( "eres", ";True neutron KE (MeV);Fractional energy resolution", 70, 0., 700., 100, -1., 1. )
    h_dict["eff"] = ROOT.TH1D( "eff", ";True neutron KE (MeV)", 70, 0., 700. )
    h_dict["denom"] = ROOT.TH1D( "denom", ";True neutron KE (MeV)", 70, 0., 700. )

    cats = ["signal", "gasbkg", "hall", "rock"]
    for c in cats:
        h_dict["reco"][c] = ROOT.TH1D( "reco_%s" % c, ";Reco neutron KE (MeV)", 70, 0., 700. )

def loop(GTree, RTree, HTree, Emin, dt, muRock, muHall, Plot_dict):

    rando = ROOT.TRandom3(12345)

    # save indices for rock, hall
    irock = 0
    ihall = 0

    # Loop over GAr interactions
    Ngas = GTree.GetEntries()
    Nrock = RTree.GetEntries()
    Nhall = HTree.GetEntries()
    for igas in range(Ngas):
        GasEvent = Event(GTree, igas)

        vtx = ROOT.TVector3(GasEvent.vtxX, GasEvent.vtxY, GasEvent.vtxZ)
        # impose a fiducial volume cut; for simplicity
        R = sqrt((vtx.y()+217.)**2 + (vtx.z()-585.)**2) # YZ circle centered at (-217, 585)
        if R > 200. or abs(vtx.x()) > 200.: continue # x is centered -250, 250

        # Get true neutrons
        tid_KE = {}
        for pn in GasEvent.primaries:
            tid_KE[pn.pTID] = pn.pKE
            h_dict["denom"].Fill( pn.pKE )

        # cone-filter the neutron candidates
        gas_candidates = neutron_cut(GasEvent.candidates, vtx)

        reco_tids = []
        # pick a time from a flat distribution from 0, 10 microseconds
        t0 = random.uniform(0.,10000.) # ns
        reco_t0 = random.normalvariate(t0, dt)
        for cand in gas_candidates:
            cand.nPosT += t0
            cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

            np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
            reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
            if reco_KE is not None: # physical neutron with 0 < beta < 1
                tid = cand.nParTID
                if tid in tid_KE: # signal neutron!!
                    trueKE = tid_KE[tid]
                    if trueKE < 1.E-6: continue # is this possible? it shouldn't be
                    if tid in reco_tids:
                        print "Duplicate reconstruction!!"
                        h_dict["reco"]["gasbkg"].Fill( reco_KE )
                    else:
                        h_dict["reco"]["signal"].Fill( reco_KE )
                        h_dict["eff"].Fill( trueKE ) # eff numerator fill with true KE
                        h_dict["eres"].Fill( trueKE, (reco_KE-trueKE)/trueKE )
                    reco_tids.append(tid)
                else: # not true primary neutron, i.e. gas background
                    h_dict["reco"]["gasbkg"].Fill( reco_KE )

        # simulate background events from the start of the spill to 200 ns after the signal
        # de-excitations can produce n-like signals microseconds after the interaction, but we don't 
        # care about stuff that happens long after the search window
        bkg_end = 10000. #t0 + 200.

        # add rock background
        n_rock = irock + rando.Poisson(muRock*bkg_end/10000.)
        while irock < n_rock:
            RockEvent = Event(RTree, irock)

            rock_candidates = neutron_cut(RockEvent.candidates, vtx)

            # pick a time
            trock = random.uniform(0.,bkg_end) # ns
            for cand in rock_candidates:
                cand.nPosT += trock
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

                np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
                reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
                if reco_KE is not None: # physical neutron
                    h_dict["reco"]["rock"].Fill( reco_KE )

            irock += 1
            if irock == Nrock: # we've run out of rock events; loop again
                print "We've run out of rock events; resetting"
                n_rock -= Nrock
                irock = 0

        # add hall background
        n_hall = ihall + rando.Poisson(muHall*bkg_end/10000.)
        while ihall < n_hall:
            HallEvent = Event(GTree, ihall)

            hall_candidates = neutron_cut(HallEvent.candidates, vtx)

            # pick a time
            thall = random.uniform(0.,bkg_end) # ns
            for cand in hall_candidates:
                cand.nPosT += thall
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

                np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
                reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
                if reco_KE is not None: # physical neutron
                    h_dict["reco"]["hall"].Fill( reco_KE )


            ihall += 1
            if ihall == Nhall: # we've run out of hall events; loop again
                print "We've run out of hall events; resetting"
                n_hall -= Nhall
                ihall = 0


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--outfile', help='Output File Name', default='FSout.root')

    parser.add_option('--topdir', help='Directory containing Input', default='/dune/app/users/marshalc/ND_neutron/ND_ECAL_neutrons')
    parser.add_option('--Gfile_str', help='Input Files for Gas Argon File String', default='outgas.root')
    parser.add_option('--Rfile_str', help='Input Files for Rock File String', default='outrock.root')
    parser.add_option('--Hfile_str', help='Input Files for Hall File String', default='outenc.root')
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

    Gfile = ROOT.TFile('%s/%s' % (topdir, Gfile_str))
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
    for entry in rock_meta: rock_pot = entry.pot # only one entry
    muRock = RTree.GetEntries() * spill_pot / rock_pot

    enc_meta = Hfile.Get('potTree')
    enc_pot = 0.
    for entry in enc_meta: enc_pot = entry.pot # only one entry
    muHall = HTree.GetEntries() * spill_pot / enc_pot

    h_dict = {}
    Initialize_Plot_Dict(h_dict)

    print "mean events rock %1.1f hall %1.1f" % (muRock, muHall)

    loop( GTree, RTree, HTree, thresh, dt, muRock, muHall, h_dict )

    tfout = ROOT.TFile( args.outfile, "RECREATE" )
    for key in h_dict:
        if key != "reco":
            h_dict[key].Write()
        else:
            for key2 in h_dict[key]:
                h_dict[key][key2].Write()


