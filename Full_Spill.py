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

def Cone_Reject(cans, vtx):
    candidates = sorted(cans, key = lambda cluster: dist(vtx, cluster.getPos()))
    output = [] # candidate list after cone rejection
    rejec_clusters = []
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
    h_dict['Reco'] = ROOT.TH1D('Reco', 'Reconstructed Neutron Energy Distribution; Reconstructed Neutron Energy (MeV);', 70, 0., 700.)
    h_dict['FReco'] = ROOT.TH1D('Fractional_Reco', 'Fractional Reconstructed Energy Residual; Fractional Residual;', 100, -1., 1.)
    h_dict['Eff'] =  ROOT.TH1D('Efficiency', 'Reconstruction Efficiency; Efficiency;', 70, 0., 700.)

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

        # cone-filter the neutron candidates
        gas_candidates = Cone_Reject(GasEvent.candidates, vtx)
        continue

        # pick a time from a flat distribution from 0, 10 microseconds
        t0 = random.uniform(0.,10000.) # ns
        for cand in GasEvent.candidates:
            cand.nPosT += t0
            cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

        reco_t0 = random.normalvariate(t0, dt)

        # simulate background events from the start of the spill to 200 ns after the signal
        # de-excitations can produce n-like signals microseconds after the interaction, but we don't 
        # care about stuff that happens long after the search window
        bkg_end = t0 + 200.

        # add rock background
        n_rock = irock + rando.Poisson(muRock*bkg_end/10000.)
        while irock < n_rock:
            RockEvent = Event(RTree, irock)

            # pick a time
            trock = random.uniform(0.,bkg_end) # ns
            for cand in RockEvent.candidates:
                cand.nPosT += trock
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

                np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
                reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
                if reco_KE is not None: # physical neutron
                    print "Rock background neutron of %1.1f MeV" % reco_KE

            irock += 1
            if irock == Nrock: # we've run out of rock events; loop again
                print "We've run out of rock events; resetting"
                n_rock -= Nrock
                irock = 0

        # add hall background
        n_hall = ihall + rando.Poisson(muHall*bkg_end/10000.)
        while ihall < n_hall:
            HallEvent = Event(GTree, ihall)

            # pick a time
            thall = random.uniform(0.,bkg_end) # ns
            for cand in HallEvent.candidates:
                cand.nPosT += thall
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty


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

    plots = {}
    Initialize_Plot_Dict(plots)

    print "mean events rock %1.1f hall %1.1f" % (muRock, muHall)

    loop( GTree, RTree, HTree, thresh, dt, muRock, muHall, plots )


