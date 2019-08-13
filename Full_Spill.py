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
import time
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

        if d < 50.: # neutron candidate is in line with charged track from gas TPC
            #print "Rejecting neutron %1.1f cm from axis of %d" % (d, cylinder[0])
            return True
    return False

def mind( cylinders, cand ):
    p = cand.getPos()
    min_d = 99999999.9
    for cylinder in cylinders:
        a = cylinder[1]
        n = cylinder[2].Unit()
        d = distFromLine(a, n, p)

        if d < min_d: # neutron candidate is in line with charged track from gas TPC
            min_d = d
    return min_d 

def isVeto( veto_time_vtx, cand ):
    p = cand.getPos()
    nt = cand.nPosT
    for idx,vtxt,vtx in veto_time_vtx:
        dt = nt - vtxt
        d = (p-vtx).Mag()

        if dt > -5. and dt < 20. and d < 200.:
            return True

    return False

def vetoParams( veto_time_vtx, cand ):
    p = cand.getPos()
    nt = cand.nPosT

    ret = []
    for idx,vtxt,vtx in veto_time_vtx:
        dt = nt - vtxt
        d = (p-vtx).Mag()
        ret.append( (idx, dt, d) )
    return ret

def Initialize_Plot_Dict(h_dict):
    h_dict["reco"] = {}
    h_dict["cyl_min"] = {}
    h_dict["iso"] = {}
    h_dict["iso_cyl"] = {}
    h_dict["dt_d"] = {}
    h_dict["d_dt_min"] = {}
    h_dict["eres"] = ROOT.TH2D( "eres", ";True neutron KE (MeV);Fractional energy resolution", 70, 0., 700., 100, -1., 1. )
    h_dict["eff"] = {}
    h_dict["denom"] = ROOT.TH1D( "denom", ";True neutron KE (MeV)", 70, 0., 700. )

    cats = ["signal", "duplicate", "gasn", "gasg", "halln", "hallg", "rockn", "rockg"]
    cuts = ["none", "gamma", "cyl", "veto", "iso"]
    for cat in cats:
        h_dict["reco"][cat] = {}
        for cut in cuts:
            h_dict["reco"][cat][cut] = ROOT.TH1D( "reco_%s_%s" % (cat,cut), ";Reco neutron KE (MeV)", 70, 0., 700. )

        h_dict["cyl_min"][cat] = ROOT.TH1D( "cyl_min_%s" % cat, ";Min dist to charged traj (cm)", 100, 0., 500. )
        h_dict["iso"][cat] = ROOT.TH1D( "iso_%s" % cat, ";Isolation (cm)", 100, 0., 500. )
        h_dict["iso_cyl"][cat] = ROOT.TH2D( "isocyl_%s" % cat, ";Isolation (cm);Cylinder (cm)", 50, 0., 500., 50, 0., 500. )

        h_dict["dt_d"][cat] = ROOT.TH2D( "dt_d_%s" % cat, ";#Deltat to ECAL activity (ns);Dist to ECAL activity (m)", 250, -150., 200., 100, 0., 10. )
        h_dict["d_dt_min"][cat] = ROOT.TH1D( "d_dt_min_%s" % cat, ";Dist to ECAL activity (m)", 100, 0., 10. )

    for cut in cuts:
        h_dict["eff"][cut] = ROOT.TH1D( "eff_%s" % cut, ";True neutron KE (MeV)", 70, 0., 700. )

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

    stime = time.time()
    for igas in range(Ngas):
        if igas and igas % 1000 == 0:
            elapsed = time.time() - stime
            frac = float(igas)/Ngas
            remaining = elapsed/frac - elapsed

            info = ROOT.ProcInfo_t()
            ROOT.gSystem.GetProcInfo(info)
            vm = info.fMemVirtual
            rm = info.fMemResident
            print "Spill %d of %d...VM %1.1fMB RM %1.1fMB, %1.0f minutes remaining" % (igas, Ngas, vm/1000., rm/1000., remaining/60.)

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
                RockEvent = Event(RTree, idx)
                if RockEvent.ECAL_visE > 10. and (RockEvent.ECAL_vetoT+t < reco_t0-5 or RockEvent.ECAL_vetoT+t > reco_t0 + 10):
                    veto_time_vtx.append( (idx, RockEvent.ECAL_vetoT+t, RockEvent.ECAL_vetoP) )

        for i,t in enumerate(hall_times):
            # only events just prior to signal neutron candidates matter, i.e. 200 ns around is sufficient
            if abs(t-t0) < 300.:
                idx = i + ihall
                if idx > Nhall: idx -= Nhall # reset the hall counter
                HallEvent = Event(HTree, idx)
                if HallEvent.ECAL_visE > 10. and (HallEvent.ECAL_vetoT+t < reco_t0-5 or HallEvent.ECAL_vetoT+t > reco_t0 + 10):
                    veto_time_vtx.append( (idx, HallEvent.ECAL_vetoT+t, HallEvent.ECAL_vetoP) )

        # cone-filter the signal neutron candidates
        gas_candidates = neutron_cut(GasEvent.candidates, vtx)

        # Get true neutrons, and other primaries that might produce backgrounds
        tid_KE = {}
        cylinders = []
        tid_pdg = {}
        for pn in GasEvent.primaries:
            tid_pdg[pn.pTID] = pn.pPDG
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

            cat = None

            tid = cand.nPrimTID
            if tid in tid_KE: # signal neutron!!
                trueKE = tid_KE[tid]
                if tid in reco_tids: # duplicate neutron; count the second one as background
                    cat = "duplicate"
                else: # Good signal neutron
                    cat = "signal"
                    h_dict["eff"]["gamma"].Fill( trueKE )
                    if passCyl: h_dict["eff"]["cyl"].Fill( trueKE )
                    if passCyl and passVeto: 
                        h_dict["eff"]["veto"].Fill( trueKE ) # eff numerator fill with true KE
                        h_dict["eres"].Fill( trueKE, (reco_KE-trueKE)/trueKE )
                    if passCyl and passVeto and cand.nIso>50.: h_dict["eff"]["iso"].Fill( trueKE )
                    reco_tids.append(tid)
            else: # not true primary neutron, i.e. gas background
                cat = "gasn" if true_neutron else "gasg"

            h_dict["cyl_min"][cat].Fill( mind(cylinders,cand) )
            h_dict["iso_cyl"][cat].Fill( cand.nIso, mind(cylinders,cand) )

            h_dict["reco"][cat]["gamma"].Fill( reco_KE )
            if passCyl: 
                h_dict["reco"][cat]["cyl"].Fill( reco_KE )
                dtd = vetoParams( veto_time_vtx, cand )
                dt_min = 9999999999.9
                d_min = None
                for idx, deltat, d in dtd:
                    h_dict["dt_d"][cat].Fill( deltat, d/100. )
                    if deltat > -5 and deltat < dt_min:
                        dt_min = deltat
                        d_min = d/100.
                if d_min is not None: h_dict["d_dt_min"][cat].Fill( d_min )
                else: h_dict["d_dt_min"][cat].Fill( 9.99 )
            if passCyl and passVeto: 
                h_dict["reco"][cat]["veto"].Fill( reco_KE )
                h_dict["iso"][cat].Fill( cand.nIso )
            if passCyl and passVeto and cand.nIso > 50.:
                h_dict["reco"][cat]["iso"].Fill( reco_KE )
                
        # add hall background
        for thall in hall_times:
            HallEvent = Event(HTree, ihall)

            hall_candidates = neutron_cut(HallEvent.candidates, vtx)

            for cand in hall_candidates:
                true_neutron = (cand.nTruePDG == 2212 or cand.nTruePDG == 2112)
                cand.nPosT += thall
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

                passCyl = (not inChargedCylinder(cylinders, cand))
    
                # Veto the event if it could have come from background
                passVeto = (not isVeto(veto_time_vtx, cand))

                np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
                reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
                if reco_KE is None or reco_KE < 5.: continue # unphysical neutron

                cat = "halln" if true_neutron else "hallg"

                h_dict["cyl_min"][cat].Fill( mind(cylinders,cand) )
                h_dict["iso_cyl"][cat].Fill( cand.nIso, mind(cylinders,cand) )

                h_dict["reco"][cat]["gamma"].Fill( reco_KE )
                if passCyl: 
                    h_dict["reco"][cat]["cyl"].Fill( reco_KE )
                    dtd = vetoParams( veto_time_vtx, cand )
                    dt_min = 9999999999.9
                    d_min = None
                    for idx, deltat, d in dtd:
                        h_dict["dt_d"][cat].Fill( deltat, d/100. )
                        if deltat > -5 and deltat < dt_min:
                            dt_min = deltat
                            d_min = d/100.
                    if d_min is not None: h_dict["d_dt_min"][cat].Fill( d_min )
                    else: h_dict["d_dt_min"][cat].Fill( 9.99 )
                if passCyl and passVeto: 
                    h_dict["reco"][cat]["veto"].Fill( reco_KE )
                    h_dict["iso"][cat].Fill( cand.nIso )
                if passCyl and passVeto and cand.nIso > 50.:
                    h_dict["reco"][cat]["iso"].Fill( reco_KE )

            ihall += 1
            if ihall == Nhall: # we've run out of hall events; loop again
                print "We've run through %d hall events; resetting" % Nhall
                n_hall -= Nhall
                ihall = 0

        # add rock background
        for trock in rock_times:
            RockEvent = Event(RTree, irock)

            rock_candidates = neutron_cut(RockEvent.candidates, vtx)

            for cand in rock_candidates:
                true_neutron = (cand.nTruePDG == 2212 or cand.nTruePDG == 2112)
                cand.nPosT += trock
                cand.nPosT += random.normalvariate(0., dt) # smear by timing uncertainty

                passCyl = (not inChargedCylinder(cylinders, cand))
    
                # Veto the event if it could have come from background
                passVeto = (not isVeto(veto_time_vtx, cand))


                np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
                reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
                if reco_KE is None or reco_KE < 5.: continue # unphysical neutron

                cat = "rockn" if true_neutron else "rockg"

                h_dict["cyl_min"][cat].Fill( mind(cylinders,cand) )
                h_dict["iso_cyl"][cat].Fill( cand.nIso, mind(cylinders,cand) )

                h_dict["reco"][cat]["gamma"].Fill( reco_KE )
                if passCyl: 
                    h_dict["reco"][cat]["cyl"].Fill( reco_KE )
                    dtd = vetoParams( veto_time_vtx, cand )
                    dt_min = 9999999999.9
                    d_min = None
                    for idx, deltat, d in dtd:
                        h_dict["dt_d"][cat].Fill( deltat, d/100. )
                        if deltat > -5 and deltat < dt_min:
                            dt_min = deltat
                            d_min = d/100.
                    if d_min is not None: h_dict["d_dt_min"][cat].Fill( d_min )
                if passCyl and passVeto: 
                    h_dict["reco"][cat]["veto"].Fill( reco_KE )
                    h_dict["iso"][cat].Fill( cand.nIso )
                if passCyl and passVeto and cand.nIso > 50.:
                    h_dict["reco"][cat]["iso"].Fill( reco_KE )

            irock += 1
            if irock == Nrock: # we've run out of rock events; loop again
                print "We've run through %d rock events; resetting" % Nrock
                n_rock -= Nrock
                irock = 0


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--outfile', help='Output File Name', default='FSout.root')
    parser.add_option('--topdir', help='Directory containing Input', default='/pnfs/dune/persistent/users/marshalc/neutronSim/ntuple')
    parser.add_option('--Gfile_str', help='Input Files for Gas Argon File String', default='outGArTPCv4.root')
    parser.add_option('--Rfile_str', help='Input Files for Rock File String', default='outRockv2.root')
    parser.add_option('--Hfile_str', help='Input Files for Hall File String', default='outDetEnclosurev2.root')
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

    random.seed(12345)

    Gfile = ROOT.TFile('%s/%s' % (topdir, Gfile_str))
    Rfile = ROOT.TFile('%s/%s' % (topdir, Rfile_str))
    Hfile = ROOT.TFile('%s/%s' % (topdir, Hfile_str))

    GTree = Gfile.Get('tree')
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

    loop( GTree, RTree, HTree, thresh, dt, muRock, muHall, h_dict, 60000 )

    tfout = ROOT.TFile( args.outfile, "RECREATE" )
    for key in h_dict:
        try: 
            h_dict[key].Write()
        except AttributeError:
            for key2 in h_dict[key]:
                try:
                    h_dict[key][key2].Write()
                except AttributeError:
                    for key3 in h_dict[key][key2]:
                        h_dict[key][key2][key3].Write()
