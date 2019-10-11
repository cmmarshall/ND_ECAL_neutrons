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

beam_dir = ROOT.TVector3( 0., -0.10418161231515577, 0.9945582897223343 )
cats = ["signal", "duplicate", "gasn", "gasg", "halln", "hallg", "rockn", "rockg"]
cuts = ["none", "gamma", "cyl", "veto", "iso"]

def dist(TVec1, TVec2):
    return (TVec1-TVec2).Mag()

# distance from point p to line specified by point a and unit vector n
def distFromLine( a, n, p ):
    return ((a-p) - n*((a-p).Dot(n))).Mag()

def obvious_photon( cand ):
    return cand.nNcell >= 6
    
def obvious_neutron( cand ):
    return cand.nMaxCell > 6.


def filter_cut(cans, vtx, vtx_t0):

    # order the neutral candidates by distance to the vertex, so that closer hits will be analyzed first
    # the first interaction will typically be the closest to the vertex
    candidates = sorted(cans, key = lambda cluster: dist(vtx, cluster.getPos()))

    # Remove non-physical neutrons, and really low energy blips
    i = 0
    while i < len(candidates):
        can = candidates[i]
        #nCut = (can.nMaxCell > 10. or (can.nMaxCell>3. and can.nNcell <= 4))
        reco_KE = GetRecoE(vtx, can.getPos(), can.nPosT-vtx_t0)
        if can.nE > 3. and reco_KE is not None and reco_KE > 5.:
            i += 1
        else:
            candidates.pop(i)

    # Remove candidates that are spatially close to other candidates, which often are due to multi-scatters
    i = 0
    while i < len(candidates):
        a = candidates[i].getPos()
        n = (vtx-a).Unit()
        j = i+1

        if obvious_neutron(candidates[i]):
            candidates[i].SetRecoNeutron()
        elif obvious_photon(candidates[i]):
            candidates[i].SetRecoGamma()

        while j < len(candidates):
            distance = distFromLine( a, n, candidates[j].getPos() )
            distance2 = dist(a,candidates[j].getPos())
            if distance < 50. or distance2 < 50:
                candidates[i].UpCylinder(candidates[j].nNcell, candidates[j].nE)
                if obvious_photon(candidates[j]) and not candidates[i].RecoNeutron:
                    candidates[i].SetRecoGamma()
                elif obvious_neutron(candidates[j]) and not candidates[i].RecoGamma:
                    candidates[i].SetRecoNeutron()
                elif obvious_photon(candidates[j]) or obvious_neutron(candidates[j]):
                    candidates[i].SetRecoWTF()
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

def inChargedCylinder( cylinders, cand ):
    p = cand.getPos()
    for cylinder in cylinders:
        a = cylinder[1]
        n = cylinder[2].Unit()
        d = distFromLine(a, n, p)

        if d < 70.:
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
    for vtxt,vtx in veto_time_vtx:
        dt = nt - vtxt
        d = (p-vtx).Mag()

        if dt > -5. and dt < 20. and d < 200.:
            return True

        # velocity between 10 and 30 cm/ns, i.e. beta between 0.3 and 1
        #if dt > 0 and d/dt > 10. and d/dt < 30.:
        #    return True

    return False

def vetoParams( veto_time_vtx, cand ):
    p = cand.getPos()
    nt = cand.nPosT

    ret = []
    for vtxt,vtx in veto_time_vtx:
        dt = nt - vtxt
        d = (p-vtx).Mag()
        ret.append( (dt, d) )
    return ret

def Initialize_Plot_Dict(h_dict):
    h_dict["reco"] = {}
    h_dict["reco_leading"] = {}
    h_dict["cyl_min"] = {}
    h_dict["iso"] = {}
    h_dict["iso_cyl"] = {}
    h_dict["dt_d"] = {}
    h_dict["d_dt_min"] = {}
    h_dict["visE"] = {}
    h_dict["tot_can"] = {}
    h_dict["eresgood"] = ROOT.TH2D( "eresgood", ";True neutron KE (MeV);Fractional energy resolution", 70, 0., 700., 150, -1., 2. )
    h_dict["eresbad"] = ROOT.TH2D( "eresbad", ";True neutron KE (MeV);Fractional energy resolution", 70, 0., 700., 150, -1., 2. )
    h_dict["eff"] = {}
    h_dict["eff_leading"] = {}
    h_dict["denom"] = ROOT.TH1D( "denom", ";True neutron KE (MeV)", 70, 0., 700. )
    h_dict["denom_leading"] = ROOT.TH1D( "denom_leading", ";True neutron KE (MeV)", 70, 0., 700. )

    h_dict["depth"] = {}
    h_dict["kink"] = {}
    h_dict["cell_max"] = {}
    h_dict["cell_n"] = {}
    h_dict["max_n"] = {}
    h_dict["cell_gaps"] = {}
    h_dict["cell_sigmas"] = {}
    h_dict["n_cylinder"] = {}
    h_dict["cell_layers"] = {}
    h_dict["odvies"] = {}

    for cat in cats:
        h_dict["reco"][cat] = {}
        h_dict["reco_leading"][cat] = {}
        for cut in cuts:
            h_dict["reco"][cat][cut] = ROOT.TH2D( "reco_%s_%s" % (cat,cut), ";Reco angle (deg);Reco neutron KE (MeV)", 36, 0., 180., 70, 0., 700. )
            h_dict["reco_leading"][cat][cut] = ROOT.TH2D( "reco_leading_%s_%s" % (cat,cut), ";Reco angle (deg);Reco neutron KE (MeV)", 36, 0., 180., 70, 0., 700. )

        h_dict["cyl_min"][cat] = ROOT.TH1D( "cyl_min_%s" % cat, ";Min dist to charged traj (cm)", 100, 0., 500. )
        h_dict["iso"][cat] = ROOT.TH1D( "iso_%s" % cat, ";Isolation (cm)", 100, 0., 500. )
        h_dict["iso_cyl"][cat] = ROOT.TH2D( "isocyl_%s" % cat, ";Isolation (cm);Cylinder (cm)", 50, 0., 500., 50, 0., 500. )
        h_dict["dt_d"][cat] = ROOT.TH2D( "dt_d_%s" % cat, ";#Deltat to ECAL activity (ns);Dist to ECAL activity (m)", 250, -150., 200., 100, 0., 10. )
        h_dict["d_dt_min"][cat] = ROOT.TH1D( "d_dt_min_%s" % cat, ";Dist to ECAL activity (m)", 100, 0., 10. )
        h_dict["visE"][cat] = ROOT.TH1D( "visE_%s" % cat, ";ECAL visible energy (MeV)", 100, 0., 3000. )
        h_dict["tot_can"][cat] = ROOT.TH1D( "tot_can_%s" % cat, ";Total neutron candidates", 100, 0., 200. )
        h_dict["depth"][cat] = ROOT.TH1D( "depth_%s" % cat, ";Depth (cm)", 100, 0., 100. )
        h_dict["kink"][cat] = ROOT.TH1D( "kink_%s" % cat, ";Kink angle (deg)", 90, 0., 180. )

        h_dict["cell_max"][cat] = ROOT.TH2D( "cell_max_%s" % cat, ";Cluster energy (MeV);Max cell (MeV)", 100, 0., 100., 100, 0., 100. )
        h_dict["cell_n"][cat] = ROOT.TH2D( "cell_n_%s" % cat, ";Cluster energy (MeV);N cells", 100, 0., 100., 100, 0., 100. )
        h_dict["max_n"][cat] = ROOT.TH2D( "max_n_%s" % cat, ";Max cell (MeV);N cells", 100, 0., 20., 20, 0., 20. )
        h_dict["cell_gaps"][cat] = ROOT.TH2D( "cell_gaps_%s" % cat, ";Cluster energy (MeV);Gaps", 100, 0., 100., 30, 0., 30. )
        h_dict["cell_sigmas"][cat] = ROOT.TH2D( "cell_sigmas_%s" % cat, ";Cluster sigmaX;Cluster sigmaY", 100, 0., 1., 100, 0., 1. )
        h_dict["n_cylinder"][cat] = ROOT.TH2D( "n_cylinder_%s" % cat, ";Cylinder energy (MeV);Number in cylinder", 100, 0., 100., 100, 0., 100. )
        h_dict["cell_layers"][cat] = ROOT.TH2D( "cell_layers_%s" % cat, ";N cell;N layers", 20, 0., 20., 20, 0., 20. )
        h_dict["odvies"][cat] = ROOT.TH1D( "odvies_%s" % cat, ";0=n 1=g 2=neither", 4, 0., 4. )

    for cut in cuts:
        h_dict["eff"][cut] = ROOT.TH2D( "eff_%s" % cut, ";Reco angle (deg);True neutron KE (MeV)", 36, 0., 180., 70, 0., 700. )
        h_dict["eff_leading"][cut] = ROOT.TH2D( "eff_leading_%s" % cut, ";Reco angle (deg);True neutron KE (MeV)", 36, 0., 180., 70, 0., 700. )

def writeHistos(outfile):
    tfout = ROOT.TFile( outfile, "RECREATE" )
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
    tfout.Close()

def candidateLoop( h_dict, candidates, source, spill_info ):

    vtx = spill_info["vtx"]
    reco_t0 = spill_info["reco_t0"]
    maxKink = spill_info["maxKink"]
    total_ECAL_visE = spill_info["total_ECAL_visE"]
    total_nCandidates = spill_info["total_nCandidates"]
    veto_time_vtx = spill_info["veto_time_vtx"] 
    cylinders = spill_info["cylinders"]
    tid_KE = spill_info["tid_KE"]
    leading_tid = spill_info["leading_tid"] 

    quiet = total_ECAL_visE < 600. and total_nCandidates < 100

    reco_tids = [] # save track IDs of reco neutrons to check for duplicates, i.e. 1 neutron --> 2 or more candidates
    leading_reco_KE = [0. for c in cuts]
    leading_reco_cat = [None for c in cuts]
    leading_reco_angle = [None for c in cuts]

    for cand in candidates:

        np = ROOT.TVector3( cand.nPosX, cand.nPosY, cand.nPosZ )
        reco_KE = GetRecoE(vtx,np,cand.nPosT-reco_t0)
        true_neutron = (cand.nTruePDG == 2212 or cand.nTruePDG == 2112)

        # reco angle in degrees
        reco_angle = 57.29577951308232 * (np-vtx).Angle(beam_dir)

        # remove super-luminal neutrons, or neutrons so delayed w.r.t. neutrino that they have very low reco energy
        if reco_KE is None or reco_KE < 5.: continue

        # Cuts
        passGamma = cand.RecoNeutron or (cand.nNcell < 5.*(cand.nMaxCell-4.))
        passCyl = (not inChargedCylinder(cylinders, cand))
        passVeto = (not isVeto(veto_time_vtx, cand))

        # Determine category (signal/bkg, n/gamma)
        cat = None
        if source == "gas":
            tid = cand.nPrimTID
            if tid in tid_KE: # signal neutron!!
                trueKE = tid_KE[tid]
                if tid in reco_tids: # duplicate neutron; count the second one as background
                    cat = "duplicate"
                else: # Good signal neutronm fill efficiencies
                    cat = "signal"

                    h_dict["eff"]["none"].Fill( reco_angle, trueKE )
                    if tid == leading_tid: h_dict["eff_leading"]["none"].Fill( reco_angle, trueKE )
                    if passGamma: 
                        h_dict["eff"]["gamma"].Fill( reco_angle, trueKE )
                        if tid == leading_tid: h_dict["eff_leading"]["gamma"].Fill( reco_angle, trueKE )
                        if passCyl: 
                            h_dict["eff"]["cyl"].Fill( reco_angle, trueKE )
                            if tid == leading_tid: h_dict["eff_leading"]["cyl"].Fill( reco_angle, trueKE )
                            if passVeto: 
                                h_dict["eff"]["veto"].Fill( reco_angle, trueKE ) # eff numerator fill with true KE
                                if tid == leading_tid: h_dict["eff_leading"]["veto"].Fill( reco_angle, trueKE )
                                if tid == cand.nParTID:
                                    h_dict["eresgood"].Fill( trueKE, (reco_KE-trueKE)/trueKE )
                                else:
                                    h_dict["eresbad"].Fill( trueKE, (reco_KE-trueKE)/trueKE )
                                if cand.nIso > 70.: 
                                    h_dict["eff"]["iso"].Fill( reco_angle, trueKE )
                                    if tid == leading_tid: h_dict["eff_leading"]["iso"].Fill( reco_angle, trueKE )
                    reco_tids.append(tid)
            else: # not true primary neutron, i.e. gas background
                cat = "gasn" if true_neutron else "gasg"
        else:
            cat = "%sn" % source
            if not true_neutron:
                cat = "%sg" % source

        h_dict["cell_max"][cat].Fill( cand.nE, cand.nMaxCell )
        h_dict["cell_n"][cat].Fill( cand.nE, cand.nNcell )
        h_dict["max_n"][cat].Fill( cand.nMaxCell, cand.nNcell )
        h_dict["cell_gaps"][cat].Fill( cand.nE, cand.nGaps )
        h_dict["cell_sigmas"][cat].Fill( cand.nSigmaX, cand.nSigmaY )
        h_dict["n_cylinder"][cat].Fill( cand.eCylinder, cand.nCylinder )
        h_dict["cell_layers"][cat].Fill( cand.nNcell, cand.nNlayers )

        if cand.RecoNeutron and not cand.RecoGamma: h_dict["odvies"][cat].Fill(0)
        elif cand.RecoGamma and not cand.RecoNeutron: h_dict["odvies"][cat].Fill(1)
        elif not cand.RecoGamma and not cand.RecoNeutron: h_dict["odvies"][cat].Fill(2)
        else: h_dict["odvies"][cat].Fill(3)

        # Fill reconstruction histograms
        h_dict["reco"][cat]["none"].Fill( reco_angle, reco_KE )
        if reco_KE > leading_reco_KE[0]:
            leading_reco_KE[0] = reco_KE
            leading_reco_cat[0] = cat
            leading_reco_angle[0] = reco_angle
        if passGamma:
            h_dict["cyl_min"][cat].Fill( mind(cylinders,cand) )
            h_dict["iso_cyl"][cat].Fill( cand.nIso, mind(cylinders,cand) )

            h_dict["reco"][cat]["gamma"].Fill( reco_angle, reco_KE )
            if reco_KE > leading_reco_KE[1]:
                leading_reco_KE[1] = reco_KE
                leading_reco_cat[1] = cat
                leading_reco_angle[1] = reco_angle
            if passCyl: 
                h_dict["reco"][cat]["cyl"].Fill( reco_angle, reco_KE )
                if reco_KE > leading_reco_KE[2]:
                    leading_reco_KE[2] = reco_KE
                    leading_reco_cat[2] = cat
                    leading_reco_angle[2] = reco_angle
                dtd = vetoParams( veto_time_vtx, cand )
                dt_min = 9999999999.9
                d_min = None
                for deltat, d in dtd:
                    h_dict["dt_d"][cat].Fill( deltat, d/100. )
                    if deltat > -5 and deltat < dt_min:
                        dt_min = deltat
                        d_min = d/100.
                if d_min is not None: h_dict["d_dt_min"][cat].Fill( d_min )
                else: h_dict["d_dt_min"][cat].Fill( 9.99 )
                if passVeto: 
                    h_dict["reco"][cat]["veto"].Fill( reco_angle, reco_KE )
                    if reco_KE > leading_reco_KE[3]:
                        leading_reco_KE[3] = reco_KE
                        leading_reco_cat[3] = cat
                        leading_reco_angle[3] = reco_angle

                    h_dict["visE"][cat].Fill( total_ECAL_visE )
                    h_dict["tot_can"][cat].Fill( total_nCandidates )    
                    h_dict["iso"][cat].Fill( cand.nIso )
                    if abs(cand.nPosX) < 355.: h_dict["depth"][cat].Fill( cand.getDepth() )
                    h_dict["kink"][cat].Fill( maxKink )

                    if cand.nIso > 70.:
                        h_dict["reco"][cat]["iso"].Fill( reco_angle, reco_KE )
                        if reco_KE > leading_reco_KE[4]:
                            leading_reco_KE[4] = reco_KE
                            leading_reco_cat[4] = cat
                            leading_reco_angle[4] = reco_angle

    return leading_reco_KE, leading_reco_cat, leading_reco_angle

def loop(GTree, RTree, HTree, Emin, dt, muRock, muHall, Plot_dict, Ngas=None, write=None, outfile=None):

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
            print "Spill %d of %d...VM %1.1fMB RM %1.1fMB, %1.1f minutes remaining" % (igas, Ngas, vm/1000., rm/1000., remaining/60.)

        if write is not None and outfile is not None and igas and not igas % write:
            writeHistos(outfile)

        GasEvent = Event(GTree, igas)

        vtx = ROOT.TVector3(GasEvent.vtxX, GasEvent.vtxY, GasEvent.vtxZ)
        # impose a fiducial volume cut; for simplicity
        R = sqrt((vtx.y()+217.)**2 + (vtx.z()-585.)**2) # YZ circle centered at (-217, 585)
        if R > 200. or abs(vtx.x()) > 200.: 
            continue # x is centered -250, 250

        # Signal time
        t0 = random.uniform(0.,10000.) # ns
        reco_t0 = random.normalvariate(t0, dt)

        total_nCandidates = 0
        for c in GasEvent.candidates:
            c.nPosT += random.normalvariate(t0, dt)
            if c.nPosT > reco_t0 and c.nPosT < reco_t0+50.:
                total_nCandidates += 1

        gas_candidates = filter_cut(GasEvent.candidates, vtx, reco_t0)

        n_rock = rando.Poisson(muRock)
        n_hall = rando.Poisson(muHall)

        # generate potential veto vertex times
        veto_time_vtx = []
        rock_times = [ random.uniform(0., 10000.) for i in range(n_rock) ]
        hall_times = [ random.uniform(0., 10000.) for i in range(n_hall) ]

        total_ECAL_visE = GasEvent.ECAL_visE

        hall_candidates = []
        rock_candidates = []

        # Generate vertex time stamps and add to the veto list, if applicable
        for i,t in enumerate(rock_times):
            # 5 MeV neutron takes ~600 ns to go from upstream rock to back of ECAL
            if t-t0 > -700. and t-t0 < 300:
                RockEvent = Event(RTree, irock)
                for c in RockEvent.candidates:
                    c.nPosT += random.normalvariate(t, dt)
                    if c.nPosT > reco_t0 and c.nPosT < reco_t0+50.:
                        total_nCandidates += 1
                rock_candidates += filter_cut(RockEvent.candidates, vtx, reco_t0)

                if t-t0> -300:
                    if RockEvent.ECAL_vetoT+t > reco_t0-5. and RockEvent.ECAL_vetoT+t < reco_t0+10.:
                        total_ECAL_visE += RockEvent.ECAL_visE
                    elif RockEvent.ECAL_visE > 10.:
                        veto_time_vtx.append( (RockEvent.ECAL_vetoT+t, RockEvent.ECAL_vetoP) )

                irock += 1
                if irock == Nrock:
                    print "We've run through %d rock events; resetting" % Nrock
                    irock = 0

        for i,t in enumerate(hall_times):
            if t-t0 > -700. and t-t0 < 300:
                HallEvent = Event(HTree, ihall)
                for c in HallEvent.candidates:
                    c.nPosT += random.normalvariate(t, dt)
                    if c.nPosT > reco_t0 and c.nPosT < reco_t0+50.:
                        total_nCandidates += 1
                hall_candidates += filter_cut(HallEvent.candidates, vtx, reco_t0)

                if t-t0 > -300.:
                    if HallEvent.ECAL_vetoT+t > reco_t0-5. and HallEvent.ECAL_vetoT+t < reco_t0+10.:
                        total_ECAL_visE += HallEvent.ECAL_visE
                    elif HallEvent.ECAL_visE > 10.:
                        veto_time_vtx.append( (HallEvent.ECAL_vetoT+t, HallEvent.ECAL_vetoP) )

                ihall += 1
                if ihall == Nhall:
                    print "We've run through %d hall events; resetting" % Nhall
                    ihall = 0

        # Get true neutrons, and other primaries that might produce backgrounds
        tid_KE = {}
        cylinders = []
        leading_KE = None
        leading_tid = None
        maxKink = 0.
        for pn in GasEvent.primaries:
            if pn.pPDG == 2112:
                #print "True neutron KE %1.1f" % (pn.pKE)
                tid_KE[pn.pTID] = pn.pKE
                h_dict["denom"].Fill( pn.pKE )
                if leading_KE is None or pn.pKE > leading_KE:
                    leading_KE = pn.pKE
                    leading_tid = pn.pTID

            if pn.pPos is not None:
                cylinders.append( (pn.pPDG, pn.pPos, pn.pDir) )

            if pn.pKinkAngle > maxKink:
                maxKink = pn.pKinkAngle

        if leading_KE is not None:
            h_dict["denom_leading"].Fill( leading_KE )

        spill_info = {}
        spill_info["vtx"] = vtx
        spill_info["reco_t0"] = reco_t0
        spill_info["maxKink"] = maxKink
        spill_info["total_ECAL_visE"] = total_ECAL_visE
        spill_info["total_nCandidates"] = total_nCandidates
        spill_info["veto_time_vtx"] = veto_time_vtx
        spill_info["cylinders"] = cylinders
        spill_info["tid_KE"] = tid_KE
        spill_info["leading_tid"] = leading_tid

        gas_KE,  gas_cat,  gas_angle  = candidateLoop( h_dict, gas_candidates,  "gas",  spill_info )
        hall_KE, hall_cat, hall_angle = candidateLoop( h_dict, hall_candidates, "hall", spill_info )
        rock_KE, rock_cat, rock_angle = candidateLoop( h_dict, rock_candidates, "rock", spill_info )

        # Fill leading neutron plots
        for icut in range(len(gas_KE)):
            angle = gas_angle[icut]
            cat = gas_cat[icut]
            leading_KE = gas_KE[icut]
            if hall_KE[icut] > gas_KE[icut] and hall_KE[icut] > rock_KE[icut]:
                leading_KE = hall_KE[icut]
                angle = hall_angle[icut]
                cat = hall_cat[icut]
            elif rock_KE[icut] > hall_KE[icut] and rock_KE[icut] > gas_KE[icut]:
                leading_KE = rock_KE[icut]
                angle = rock_angle[icut]
                cat = rock_cat[icut]
            # else gas is leading
            if leading_KE > 0.: # there is at least one neutron
                h_dict["reco_leading"][cat][cuts[icut]].Fill( angle, leading_KE )

    print "Did %d gas events in %1.1f minutes" % (Ngas, (time.time()-stime)/60.)


if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option('--outfile', help='Output File Name', default='FSout.root')
    parser.add_option('--topdir', help='Directory containing Input', default='/pnfs/dune/persistent/users/marshalc/neutronSim/ntuple')
    parser.add_option('--Gfile_str', help='Input Files for Gas Argon File String', default='outGArTPCv7.root')
    parser.add_option('--Rfile_str', help='Input Files for Rock File String', default='outRockv5.root')
    parser.add_option('--Hfile_str', help='Input Files for Hall File String', default='outDetEnclosurev5.root')
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

    #loop( GTree, RTree, HTree, thresh, dt, muRock, muHall, h_dict, 100000, 10000, args.outfile )
    loop( GTree, RTree, HTree, thresh, dt, muRock, muHall, h_dict, GTree.GetEntries(), 10000, args.outfile )
    writeHistos( args.outfile )


