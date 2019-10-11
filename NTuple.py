from array import array
from ROOT import TVector3
from math import sqrt

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
t_nE = array( 'd', MAXCANDIDATES*[0.] )
t_nIso = array( 'd', MAXCANDIDATES*[0.] )
t_nNcell = array( 'i', MAXCANDIDATES*[0] )
t_nNlayers = array( 'i', MAXCANDIDATES*[0] )
t_nGaps = array( 'i', MAXCANDIDATES*[0] )
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
t_pKinkAngle = array( 'd', MAXNEUTRONS*[0.] )

def setBranches( tree ):
    tree.SetBranchAddress( "run", t_run )
    tree.SetBranchAddress( "event", t_event )
    tree.SetBranchAddress( "vtxX", t_vtxX )
    tree.SetBranchAddress( "vtxY", t_vtxY )
    tree.SetBranchAddress( "vtxZ", t_vtxZ )
    tree.SetBranchAddress( "vtxA", t_vtxA )
    tree.SetBranchAddress( "ECAL_visE", t_ECAL_visE )
    tree.SetBranchAddress( "ECAL_vetoX", t_ECAL_vetoX )
    tree.SetBranchAddress( "ECAL_vetoY", t_ECAL_vetoY )
    tree.SetBranchAddress( "ECAL_vetoZ", t_ECAL_vetoZ )
    tree.SetBranchAddress( "ECAL_vetoT", t_ECAL_vetoT )
    tree.SetBranchAddress( "nCandidates", t_nCandidates )
    tree.SetBranchAddress( "nPosX", t_nPosX )
    tree.SetBranchAddress( "nPosY", t_nPosY )
    tree.SetBranchAddress( "nPosZ", t_nPosZ )
    tree.SetBranchAddress( "nPosT", t_nPosT )
    tree.SetBranchAddress( "nSigmaX", t_nSigmaX )
    tree.SetBranchAddress( "nSigmaY", t_nSigmaY )
    tree.SetBranchAddress( "nE", t_nE )
    tree.SetBranchAddress( "nIso", t_nIso )
    tree.SetBranchAddress( "nNcell", t_nNcell )
    tree.SetBranchAddress( "nMaxCell", t_nMaxCell )
    tree.SetBranchAddress( "nNlayers", t_nNlayers )
    tree.SetBranchAddress( "nGaps", t_nGaps )
    tree.SetBranchAddress( "nTruePDG", t_nTruePDG )
    tree.SetBranchAddress( "nTrueKE", t_nTrueKE )
    tree.SetBranchAddress( "nParTID", t_nParTID )
    tree.SetBranchAddress( "nPrimTID", t_nPrimTID )
    tree.SetBranchAddress( "pParticles", t_pParticles )
    tree.SetBranchAddress( "pTID", t_pTID )
    tree.SetBranchAddress( "pPDG", t_pPDG )
    tree.SetBranchAddress( "pKE", t_pKE )
    tree.SetBranchAddress( "pExitX", t_pExitX )
    tree.SetBranchAddress( "pExitY", t_pExitY )
    tree.SetBranchAddress( "pExitZ", t_pExitZ )
    tree.SetBranchAddress( "pExitdX", t_pExitdX )
    tree.SetBranchAddress( "pExitdY", t_pExitdY )
    tree.SetBranchAddress( "pExitdZ", t_pExitdZ )
    tree.SetBranchAddress( "pKinkAngle", t_pKinkAngle )

class Candidate:
    def __init__(self, idx):
        self.nPosX = t_nPosX[idx]
        self.nPosY = t_nPosY[idx]
        self.nPosZ = t_nPosZ[idx]
        self.nPosT = t_nPosT[idx]
        self.nSigmaX = t_nSigmaX[idx]
        self.nSigmaY = t_nSigmaY[idx]
        self.nE = t_nE[idx]
        self.nIso = t_nIso[idx]
        self.nNcell = t_nNcell[idx]
        self.nMaxCell = t_nMaxCell[idx]
        self.nNlayers = t_nNlayers[idx]
        self.nGaps = t_nGaps[idx]
        self.nTruePDG = t_nTruePDG[idx]
        self.nTrueKE = t_nTrueKE[idx]
        self.nParTID = t_nParTID[idx]
        self.nPrimTID = t_nPrimTID[idx]

        # additional externally-set cut variables
        self.RecoGamma = False
        self.RecoNeutron = False
        self.nCylinder = 0
        self.eCylinder = 0.
        self.inNeutronCylinder = False
        self.inGammaCylinder = False

    def getPos(self):
        pos = TVector3( self.nPosX, self.nPosY, self.nPosZ )
        return pos

    def getDepth(self):
        ret = None
        if abs(self.nPosX) > 356.: # endcap
            ret = min(abs(self.nPosX+356.), abs(self.nPosX-356.))
        elif sqrt( (self.nPosY+217.)**2 + (self.nPosZ-585.)**2 ) > 250.: # barrel
            ret = sqrt( (self.nPosY+217.)**2 + (self.nPosZ-585.)**2 ) - 250.
        else:
            print "N pos X = (%1.0f, %1.0f, %1.0f) R = %1.0f" % (self.nPosX, self.nPosY, self.nPosZ, sqrt( (self.nPosY+217.)**2 + (self.nPosZ-585.)**2 ))
            ret = -1.

        return ret

    def UpCylinder(self, n, e):
        self.nCylinder += n
        self.eCylinder += e

    def SetRecoGamma(self):
        self.RecoGamma = True
        self.RecoNeutron = False
    def SetRecoNeutron(self):
        self.RecoGamma = False
        self.RecoNeutron = True
    def SetRecoWTF(self):
        self.RecoGamma = True
        self.RecoNeutron = True


class PrimaryParticle:
    def __init__(self, idx):
        self.pTID = t_pTID[idx]
        self.pPDG = t_pPDG[idx]
        self.pKE = t_pKE[idx]
        self.pKinkAngle = t_pKinkAngle[idx] * 180. / 3.1416

        self.pPos = None
        self.pDir = None
        if t_pExitdX[idx] > -9999.:
            self.pPos = TVector3( t_pExitX[idx], t_pExitY[idx], t_pExitZ[idx] )
            self.pDir = TVector3( t_pExitdX[idx], t_pExitdY[idx], t_pExitdZ[idx] )

class Event:
    def __init__(self,tree,idx):
        tree.GetEntry(idx)

        self.run = t_run[0]
        self.event = t_event[0]
        self.vtxX = t_vtxX[0]
        self.vtxY = t_vtxY[0]
        self.vtxZ = t_vtxZ[0]
        self.vtxA = t_vtxA[0]
        self.ECAL_visE = t_ECAL_visE[0]
        self.ECAL_vetoT = t_ECAL_vetoT[0]
        self.ECAL_vetoP = TVector3( t_ECAL_vetoX[0], t_ECAL_vetoY[0], t_ECAL_vetoZ[0] )
        self.nCandidates = t_nCandidates[0]
        self.candidates = [Candidate(i) for i in range(self.nCandidates)]
        self.pParticles = t_pParticles[0]
        self.primaries = [PrimaryParticle(i) for i in range(self.pParticles)]

