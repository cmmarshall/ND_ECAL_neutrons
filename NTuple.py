from array import array

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
t_nCandidates = array( 'i', [0] )
t_nPosX = array( 'd', MAXCANDIDATES*[0.] )
t_nPosY = array( 'd', MAXCANDIDATES*[0.] )
t_nPosZ = array( 'd', MAXCANDIDATES*[0.] )
t_nPosT = array( 'd', MAXCANDIDATES*[0.] )
t_nE = array( 'd', MAXCANDIDATES*[0.] )
t_nNcell = array( 'i', MAXCANDIDATES*[0] )
t_nMaxCell = array( 'd', MAXCANDIDATES*[0.] )
t_nTruePDG = array( 'i', MAXCANDIDATES*[0] )
t_nTrueKE = array( 'd', MAXCANDIDATES*[0.] )
t_nParTID = array('i', MAXCANDIDATES*[0])
t_pNeutrons = array( 'i', [0] )
t_pTID = array( 'i', MAXNEUTRONS*[0] )
t_pKE = array( 'd', MAXNEUTRONS*[0.] )

def setBranches( tree ):
    tree.SetBranchAddress( "run", t_run )
    tree.SetBranchAddress( "event", t_event )
    tree.SetBranchAddress( "vtxX", t_vtxX )
    tree.SetBranchAddress( "vtxY", t_vtxY )
    tree.SetBranchAddress( "vtxZ", t_vtxZ )
    tree.SetBranchAddress( "vtxA", t_vtxA )
    tree.SetBranchAddress( "ECAL_visE", t_ECAL_visE )
    tree.SetBranchAddress( "nCandidates", t_nCandidates )
    tree.SetBranchAddress( "nPosX", t_nPosX )
    tree.SetBranchAddress( "nPosY", t_nPosY )
    tree.SetBranchAddress( "nPosZ", t_nPosZ )
    tree.SetBranchAddress( "nPosT", t_nPosT )
    tree.SetBranchAddress( "nE", t_nE )
    tree.SetBranchAddress( "nNcell", t_nNcell )
    tree.SetBranchAddress( "nMaxCell", t_nMaxCell )
    tree.SetBranchAddress( "nTruePDG", t_nTruePDG )
    tree.SetBranchAddress( "nTrueKE", t_nTrueKE )
    tree.SetBranchAddress( "nParTID", t_nParTID )
    tree.SetBranchAddress( "pNeutrons", t_pNeutrons )
    tree.SetBranchAddress( "pTID", t_pTID )
    tree.SetBranchAddress( "pKE", t_pKE )

class Candidate:
    def __init__(self, idx):
        self.nPosX = t_nPosX[idx]
        self.nPosY = t_nPosY[idx]
        self.nPosZ = t_nPosZ[idx]
        self.nPosT = t_nPosT[idx]
        self.nE = t_nE[idx]
        self.nNcell = t_nNcell[idx]
        self.nMaxCell = t_nMaxCell[idx]
        self.nTruePDG = t_nTruePDG[idx]
        self.nTrueKE = t_nTrueKE[idx]
        self.nParTID = t_nParTID[idx]

class PrimaryNeutron:
    def __init__(self, idx):
        self.pTID = t_pTID[idx]
        self.pKE = t_pKE[idx]

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
        self.nCandidates = t_nCandidates[0]
        self.candidates = [Candidate(i) for i in range(self.nCandidates)]
        self.pNeutrons = t_pNeutrons[0]
        self.primaries = [PrimaryNeutron(i) for i in range(self.pNeutrons)]

