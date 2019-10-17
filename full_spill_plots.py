import ROOT
ROOT.gROOT.SetBatch(1)

cats = ["signal", "duplicate", "gasn", "gasg", "halln", "hallg", "rockn", "rockg"]
cuts = ["none", "gamma", "cyl", "veto", "iso"]
cutnames = [ "None", "Photon", "Charge trk", "ECAL veto", "Isolation" ]
names = ["Signal", "Duplicate", "Gas TPC n","Gas TPC #gamma", "Hall n", "Hall #gamma", "Rock n", "Rock #gamma"]
regions = { "forward":[1,12], "mid":[13,24], "backward":[25,36], "all":[1,36] }

colorsl = [ ROOT.kBlack, ROOT.kGray, ROOT.kRed, ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan, ROOT.kGreen, ROOT.kGreen+3 ]
colorsf = [ ROOT.kWhite, ROOT.kGray, ROOT.kRed, ROOT.kMagenta, ROOT.kBlue, ROOT.kCyan, ROOT.kGreen, ROOT.kGreen+3 ]

ymax = [ 0. for r in regions ]

def plotRes( h2good, h2bad, b0, b1, name ):
    can = ROOT.TCanvas()
    h1good = h2good.ProjectionY( "tmpg_%d_%d" % (b0,b1), b0, b1 )
    h1bad = h2bad.ProjectionY( "tmpb_%d_%d" % (b0,b1), b0, b1 )
    h1good.Rebin(2)
    h1bad.Rebin(2)
    stk = ROOT.THStack( "stk", "%1.0f < T_{n} < %1.0f MeV;(Reco - True) / True KE" % (h2good.GetXaxis().GetBinLowEdge(b0), h2good.GetXaxis().GetBinLowEdge(b1+1)) )

    h1good.SetLineColor(2)
    h1good.SetFillColor(2)
    h1bad.SetLineColor(4)
    h1bad.SetFillColor(4)
    stk.Add( h1bad )
    stk.Add( h1good )
    leg = ROOT.TLegend( 0.5, 0.5, 0.846, 0.846 )
    leg.SetFillStyle(0)
    leg.AddEntry( h1good, "#splitline{1st scatter (%1.0f%%)}{#mu = %1.2f #sigma = %1.2f}" % (100.*h1good.Integral()/(h1good.Integral()+h1bad.Integral()), h1good.GetMean(), h1good.GetRMS()), "f" )
    leg.AddEntry( h1bad, "#splitline{Rescatter (%1.0f%%)}{#mu = %1.2f #sigma = %1.2f}" % (100.*h1bad.Integral()/(h1good.Integral()+h1bad.Integral()), h1bad.GetMean(), h1bad.GetRMS()), "f" )

    stk.Draw("hist")
    leg.Draw()
    can.Print( "plots/%s.png" % name )

def plotEff( hEff, hDenom, rebin, name ):
    can = ROOT.TCanvas()

    colors = [ ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan ]
    leg = ROOT.TLegend( 0.2, 0.65, 0.8, 0.846 )
    leg.SetNColumns(2)
    for i,cut in enumerate(cuts):
        hEffAll = hEff[i].ProjectionY("effall_%s" % cut)
        hEffForward = hEff[i].ProjectionY("effforward_%s" % cut, 1, 12)
        hEffAll.Rebin(rebin)
        hEffForward.Rebin(rebin)
        hClone = hDenom.Clone( "tmp_%s" % cut )
        hClone.Rebin(rebin)
        hEffAll.Divide( hClone )
        hEffAll.SetLineColor( colors[i] )
        hEffAll.SetMinimum(0.)
        hEffAll.SetMaximum(1.)
        opt = "hist" if not i else "hist same"
        hEffAll.Draw(opt)
        leg.AddEntry( hEffAll, cutnames[i], "l" )
        if i == len(cuts)-1:
            hEffForward.Divide(hClone)
            hEffForward.SetLineColor( colors[len(cuts)] )
            hEffForward.Draw("hist same")
            leg.AddEntry( hEffForward, "Forward", "l" )

    leg.Draw()
    can.Print( "plots/%s.png" % name )

def plotRecoRegions( hReco, low, high, name, title, updateMax = False ):
    can = ROOT.TCanvas()

    for a,r in enumerate(regions):

        leg = ROOT.TLegend( 0.5, 0.6, 0.846, 0.846 )
        leg.SetFillStyle(0)
        hPurity = [None for h in hReco]
        for i in range(len(hReco)):
            t = "f" if i else "l"
            hPurity[i] = ROOT.TH1D( "hpurity_%d_%s_%s" % (i,r,name), "", hReco[0].GetNbinsY(), 0., hReco[0].GetYaxis().GetBinLowEdge(hReco[0].GetNbinsY()+1) )
            hPurity[i].SetLineColor(colorsl[i])
            hPurity[i].SetFillColor(colorsf[i])
            leg.AddEntry( hPurity[i], names[i], t )

        stk = ROOT.THStack( "stk_%s" % r, title )

        blow = regions[r][0]
        bhigh = regions[r][1]
        for i in range(len(cats)-1, -1, -1):
            h = hReco[i].ProjectionY("%s_%s" % (cats[i],r), blow, bhigh)
            h.SetLineColor(colorsl[i])
            h.SetFillColor(colorsf[i])
            h.GetXaxis().SetRangeUser(low+0.001, high-0.001)
            stk.Add( h )
        stk.Draw()
        if updateMax:
            ymax[a] = 1.1*stk.GetStack().Last().GetBinContent( hReco[0].GetYaxis().FindBin(50.001) )
        #stk.SetMaximum( ymax[a] )
        stk.Draw()
        leg.Draw()
        can.Print( "plots/%s_%s.png" % (r,name) )

        stksum = stk.GetStack().Last()

        # purity stacks
        puritystk = ROOT.THStack( "pstk_%s" % r, "%s;Fraction" % title )
        for b in range( 1, hReco[i].GetNbinsY()+1 ):
            for i,h in enumerate(hReco):
                if stksum.GetBinContent(b) > 0.:
                    hPurity[i].SetBinContent( b, h.ProjectionY("tmp", blow, bhigh).GetBinContent(b) / stksum.GetBinContent(b) )
                else:
                    hPurity[i].SetBinContent( b, 0. )

        for h in hPurity:
            puritystk.Add( h )

        puritystk.Draw()
        can.Print( "plots/purity_%s_%s.png" % (r,name) )        

def plotReco( hReco, low, high, name, title, accLow = True ):
    can = ROOT.TCanvas()

    leg = ROOT.TLegend( 0.5, 0.6, 0.846, 0.846 )
    leg.SetFillStyle(0)
    htmp = [None for h in hReco]
    for i in range(len(hReco)):
        t = "f" if i else "l"
        htmp[i] = ROOT.TH1D( "htmp_%d" % i, "", 1, 0, 1 )
        htmp[i].SetLineColor(colorsl[i])
        htmp[i].SetFillColor(colorsf[i])
        leg.AddEntry( htmp[i], names[i], t )

    stk = ROOT.THStack( "stk", title )

    for i in range(len(cats)-1, -1, -1):
        hReco[i].SetLineColor(colorsl[i])
        hReco[i].SetFillColor(colorsf[i])
        hReco[i].GetXaxis().SetRangeUser(low+0.001, high-0.001)
        stk.Add( hReco[i] )
    stk.Draw()
    stk.GetXaxis().SetRangeUser(low+0.001, high-0.001)
    stk.Draw()
    leg.Draw()
    can.Print( "plots/%s.png" % name )

    eff = ROOT.TGraph()
    pur = ROOT.TGraph()
    exp = ROOT.TGraph()
    for b in range( 1, hReco[0].GetNbinsX()+1 ):
        cut = hReco[0].GetBinLowEdge(b+1)
        sig = hReco[0].Integral(1,b) if accLow else hReco[0].Integral(b+1, hReco[0].GetNbinsX()+1)
        bkg = 0.
        totbkg = 0.
        for i in range(1, len(hReco)):
            addbkg = hReco[i].Integral(1,b) if accLow else hReco[i].Integral(b+1, hReco[i].GetNbinsX()+1)
            bkg += addbkg
            totbkg += hReco[i].Integral(1, hReco[i].GetNbinsX()+1)
        totsig = hReco[0].Integral(1,hReco[0].GetNbinsX()+1)

        eff.SetPoint( b-1, cut, sig/totsig )
        if sig+bkg > 0.:
            pur.SetPoint( b-1, cut, sig/(sig+bkg) )
            exp.SetPoint( b-1, cut, (sig/totsig)*sig/(sig+bkg) )
        else:
            pur.SetPoint( b-1, cut, 0. )
            exp.SetPoint( b-1, cut, 0. )

    eff.SetLineColor(2)
    pur.SetLineColor(4)
    eff.SetMarkerColor(2)
    pur.SetMarkerColor(4)
    eff.SetMarkerSize(0.2)
    pur.SetMarkerSize(0.2)
    exp.SetMarkerSize(0.2)

    exp.Draw("APC")
    exp.GetYaxis().SetRangeUser(0., 1.0)
    exp.Draw("APC")
    eff.Draw("PC same")
    pur.Draw("PC same")
    can.Print( "plots/cut_%s.png" % name )

def plot2Dall( h, name, title ):
    can = ROOT.TCanvas()
    for i,cat in enumerate(cats):
        h[i].SetTitle( title )
        h[i].Draw("colz")
        can.Print( "plots/%s_%s.png" % (name,cat) )

def plot2Done( h, name, title ):
    can = ROOT.TCanvas()
    h.SetTitle( title )
    h.Draw("colz")
    can.Print( "plots/%s.png" % name )

if __name__ == "__main__":

    tf = ROOT.TFile( "FSoutRHC.root" )

    hRes2Dgood = tf.Get( "eresgood" )
    hRes2Dbad = tf.Get( "eresbad" )
    hEff = [None for cut in cuts]
    hDenom = tf.Get( "denom" )
    hEffLeading = [None for cut in cuts]
    hDenomLeading = tf.Get( "denom_leading" )

    hReco = [[None for cat in cats] for cut in cuts]
    hRecoLeading = [[None for cat in cats] for cut in cuts]
    hMinCyl = [None for cat in cats]
    hIso = [None for cat in cats]
    hDist = [None for cat in cats]
    hIsoCyl = [None for cat in cats]
    hTimeDist = [None for cat in cats]
    hVisE = [None for cat in cats]
    hTotCan = [None for cat in cats]
    hKink = [None for cat in cats]
    hDepth = [None for cat in cats]
    hMax_n = [None for cat in cats]
    for i,cat in enumerate(cats):
        for j,cut in enumerate(cuts):
            hReco[j][i] = tf.Get( "reco_%s_%s" % (cat,cut) )
            hRecoLeading[j][i] = tf.Get( "reco_leading_%s_%s" % (cat,cut) )
        hMinCyl[i] = tf.Get( "cyl_min_%s" % cat ) 
        hIsoCyl[i] = tf.Get( "isocyl_%s" % cat )  
        hTimeDist[i] = tf.Get( "dt_d_%s" % cat ) 
        hIso[i] = tf.Get( "iso_%s" % cat )
        hDist[i] = tf.Get( "d_dt_min_%s" % cat )
        hVisE[i] = tf.Get( "visE_%s" % cat )
        hTotCan[i] = tf.Get( "tot_can_%s" % cat )
        hKink[i] = tf.Get( "kink_%s" % cat )
        hDepth[i] = tf.Get( "depth_%s" % cat )
        hMax_n[i] = tf.Get( "max_n_%s" % cat )
        

    for j,cut in enumerate(cuts):
        hEff[j] = tf.Get( "eff_%s" % cut )
        hEffLeading[j] = tf.Get( "eff_leading_%s" % cut )

    for b in range( 1, 66, 5 ):
        plotRes( hRes2Dgood, hRes2Dbad, b, b+4, "KE_res_%03d_%03d" % ((b-1)*10, (b+4)*10) )

    plotEff( hEff, hDenom, 5, "eff" )
    plotEff( hEffLeading, hDenomLeading, 5, "eff_leading" )

    plotRecoRegions( hReco[0], 0., 700., "KE_reco_0none", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hReco[1], 0., 700., "KE_reco_1gamma", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hReco[2], 0., 700., "KE_reco_2cyl", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hReco[3], 0., 700., "KE_reco_3veto", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hReco[4], 0., 700., "KE_reco_4iso", ";Reconstructed KE (MeV)" )

    plotRecoRegions( hRecoLeading[0], 0., 700., "leading_KE_reco_0none", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hRecoLeading[1], 0., 700., "leading_KE_reco_1gamma", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hRecoLeading[2], 0., 700., "leading_KE_reco_2cyl", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hRecoLeading[3], 0., 700., "leading_KE_reco_3veto", ";Reconstructed KE (MeV)" )
    plotRecoRegions( hRecoLeading[4], 0., 700., "leading_KE_reco_4quiet", ";Reconstructed KE (MeV)" )

    plotReco( hMinCyl, 0., 500., "min_cyl", ";Distance to charged particle (cm)", False )
    plotReco( hIso, 0., 500., "iso", ";Distance to nearest cluster (cm)", False )
    plotReco( hDist, 0., 10., "dist", ";Distance to ECAL veto (m)", False )
    plotReco( hVisE, 0., 1500., "ecal_visE", ";ECAL visible energy (MeV)", True )
    plotReco( hTotCan, 0., 200., "tot_can", ";Total neutron candidates", True )
    plotReco( hKink, 0., 180., "kink", ";Max kink angle (degrees)", True )
    plotReco( hDepth, 0., 100., "depth", ";Neutron depth (cm)", True )

    plot2Dall( hTimeDist, "veto_time_dist", ";#Deltat to ECAL activity (ns);Distance to ECAL activity (m)" )
    plot2Dall( hIsoCyl, "iso_cyl", ";Dist to nearest ECAL cluster (cm);Distance to nearest TPC track (cm)" )
    plot2Done( hRes2Dgood, "eres2d", ";Neutron KE (MeV);Fracional energy residual" )
    plot2Done( hRes2Dbad, "eres2dbad", ";Neutron KE (MeV);Fracional energy residual" )

    c2 = ROOT.TCanvas()
    c2.SetLogz()
    hMax_n[0].Draw("colz")
    c2.Print( "plots/max_n_neutron_log.png" )
    hMax_n[3].Draw("colz")
    c2.Print( "plots/max_n_gamma_log.png" )

