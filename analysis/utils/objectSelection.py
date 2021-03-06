import awkward as ak

def selectMuons(muons):
    muonSelectTight = ((muons.pt>=30) & 
                       (abs(muons.eta)<=2.4) & 
                       (muons.tightId) & 
                       (muons.pfRelIso04_all < 0.15)
                      )

    muonSelectLoose = ((muons.pt>=15) & 
                       (abs(muons.eta)<=2.4) & 
                       ((muons.isPFcand) & (muons.isTracker | muons.isGlobal)) & 
                       (muons.pfRelIso04_all < 0.25) &
                       ~muonSelectTight
                      )

    return muons[muonSelectTight], muons[muonSelectLoose]        


def selectMuonsMiniIso(muons):
    muonSelectTight = ((muons.pt>=55) & 
                       (abs(muons.eta)<=2.4) & 
                       (muons.tightId) & 
                       (muons.miniPFRelIso_all < 0.1)
                      )

    muonSelectLoose = ((muons.pt>=30) & 
                       (abs(muons.eta)<=2.4) & 
                       ((muons.isPFcand) & (muons.isTracker | muons.isGlobal)) & 
                       (muons.miniPFRelIso_all < 0.4) &
                       ~muonSelectTight
                      )

    return muons[muonSelectTight], muons[muonSelectLoose]        



def selectElectrons(electrons, photons):
    #define electron cuts
    scEta = electrons.eta + electrons.deltaEtaSC
    eleEtaGap = (abs(scEta) < 1.4442) | (abs(scEta) > 1.566)
    elePassDXY = ((abs(electrons.eta) < 1.479) & (abs(electrons.dxy) < 0.05) |
                 (abs(electrons.eta) > 1.479)  & (abs(electrons.dxy) < 0.1)
                )
    elePassDZ = ((abs(electrons.eta) < 1.479) & (abs(electrons.dz) < 0.1) |
                 (abs(electrons.eta) > 1.479)  & (abs(electrons.dz) < 0.2)
                )

    electronSelectTight = ((electrons.pt>=35) & 
                           (abs(electrons.eta)<=2.4) & 
                           eleEtaGap &      
                           electrons.isTight &
                           elePassDXY & 
                           elePassDZ
                          )

    vetoID = (electrons.vidNestedWPBitmap & 7) >= 1
    for i_cut in [1,2,3,4,5,6,8,9]:
        vetoID = vetoID & (((electrons.vidNestedWPBitmap>>(i_cut*3)) & 7) >= 1)

    vetoIsoCut = ak.where(scEta<1.479, 0.198 + 0.506/electrons.pt, 0.203 + 0.963/electrons.pt)
    phoMatch = ak.unflatten(ak.flatten(photons)[ak.flatten(electrons.photonIdxG)],ak.num(electrons))
    phoRelIso = ak.where(electrons.photonIdxG>-1, phoMatch.pfRelIso03_all, 999)
    looseIso = ak.where(electrons.pfRelIso03_all < phoRelIso, electrons.pfRelIso03_all, phoRelIso)
    vetoIDIso = looseIso<vetoIsoCut
    vetoID = vetoID & vetoIDIso
    
    electronSelectLoose = ((electrons.pt>=15) & 
                           (abs(electrons.eta)<=2.4) & 
                           eleEtaGap &      
                           vetoID &
                           elePassDXY & 
                           elePassDZ & 
                           ~electronSelectTight
                          )

    return electrons[electronSelectTight], electrons[electronSelectLoose]

def selectElectronsPFIso(electrons, photons):
    #define electron cuts
    scEta = electrons.eta + electrons.deltaEtaSC
    eleEtaGap = (abs(scEta) < 1.4442) | (abs(scEta) > 1.566)
    elePassDXY = ((abs(electrons.eta) < 1.479) & (abs(electrons.dxy) < 0.05) |
                 (abs(electrons.eta) > 1.479)  & (abs(electrons.dxy) < 0.1)
                )
    elePassDZ = ((abs(electrons.eta) < 1.479) & (abs(electrons.dz) < 0.1) |
                 (abs(electrons.eta) > 1.479)  & (abs(electrons.dz) < 0.2)
                )

    electronSelectTight = ((electrons.pt>=35) & 
                           (abs(electrons.eta)<=2.4) & 
                           eleEtaGap &      
                           electrons.isTight &
                           elePassDXY & 
                           elePassDZ
                          )

    electronSelectLoose = ((electrons.pt>=15) & 
                           (abs(electrons.eta)<=2.4) & 
                           eleEtaGap &      
                           electrons.isVeto &
                           elePassDXY & 
                           elePassDZ & 
                           ~electronSelectTight
                          )

    return electrons[electronSelectTight], electrons[electronSelectLoose]

def selectElectronsMiniIso(electrons, photons):
    #define electron cuts
    scEta = electrons.eta + electrons.deltaEtaSC
    eleEtaGap = (abs(scEta) < 1.4442) | (abs(scEta) > 1.566)
    elePassDXY = ((abs(electrons.eta) < 1.479) & (abs(electrons.dxy) < 0.05) |
                 (abs(electrons.eta) > 1.479)  & (abs(electrons.dxy) < 0.1)
                )
    elePassDZ = ((abs(electrons.eta) < 1.479) & (abs(electrons.dz) < 0.1) |
                 (abs(electrons.eta) > 1.479)  & (abs(electrons.dz) < 0.2)
                )

    electronSelectTight = ((electrons.pt>=50) & 
                           (abs(electrons.eta)<=2.2) & 
                           eleEtaGap &      
                           electrons.mvaFall17V2noIso_WP90 & 
                           (electrons.miniPFRelIso_all < 0.1)
                          )

    electronSelectLoose = ((electrons.pt>=35) & 
                           (abs(electrons.eta)<=2.2) & 
                           eleEtaGap &      
                           electrons.mvaFall17V2noIso_WPL &
                           ~electrons.mvaFall17V2noIso_WPL & 
                           (electrons.miniPFRelIso_all < 0.4)

                          )

    return electrons[electronSelectTight], electrons[electronSelectLoose]

def selectPhotons(photons, muons, electrons):
    phoMu, phoMuDR  = photons.nearest(muons,return_metric=True)
    phoMuMask = ak.fill_none(phoMuDR > 0.4, True)

    phoEle, phoEleDR = photons.nearest(electrons, return_metric=True)
    phoEleMask = ak.fill_none(phoEleDR > 0.4, True)
    photonSelect = ((photons.pt>20) & 
                    (abs(photons.eta) < 1.4442) &
#                     (photons.isScEtaEE | photons.isScEtaEB) &
                    (photons.electronVeto) & 
                    ~photons.pixelSeed & 
                    phoMuMask & phoEleMask &
                    (photons.cutBased>=2)
                   )

    return photons[photonSelect]


def selectJets(jets, muons, electrons, photons, bTagWP):
    # ##check dR jet,lepton & jet,photon
    jetMu, jetMuDR = jets.nearest(muons, return_metric=True)
    jetMuMask = ak.fill_none(jetMuDR > 0.4, True)

    jetEle, jetEleDR = jets.nearest(electrons, return_metric=True)
    jetEleMask = ak.fill_none(jetEleDR > 0.4, True)

    jetPho, jetPhoDR = jets.nearest(photons, return_metric=True)
    jetPhoMask = ak.fill_none(jetPhoDR > 0.4, True)

    jetSel = ((jets.pt>30) & 
              (abs(jets.eta)<2.4)&
              jetMuMask & jetEleMask & jetPhoMask & 
              (jets.isLoose) 
             )

    btagged = jets[jetSel].btagDeepB>bTagWP  
    
    return jets[jetSel], btagged



def selectFatJets(fatJets):
    fatJetSel = ((fatJets.pt>300)&
                 (abs(fatJets.eta)<2.4) & 
                 fatJets.isLoose)
    return fatJets[fatJetSel]
