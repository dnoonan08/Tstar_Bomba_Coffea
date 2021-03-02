import awkward as ak
import numpy as np
from .jetPermutations import permutations,argpermutations

def chi2_Top(leptons, neutrinos, jets, N=4, hists=None, lepFlavor=None, weights=None, dataset=None, btagCut=0.6321):
    # tt chi2

    jetPerms = permutations(jets,N)
    bL = jetPerms[:,:,0]
    bH = jetPerms[:,:,1]

    nTags = ak.sum(jets.btagDeepB>btagCut,axis=1)
    nTags = ak.where(nTags>2,2,nTags)

    nBCandTags = ((bL.btagDeepB>btagCut)*1.+(bH.btagDeepB>btagCut)*1.)
    passBtags = nTags==nBCandTags

    lepList, nuList, jList = ak.unzip(ak.cartesian([leptons, neutrinos, jetPerms[passBtags]]))
            
    topLep = (lepList + nuList + bL)
    W = (jList[:,:,2] + jList[:,:,3])
    topHad = (bH + W)

    chiTopLep = (topLep.mass-172.5)/30.
    chiTopHad = (topHad.mass - 172.5)/34.
    chiW = (W.mass - 80.4)/24.

    chi2 = chiTopHad**2 + chiTopLep**2 + chiW**2 
    minCombo = ak.argmin(chi2,axis=1,keepdims=True)
    return chi2[minCombo], topHad[minCombo].mass, topLep[minCombo].mass

#     if not hists is None:
#         if lepFlavor is None:
#             lepFlavor=np.zeros(len(leptons))-1
#         if weights is None:
#             weights=np.ones(len(leptons))
#         if dataset is None:
#             dataset='None'
#         if 'tt_WMass' in hists:
#             hists['tt_WMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                    mass=ak.flatten(W.mass[minCombo]))
#         if 'tt_topHadMass' in hists:
#             hists['tt_topHadMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                         mass=ak.flatten(topHad.mass[minCombo]))
#         if 'tt_topLepMass' in hists:
#             hists['tt_topLepMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                         mass=ak.flatten(topLep.mass[minCombo]))
#         if 'tt_chi2' in hists:
#             hists['tt_chi2'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                   chi2=ak.flatten(chi2[minCombo]))
#     return chi2[minCombo], topHad[minCombo], topLep[minCombo]




def chi2_TT_tgtg(leptons, neutrinos, jets, N=6, hists=None, lepFlavor=None, weights=None, dataset=None, btagCut=0.6321):
    # tt chi2

    jetPerms = permutations(jets,N)
    bL = jetPerms[:,:,0]
    bH = jetPerms[:,:,1]

    nTags = ak.sum(jets.btagDeepB>btagCut,axis=1)
    nTags = ak.where(nTags>2,2,nTags)

    nBCandTags = ((bL.btagDeepB>btagCut)*1.+(bH.btagDeepB>btagCut)*1.)
    passBtags = nTags==nBCandTags

    lepList, nuList, jList = ak.unzip(ak.cartesian([leptons, neutrinos, jetPerms[passBtags]]))

    topLep = (lepList + nuList + bL)
    W = (jList[:,:,2] + jList[:,:,3])
    topHad = (bH + W)
    tStarLep = topLep+jList[:,:,4]
    tStarHad = topHad+jList[:,:,5]

    chiTopLep = (topLep.mass-172.5)/30.
    chiTopHad = (topHad.mass - 172.5)/34.
    chiW = (W.mass - 80.4)/24.
    chiTstar_h = (tStarHad.mass - tStarLep.mass)/230.

    chi2 = chiTopHad**2 + chiTopLep**2 + chiW**2 + chiTstar**2 

    tStarMass = (tStarLep.mass + tStarHad.mass)/2.

    minCombo = ak.argmin(chi2,axis=1,keepdims=True)
    return chi2[minCombo], tStarMass[minCombo], tStarHad.mass[minCombo], tStarLep.mass[minCombo]

#     if not hists is None:
#         if lepFlavor is None:
#             lepFlavor=np.zeros(len(leptons))-1
#         if weights is None:
#             weights=np.ones(len(leptons))
#         if dataset is None:
#             dataset='None'
            
#         if 'tgtg_WMass' in hists:
#             hists['tgtg_WMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                      mass=ak.flatten(W.mass[minCombo]))
#         if 'tgtg_topHadMass' in hists:
#             hists['tgtg_topHadMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                           mass=ak.flatten(topHad.mass[minCombo]))
#         if 'tgtg_topLepMass' in hists:
#             hists['tgtg_topLepMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                           mass=ak.flatten(topLep.mass[minCombo]))
#         if 'tgtg_tstarMass' in hists:            
#             hists['tgtg_tstarMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                          mass=ak.flatten(tStarMass[minCombo]),
#                                          phopt = ak.flatten(phoList[minCombo].pt))
#         if 'tgtg_tstarHadMass' in hists:            
#             hists['tgtg_tstarHadMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                             mass=ak.flatten(tStarHad.mass[minCombo]),
#                                             phopt = ak.flatten(phoList[minCombo].pt))
#         if 'tgtg_tstarLepMass' in hists:            
#             hists['tgtg_tstarLepMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                             mass=ak.flatten(tStarLep.mass[minCombo]),
#                                             phopt = ak.flatten(phoList[minCombo].pt))
#         if 'tgtg_chi2' in hists:
#             hists['tgtg_chi2'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
#                                    chi2=ak.flatten(chi2[minCombo]))

#     return chi2[minCombo], tStarMass[minCombo], tStarHad.mass[minCombo], tStarLep.mass[minCombo]

# def chi2_TT_tgty(leptons, neutrinos, jets, photons, N=5, hists=None, lepFlavor=None, weights=None, dataset=None, btagCut=0.6321):

#     #jetPerms = permutations(jets,N)
#     jetPerms = argpermutations(jets,N)

#     bL = jets[jetPerms[:,:,0]]
#     bH = jets[jetPerms[:,:,1]]

#     nTags = ak.sum(jets.btagDeepB>btagCut,axis=1)
#     nTags = ak.where(nTags>2,2,nTags)

#     nBCandTags = ((bL.btagDeepB>btagCut)*1.+(bH.btagDeepB>btagCut)*1.)
#     passBtags = nTags==nBCandTags
#     del bL,bH
#     phoList, lepList, nuList, jList = ak.unzip(ak.cartesian([photons, leptons, neutrinos, jetPerms[passBtags]]))

#     del jetPerms

#     bL = jets[jList[:,:,0]]
#     bH = jets[jList[:,:,1]]
    
#     topLep = (lepList + nuList + bL)
#     W = (jets[jList[:,:,2]] + jets[jList[:,:,3]])
#     topHad = (bH + W)
#     tStarLep_l = topLep + phoList
#     tStarHad_l = topHad + jets[jList[:,:,4]]

#     tStarLep_h = topLep + jets[jList[:,:,4]]
#     tStarHad_h = topHad + phoList

#     del bL, bH, jList

#     chiTopLep = (topLep.mass-172.5)/30.
#     chiTopHad = (topHad.mass - 172.5)/34.
#     chiW = (W.mass - 80.4)/24.
#     chiTstar_l = (tStarHad_l.mass - tStarLep_l.mass)/230.
#     chiTstar_h = (tStarHad_h.mass - tStarLep_h.mass)/230.

#     chi2_h = chiTopHad**2 + chiTopLep**2 + chiW**2 + chiTstar_h**2
#     chi2_l = chiTopHad**2 + chiTopLep**2 + chiW**2 + chiTstar_l**2

#     phoIsHadSide = ak.min(chi2_h,axis=1)<ak.min(chi2_l,axis=1)

#     chi2 = ak.where(phoIsHadSide,chi2_h, chi2_l)
#     tStarLepMass = ak.where(phoIsHadSide,tStarLep_h.mass, tStarLep_l.mass)
#     tStarHadMass = ak.where(phoIsHadSide,tStarHad_h.mass, tStarHad_l.mass)
#     tStarMass = (tStarLepMass + tStarHadMass)/2.

#     minCombo = ak.argmin(chi2,axis=1,keepdims=True)
    
#     return chi2[minCombo], tStarMass[minCombo], tStarHadMass[minCombo], tStarLepMass[minCombo], phoIsHadSide, phoList[minCombo].pt

# #     if not hists is None:
# #         print('filling')
# #         if lepFlavor is None:
# #             lepFlavor=np.zeros(len(leptons))-1
# #         if weights is None:
# #             weights=np.ones(len(leptons))
# #         if dataset is None:
# #             dataset='None'
            
# #         if 'tgty_WMass' in hists:
# #             hists['tgty_WMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
# #                                      mass=ak.flatten(W.mass[minCombo]))
# #         if 'tgty_topHadMass' in hists:
# #             hists['tgty_topHadMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
# #                                           mass=ak.flatten(topHad.mass[minCombo]))
# #         if 'tgty_topLepMass' in hists:
# #             hists['tgty_topLepMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
# #                                           mass=ak.flatten(topLep.mass[minCombo]))
# #         if 'tgty_tstarMass' in hists:            
# #             hists['tgty_tstarMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
# #                                          mass=ak.flatten(tStarMass[minCombo]),
# #                                          phopt = ak.flatten(phoList[minCombo].pt))
# #         if 'tgty_tstarHadMass' in hists:            
# #             hists['tgty_tstarHadMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
# #                                             mass=ak.flatten(tStarHadMass[minCombo]),
# #                                             phopt = ak.flatten(phoList[minCombo].pt))
# #         if 'tgty_tstarLepMass' in hists:            
# #             hists['tgty_tstarLepMass'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
# #                                             mass=ak.flatten(tStarLepMass[minCombo]),
# #                                             phopt = ak.flatten(phoList[minCombo].pt))
# #         if 'tgty_chi2' in hists:
# #             hists['tgty_chi2'].fill(dataset=dataset, lepFlavor=lepFlavor, weight=weights,
# #                                    chi2=ak.flatten(chi2[minCombo]))

# # #    return chi2[minCombo], tStarMass[minCombo], tStarHadMass[minCombo], tStarLepMass[minCombo], phoIsHadSide
# #     return hists



def chi2_TT_tgty(leptons, neutrinos, jets, photons, N=5, hists=None, lepFlavor=None, weights=None, dataset=None, btagCut=0.6321):
    #jetPerms = permutations(jets,N)
    jetPerms = argpermutations(jets,N)

    bL = jets[jetPerms[:,:,0]]
    bH = jets[jetPerms[:,:,1]]

    nTags = ak.sum(jets.btagDeepB>btagCut,axis=1)
    nTags = ak.where(nTags>2,2,nTags)

    nBCandTags = ((bL.btagDeepB>btagCut)*1.+(bH.btagDeepB>btagCut)*1.)
    passBtags = nTags==nBCandTags
    del bL,bH, nTags, nBCandTags

    phoList, lepList, nuList, jList = ak.unzip(ak.cartesian([photons, leptons, neutrinos, jetPerms[passBtags]]))

    del jetPerms

    bL = jets[jList[:,:,0]]
    bH = jets[jList[:,:,1]]
    
    topLep = (lepList + nuList + bL)
    W = (jets[jList[:,:,2]] + jets[jList[:,:,3]])
    topHad = (bH + W)
    tStarLep_l = topLep + phoList
    tStarHad_l = topHad + jets[jList[:,:,4]]

    tStarLep_h = topLep + jets[jList[:,:,4]]
    tStarHad_h = topHad + phoList

    del bL, bH, jList

    chiTopLep = (topLep.mass-172.5)/30.
    chiTopHad = (topHad.mass - 172.5)/34.
    chiW = (W.mass - 80.4)/24.

    del W, topLep, topHad

    chiTstar_l = (tStarHad_l.mass - tStarLep_l.mass)/230.
    chiTstar_h = (tStarHad_h.mass - tStarLep_h.mass)/230.

    chi2_h = chiTopHad**2 + chiTopLep**2 + chiW**2 + chiTstar_h**2
    chi2_l = chiTopHad**2 + chiTopLep**2 + chiW**2 + chiTstar_l**2

    del chiTopHad, chiTopLep, chiW, chiTstar_h, chiTstar_l

    phoIsHadSide = ak.min(chi2_h,axis=1)<ak.min(chi2_l,axis=1)

    chi2 = ak.where(phoIsHadSide,chi2_h, chi2_l)
    tStarLepMass = ak.where(phoIsHadSide,tStarLep_h.mass, tStarLep_l.mass)
    tStarHadMass = ak.where(phoIsHadSide,tStarHad_h.mass, tStarHad_l.mass)
    tStarMass = (tStarLepMass + tStarHadMass)/2.

    minCombo = ak.argmin(chi2,axis=1,keepdims=True)

    return chi2[minCombo], tStarMass[minCombo], tStarHadMass[minCombo], tStarLepMass[minCombo], phoIsHadSide, phoList[minCombo].pt

