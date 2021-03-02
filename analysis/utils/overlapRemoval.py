#################
# OVERLAP REMOVAL
#################
# Overlap removal between related samples
# TTGamma and TTbar
# WGamma and WJets
# ZGamma and ZJets
# We need to remove events from TTbar which are already counted in the phase space in which the TTGamma sample is produced
# photon with pT> 10 GeV, eta<5, and at least dR>0.1 from other gen objects 

import numpy as np
import awkward as ak
import re



def overlapFilter(dataset, events):
    genPart = events.GenPart
    passOverlap=np.ones(len(genPart))==1

    if 'TTbar' in dataset:
        passOverlap = overlap(genPart, overlapPt=10., overlapEta=5., overlapDR=0.1)
    if 'TTGamma' in dataset:
        passOverlap = ~overlap(genPart, overlapPt=10., overlapEta=5., overlapDR=0.1)

        ## for ttgamma 1l and 2l, also remove photons with pt>100
        if dataset in ['TTGamma_SingleLept','TTGamma_Dilepton']:
            lhePart = events.LHEPart
            lhePhotons = lhePart[lhePart.pdgId==22]
            passOverlap = passOverlap & ak.flatten(lhePhotons.pt<=100)

    if 'ST_t_channel' in dataset or 'ST_tbar_channel' in dataset:
        passOverlap = overlap(genPart, overlapPt=10., overlapEta=2.6, overlapDR=0.05)
    if 'TGJets' in dataset:
        passOverlap = overlap(genPart, overlapPt=10., overlapEta=2.6, overlapDR=0.05)

            
    if re.search("^W[1234]jets$", dataset):
        passOverlap = overlap(genPart, overlapPt=10., overlapEta=2.5, overlapDR=0.5)
    if re.search("WGamma", dataset):
        passOverlap = ~overlap(genPart, overlapPt=10., overlapEta=2.5, overlapDR=0.5)

    if 'DYjetsM' in dataset:
        passOverlap = overlap(genPart, overlapPt=16., overlapEta=2.6, overlapDR=0.05)
    if 'ZGamma' in dataset:
        passOverlap = ~overlap(genPart, overlapPt=16., overlapEta=2.6, overlapDR=0.05)

    return passOverlap


def overlap(genPart, overlapPt=10., overlapEta=5., overlapDR=0.1):

    
    genmotherIdx = genPart.genPartIdxMother
    genpdgid = genPart.pdgId

    #potential overlap photons are only those passing the kinematic cuts 
    #if the overlap photon is actually from a non prompt decay (maxParent > 37), it's not part of the phase space of the separate sample 
    overlapPhoSelect = ((genPart.pt>=overlapPt) & 
                        (abs(genPart.eta) < overlapEta) & 
                        (genPart.pdgId==22) & 
                        (genPart.status==1) & 
                        (genPart.maxParent < 37)
    )
    overlapPhotons = genPart[overlapPhoSelect] 

    #also require that photons are separate from all other gen particles
    #don't consider neutrinos and don't calculate the dR between the overlapPhoton and itself
    finalGen = genPart[((genPart.status==1)|(genPart.status==71)) & (genPart.pt > 0.01) &
                       ~((abs(genPart.pdgId)==12) | (abs(genPart.pdgId)==14) | (abs(genPart.pdgId)==16)) &
                       ~overlapPhoSelect]

    #calculate dR between overlap photons and nearest gen particle
    phoGen, phoGenDR = overlapPhotons.nearest(finalGen, return_metric = True)
    phoGenMask = ak.fill_none(phoGenDR > overlapDR, True)

    #the event is overlapping with the separate sample if there is an overlap photon passing the dR cut, kinematic cuts, and not coming from hadronic activity
    isOverlap = ak.any(phoGenMask, axis=-1)
    passOverlapRemoval = ~isOverlap

    return passOverlapRemoval
