import awkward as ak
import numpy as np

from coffea import hist, util, processor
from coffea.analysis_tools import PackedSelection

from coffea.btag_tools import BTagScaleFactor
from coffea.lookup_tools import extractor, dense_lookup
from coffea.analysis_tools import Weights
from coffea.jetmet_tools import CorrectedJetsFactory, JECStack


from .utils.objectSelection import selectMuons, selectElectrons, selectPhotons, selectJets, selectFatJets
from .utils.MetZCalculator import MetPz, getNeutrinoArrays
from .utils.overlapRemoval import overlapFilter
from .utils.genParentage import maxHistoryPDGID
from .utils.crossSections import getLumiWeight
from .utils.chi2Calculator import chi2_TT_tgty,chi2_TT_tgtg,chi2_Top

from .utils.PUReweight import getPileupReweight
from .utils.ElectronScaleFactors import ElectronScaleFactors
from .utils.MuonScaleFactors import MuonScaleFactors
from .utils.BtagWeights import BtagWeights

import time
import psutil

class TstarSelector(processor.ProcessorABC):
    def __init__(self,isMC=True, topChi2=False, tgtgChi2=False, tgtyChi2=True, mcTotals=None, bJetCut=None, mcPUDist=None):
        self.isMC = isMC
        self.topChi2=topChi2
        self.tgtgChi2=tgtgChi2
        self.tgtyChi2=tgtyChi2
        self.mcTotals=mcTotals
        self.bJetCut=bJetCut
        self.mcPUDist=mcPUDist
        
        dataset_axis = hist.Cat("dataset","Dataset")

        mTop_axis = hist.Bin("mass","$m_{T}$ [GeV]", 400, 0, 400)
        mW_axis = hist.Bin("mass","$m_{W}$ [GeV]", 200, 0, 200)
        mTstar_axis = hist.Bin("mass","$m_{T*}$ [GeV]", 500, 0, 2500)
        chi2_axis = hist.Bin('chi2','$\chi^2$',1000,0,1000)
        njet_axis = hist.Bin('njet','Number of Jets',20,0,20)
        nbjet_axis = hist.Bin('nbjet','Number of B-tagged Jets',20,0,20)
        npho_axis = hist.Bin('npho','Number of photons',20,0,20)
        phoPt_axis = hist.Bin('pt','$p_{T}(\gamma)$',300,0,1500)
        phoEta_axis = hist.Bin('eta','$\eta(\gamma)$',50,-2.5,2.5)


        lep_axis = hist.Bin('lepFlavor','Lepton Flavor',2,0,2)
#         phoHighPt_axis = hist.Bin('phopt','$p_{T}(\gamma)$',np.array([20,50,100,150,200,250,300,500,1000,1500]))

        cut_axis = hist.Cat('cut','Cut')

        self._accumulator = processor.dict_accumulator({
#             'tt_WMass': hist.Hist("Counts", dataset_axis, lep_axis, mW_axis),
            'tt_topHadMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis),
            'tt_topLepMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis),
            'tt_chi2': hist.Hist("Counts", dataset_axis, lep_axis, chi2_axis),

#             'tgtg_topHadMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis),
#             'tgtg_topLepMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis),
#             'tgtg_WMass': hist.Hist("Counts", dataset_axis, lep_axis, mW_axis),
#             'tgtg_tstarHadMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis),
#             'tgtg_tstarLepMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis),
#             'tgtg_chi2': hist.Hist("Counts", dataset_axis, lep_axis, chi2_axis),

#             'tgty_WMass': hist.Hist("Counts", dataset_axis, lep_axis, mW_axis),
#             'tgty_topHadMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis),
#             'tgty_topLepMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis),
            'tgty_tstarMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis),
            'tgty_tstarHadMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis),
            'tgty_tstarLepMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis),
            'tgty_chi2': hist.Hist("Counts", dataset_axis, lep_axis, chi2_axis),

            'phoPt': hist.Hist("Counts", dataset_axis, lep_axis, phoPt_axis),
            'phoEta': hist.Hist("Counts", dataset_axis, lep_axis, phoEta_axis),
            
            'nJet':hist.Hist("Counts",dataset_axis, lep_axis, njet_axis),
            'nBJet':hist.Hist("Counts",dataset_axis, lep_axis, nbjet_axis),
            'nPho':hist.Hist("Counts",dataset_axis, lep_axis, npho_axis),
            'cutflow':hist.Hist("Counts",dataset_axis, lep_axis, cut_axis),
            })
        
        self.bCalc = BtagWeights()
        self.eleSF = ElectronScaleFactors()
        self.muSF  = MuonScaleFactors()
        
    @property
    def accumulator(self):
        return self._accumulator

        
    def process(self, events):
        ## To access variables from the ntuples, use the "events" object
        ## The dataset name is part of events.metadata
        dataset, year = events.metadata['dataset']

        # print('processor start :',psutil.Process().memory_info())
        output = self.accumulator.identity()

        bjetCut='1b'
        cutList = ['overlap','lepton', 'met', '3jet', '4jet', '5jet', '1b', 'photon', 'photon100', 'photon200', 'photon300', 'photon500']
        if not self.bJetCut is None:
            goodCuts=['','1b','2b','ex2b']
            if self.bJetCut in goodCuts:
                idx = goodCuts.index(self.bJetCut)+1
                bjetCut = self.bJetCut
                
                cutList = ['overlap','lepton', 'met', '3jet', '4jet', '5jet']+ goodCuts[1:idx] + ['photon', 'photon100', 'photon200', 'photon300', 'photon500']
            
        
        selection = PackedSelection()

        #Calculate the maximum pdgID of any of the particles in the GenPart history
        if self.isMC:
            idx = ak.to_numpy(ak.flatten(abs(events.GenPart.pdgId)))
            par = ak.to_numpy(ak.flatten(events.GenPart.genPartIdxMother))
            num = ak.to_numpy(ak.num(events.GenPart.pdgId))
            maxParentFlatten = maxHistoryPDGID(idx,par,num)
            # events["GenPart","maxParent"] = ak.broadcast_arrays(events.GenPart.pdgId),1.)*maxParentFlatten
            events["GenPart","maxParent"] = ak.unflatten(maxParentFlatten, num)

        # print('after gen history :',psutil.Process().memory_info())
        #object selection
        tightMuons, looseMuons = selectMuons(events.Muon)
        tightElectrons, looseElectrons = selectElectrons(events.Electron)
        tightPhotons = selectPhotons(events.Photon, tightMuons, tightElectrons)
        tightJets, btag = selectJets(events.Jet, tightMuons, tightElectrons, tightPhotons)
        met = events.MET        
        
        # print('after object selection :',psutil.Process().memory_info())
        #event selection
        selection.add('overlap', overlapFilter(dataset, events))
        selection.add('muTrigger', events.HLT.IsoMu24 | events.HLT.IsoTkMu24)
        selection.add('eleTrigger', events.HLT.Ele27_WPTight_Gsf)
        
        selection.add('oneMuon', (ak.num(tightMuons) == 1))
        selection.add('oneEle', (ak.num(tightElectrons) == 1))
        selection.add('muVeto', (ak.num(tightMuons) == 0))
        selection.add('eleVeto', (ak.num(tightElectrons) == 0))

        selection.add('looseVeto', (ak.num(looseMuons) == 0) & (ak.num(looseElectrons) == 0))

        
        selection.add('singleMu',selection.all(*('oneMuon','eleVeto', 'looseVeto')))
        selection.add('singleEle',selection.all(*('oneEle','muVeto', 'looseVeto')))
        
        muon_eventSelection = selection.all(*('muTrigger','singleMu'))
        ele_eventSelection = selection.all(*('eleTrigger','singleEle'))

        selection.add('lepton', muon_eventSelection | ele_eventSelection)
        
        selection.add('3jet',(ak.num(tightJets)>=3))
        selection.add('4jet',(ak.num(tightJets)>=4))
        selection.add('5jet',(ak.num(tightJets)>=5))
        selection.add('6jet',(ak.num(tightJets)>=6))
        selection.add('',np.ones(len(events))==1)
        selection.add('1b',(ak.num(tightJets[btag])>=1))
        selection.add('2b',(ak.num(tightJets[btag])>=2))
        selection.add('ex2b',(ak.num(tightJets[btag])==2))
        selection.add('met',(met.pt>20))
        selection.add('photon',(ak.num(tightPhotons)>=1))
        selection.add('photon100',(ak.any(tightPhotons.pt>100,axis=1)))
        selection.add('photon200',(ak.any(tightPhotons.pt>200,axis=1)))
        selection.add('photon300',(ak.any(tightPhotons.pt>300,axis=1)))
        selection.add('photon500',(ak.any(tightPhotons.pt>500,axis=1)))
        
        # print('after selectors :',psutil.Process().memory_info())
        evtSel = selection.all(*('overlap','lepton','3jet',bjetCut,'met'))
        lepFlavor = np.zeros_like(evtSel)-1.
        lepFlavor[ele_eventSelection] = 0
        lepFlavor[muon_eventSelection] = 1
            
        weightCollection = Weights(len(events))
        weightCollection.add('lumiWeight',events.Generator.weight*getLumiWeight(dataset,self.mcTotals, year))


        puWeight = getPileupReweight(events.Pileup.nTrueInt, self.mcPUDist, dataset, year)
        weightCollection.add('puWeight', weight=puWeight[0], weightUp=puWeight[1], weightDown=puWeight[2])


        muWeight = self.muSF.getMuonScaleFactors(tightMuons, year, 'total')
        weightCollection.add('muSF',weight=muWeight[0], weightUp=muWeight[1], weightDown=muWeight[2])
            
        eleWeight = self.eleSF.getEleScaleFactors(tightElectrons, year, 'total')
        weightCollection.add('eleSF',weight=eleWeight[0], weightUp=eleWeight[1], weightDown=eleWeight[2])
            
        btagsWeight = self.bCalc.getBtagWeights(tightJets, btag, year, dataset)
        weightCollection.add('btagSF',weight=btagsWeight[0], weightUp=btagsWeight[1], weightDown=btagsWeight[2])


        # print('before met :',psutil.Process().memory_info())
        ##construct awkward array four vectors for the neutrino, with each of the hypothesis pz
        nu_pz = MetPz(tightElectrons, tightMuons, met, solutionType=-1)
        # print('after met :',psutil.Process().memory_info())

        neutrinos = getNeutrinoArrays(nu_pz, met)

        tightLep = ak.with_name(ak.concatenate([tightElectrons, tightMuons], axis=1), "PtEtaPhiMCandidate")

        weights=weightCollection.weight(None)
        
        # print('before chi2 :',psutil.Process().memory_info())
        if (ak.sum(evtSel)>0) and self.topChi2:
            sel_jets = tightJets[evtSel]
            sel_leptons = tightLep[evtSel]
            sel_neutrinos = neutrinos[evtSel]

            chi2, topHadMass, topLepMass = chi2_Top(sel_leptons, sel_neutrinos, sel_jets, 4, btagCut=0.6321)

            output['tt_topHadMass'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[evtSel],
                                        mass=ak.flatten(topHadMass))

            output['tt_topLepMass'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[evtSel],
                                        mass=ak.flatten(topLepMass))

            output['tt_chi2'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[evtSel],
                                  chi2=ak.flatten(chi2))
            
        evtSelPho = selection.all(*('overlap','lepton','5jet',bjetCut,'met','photon'))

        if ak.sum(evtSelPho)>0 and self.tgtyChi2:
            N=500
            for n in range(0,ak.sum(evtSelPho),N):
            
                sel_jets = tightJets[evtSelPho,:5][n:n+N]
                sel_leptons = tightLep[evtSelPho][n:n+N]
                sel_neutrinos = neutrinos[evtSelPho][n:n+N]
                sel_photons = tightPhotons[evtSelPho][n:n+N]
                sel_weights = weights[evtSelPho][n:n+N]
                sel_lepFlav = lepFlavor[evtSelPho][n:n+N]
            
                chi2, tStarMass, tStarHadMass, tStarLepMass, phoIsHadSide, phoPt = chi2_TT_tgty(sel_leptons, sel_neutrinos, sel_jets, sel_photons, btagCut=0.6321)

                output['tgty_tstarMass'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                              mass=ak.flatten(tStarMass))

                output['tgty_tstarHadMass'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                                 mass=ak.flatten(tStarHadMass))

                output['tgty_tstarLepMass'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                                 mass=ak.flatten(tStarLepMass))

                output['tgty_chi2'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                         chi2=ak.flatten(chi2))

                del chi2, tStarMass, tStarHadMass, tStarLepMass, phoIsHadSide, phoPt

        # print('after chi2 :',psutil.Process().memory_info())
        evtPreSel = selection.all(*('overlap','lepton','met'))
        output['nJet'].fill(dataset=dataset, lepFlavor=lepFlavor[evtPreSel], weight=weights[evtPreSel],
                            njet = (ak.num(tightJets))[evtPreSel])
        output['nBJet'].fill(dataset=dataset, lepFlavor=lepFlavor[evtPreSel], weight=weights[evtPreSel],
                             nbjet = (ak.num(tightJets[btag]))[evtPreSel])

        evtSel = selection.all(*('overlap','lepton','met','5jet',bjetCut))
        output['nPho'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[evtSel],
                            npho = (ak.num(tightPhotons))[evtSel])

        evtSel = selection.all(*('overlap','lepton','met','5jet',bjetCut,'photon'))
        output['phoPt'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[evtSel],
                             pt = ak.flatten(tightPhotons[evtSel,:1].pt))
        output['phoEta'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[evtSel],
                             eta = ak.flatten(tightPhotons[evtSel,:1].eta))

        for i in range(len(cutList)):
            cut = cutList[i]
            cutMask = selection.all(*cutList[:i+1])
            
            output['cutflow'].fill(dataset=dataset,
                                   lepFlavor=np.zeros_like(cutMask)[cutMask & ele_eventSelection],
                                   cut=cut,
                                   weight=weights[cutMask & ele_eventSelection],
                                  )

            output['cutflow'].fill(dataset=dataset,
                                   lepFlavor=np.ones_like(cutMask)[cutMask & muon_eventSelection],
                                   cut=cut,
                                   weight=weights[cutMask & muon_eventSelection],
                                  )
        del selection
        # print('processor end :',psutil.Process().memory_info())
        return output

    def postprocess(self, accumulator):
        return accumulator
        
