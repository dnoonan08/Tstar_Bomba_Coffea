import awkward as ak
import numpy as np

import pickle

from coffea import hist, util, processor
from coffea.analysis_tools import PackedSelection

from coffea.analysis_tools import Weights
#from coffea.btag_tools import BTagScaleFactor
#from coffea.lookup_tools import extractor, dense_lookup
#from coffea.jetmet_tools import CorrectedJetsFactory, JECStack


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
from .utils.PhotonScaleFactors import PhotonScaleFactors

from .utils.JetCorrections import JetCorrections

import time
import psutil

def bjetCut_none(array):
    return ak.num(array)>=0
def bjetCut_1b(array):
    return ak.num(array)>=1
def bjetCut_2b(array):
    return ak.num(array)>=2
def bjetCut_ex2b(array):
    return ak.num(array)==2

class TstarSelector(processor.ProcessorABC):
    def __init__(self,topChi2=False, tgtgChi2=False, tgtyChi2=True, mcTotals=None, bJetCut=None, mcPUDist=None):

        self.topChi2=topChi2
        self.tgtgChi2=tgtgChi2
        self.tgtyChi2=tgtyChi2
        self.mcTotals=mcTotals
        self.bJetCut=bJetCut
        self.mcPUDist=mcPUDist
        
        dataset_axis = hist.Cat("dataset","Dataset")
        systematic_axis = hist.Cat("systematic","Systematic")

        mTop_axis = hist.Bin("mass","$m_{T}$ [GeV]", 400, 0, 400)
        mW_axis = hist.Bin("mass","$m_{W}$ [GeV]", 200, 0, 200)
        mTstar_axis = hist.Bin("mass","$m_{T*}$ [GeV]", 500, 0, 2500)
        chi2_axis = hist.Bin('chi2','$\chi^2$',1000,0,1000)
        njet_axis = hist.Bin('njet','Number of Jets',20,0,20)
        nbjet_axis = hist.Bin('nbjet','Number of B-tagged Jets',20,0,20)
        npho_axis = hist.Bin('npho','Number of photons',20,0,20)
        phoPt_axis = hist.Bin('pt','$p_{T}(\gamma)$',300,0,1500)
        phoEta_axis = hist.Bin('eta','$\eta(\gamma)$',50,-2.5,2.5)

        leptonPt_axis = hist.Bin('pt','$p_{T}$',60,0,300)
        leptonEta_axis = hist.Bin('eta','$\eta$',50,-2.5,2.5)

        jetPt_axis = hist.Bin('pt','$p_{T}$',100,0,1000)
        jetEta_axis = hist.Bin('eta','$\eta$',50,-2.5,2.5)

        lep_axis = hist.Bin('lepFlavor','Lepton Flavor',2,0,2)
        phoHighPt_axis = hist.Bin('phopt','$p_{T}(\gamma)$',np.array([20,50,100,150,200,250,300,500,1000,1500]))

        lepPhoMass_axis = hist.Bin('mass','$m(\ell, \gamma)$',200,0,200)
        jetCat_axis = hist.Cat('jetCat','jetCategory')

        cut_axis = hist.Cat('cut','Cut')

        self._accumulator = processor.dict_accumulator({

            'leptonPt'  : hist.Hist("Counts", dataset_axis, lep_axis, leptonPt_axis, systematic_axis),
            'leptonEta' : hist.Hist("Counts", dataset_axis, lep_axis, leptonEta_axis, systematic_axis),

            'jetPt'  : hist.Hist("Counts", dataset_axis, lep_axis, jetPt_axis, systematic_axis),
            'jetEta' : hist.Hist("Counts", dataset_axis, lep_axis, jetEta_axis, systematic_axis),

            'leptonPt_PhoSel'  : hist.Hist("Counts", dataset_axis, lep_axis, leptonPt_axis, systematic_axis),
            'leptonEta_PhoSel' : hist.Hist("Counts", dataset_axis, lep_axis, leptonEta_axis, systematic_axis),

            'jetPt_PhoSel'  : hist.Hist("Counts", dataset_axis, lep_axis, jetPt_axis, systematic_axis),
            'jetEta_PhoSel' : hist.Hist("Counts", dataset_axis, lep_axis, jetEta_axis, systematic_axis),

            'lepPhoMass'  : hist.Hist("Counts", dataset_axis, lep_axis, lepPhoMass_axis, jetCat_axis, phoHighPt_axis, systematic_axis),

            'tt_topHadMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis, systematic_axis),
            'tt_topLepMass': hist.Hist("Counts", dataset_axis, lep_axis, mTop_axis, systematic_axis),
            'tt_chi2': hist.Hist("Counts", dataset_axis, lep_axis, chi2_axis, systematic_axis),

            'tgty_tstarMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis, systematic_axis, phoHighPt_axis),
            'tgty_tstarHadMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis, systematic_axis, phoHighPt_axis),
            'tgty_tstarLepMass': hist.Hist("Counts", dataset_axis, lep_axis, mTstar_axis, systematic_axis, phoHighPt_axis),
            'tgty_chi2': hist.Hist("Counts", dataset_axis, lep_axis, chi2_axis, systematic_axis),

            'phoPt': hist.Hist("Counts", dataset_axis, lep_axis, phoPt_axis, systematic_axis),
            'phoEta': hist.Hist("Counts", dataset_axis, lep_axis, phoEta_axis, systematic_axis),
            
            'phoPt_3j4j': hist.Hist("Counts", dataset_axis, lep_axis, phoPt_axis, systematic_axis),
            'phoEta_3j4j': hist.Hist("Counts", dataset_axis, lep_axis, phoEta_axis, systematic_axis),
            
            'nJet':hist.Hist("Counts",dataset_axis, lep_axis, njet_axis, systematic_axis),
            'nBJet':hist.Hist("Counts",dataset_axis, lep_axis, nbjet_axis, systematic_axis),
            'nPho':hist.Hist("Counts",dataset_axis, lep_axis, npho_axis, systematic_axis),

            'nJet_PhoSel':hist.Hist("Counts",dataset_axis, lep_axis, njet_axis, systematic_axis),
            'nBJet_PhoSel':hist.Hist("Counts",dataset_axis, lep_axis, nbjet_axis, systematic_axis),

            'cutflow':hist.Hist("Counts",dataset_axis, lep_axis, cut_axis, systematic_axis),
            })
        
        self.bCalc = BtagWeights()
        self.eleSF = ElectronScaleFactors()
        self.muSF  = MuonScaleFactors()
        self.phoSF = PhotonScaleFactors()

        self.jetCorr = JetCorrections()
        

    @property
    def accumulator(self):
        return self._accumulator

        
    def process(self, events):
        ## To access variables from the ntuples, use the "events" object
        ## The dataset name is part of events.metadata
        dataset, year = events.metadata['dataset']

        isMC = True
        if 'Data_Single' in dataset:
            isMC = False

        output = self.accumulator.identity()

        bjetCut='1b'
        bjetCutFunc=bjetCut_1b
        cutList = ['overlap','lepton', 'met', '3jet', '4jet', '5jet', '1b', 'photon', 'photon100', 'photon200', 'photon300', 'photon500']
        if not self.bJetCut is None:
            goodCuts=['','1b','2b','ex2b']
            if self.bJetCut in goodCuts:
                idx = goodCuts.index(self.bJetCut)+1
                bjetCut = self.bJetCut
                if bjetCut=='':     bjetCutFunc = bjetCut_none
                if bjetCut=='1b':   bjetCutFunc = bjetCut_1b
                if bjetCut=='2b':   bjetCutFunc = bjetCut_2b
                if bjetCut=='ex2b': bjetCutFunc = bjetCut_ex2b
                
                cutList = ['overlap','lepton', 'met', '3jet', '4jet', '5jet']+ goodCuts[1:idx] + ['photon', 'photon100', 'photon200', 'photon300', 'photon500']
            
        metCut = -1

        selection = PackedSelection()

        #Calculate the maximum pdgID of any of the particles in the GenPart history
        if isMC:
            idx = ak.to_numpy(ak.flatten(abs(events.GenPart.pdgId)))
            par = ak.to_numpy(ak.flatten(events.GenPart.genPartIdxMother))
            num = ak.to_numpy(ak.num(events.GenPart.pdgId))
            maxParentFlatten = maxHistoryPDGID(idx,par,num)
            # events["GenPart","maxParent"] = ak.broadcast_arrays(events.GenPart.pdgId),1.)*maxParentFlatten
            events["GenPart","maxParent"] = ak.unflatten(maxParentFlatten, num)
            overlap = overlapFilter(dataset, events)

        #object selection
        tightMuons, looseMuons = selectMuons(events.Muon)
        tightElectrons, looseElectrons = selectElectrons(events.Electron, events.Photon)
        tightPhotons = selectPhotons(events.Photon, tightMuons, tightElectrons)

        if isMC:
            corr_jets, metVars = self.jetCorr.GetCorrectedJetMet(events, year)
        else:
            corr_jets = events.Jet
            metVars = {'nominal':events.MET}
        
        tightJetVars = {}
        btagVars = {}

        btagWP = {2016: 0.6321, 2017: 0.4941, 2018: 0.4184}

        tightJetVars['nominal'], btagVars['nominal'] = selectJets(corr_jets, tightMuons, tightElectrons, tightPhotons, btagWP[year])
        jet_variations = [x for x in corr_jets.fields if x.startswith('JE')]
        for k in jet_variations:
            tightJetVars[f'{k}_up'], btagVars[f'{k}_up'] = selectJets(corr_jets[k].up, tightMuons, tightElectrons, tightPhotons, btagWP[year])
            tightJetVars[f'{k}_down'], btagVars[f'{k}_down'] = selectJets(corr_jets[k].down, tightMuons, tightElectrons, tightPhotons, btagWP[year])



        #for now, just keep the nominal ones:
        tightJets = tightJetVars['nominal']
        btag = btagVars['nominal']
        met = metVars['nominal']

        # tightJets, btag = selectJets(events.Jet, tightMuons, tightElectrons, tightPhotons)
        # met = events.MET

        if not isMC:
            overlap = np.ones(len(events),dtype=bool)
        

        filters = (events.Flag.BadPFMuonFilter & 
                   events.Flag.EcalDeadCellTriggerPrimitiveFilter &
                   events.Flag.HBHENoiseFilter & 
                   events.Flag.HBHENoiseIsoFilter & 
                   events.Flag.globalSuperTightHalo2016Filter & 
                   events.Flag.goodVertices)
        pv = events.PV
        goodPV = ((pv.ndof>4) & ((pv.x**2 + pv.y**2)<=2) & abs(pv.z<=24))

        #event selection
        selection.add('overlap', overlap & filters & goodPV)
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

        if not isMC:
            if 'SingleEle' in dataset:
                muon_eventSelection = np.zeros_like(muon_eventSelection,dtype=bool)
            if 'SingleMu' in dataset:
                ele_eventSelection = np.zeros_like(ele_eventSelection,dtype=bool)

        selection.add('lepton', muon_eventSelection | ele_eventSelection)
        
        selection.add('2jet',(ak.num(tightJets)>=2))
        selection.add('3jet',(ak.num(tightJets)>=3))
        selection.add('4jet',(ak.num(tightJets)>=4))
        selection.add('5jet',(ak.num(tightJets)>=5))
        selection.add('6jet',(ak.num(tightJets)>=6))
        selection.add('3j4jCR',(ak.num(tightJets)==3) | (ak.num(tightJets)==4))
        selection.add('',np.ones(len(events),dtype=bool))
        selection.add('1b',(ak.num(tightJets[btag])>=1))
        selection.add('2b',(ak.num(tightJets[btag])>=2))
        selection.add('ex2b',(ak.num(tightJets[btag])==2))
        selection.add('met',(met.pt>metCut))
        if isMC:
            selection.add('photon',   (ak.num(tightPhotons)==1))
            selection.add('photon100',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>100,axis=1)))
            selection.add('photon200',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>200,axis=1)))
            selection.add('photon300',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>300,axis=1)))
            selection.add('photon500',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>500,axis=1)))
        else:
            #blind to photon region in data for now
            blindData = ak.any(tightPhotons.pt>50,axis=1) & (ak.num(tightJets)>=5)
            selection.add('photon',(ak.num(tightPhotons)==1) & ~blindData)
            selection.add('photon100',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>100,axis=1) & ~blindData))
            selection.add('photon200',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>200,axis=1) & ~blindData))
            selection.add('photon300',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>300,axis=1) & ~blindData))
            selection.add('photon500',(ak.num(tightPhotons)==1) & (ak.any(tightPhotons.pt>500,axis=1) & ~blindData))

        
        
        lepFlavor = np.zeros(len(events))-1.
        lepFlavor[ele_eventSelection] = 0
        lepFlavor[muon_eventSelection] = 1

        weightCollection = Weights(len(events))
        if isMC:
            lumiWeight =events.Generator.weight*getLumiWeight(dataset,self.mcTotals, year)
            weightCollection.add('lumiWeight',events.Generator.weight*getLumiWeight(dataset,self.mcTotals, year))

            puWeight = getPileupReweight(events.Pileup.nTrueInt, self.mcPUDist, dataset, year)
            weightCollection.add('puWeight', weight=puWeight[0], weightUp=puWeight[1], weightDown=puWeight[2])

            btagsWeight = self.bCalc.getBtagWeights(tightJets, btag, year, dataset)
            weightCollection.add('btagSF',weight=btagsWeight[0], weightUp=btagsWeight[1], weightDown=btagsWeight[2])

            muWeight = self.muSF.getMuonScaleFactors(tightMuons, year, 'total')
            weightCollection.add('muSF',weight=muWeight[0], weightUp=muWeight[1], weightDown=muWeight[2])

            eleWeight = self.eleSF.getEleScaleFactors(tightElectrons, year, 'total')
            weightCollection.add('eleSF',weight=eleWeight[0], weightUp=eleWeight[1], weightDown=eleWeight[2])

            phoWeight = self.phoSF.getPhoScaleFactors(tightPhotons, year, 'total')
            weightCollection.add('phoSF',weight=phoWeight[0], weightUp=phoWeight[1], weightDown=phoWeight[2])

            weightCollection.add('btagSF',weight=np.ones(len(events)),weightUp=np.ones(len(events)),weightDown=np.ones(len(events)))
            weightCollection.add('eleSF',weight=np.ones(len(events)),weightUp=np.ones(len(events)),weightDown=np.ones(len(events)))
            weightCollection.add('muSF',weight=np.ones(len(events)),weightUp=np.ones(len(events)),weightDown=np.ones(len(events)))
            weightCollection.add('phoSF',weight=np.ones(len(events)),weightUp=np.ones(len(events)),weightDown=np.ones(len(events)))


            if ak.all(ak.num(events.PSWeight)==4):
                renormWeight = ak.where( (events.Generator.weight==0) | (events.LHEWeight.originalXWGTUP==0), 
                                         ak.sum(events.PSWeight)/4., 
                                         events.Generator.weight/events.LHEWeight.originalXWGTUP)

                renormWeight = ak.where( renormWeight==0, 1, renormWeight)

                weightCollection.add('ISR', weight=np.ones(len(events)), weightUp=events.PSWeight[:,2]/renormWeight, weightDown=events.PSWeight[:,0]/renormWeight)
                weightCollection.add('FSR', weight=np.ones(len(events)), weightUp=events.PSWeight[:,3]/renormWeight, weightDown=events.PSWeight[:,1]/renormWeight)
            else:
                weightCollection.add('ISR', weight=np.ones(len(events)), weightUp=np.ones(len(events)), weightDown=np.ones(len(events)))
                weightCollection.add('FSR', weight=np.ones(len(events)), weightUp=np.ones(len(events)), weightDown=np.ones(len(events)))


            if ak.all(ak.num(events.LHEScaleWeight)==9):
                nom = events.LHEScaleWeight[:,4]
                scales = events.LHEScaleWeight[:,[0,1,3,5,7,8]]
                q2Up = ak.max(scales,axis=1)/nom
                q2Down = ak.max(scales,axis=1)/nom    
            elif ak.all(ak.num(events.LHEScaleWeight)==9):
                scales = events.LHEScaleWeight[:,[0,5,15,24,34,39]]
                q2Up = ak.max(scales,axis=1)
                q2Down = ak.max(scales,axis=1)
            else:
                q2Up = np.ones(len(events))
                q2Down = np.ones(len(events))
            weightCollection.add('Q2', weight=np.ones(len(events)), weightUp=q2Up, weightDown=q2Down)

            pdfUnc = ak.std(events.LHEPdfWeight,axis=1)/ak.mean(events.LHEPdfWeight,axis=1)
            weightCollection.add('PDF', weight=np.ones(len(events)), weightUp=np.ones(len(events)) + pdfUnc, weightDown=np.ones(len(events)) - pdfUnc)

            systList = ['noWeight','nominal',
                        'puWeightUp',  'muSFUp',  'eleSFUp',  'btagSFUp',  'ISRUp',  'FSRUp',  'Q2Up',  'PDFUp',
                        'puWeightDown','muSFDown','eleSFDown','btagSFDown','ISRDown','FSRDown','Q2Down','PDFDown']
        else:
            systList=['nominal']


        tightLep = ak.with_name(ak.concatenate([tightElectrons, tightMuons], axis=1), "PtEtaPhiMCandidate")

#        weights=weightCollection.weight(None)
        weights = {}
        for s in systList:
            if s=='noWeight':
                weights[s] = np.ones(len(events))
            elif s in ['nominal']:
                weights[s] = weightCollection.weight(None)
            else:
                weights[s] = weightCollection.weight(s)

        # if (ak.sum(evtSel)>0) and self.topChi2:
        if self.topChi2:
            for jetK in tightJetVars:
                _jet = tightJetVars[jetK]
                _tag = btagVars[jetK]
                _met = metVars[jetK]
                evtSel = selection.all(*('overlap','lepton')) & (ak.num(_jet)>=3) & bjetCutFunc(_jet[_tag]) & (_met.pt>metCut)

                if ak.sum(evtSel)==0: continue

                ##construct awkward array four vectors for the neutrino, with each of the hypothesis pz
                nu_pz = MetPz(tightElectrons, tightMuons, _met, solutionType=-1)
                neutrinos = getNeutrinoArrays(nu_pz, _met)

                sel_jets = _jet[evtSel]
                sel_leptons = tightLep[evtSel]
                sel_neutrinos = neutrinos[evtSel]
                sel_lepFlav = lepFlavor[evtSel]

                chi2, topHadMass, topLepMass = chi2_Top(sel_leptons, sel_neutrinos, sel_jets, 4, btagCut=btagWP[year])

                #if using nominal jet correction, loop over other systematics 
                jetK_systList = systList if jetK is 'nominal' else ['nominal']
                for syst in systList:
                    sel_weights = weights[syst][evtSel]
                    output['tt_topHadMass'].fill(dataset=dataset, lepFlavor=sel_lepFlavor, weight=sel_weights,
                                                 mass=ak.flatten(topHadMass),
                                                 systematic=syst if jetK is 'nominal' else jetK)

                    output['tt_topLepMass'].fill(dataset=dataset, lepFlavor=sel_lepFlavor, weight=sel_weights,
                                                 mass=ak.flatten(topLepMass),
                                                 systematic=syst if jetK is 'nominal' else jetK)
                    output['tt_chi2'].fill(dataset=dataset, lepFlavor=sel_lepFlavor, weight=sel_weights,
                                           chi2=ak.flatten(chi2),
                                           systematic=syst if jetK is 'nominal' else jetK)



        # evtSelPho = selection.all(*('overlap','lepton','5jet',bjetCut,'met','photon'))
        # if ak.sum(evtSelPho)>0 and self.tgtyChi2:

        if self.tgtyChi2:
            for jetK in tightJetVars:
                _jet = tightJetVars[jetK]
                _tag = btagVars[jetK]
                _met = metVars[jetK]
                evtSelPho = selection.all(*('overlap','lepton','photon')) & (ak.num(_jet)>=5) & bjetCutFunc(_jet[_tag]) & (_met.pt>metCut)
                
                if ak.sum(evtSelPho)==0: continue

                ##construct awkward array four vectors for the neutrino, with each of the hypothesis pz
                nu_pz = MetPz(tightElectrons, tightMuons, _met, solutionType=-1)
                neutrinos = getNeutrinoArrays(nu_pz, _met)

                N=500
                for n in range(0,ak.sum(evtSelPho),N):
            
                    sel_jets = _jet[evtSelPho,:5][n:n+N]
                    sel_leptons = tightLep[evtSelPho][n:n+N]
                    sel_neutrinos = neutrinos[evtSelPho][n:n+N]
                    sel_photons = tightPhotons[evtSelPho][n:n+N]
                    sel_lepFlav = lepFlavor[evtSelPho][n:n+N]
            
                    chi2, tStarMass, tStarHadMass, tStarLepMass, phoIsHadSide, phoPt = chi2_TT_tgty(sel_leptons, sel_neutrinos, sel_jets, sel_photons, btagCut=btagWP[year])

                    #if using nominal jet correction, loop over other systematics 
                    jetK_systList = systList if jetK is 'nominal' else ['nominal']
                    for syst in systList:
                        sel_weights = weights[syst][evtSelPho][n:n+N]
                        output['tgty_tstarMass'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                                      mass=ak.flatten(tStarMass),
                                                      phopt=ak.flatten(phoPt),
                                                      systematic=syst if jetK is 'nominal' else jetK)

                        output['tgty_tstarHadMass'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                                         mass=ak.flatten(tStarHadMass),
                                                         phopt=ak.flatten(phoPt),
                                                         systematic=syst if jetK is 'nominal' else jetK)

                        output['tgty_tstarLepMass'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                                         mass=ak.flatten(tStarLepMass),
                                                         phopt=ak.flatten(phoPt),
                                                         systematic=syst if jetK is 'nominal' else jetK)

                        output['tgty_chi2'].fill(dataset=dataset, lepFlavor=sel_lepFlav, weight=sel_weights,
                                                 chi2=ak.flatten(chi2),
                                                 systematic=syst if jetK is 'nominal' else jetK)

                    del chi2, tStarMass, tStarHadMass, tStarLepMass, phoIsHadSide, phoPt

        for syst in systList:
            evtPreSel = selection.all(*('overlap','lepton','met'))
            output['nJet'].fill(dataset=dataset, lepFlavor=lepFlavor[evtPreSel], weight=weights[syst][evtPreSel],
                                njet = (ak.num(tightJets))[evtPreSel],
                                systematic=syst)
            output['nBJet'].fill(dataset=dataset, lepFlavor=lepFlavor[evtPreSel], weight=weights[syst][evtPreSel],
                                 nbjet = (ak.num(tightJets[btag]))[evtPreSel],
                                 systematic=syst)

            evtPreSel = selection.all(*('overlap','lepton','met','photon'))
            output['nJet_PhoSel'].fill(dataset=dataset, lepFlavor=lepFlavor[evtPreSel], weight=weights[syst][evtPreSel],
                                       njet = (ak.num(tightJets))[evtPreSel],
                                       systematic=syst)
            output['nBJet_PhoSel'].fill(dataset=dataset, lepFlavor=lepFlavor[evtPreSel], weight=weights[syst][evtPreSel],
                                        nbjet = (ak.num(tightJets[btag]))[evtPreSel],
                                        systematic=syst)

            evtSel = selection.all(*('overlap','lepton','met','5jet',bjetCut))
            output['nPho'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                npho = (ak.num(tightPhotons))[evtSel],
                                systematic=syst)

            output['leptonPt'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                    pt = ak.flatten(tightLep[evtSel].pt),
                                    systematic=syst)
            output['leptonEta'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                     eta = ak.flatten(tightLep[evtSel].eta),
                                     systematic=syst)

            output['jetPt'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                 pt = ak.flatten(tightJets[evtSel,:1].pt),
                                 systematic=syst)
            output['jetEta'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                  eta = ak.flatten(tightJets[evtSel,:1].eta),
                                  systematic=syst)


            evtSel = selection.all(*('overlap','lepton','met','5jet',bjetCut,'photon'))
            output['phoPt'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                 pt = ak.flatten(tightPhotons[evtSel,:1].pt),
                                 systematic=syst)
            output['phoEta'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                  eta = ak.flatten(tightPhotons[evtSel,:1].eta),
                                  systematic=syst)
            
            output['leptonPt_PhoSel'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                           pt = ak.flatten(tightLep[evtSel].pt),
                                           systematic=syst)
            output['leptonEta_PhoSel'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                            eta = ak.flatten(tightLep[evtSel].eta),
                                            systematic=syst)

            output['jetPt_PhoSel'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                        pt = ak.flatten(tightJets[evtSel,:1].pt),
                                        systematic=syst)
            output['jetEta_PhoSel'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                         eta = ak.flatten(tightJets[evtSel,:1].eta),
                                         systematic=syst)

            evtSel = selection.all(*('overlap','lepton','met','3j4jCR',bjetCut,'photon'))
            output['phoPt_3j4j'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                      pt = ak.flatten(tightPhotons[evtSel,:1].pt),
                                      systematic=syst)
            output['phoEta_3j4j'].fill(dataset=dataset, lepFlavor=lepFlavor[evtSel], weight=weights[syst][evtSel],
                                       eta = ak.flatten(tightPhotons[evtSel,:1].eta),
                                       systematic=syst)

            #separate jet categories
            jetCatCuts = {'2j0t':(ak.num(tightJets)==2) & (ak.num(tightJets[btag])==0),
                          '3j0t':(ak.num(tightJets)==3) & (ak.num(tightJets[btag])==0),
                          '4j0t':(ak.num(tightJets)==4) & (ak.num(tightJets[btag])==0),
                          '5j0t':(ak.num(tightJets)>=5) & (ak.num(tightJets[btag])==0),
                          '2j1t':(ak.num(tightJets)==2) & (ak.num(tightJets[btag])>=1),
                          '3j1t':(ak.num(tightJets)==3) & (ak.num(tightJets[btag])>=1),
                          '4j1t':(ak.num(tightJets)==4) & (ak.num(tightJets[btag])>=1),
                          '5j1t':(ak.num(tightJets)>=5) & (ak.num(tightJets[btag])>=1)}
            for jetCat in jetCatCuts:
                evtSel = selection.all(*('overlap','lepton','met','photon')) & jetCatCuts[jetCat]
                if ak.sum(evtSel)==0: continue
                output['lepPhoMass'].fill(dataset=dataset, lepFlavor = lepFlavor[evtSel], weight=weights[syst][evtSel], systematic=syst, 
                                          jetCat=jetCat,
                                          mass = ak.flatten( (tightLep[evtSel] + tightPhotons[evtSel]).mass), 
                                          phopt = ak.flatten( tightPhotons[evtSel].pt))

            
        for i in range(len(cutList)):
            cut = cutList[i]
            cutMask = selection.all(*cutList[:i+1])
            for syst in systList:
                output['cutflow'].fill(dataset=dataset,
                                       lepFlavor=np.zeros_like(cutMask)[cutMask & ele_eventSelection],
                                       cut=cut,
                                       weight=weights[syst][cutMask & ele_eventSelection],
                                       systematic=syst
                                   )

                output['cutflow'].fill(dataset=dataset,
                                       lepFlavor=np.ones_like(cutMask)[cutMask & muon_eventSelection],
                                       cut=cut,
                                       weight=weights[syst][cutMask & muon_eventSelection],
                                       systematic=syst
                                   )

        # for i in range(len(cutList)):
        #     cut = cutList[i]
        #     cutMask = selection.all(*cutList[:i+1]) & muon_eventSelection
        #     print(cut, ak.sum(cutMask))
            
        return output

    def postprocess(self, accumulator):
        return accumulator
        
