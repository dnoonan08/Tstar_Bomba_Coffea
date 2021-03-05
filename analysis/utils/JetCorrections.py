import awkward as ak
import numpy as np

import os.path
cwd = os.path.dirname(__file__)

from coffea.lookup_tools import extractor, dense_lookup
from coffea.jetmet_tools import CorrectedJetsFactory, JECStack

from .objectSelection import selectJets



class JetCorrections:
    def __init__(self):

        Jetext = extractor()
        Jetext.add_weight_sets([
            f"* * {cwd}/ScaleFactors/JEC/Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.jec.txt",
            f"* * {cwd}/ScaleFactors/JEC/Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.jec.txt",
            f"* * {cwd}/ScaleFactors/JEC/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.junc.txt",
            f"* * {cwd}/ScaleFactors/JEC/Summer16_25nsV1_MC_PtResolution_AK4PFchs.jr.txt",
            f"* * {cwd}/ScaleFactors/JEC/Summer16_25nsV1_MC_SF_AK4PFchs.jersf.txt",
            
            f"* * {cwd}/ScaleFactors/JEC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.jec.txt",
            f"* * {cwd}/ScaleFactors/JEC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.jec.txt",
            f"* * {cwd}/ScaleFactors/JEC/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.junc.txt",
            f"* * {cwd}/ScaleFactors/JEC/Fall17_V3_MC_PtResolution_AK4PFchs.jr.txt",
            f"* * {cwd}/ScaleFactors/JEC/Fall17_V3_MC_SF_AK4PFchs.jersf.txt",

            f"* * {cwd}/ScaleFactors/JEC/Autumn18_V19_MC_L1FastJet_AK4PFchs.jec.txt",
            f"* * {cwd}/ScaleFactors/JEC/Autumn18_V19_MC_L2Relative_AK4PFchs.jec.txt",
            f"* * {cwd}/ScaleFactors/JEC/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.junc.txt",
            f"* * {cwd}/ScaleFactors/JEC/Autumn18_V7b_MC_PtResolution_AK4PFchs.jr.txt",
            f"* * {cwd}/ScaleFactors/JEC/Autumn18_V7b_MC_SF_AK4PFchs.jr.txt",
        ])
        Jetext.finalize()
        Jetevaluator = Jetext.make_evaluator()

        jec_names = ['Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs',
                     'Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs', 
                     'RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs_Total',
                     'Summer16_25nsV1_MC_PtResolution_AK4PFchs',
                     'Summer16_25nsV1_MC_SF_AK4PFchs']
        jec_inputs = {name: Jetevaluator[name] for name in jec_names}
        jec_stack = JECStack(jec_inputs)

        name_map = jec_stack.blank_name_map
        name_map['JetPt'] = 'pt'
        name_map['JetMass'] = 'mass'
        name_map['JetEta'] = 'eta'
        name_map['JetA'] = 'area'

        name_map['ptGenJet'] = 'pt_gen'
        name_map['ptRaw'] = 'pt_raw'
        name_map['massRaw'] = 'mass_raw'
        name_map['Rho'] = 'rho'

        self.jet_factory={}
        self.jet_factory[2016] = CorrectedJetsFactory(name_map, jec_stack)

        # jec_names = ['Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs',
        #              'Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs',
        #              'RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs_Total',
        #              'Fall17_V3_MC_PtResolution_AK4PFchs',
        #              'Fall17_V3_MC_SF_AK4PFchs']
        # jec_inputs = {name: Jetevaluator[name] for name in jec_names}
        # jec_stack = JECStack(jec_inputs)

        # name_map = jec_stack.blank_name_map
        # name_map['JetPt'] = 'pt'
        # name_map['JetMass'] = 'mass'
        # name_map['JetEta'] = 'eta'
        # name_map['JetA'] = 'area'

        # name_map['ptGenJet'] = 'pt_gen'
        # name_map['ptRaw'] = 'pt_raw'
        # name_map['massRaw'] = 'mass_raw'
        # name_map['Rho'] = 'rho'

        # self.jet_factory[2017] = CorrectedJetsFactory(name_map, jec_stack)


        # jec_names = ['Autumn18_V19_MC_L1FastJet_AK4PFchs',
        #              'Autumn18_V19_MC_L2Relative_AK4PFchs',
        #              'RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs_Total',
        #              'Autumn18_V7b_MC_PtResolution_AK4PFchs',
        #              'Autumn18_V7b_MC_SF_AK4PFchs']
        # jec_inputs = {name: Jetevaluator[name] for name in jec_names}
        # jec_stack = JECStack(jec_inputs)

        # name_map = jec_stack.blank_name_map
        # name_map['JetPt'] = 'pt'
        # name_map['JetMass'] = 'mass'
        # name_map['JetEta'] = 'eta'
        # name_map['JetA'] = 'area'

        # name_map['ptGenJet'] = 'pt_gen'
        # name_map['ptRaw'] = 'pt_raw'
        # name_map['massRaw'] = 'mass_raw'
        # name_map['Rho'] = 'rho'

        # self.jet_factory[2018] = CorrectedJetsFactory(name_map, jec_stack)


    def selectCorrectedJetMet(self, events, muons, electrons, photons, year):
        events["Jet","pt_raw"]=(1 - events.Jet.rawFactor)*events.Jet.pt
        events["Jet","mass_raw"]=(1 - events.Jet.rawFactor)*events.Jet.mass
        events["Jet","pt_gen"]=ak.values_astype(ak.fill_none(events.Jet.matched_gen.pt, 0), np.float32)
        events["Jet","rho"]= ak.broadcast_arrays(events.fixedGridRhoFastjetAll, events.Jet.pt)[0]

        events_cache = events.caches[0]
        corrected_jets = self.jet_factory[year].build(events.Jet, lazy_cache=events_cache)

        #make dictionaries to store jet variations, btags, and met differences
        tightJets = {}
        btag = {}
        met = {}

        #get JER and JES uncertainties types
        variations = [x for x in corrected_jets.fields if x.startswith('JE')]
        
        tightJets['nominal'], btag['nominal'] = selectJets(corrected_jets, muons, electrons, photons)
        met['nominal'] = events.MET

        nomJetSum = self.jetSum(corrected_jets)
        
        for k in variations:
            tightJets[f'{k}_up'], btag[f'{k}_up'] = selectJets(corrected_jets[k].up, muons, electrons, photons)
            met[f'{k}_up'] = self.metCorrection(corrected_jets[k].up, nomJetSum, events.MET)
            tightJets[f'{k}_down'], btag[f'{k}_down'] = selectJets(corrected_jets[k].down, muons, electrons, photons)
            met[f'{k}_down'] = self.metCorrection(corrected_jets[k].down, nomJetSum, events.MET)
        
        return tightJets, btag, met
        
    #calculate the sum of jets, to be used for updating met pre/post jet systematic shifts
    def jetSum(self, jet):
        px = ak.sum(jet.px,axis=-1)
        py = ak.sum(jet.py,axis=-1)
        pt = (px**2 + py**2)**0.5
        phi = np.arccos(px/pt)*np.sign(py)
        metChange = ak.zip(
            {
                "pt": pt,
                "eta": pt*0,
                "phi": phi,
                "mass": pt*0,
            },
            with_name="PtEtaPhiMLorentzVector",
        )
        return metChange

    # calculate a new MET vector, based on change to the jet energies
    def metCorrection(self, jet, nom, met):
        jetVar = self.jetSum(jet)
        
        return met + nom + -1*jetVar
