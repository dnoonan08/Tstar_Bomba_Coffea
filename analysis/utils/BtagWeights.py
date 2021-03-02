import os.path
cwd = os.path.dirname(__file__)

from coffea.btag_tools import BTagScaleFactor

from coffea.lookup_tools import extractor, dense_lookup

import awkward as ak

class BtagWeights:
    def __init__(self):

        ext = extractor()
        for y in [2016,2017,2018]:
            ext.add_weight_sets([f"Eff{y}_top_l Top_l_efficiency {cwd}/ScaleFactors/Btag/btag_efficiencies_{y}.root"])
            ext.add_weight_sets([f"Eff{y}_top_c Top_c_efficiency {cwd}/ScaleFactors/Btag/btag_efficiencies_{y}.root"])
            ext.add_weight_sets([f"Eff{y}_top_b Top_b_efficiency {cwd}/ScaleFactors/Btag/btag_efficiencies_{y}.root"])
            ext.add_weight_sets([f"Eff{y}_other_l Other_l_efficiency {cwd}/ScaleFactors/Btag/btag_efficiencies_{y}.root"])
            ext.add_weight_sets([f"Eff{y}_other_c Other_c_efficiency {cwd}/ScaleFactors/Btag/btag_efficiencies_{y}.root"])
            ext.add_weight_sets([f"Eff{y}_other_b Other_b_efficiency {cwd}/ScaleFactors/Btag/btag_efficiencies_{y}.root"])
        ext.finalize()
        self.efficiencies = ext.make_evaluator()
        self.bJetSF = {2016 : BTagScaleFactor(f"{cwd}/ScaleFactors/Btag/DeepCSV_2016LegacySF_V1.csv","MEDIUM"),
                       2017 : BTagScaleFactor(f"{cwd}/ScaleFactors/Btag/DeepCSV_94XSF_V3_B_F.csv","MEDIUM"),
                       2018 : BTagScaleFactor(f"{cwd}/ScaleFactors/Btag/DeepCSV_102XSF_V1.csv","MEDIUM")
                      }
        
    def getBtagWeights(self, jets, btagged, year, dataset):

        bJetSF      = self.bJetSF[year]('central', jets.hadronFlavour, abs(jets.eta), jets.pt)
        bJetSF_up   = self.bJetSF[year]('up',      jets.hadronFlavour, abs(jets.eta), jets.pt)
        bJetSF_down = self.bJetSF[year]('down',    jets.hadronFlavour, abs(jets.eta), jets.pt)

        ## mc efficiency lookup, data efficiency is eff* scale factor
        if 'TTGamma' in dataset or 'TTbar' in dataset or 'TstarTstar' in dataset or 'ST_' in dataset:
            effType='top'
        else:
            effType='other'

        btagEfficiencies = ak.where(jets.hadronFlavour==5, 
                                    self.efficiencies[f'Eff{year}_{effType}_b'](jets.pt,abs(jets.eta)),
                                    ak.where(jets.hadronFlavour==4, 
                                             self.efficiencies[f'Eff{year}_{effType}_c'](jets.pt,abs(jets.eta)),
                                             self.efficiencies[f'Eff{year}_{effType}_l'](jets.pt,abs(jets.eta))
                                         )
                                )

        btagEfficienciesData = btagEfficiencies*bJetSF
        btagEfficienciesData_up   = btagEfficiencies*bJetSF_up
        btagEfficienciesData_down = btagEfficiencies*bJetSF_down

        ##probability is the product of all efficiencies of tagged jets, times product of 1-eff for all untagged jets
        ## https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
        pMC          = ak.prod(btagEfficiencies[btagged], axis=-1)           * ak.prod((1.-btagEfficiencies[~(btagged)]), axis=-1) 
        pData        = ak.prod(btagEfficienciesData[btagged], axis=-1)       * ak.prod((1.-btagEfficienciesData[~(btagged)]),axis=-1)
        pData_up   = ak.prod(btagEfficienciesData_up[btagged], axis=-1)  * ak.prod((1.-btagEfficienciesData_up[~(btagged)]),axis=-1)
        pData_down = ak.prod(btagEfficienciesData_down[btagged],axis=-1) * ak.prod((1.-btagEfficienciesData_down[~(btagged)]),axis=-1)
    
        pMC = ak.where(pMC==0,1,pMC)
        btagWeight = pData/pMC
        btagWeight_up = pData_up/pMC
        btagWeight_down = pData_down/pMC

        return btagWeight, btagWeight_up, btagWeight_down
