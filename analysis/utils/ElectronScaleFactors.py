import os.path
cwd = os.path.dirname(__file__)


from coffea.lookup_tools import extractor, dense_lookup

import awkward as ak

class ElectronScaleFactors:
    def __init__(self):

        ext = extractor()

        ext.add_weight_sets([f"id_2016 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_id_2016.root",
                             f"reco_2016 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_reco_2016.root",
                             f"trig_2016 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_trig_2016.root",
                             f"id_2017 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_id_2017.root",
                             f"reco_2017 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_reco_2017.root",
                             f"trig_2017 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_trig_2017.root",
                             f"id_2018 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_id_2018.root",
                             f"reco_2018 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_reco_2018.root",
                             f"trig_2018 EGamma_SF2D {cwd}/ScaleFactors/Electrons/ele_trig_2018.root",
                         ])

        ext.add_weight_sets([f"id_2016_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_id_2016.root",
                             f"reco_2016_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_reco_2016.root",
                             f"trig_2016_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_trig_2016.root",
                             f"id_2017_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_id_2017.root",
                             f"reco_2017_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_reco_2017.root",
                             f"trig_2017_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_trig_2017.root",
                             f"id_2018_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_id_2018.root",
                             f"reco_2018_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_reco_2018.root",
                             f"trig_2018_err EGamma_SF2D_error {cwd}/ScaleFactors/Electrons/ele_trig_2018.root",
                         ])

        ext.finalize()
        self.lookups = ext.make_evaluator()

    def getEleScaleFactors(self, electrons, year, syst='Total'):
        
        id_SF = self.lookups[f'id_{year}'](electrons.eta, electrons.pt)
        reco_SF = self.lookups[f'id_{year}'](electrons.eta, electrons.pt)
        trig_SF = self.lookups[f'id_{year}'](electrons.eta, electrons.pt)

        id_SF_err = self.lookups[f'id_{year}_err'](electrons.eta, electrons.pt)
        reco_SF_err = self.lookups[f'id_{year}_err'](electrons.eta, electrons.pt)
        trig_SF_err = self.lookups[f'id_{year}_err'](electrons.eta, electrons.pt)

        sf = id_SF * reco_SF * trig_SF
        if syst.lower()=='total':
            sfUp = (id_SF + id_SF_err) * (reco_SF + reco_SF_err) * (trig_SF + trig_SF_err)
            sfDo = (id_SF - id_SF_err) * (reco_SF - reco_SF_err) * (trig_SF - trig_SF_err)
        elif syst.lower()=='id':
            sfUp = (id_SF + id_SF_err) * (reco_SF + reco_SF_err) * (trig_SF)
            sfDo = (id_SF - id_SF_err) * (reco_SF - reco_SF_err) * (trig_SF)
        elif syst.lower()=='trig':
            sfUp = (id_SF) * (reco_SF) * (trig_SF + trig_SF_err)
            sfDo = (id_SF) * (reco_SF) * (trig_SF - trig_SF_err)
        else:
            print('Unrecognized syst level, returning no variation')
            sfUp = sf
            sfDo = sf

        return ak.prod(sf,axis=-1), ak.prod(sfUp,axis=-1), ak.prod(sfDo,axis=-1)

            
