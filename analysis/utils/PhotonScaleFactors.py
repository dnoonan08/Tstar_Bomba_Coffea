import os.path
cwd = os.path.dirname(__file__)


from coffea.lookup_tools import extractor, dense_lookup

import awkward as ak

class PhotonScaleFactors:
    def __init__(self):

        ext = extractor()

        directory=f"{cwd}/ScaleFactors/Photons"
        inputs = [f"phoId_2016 phoID_EtaPt_2016 {directory}/photonSF_Full.root",
                  f"phoId_err_2016 phoID_EtaPt_2016_error {directory}/photonSF_Full.root",
                  f"phoVeto_2016 phoVeto_AbsEtaPt_2016 {directory}/photonSF_Full.root",
                  f"phoVeto_err_2016 phoVeto_AbsEtaPt_2016_error {directory}/photonSF_Full.root",

                  f"phoId_2017 phoID_EtaPt_2017 {directory}/photonSF_Full.root", 
                  f"phoId_err_2017 phoID_EtaPt_2017_error {directory}/photonSF_Full.root", 

                  f"phoId_2018 phoID_EtaPt_2018 {directory}/photonSF_Full.root", 
                  f"phoId_err_2018 phoID_EtaPt_2018_error {directory}/photonSF_Full.root", 
                  f"phoVeto_2018 phoVeto_PtAbsEta_2018 {directory}/photonSF_Full.root", 
                  f"phoVeto_err_2018 phoVeto_unc_PtAbsEta_2018 {directory}/photonSF_Full.root"] 


        ext.add_weight_sets(inputs)

        ext.finalize()
        self.lookups = ext.make_evaluator()

    def getPhoScaleFactors(self, photons, year, syst='Total'):
        
        if year==2016:
            id_SF = self.lookups['phoId_2016'](photons.eta, photons.pt)
            id_SF_err = self.lookups['phoId_err_2016'](photons.eta, photons.pt)
            veto_SF = self.lookups['phoVeto_2016'](abs(photons.eta), photons.pt)
            veto_SF_err = self.lookups['phoVeto_err_2016'](abs(photons.eta), photons.pt)
        if year==2017:
            id_SF = self.lookups['phoId_2017'](photons.eta, photons.pt)
            id_SF_err = self.lookups['phoId_err_2017'](photons.eta, photons.pt)
            veto_SF = ak.where(abs(photons.eta)<1.5, 0.967260, 0.914782)
            veto_SF_err = ak.where(abs(photons.eta)<1.5, 0.004449, 0.012542)
        if year==2018:
            id_SF = self.lookups['phoId_2018'](photons.eta, photons.pt)
            id_SF_err = self.lookups['phoId_err_2018'](photons.eta, photons.pt)
            veto_SF = self.lookups['phoVeto_2018'](photons.pt, abs(photons.eta))
            veto_SF_err = self.lookups['phoVeto_err_2018'](photons.pt, abs(photons.eta))

        sf = id_SF * veto_SF
        if syst.lower()=='total':
            sfUp = (id_SF + id_SF_err) * (veto_SF + veto_SF_err)
            sfDo = (id_SF - id_SF_err) * (veto_SF - veto_SF_err)
        elif syst.lower()=='id':
            sfUp = (id_SF + id_SF_err) * (veto_SF)
            sfDo = (id_SF - id_SF_err) * (veto_SF)
        elif syst.lower()=='veto':
            sfUp = (id_SF) * (veto_SF + veto_SF_err)
            sfDo = (id_SF) * (veto_SF - veto_SF_err)
        else:
            print('Unrecognized syst level, returning no variation')
            sfUp = sf
            sfDo = sf

        return ak.prod(sf,axis=-1), ak.prod(sfUp,axis=-1), ak.prod(sfDo,axis=-1)

            
