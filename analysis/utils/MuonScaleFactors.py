import os.path
cwd = os.path.dirname(__file__)


from coffea.lookup_tools import extractor, dense_lookup

import awkward as ak

class MuonScaleFactors:
    def __init__(self):

        ext = extractor()

        directory=f"{cwd}/ScaleFactors/Muons"
        inputs = [f"id_2016BtoF NUM_TightID_DEN_genTracks_eta_pt {directory}/muon_id_2016BtoF.root",
                  f"iso_2016BtoF NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt {directory}/muon_iso_2016BtoF.root",
                  f"trig_2016BtoF abseta_pt_ratio {directory}/muon_trig_2016BtoF.root",
                  f"id_2016GH NUM_TightID_DEN_genTracks_eta_pt {directory}/muon_id_2016GH.root",
                  f"iso_2016GH NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt {directory}/muon_iso_2016GH.root",
                  f"trig_2016GH abseta_pt_ratio {directory}/muon_trig_2016GH.root",
                  f"id_2017 NUM_TightID_DEN_genTracks_pt_abseta {directory}/muon_id_2017.root",
                  f"iso_2017 NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta {directory}/muon_iso_2017.root",
                  f"trig_2017 abseta_pt_ratio {directory}/muon_trig_2017.root",
                  f"id_2018 NUM_TightID_DEN_TrackerMuons_pt_abseta {directory}/muon_id_2018.root",
                  f"iso_2018 NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta {directory}/muon_iso_2018.root",
                  f"trig_2018Pre abseta_pt_ratio {directory}/muon_trig_2018Pre.root",
                  f"trig_2018Post abseta_pt_ratio {directory}/muon_trig_2018Post.root"]

        ext.add_weight_sets(inputs)

        inputs = [f"id_2016BtoF_err NUM_TightID_DEN_genTracks_eta_pt_error {directory}/muon_id_2016BtoF.root",
                  f"iso_2016BtoF_err NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_error {directory}/muon_iso_2016BtoF.root",
                  f"trig_2016BtoF_err abseta_pt_ratio_error {directory}/muon_trig_2016BtoF.root",
                  f"id_2016GH_err NUM_TightID_DEN_genTracks_eta_pt_error {directory}/muon_id_2016GH.root",
                  f"iso_2016GH_err NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_error {directory}/muon_iso_2016GH.root",
                  f"trig_2016GH_err abseta_pt_ratio_error {directory}/muon_trig_2016GH.root",
                  f"id_2017_err NUM_TightID_DEN_genTracks_pt_abseta_error {directory}/muon_id_2017.root",
                  f"iso_2017_err NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_error {directory}/muon_iso_2017.root",
                  f"trig_2017_err abseta_pt_ratio_error {directory}/muon_trig_2017.root",
                  f"id_2018_err NUM_TightID_DEN_TrackerMuons_pt_abseta_error {directory}/muon_id_2018.root",
                  f"iso_2018_err NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_error {directory}/muon_iso_2018.root",
                  f"trig_2018Pre_err abseta_pt_ratio_error {directory}/muon_trig_2018Pre.root",
                  f"trig_2018Post_err abseta_pt_ratio_error {directory}/muon_trig_2018Post.root"]

        ext.add_weight_sets(inputs)

        ext.finalize()
        self.lookups = ext.make_evaluator()

    def getMuonScaleFactors(self, muons, year, syst='Total'):
        if year==2016:
            id_SF_A = self.lookups["id_2016BtoF"](muons.eta,muons.pt)*19.695422959/35.921875595
            id_SF_B = self.lookups["id_2016GH"](muons.eta,muons.pt)*16.226452636/35.921875595
            id_SF = id_SF_A + id_SF_B

            iso_SF_A = self.lookups["iso_2016BtoF"](muons.eta,muons.pt)*19.695422959/35.921875595
            iso_SF_B = self.lookups["iso_2016GH"](muons.eta,muons.pt)*16.226452636/35.921875595
            iso_SF = iso_SF_A + iso_SF_B

            trig_SF_A = self.lookups["trig_2016BtoF"](abs(muons.eta),muons.pt)*19.695422959/35.921875595
            trig_SF_B = self.lookups["trig_2016GH"](abs(muons.eta),muons.pt)*16.226452636/35.921875595
            trig_SF = trig_SF_A + trig_SF_B

            id_SF_A_err = self.lookups["id_2016BtoF_err"](muons.eta,muons.pt)*19.695422959/35.921875595
            id_SF_B_err = self.lookups["id_2016GH_err"](muons.eta,muons.pt)*16.226452636/35.921875595
            id_SF_err = id_SF_A_err + id_SF_B_err

            iso_SF_A_err = self.lookups["iso_2016BtoF_err"](muons.eta,muons.pt)*19.695422959/35.921875595
            iso_SF_B_err = self.lookups["iso_2016GH_err"](muons.eta,muons.pt)*16.226452636/35.921875595
            iso_SF_err = iso_SF_A_err + iso_SF_B_err

            trig_SF_A_err = self.lookups["trig_2016BtoF_err"](abs(muons.eta),muons.pt)*19.695422959/35.921875595
            trig_SF_B_err = self.lookups["trig_2016GH_err"](abs(muons.eta),muons.pt)*16.226452636/35.921875595
            trig_SF_err = trig_SF_A_err + trig_SF_B_err

        if year==2017:
            id_SF = self.lookups['id_2017'](muons.pt, abs(muon.eta))
            iso_SF = self.lookups['iso_2017'](muons.pt, abs(muon.eta))
            trig_SF = self.lookups['trig_2017'](abs(muon.eta), muons.pt)

            id_SF_err = self.lookups['id_2017_err'](muons.pt, abs(muon.eta))
            iso_SF_err = self.lookups['iso_2017_err'](muons.pt, abs(muon.eta))
            trig_SF_err = self.lookups['trig_2017_err'](abs(muon.eta), muons.pt)

        if year==2018:
            id_SF = self.lookups['id_2018'](muons.pt, abs(muon.eta))
            iso_SF = self.lookups['iso_2018'](muons.pt, abs(muon.eta))
            trig_SF_A = self.lookups['trig_2018Pre'](abs(muon.eta), muons.pt)*8.950818835/59.740565202
            trig_SF_B = self.lookups['trig_2018Post'](abs(muon.eta), muons.pt)*50.789746366/59.740565202
            trig_SF = trig_SF_A + trig_SF_B

            id_SF_err = self.lookups['id_2018_err'](muons.pt, abs(muon.eta))
            iso_SF_err = self.lookups['iso_2018_err'](muons.pt, abs(muon.eta))
            trig_SF_A_err = self.lookups['trig_2018Pre_err'](abs(muon.eta), muons.pt)*8.950818835/59.740565202
            trig_SF_B_err = self.lookups['trig_2018Post_err'](abs(muon.eta), muons.pt)*50.789746366/59.740565202
            trig_SF_err = trig_SF_A_err + trig_SF_B_err

        sf = id_SF * iso_SF * trig_SF
        if syst.lower()=='total':
            sfUp = (id_SF + id_SF_err) * (iso_SF + iso_SF_err) * (trig_SF + trig_SF_err)
            sfDo = (id_SF - id_SF_err) * (iso_SF - iso_SF_err) * (trig_SF - trig_SF_err)
        elif syst.lower()=='id':
            sfUp = (id_SF + id_SF_err) * (iso_SF + iso_SF_err) * (trig_SF)
            sfDo = (id_SF - id_SF_err) * (iso_SF - iso_SF_err) * (trig_SF)
        elif syst.lower()=='trig':
            sfUp = (id_SF) * (iso_SF) * (trig_SF + trig_SF_err)
            sfDo = (id_SF) * (iso_SF) * (trig_SF - trig_SF_err)
        else:
            print('Unrecognized syst level, returning no variation')
            sfUp = sf
            sfDo = sf

        return ak.prod(sf,axis=-1), ak.prod(sfUp,axis=-1), ak.prod(sfDo,axis=-1)
            
        
