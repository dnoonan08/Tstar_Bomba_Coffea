
eosDir = 'root://cmseos.fnal.gov//store/user/lpctop/TTGamma_FullRun2/BSM_Skims/2016'


fileset = {
           ('TstarTstar_TGluTGamma_M800',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M800_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M900',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M900_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M1000',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M1000_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M1100',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M1100_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M1200',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M1200_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M1300',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M1300_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M1400',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M1400_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M1500',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M1500_2016_skim.root'],
           ('TstarTstar_TGluTGamma_M1600',2016):[f'{eosDir}/TstarTstarToTgammaTgluon_M1600_2016_skim.root'],

           ('TstarTstar_TGluTGlu_M800',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M800_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M900',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M900_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M1000',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M1000_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M1100',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M1100_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M1200',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M1200_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M1300',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M1300_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M1400',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M1400_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M1500',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M1500_2016_skim.root'],
           ('TstarTstar_TGluTGlu_M1600',2016):[f'{eosDir}/TstarTstarToTgluonTgluon_M1600_2016_skim.root'],

           ('TstarTstar_TGammaTGamma_M800',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M800_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M900',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M900_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M1000',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M1000_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M1100',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M1100_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M1200',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M1200_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M1300',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M1300_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M1400',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M1400_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M1500',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M1500_2016_skim.root'],
           ('TstarTstar_TGammaTGamma_M1600',2016):[f'{eosDir}/TstarTstarToTgammaTgamma_M1600_2016_skim.root'],

           ('TTGamma_SingleLept',2016):[f'{eosDir}/TTGamma_SingleLept_2016_skim.root'],
           ('TTGamma_SingleLept_Pt100',2016):[f'{eosDir}/TTGamma_SingleLept_Pt100_2016_skim.root'],
           ('TTGamma_SingleLept_Pt200',2016):[f'{eosDir}/TTGamma_SingleLept_Pt200_2016_skim.root'],
           ('TTGamma_Dilepton',2016):[f'{eosDir}/TTGamma_Dilepton_2016_skim.root'],
           ('TTGamma_Dilepton_Pt100',2016):[f'{eosDir}/TTGamma_Dilepton_Pt100_2016_skim.root'],
           ('TTGamma_Dilepton_Pt200',2016):[f'{eosDir}/TTGamma_Dilepton_Pt200_2016_skim.root'],
           ('TTGamma_Hadronic',2016):[f'{eosDir}/TTGamma_Hadronic_2016_skim.root'],


           ('TTbarPowheg_SingleLept',2016):[f'{eosDir}/TTbarPowheg_Semilept_2016_skim_{i}of10.root' for i in range(1,11)],
           ('TTbarPowheg_Dilepton',2016):[f'{eosDir}/TTbarPowheg_Dilepton_2016_skim_{i}of10.root' for i in range(1,11)],
           ('TTbarPowheg_Hadronic',2016):[f'{eosDir}/TTbarPowheg_Hadronic_2016_skim_{i}of10.root' for i in range(1,11)],

           ('TGJets',2016):[f"{eosDir}/TGJets_2016_skim.root",f"{eosDir}/TGJets_ext1_2016_skim.root"],

           ('ST_s_channel',2016):[f"{eosDir}/ST_s_channel_2016_skim.root"],
           ('ST_tW_channel',2016):[f"{eosDir}/ST_tW_channel_2016_skim.root"],
           ('ST_tbarW_channel',2016):[f"{eosDir}/ST_tbarW_channel_2016_skim.root"],
           ('ST_t_channel',2016):[f"{eosDir}/ST_t_channel_2016_skim_{i}of5.root" for i in range(1,6)],
           ('ST_tbar_channel',2016):[f"{eosDir}/ST_tbar_channel_2016_skim_{i}of5.root" for i in range(1,6)],

    
           ('WGamma',2016):[f"{eosDir}/WGamma_2016_skim.root"],
           ('ZGamma_01J_5f_lowMass',2016):[f"{eosDir}/ZGamma_01J_5f_lowMass_2016_skim_{i}of5.root" for i in range(1,6)],

           ('TTWtoLNu',2016):[f"{eosDir}/TTWtoLNu_ext1_2016_skim.root",f"{eosDir}/TTWtoLNu_ext2_2016_skim.root"],
           ('TTWtoQQ',2016):[f"{eosDir}/TTWtoQQ_2016_skim.root"],
           ('TTZtoLL',2016):[f"{eosDir}/TTZtoLL_ext1_2016_skim.root",f"{eosDir}/TTZtoLL_ext2_2016_skim.root",f"{eosDir}/TTZtoLL_ext3_2016_skim.root"],
           ('TTZtoLL_M1to10',2016):[f"{eosDir}/TTZtoLL_M1to10_2016_skim.root"],
           ('TTZtoQQ',2016):[f"{eosDir}/TTZtoQQ_2016_skim.root"],

           ('W1jets',2016):[f"{eosDir}/W1jets_2016_skim.root"],
           ('W2jets',2016):[f"{eosDir}/W2jets_2016_skim.root",f"{eosDir}/W2jets_ext1_2016_skim.root"],
           ('W3jets',2016):[f"{eosDir}/W3jets_2016_skim.root",f"{eosDir}/W3jets_ext1_2016_skim.root"],
           ('W4jets',2016):[f"{eosDir}/W4jets_2016_skim.root",f"{eosDir}/W4jets_ext1_2016_skim.root",f"{eosDir}/W4jets_ext2_2016_skim.root"],

           ('DYjetsM10to50',2016):[f'{eosDir}/DYjetsM10to50_2016_skim.root'],
           ('DYjetsM50',2016):[f'{eosDir}/DYjetsM50_ext1_2016_skim_{i}of10.root' for i in range(1,11)] + [f'{eosDir}/DYjetsM50_ext2_2016_skim_{i}of10.root' for i in range(1,11)],

           ('VVTo2L2Nu',2016):[f"{eosDir}/VVTo2L2Nu_amcatnlo_2016_skim.root",f"{eosDir}/VVTo2L2Nu_amcatnlo_ext1_2016_skim.root"],

           ('WWTo4Q',2016):[f"{eosDir}/WWTo4Q_powheg_2016_skim.root"],
           ('WWToLNuQQ',2016):[f"{eosDir}/WWToLNuQQ_powheg_2016_skim.root",f"{eosDir}/WWToLNuQQ_powheg_ext1_2016_skim.root"],

           ('WZTo1L1Nu2Q',2016):[f"{eosDir}/WZTo1L1Nu2Q_amcatnlo_2016_skim_{i}of5.root" for i in range(1,6)],
           ('WZTo1L3Nu',2016):[f"{eosDir}/WZTo1L3Nu_amcatnlo_2016_skim.root"],
           ('WZTo2L2Q',2016):[f"{eosDir}/WZTo2L2Q_amcatnlo_2016_skim_{i}of10.root" for i in range(1,6)],
           ('WZTo3L1Nu',2016):[f"{eosDir}/WZTo3LNu_powheg_2016_skim.root",f"{eosDir}/WZTo3LNu_powheg_ext1_2016_skim_1of10.root",f"{eosDir}/WZTo3LNu_powheg_ext1_2016_skim_2of10.root",f"{eosDir}/WZTo3LNu_powheg_ext1_2016_skim_3of10.root",f"{eosDir}/WZTo3LNu_powheg_ext1_2016_skim_4of10.root",f"{eosDir}/WZTo3LNu_powheg_ext1_2016_skim_5of10.root"],

           ('ZZTo2L2Q',2016):[f"{eosDir}/ZZTo2L2Q_powheg_2016_skim.root"],
           ('ZZTo2Q2Nu',2016):[f"{eosDir}/ZZTo2Q2Nu_amcatnlo_2016_skim_1of5.root",f"{eosDir}/ZZTo2Q2Nu_amcatnlo_2016_skim_2of5.root",f"{eosDir}/ZZTo2Q2Nu_amcatnlo_2016_skim_3of5.root",f"{eosDir}/ZZTo2Q2Nu_amcatnlo_2016_skim_4of5.root",f"{eosDir}/ZZTo2Q2Nu_amcatnlo_2016_skim_5of5.root"],
           ('ZZTo4L',2016):[f"{eosDir}/ZZTo4L_powheg_2016_skim.root"],
        }


filesetData = {
    ('Data_SingleEle_b', 2016) : [f'{eosDir}/Data_SingleEle_b_2016_skim_{i}of10.root' for i in range(1,11)],
    ('Data_SingleEle_c', 2016) : [f'{eosDir}/Data_SingleEle_c_2016_skim_{i}of10.root' for i in range(1,11)],
    ('Data_SingleEle_d', 2016) : [f'{eosDir}/Data_SingleEle_d_2016_skim_{i}of10.root' for i in range(1,11)],
    ('Data_SingleEle_e', 2016) : [f'{eosDir}/Data_SingleEle_e_2016_skim_{i}of10.root' for i in range(1,11)],
    ('Data_SingleEle_f', 2016) : [f'{eosDir}/Data_SingleEle_f_2016_skim_{i}of10.root' for i in range(1,11)],
    ('Data_SingleEle_g', 2016) : [f'{eosDir}/Data_SingleEle_g_2016_skim_{i}of10.root' for i in range(1,11)],
    ('Data_SingleEle_h', 2016) : [f'{eosDir}/Data_SingleEle_h_2016_skim_{i}of10.root' for i in range(1,11)],
    ('Data_SingleMu_b', 2016)  : [f'{eosDir}/Data_SingleMu_b_2016_skim_{i}of10.root'  for i in range(1,11)],
    ('Data_SingleMu_c', 2016)  : [f'{eosDir}/Data_SingleMu_c_2016_skim_{i}of10.root'  for i in range(1,11)],
    ('Data_SingleMu_d', 2016)  : [f'{eosDir}/Data_SingleMu_d_2016_skim_{i}of10.root'  for i in range(1,11)],
    ('Data_SingleMu_e', 2016)  : [f'{eosDir}/Data_SingleMu_e_2016_skim_{i}of10.root'  for i in range(1,11)],
    ('Data_SingleMu_f', 2016)  : [f'{eosDir}/Data_SingleMu_f_2016_skim_{i}of10.root'  for i in range(1,11)],
    ('Data_SingleMu_g', 2016)  : [f'{eosDir}/Data_SingleMu_g_2016_skim_{i}of10.root'  for i in range(1,11)],
    ('Data_SingleMu_h', 2016)  : [f'{eosDir}/Data_SingleMu_h_2016_skim_{i}of10.root'  for i in range(1,11)],
}

def updateDatasetInfo():
    nEvents = {}
    mcPUDist = {}
    import uproot
    for d in fileset:
        if not d in nEvents:
            nEvents[d] = 0
        for fName in fileset[d]:
            with uproot.open(fName) as _file:
                with _file['hEvents'] as hEvents:
                    nEvents[d] += hEvents.values()[1]
                with _file['hPUTrue'] as hPUTrue:
                    if not d in mcPUDist:
                        mcPUDist[d] = hPUTrue.values()
                    else:
                        mcPUDist[d] += hPUTrue.values()

    from coffea import util
    import os.path
    cwd = os.path.dirname(__file__)

    util.save(nEvents,f'{cwd}/MCEventCounts.coffea')
    util.save(mcPUDist,f'{cwd}/MCPileupDists.coffea')


if __name__=="__main__":
    updateDatasetInfo()
