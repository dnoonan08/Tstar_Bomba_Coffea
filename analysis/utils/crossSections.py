# Everything in picobarns

def getLumiWeight(dataset, mcTotals, year=2016):

    if mcTotals is None:
        return 1.
    if dataset in crossSections:
        xsec=crossSections[dataset]
        if isinstance(xsec,dict):
            xsec=xsec[year]
        if (dataset,year) in mcTotals:
            nEvents=mcTotals[(dataset,year)]
            return xsec*lumis[year]/nEvents
        else:
            return 1.
        
    else:
        return 1.


lumis = {2016 : 35860.,
         2017 : 41529.,
         2018 : 58877.,
         }

crossSections = {
    "TstarTstar_TGluTGamma_M700"  :  4.92,
    "TstarTstar_TGluTGamma_M800"  :  1.68,
    "TstarTstar_TGluTGamma_M900"  :  0.636,
    "TstarTstar_TGluTGamma_M1000" :  0.262,
    "TstarTstar_TGluTGamma_M1100" :  0.116,
    "TstarTstar_TGluTGamma_M1200" :  0.0537,
    "TstarTstar_TGluTGamma_M1300" :  0.0261,
    "TstarTstar_TGluTGamma_M1400" :  0.0131,
    "TstarTstar_TGluTGamma_M1500" :  0.00677,
    "TstarTstar_TGluTGamma_M1600" :  0.00359,

    "TstarTstar_TGluTGlu_M700"  :  4.92,
    "TstarTstar_TGluTGlu_M800"  :  1.68,
    "TstarTstar_TGluTGlu_M900"  :  0.636,
    "TstarTstar_TGluTGlu_M1000" :  0.262,
    "TstarTstar_TGluTGlu_M1100" :  0.116,
    "TstarTstar_TGluTGlu_M1200" :  0.0537,
    "TstarTstar_TGluTGlu_M1300" :  0.0261,
    "TstarTstar_TGluTGlu_M1400" :  0.0131,
    "TstarTstar_TGluTGlu_M1500" :  0.00677,
    "TstarTstar_TGluTGlu_M1600" :  0.00359,
    
    'TTGamma_Dilepton'    : 1.495*1.4852,
    'TTGamma_SingleLept'  : 5.056*1.4852,
    'TTGamma_Hadronic'    : 4.149*1.4852,

    'TTGamma_Dilepton_Pt100'    :  0.03412*1.4852,
    'TTGamma_SingleLept_Pt100'  :  0.1309*1.4852,
    'TTGamma_Hadronic_Pt100'    :  0.1249*1.4852,

    'TTGamma_Dilepton_Pt200'   :  0.006797*1.4852,
    'TTGamma_SingleLept_Pt200' :  0.02685*1.4852,
    'TTGamma_Hadronic_Pt200'   :  0.02687*1.4852,

    'TTbar'             : 831.76,

    'TTbarPowheg_Dilepton'    :  88.29,
    'TTbarPowheg_SingleLept'  : 365.34,
    'TTbarPowheg_Hadronic'    : 377.96,

    'W1jets'            : 11775.9345,
    'W2jets'            :  3839.4345,
    'W3jets'            :  1165.8108,
    'W4jets'            :   592.9176,

    'DYjetsM50'         :  6077.22,
    'DYjetsM10to50'     : 18610.0,

    'TTWtoQQ'           :  0.40620,
    'TTWtoLNu'          :  0.2043,
    'TTZtoLL'           :  0.2728,
    'TTZtoLL_M1to10'    :  0.0493,
    'TTZtoQQ'           :  0.5297,
    
    'ZGamma'            :  {2016 :  98.3, 
                            2017 : 105.4, 
                            2018 : 105.4},
    'WGamma'            :  {2016 : 489., 
                            2017 : 463.9*1.295, 
                            2018 : 463.9*1.295},

    'TGJets'             : 2.967,

    'ST_tW_channel'      :  35.85,
    'ST_tbarW_channel'   :  35.85,
    'ST_t_channel'       :  136.02,
    'ST_tbar_channel'    :  80.95,
    'ST_s_channel'       :  3.68064,

    'WWTo4Q' : 51.723,
    'WWToLNuQQ' : 49.997,
    
    'WZTo1L1Nu2Q' : 10.71,
    'WZTo1L3Nu'   :  3.033,
    'WZTo2L2Q'    :  5.595,
    'WZTo3L1Nu'   :  4.42965,

    'ZZTo2L2Q'  : 3.28,
    'ZZTo2Q2Nu' : 4.04,
    'ZZTo4L'    : 1.3816,
    
    'VVTo2L2Nu' : 11.95,
    
    # 'WW'                :  75.8,
    # 'WZ'                :  27.6,
    # 'ZZ'                :  12.14,

    # 'QCD_Pt20to30_Mu'    : 2960198.4,
    # 'QCD_Pt30to50_Mu'    : 1652471.46,
    # 'QCD_Pt50to80_Mu'    :  437504.1,
    # 'QCD_Pt80to120_Mu'   :  106033.6648,
    # 'QCD_Pt120to170_Mu'  :   25190.5151,
    # 'QCD_Pt170to300_Mu'  :    8654.4932,
    # 'QCD_Pt300to470_Mu'  :     797.3527,
    # 'QCD_Pt470to600_Mu'  :      79.0255,
    # 'QCD_Pt600to800_Mu'  :      25.0951,
    # 'QCD_Pt800to1000_Mu' :       4.7074,
    # 'QCD_Pt1000toInf_Mu' :       1.6213,

    # 'QCD_Pt20to30_Ele'   : 5352960.,
    # 'QCD_Pt30to50_Ele'   : 9928000.,
    # 'QCD_Pt50to80_Ele'   : 2890800.,
    # 'QCD_Pt80to120_Ele'  :  350000.,
    # 'QCD_Pt120to170_Ele' :   62964.,
    # 'QCD_Pt170to300_Ele' :   18810.,
    # 'QCD_Pt300toInf_Ele' :   1350.,

    # 'GJets_HT40To100'  : 20790.,
    # 'GJets_HT100To200' : 9238.,
    # 'GJets_HT200To400' : 2305.,
    # 'GJets_HT400To600' : 274.4,
    # 'GJets_HT600ToInf' : 93.46,

}
