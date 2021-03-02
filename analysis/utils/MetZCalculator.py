import awkward as ak
import numpy as np

##calculate MET pz component solutions

def MetPz(ele, mu, met, solutionType=-1):
    lep_m = ak.num(ele)*0.000511 + ak.num(mu)*0.105658367

    #sum across, to get a 1 dim array
    # would be invalid for nLep!=1, but those cases are not valid in this calculation anyways
    lep_px = ak.sum(ele.px, axis=1) + ak.sum(mu.px, axis=1)
    lep_py = ak.sum(ele.py ,axis=1) + ak.sum(mu.py, axis=1)
    lep_pz = ak.sum(ele.pz ,axis=1) + ak.sum(mu.pz, axis=1)
    lep_e  = ak.sum(ele.energy ,axis=1) + ak.sum(mu.energy, axis=1)

    a = 80.4*80.4 - lep_m*lep_m + 2.0*lep_px*met.px + 2.0*lep_py*met.py
    A = 4.0*(lep_e*lep_e - lep_pz*lep_pz)
    B = -4.0*a*lep_pz
    C = 4.0*lep_e*lep_e*(met.px*met.px + met.py*met.py) - a*a

    tmproot = B*B - 4.0*A*C;

    #if no leptons, set A to 1 to avoid inf
    nLep = (ak.num(ele)+ak.num(mu))
    A = ak.where(nLep==0, 1, A)

    nu_pz_1 = np.where(tmproot<0, -B/(2*A), (-B + (tmproot**0.5))/(2.0*A))
    nu_pz_2 = np.where(tmproot<0, 999999999999., (-B - (tmproot**0.5))/(2.0*A))

    nu_pz = ak.Array(np.array([np.asarray(nu_pz_1),np.asarray(nu_pz_2)]).transpose().tolist())

    #keep only first solution if second is root is not real
    nu_pz = ak.where(tmproot<0,nu_pz[:,:1],nu_pz[:,:])
    #keep solutions only if exactly 1 lepton
    nu_pz = ak.where(nLep==1,nu_pz[:,:],nu_pz[:,:0])

    if solutionType==0:

        lepPz = ak.concatenate([mu.pz,ele.pz],axis=1)
        temp_nu_pz, temp_lep_pz = ak.unzip(ak.cartesian([nu_pz,lepPz]))
        delta_pz = abs(temp_nu_pz - temp_lep_pz)
        
        return np.where(ak.num(nu_pz)<2, ## if one solution, use it
                 nu_pz,
                 np.where(ak.any(abs(nu_pz)>300,axis=1),  #if one solution is >300, use most central one
                          nu_pz[ak.argmin(abs(nu_pz),axis=1,keepdims=True)],
                          nu_pz[ak.argmin(delta_pz,axis=1,keepdims=True)] #otherwise, use one closest to lepton pz
                         )
                 )
        
        twoSol = nu_pz[ak.num(nu_pz)==2]
    elif solutionType==1:
        lepPz = ak.concatenate([mu.pz,ele.pz],axis=1)
        temp_nu_pz, temp_lep_pz = ak.unzip(ak.cartesian([nu_pz,lepPz]))
        delta_pz = abs(temp_nu_pz - temp_lep_pz)
        return np.where(ak.num(nu_pz)<2, ## if one solution, use it
                        nu_pz,
                        nu_pz[ak.argmin(delta_pz,axis=1,keepdims=True)])
    elif solutionType==2:
        return np.where(ak.num(nu_pz)<2, ## if one solution, use it
                        nu_pz,
                        nu_pz[ak.argmin(abs(nu_pz),axis=1,keepdims=True)])
    else:
        return nu_pz

def getNeutrinoArrays(nu_pz, met):
    _nu_pz, _nu = ak.unzip(ak.cartesian([nu_pz,met]))
    _nu_eta = np.arcsinh(_nu_pz/_nu.pt)

    neutrinos = ak.zip(
        {
            "pt": _nu.pt,
            "eta": _nu_eta,
            "phi": _nu.phi,
            "mass": abs(_nu_eta)*0
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    return neutrinos
