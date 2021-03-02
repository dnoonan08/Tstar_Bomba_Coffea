import awkward as ak
import numpy as np
import itertools

def permutations(array, n):
    tmp = ak.concatenate(ak.unzip(ak.combinations(array, n)[:, :, None]), axis=2)
    perms = np.array(list(itertools.permutations(range(n))))
    return ak.flatten(tmp[:,:,perms],axis=2)

def argpermutations(array, n):
    tmp = ak.concatenate(ak.unzip(ak.argcombinations(array, n)[:, :, None]), axis=2)
    perms = np.array(list(itertools.permutations(range(n))))
    return ak.flatten(tmp[:,:,perms],axis=2)

def permutations0(array, n):
    tmp = ak.concatenate(ak.unzip(ak.combinations(array, n)[:, :, None]), axis=2)
    return ak.concatenate((tmp[:, :, perm] for perm in itertools.permutations(range(n))), axis=1)

def jetPerm(tightJet,N):
    jetCombos = permutations(tightJet,N)
    GoodBjets = ((jetCombos[:,:,0].btagDeepB>0.6321) & 
                 (jetCombos[:,:,1].btagDeepB>0.6321))    
    return jetCombos[GoodBjets]

def jetPermArg(tightJet,N):
    jetCombos = argpermutations(tightJet,N)
    GoodBjets = ((tightJet[jetCombos[:,:,0]].btagDeepB>0.6321) & 
                 (tightJet[jetCombos[:,:,1]].btagDeepB>0.6321))
    
    return jetCombos[GoodBjets]

def jetPermSplit(tightJet,n,chunkSize):
    x = []
    for i in np.arange(0,len(tightJet),chunkSize):
        x = ak.concatenate([x, jetPerm(tightJet[i:(i+chunkSize)],n)])
    return x

def jetPermArgSplit(tightJet,n,chunkSize=1000):
    x = []
    for i in np.arange(0,len(tightJet),chunkSize):
        x = ak.concatenate([x, jetPermArg(tightJet[i:(i+chunkSize)],n)])
    return x
