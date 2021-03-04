
import numpy as np

from coffea.nanoevents import NanoAODSchema
from coffea import util, processor

from coffea.analysis_tools import PackedSelection

import uproot
uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

from analysis.processor import TstarSelector
from analysis.utils.datasets import fileset as filesetFull
from analysis.utils.datasets import filesetData as filesetDataFull


def fix_uproot():
    import uproot
    uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource


#the NanoAODSchema needs to be adjusted, to remove cross references to FSRPhotons
class SkimmedSchema(NanoAODSchema):
    def __init__(self, base_form):
        base_form["contents"].pop("Muon_fsrPhotonIdx", None)
        base_form["contents"].pop("FatJet_subJetIdx1", None)
        base_form["contents"].pop("FatJet_subJetIdx2", None)
        super().__init__(base_form)


def runJob(dask=False, workers=4, memory="2GB", chunkSize=50000, maxChunks=-1, bJetCut='1b', dataset='small', year=2016, maxFiles=None):
    if dask:
        from dask.distributed import Client
        from lpc_dask import HTCondorCluster
        import socket
        import time
        import os

        extra = ['--worker-port 10002:10100']

        hostname = socket.gethostname()
        from glob import glob
        # make list of local package directories (absolute paths) that should be sent to jobs
        initpylist = [os.path.abspath(os.path.dirname(x)) for x in glob('*/__init__.py')]
        initpylist.append("patch.sh")
    
        job_extra = {'transfer_input_files': ','.join(initpylist)}

        cluster = HTCondorCluster(scheduler_options = {'host': f'{hostname}:10000'},
                                  cores=1, 
                                  memory=memory, 
                                  disk="2GB", 
                                  python='python',
                                  nanny=False,
                                  extra=extra,
                                  job_extra=job_extra,
                              )

        cluster.scale(jobs=workers)
        
        client = Client(cluster)
        client.run(fix_uproot)

    
    fileLists = {'glugamma'      : [f'TstarTstar_TGluTGamma_M{massPoint}00' for massPoint in range(8,17)],
                 'gluglu'        : [f'TstarTstar_TGluTGlu_M{massPoint}00' for massPoint in range(8,17)],
                 'gammagamma'    : [f'TstarTstar_TGammaTGamma_M{massPoint}00' for massPoint in range(8,17)],
                 'ttgamma_1l'    : ['TTGamma_SingleLept','TTGamma_SingleLept_Pt100','TTGamma_SingleLept_Pt200'],
                 'ttgamma_2l'    : ['TTGamma_Dilepton','TTGamma_Dilepton_Pt100','TTGamma_Dilepton_Pt200'],
                 'ttgamma_0l'    : ['TTGamma_Hadronic'],
                 'ttbar_0l'      : ['TTbarPowheg_Hadronic'],
                 'ttbar_1l'      : ['TTbarPowheg_SingleLept'],
                 'ttbar_2l'      : ['TTbarPowheg_Dilepton'],
                 'wjets'         : ['W1jets','W2jets','W3jets','W4jets'],
                 'zjets'         : ['DYjetsM10to50','DYjetsM50'],
                 'vgamma'        : ['WGamma','ZGamma_01J_5f_lowMass'],
                 'ttx'           : ['TTWtoLNu','TTWtoQQ','TTZtoLL','TTZtoLL_M1to10','TTZtoQQ'],
                 'singletop'     : ['ST_s_channel','ST_tW_channel','ST_tbarW_channel','ST_t_channel','ST_tbar_channel','TGJets'],
                 'ww'            : ['WWTo4Q','WWToLNuQQ'],
                 'wz'            : ['WZTo1L1Nu2Q','WZTo1L3Nu','WZTo2L2Q','WZTo3L1Nu'],
                 'zz'            : ['ZZTo2L2Q','ZZTo2Q2Nu','ZZTo4L']
    }

    fileLists['small'] = fileLists['glugamma'] + fileLists['ttgamma_1l']

    fileLists['signal'] = fileLists['glugamma'] + fileLists['gluglu'] + fileLists['gammagamma']
    fileLists['ttgamma'] = fileLists['ttgamma_0l'] + fileLists['ttgamma_1l'] + fileLists['ttgamma_2l']
    fileLists['ttbar'] = fileLists['ttbar_0l'] + fileLists['ttbar_1l'] + fileLists['ttbar_2l']
    fileLists['diboson'] = fileLists['ww'] + fileLists['wz'] + fileLists['zz'] + ['VVTo2L2Nu']

    fileLists['other'] = fileLists['diboson'] + fileLists['ttx'] + fileLists['vgamma']
    
    fileLists['fullmc']   = [k[0] for k in filesetFull.keys()]

    fileLists['fulldata'] = [k[0] for k in filesetDataFull.keys()]
    fileLists['dataele'] = [k[0] for k in filesetDataFull.keys() if 'SingleEle' in k[0]]
    fileLists['datamu'] = [k[0] for k in filesetDataFull.keys() if 'SingleMu' in k[0]]

    isData = 'data' in dataset.lower()
    
    if not isData:
        if not dataset.lower() in fileLists:
            if (dataset,year) in filesetFull:
                flist = [dataset]
            else:
                print('unknown dataset key')
                return
        else:
            flist = fileLists[dataset.lower()]
        fileset = {k:filesetFull[k][:maxFiles] for k in filesetFull if (k[0] in flist and k[1]==year)}
    else:
        if not dataset.lower() in fileLists:
            if (dataset,year) in filesetDataFull:
                flist = [dataset]
            else:
                print('unknown dataset key')
                return
        else:
            flist = fileLists[dataset.lower()]
        fileset = {k:filesetDataFull[k][:maxFiles] for k in filesetDataFull if (k[0] in flist and k[1]==year)}

    if isData:
        nEvents = None
        mcPUDist = None
    else:
        nEvents = {}
        mcPUDist = {}
        for d in fileset:
            if not d in nEvents:
                nEvents[d] = 0
                for fName in fileset[d]:
                    with uproot.open(fName)['hEvents'] as hEvents:
                        nEvents[d] += hEvents.values()[1]
                        with uproot.open(fName)['hPUTrue'] as hPUTrue:
                            if not d in mcPUDist:
                                mcPUDist[d] = hPUTrue.values()
                            else:
                                mcPUDist[d] += hPUTrue.values()

    proc = TstarSelector(mcTotals=nEvents, mcPUDist=mcPUDist, bJetCut=bJetCut)

    exe_args = {"schema": SkimmedSchema,'workers': workers, 'savemetrics': True}

    if dask:
        exe_args = {
            'client': client,
            'savemetrics': True,
            'schema': SkimmedSchema,
            'align_clusters': True,
            'workers': 1
        }

        while len(client.ncores()) < 1:
            print('Waiting for more cores to spin up, currently there are {0} available...'.format(len(client.ncores())))
            print('Dask client info ->', client)
            time.sleep(10)


    print(fileset)

    hists = processor.run_uproot_job(fileset,
                                     treename="Events",
                                     processor_instance=proc,
                                     executor=processor.dask_executor if dask else processor.futures_executor,
                                     executor_args=exe_args,
                                     chunksize=chunkSize,
                                     maxchunks=maxChunks,
                                 )
    return hists

if __name__ == "__main__":

    from  argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(      '--dask',       help='run w/ dask', dest='dask',              default=False, action='store_true')
    parser.add_argument('-w', '--workers',    help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=4)
    parser.add_argument(      '--dataset',    help='Keywords for subsets of samples', dest='dataset', default="small")
    parser.add_argument('-y', '--year',       help='Year to run on', dest='year', type=int, default=2016)
    parser.add_argument(      '--maxFiles',   help='Number of files to use for each data type', dest='maxFiles', type=int, default=None)
    parser.add_argument(      '--memory',     help='Amount of memory to require for condor nodes', dest='memory', default="4GB")
    parser.add_argument(      '--chunkSize',  help='Coffea processor chunk size', dest='chunkSize', type=int, default=100000)
    parser.add_argument(      '--maxChunks',  help='Coffea processor maximum chunks', dest='maxChunks', type=int, default=None)
    parser.add_argument(      '--bJetCut',    help='Cut to use for bjet selection (1b, 2b, ex2b)', dest='bJetCut', default='1b')

    args = parser.parse_args()

    hists = runJob(**vars(args))

    if not hists is None:
        util.save(hists[0],f'hists_{args.dataset}_{args.year}_{args.bJetCut}.coffea')
    
