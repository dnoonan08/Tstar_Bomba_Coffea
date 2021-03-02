
import numpy as np

from coffea.nanoevents import NanoAODSchema
from coffea import util, processor

from coffea.analysis_tools import PackedSelection

import uproot
uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

from analysis.processor import TstarSelector
from analysis.utils.datasets import fileset as filesetFull


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


def runJob(dask=False, workers=4, memory="2GB", chunkSize=50000, maxChunks=-1):
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


#    flist = ['TstarTstar_TGluTGamma_M800','TstarTstar_TGluTGamma_M1000','TstarTstar_TGluTGamma_M1200','TstarTstar_TGluTGamma_M1400','TstarTstar_TGluTGamma_M1600','TTbarPowheg_SingleLept','TTGamma_SingleLept','TTGamma_SingleLept_Pt100','TTGamma_SingleLept_Pt200','WGamma','W4jets']
    flist = ['TstarTstar_TGluTGamma_M800','TstarTstar_TGluTGamma_M1000','TstarTstar_TGluTGamma_M1200','TstarTstar_TGluTGamma_M1400','TstarTstar_TGluTGamma_M1600','TTbarPowheg_SingleLept']

    fileset = {k:filesetFull[k][:1] for k in filesetFull if k[0] in flist}


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

    proc = TstarSelector(mcTotals=nEvents, mcPUDist=mcPUDist, bJetCut='1b')

    exe_args = {"schema": SkimmedSchema,'workers': workers, 'savemetrics': True}

    if dask:
        exe_args = {
            'client': client,
            'savemetrics': True,
            'schema': SkimmedSchema,
            'align_clusters': True
        }

        while len(client.ncores()) < 2:
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
    parser.add_argument(      '--local',      help='run w/ dask', dest='dask',              default=True, action='store_false')
    parser.add_argument('-w', '--workers',    help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=4)
    parser.add_argument(      '--memory',     help='Amount of memory to require for condor nodes', dest='memory', default="4GB")
    parser.add_argument(      '--chunkSize',  help='Coffea processor chunk size', dest='chunkSize', type=int, default=100000)
    parser.add_argument(      '--maxChunks',  help='Coffea processor maximum chunks', dest='maxChunks', type=int, default=None)

    args = parser.parse_args()

    hists = runJob(**vars(args))

    util.save(hists[0],'signals_1b.coffea')
    
