from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.psetName = '/afs/cern.ch/work/s/sandhya/Physics/Upgrade/RTB/btag/git/CMSSW_10_6_0_patch2/src/TreeMaker/Ntuplzr/test/myproduceNtuples_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['maxEvents=10']
config.JobType.outputFiles = ['file.root']
config.JobType.maxMemoryMB = 1500
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/group/upgrade/RTB/Iter4/'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

            
    config.Data.inputDataset = '/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/PhaseIITDRSpring19MiniAOD-NoPU_106X_upgrade2023_realistic_v3-v2/MINIAODSIM'
    config.Data.secondaryInputDataset = '/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/PhaseIITDRSpring19DR-NoPU_106X_upgrade2023_realistic_v3-v2/AODSIM'
    config.General.requestName = 'QCD_Pt_120to170_TuneCP5_14TeV_pythia8'
    p = Process(target=submit, args=(config,)) 
    p.start() 
    p.join()
    # submit(config)
    


