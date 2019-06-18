from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.psetName = '/afs/cern.ch/work/s/sandhya/Physics/Upgrade/VBF/FullSim/CMSSW_10_4_0/src/TreeMaker/Ntuplzr/test/myproduceNtuples_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['file.root']
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/group/upgrade/RTB/Iter2/'
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

            
    config.Data.inputDataset = '/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8/PhaseIIMTDTDRAutumn18MiniAOD-PU200_103X_upgrade2023_realistic_v2-v1/MINIAODSIM'
    config.General.requestName = 'QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8'
    submit(config)
    

    config.Data.inputDataset = '/DYToLL_M-50_14TeV_TuneCP5_pythia8/PhaseIIMTDTDRAutumn18MiniAOD-PU200_103X_upgrade2023_realistic_v2-v2/MINIAODSIM'
    config.General.requestName = 'DYToLL_M-50_14TeV_TuneCP5_pythia8'
    submit(config)

