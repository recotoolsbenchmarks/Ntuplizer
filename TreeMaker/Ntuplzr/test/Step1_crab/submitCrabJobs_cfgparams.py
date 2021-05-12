#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.psetName = '/afs/cern.ch/work/s/sandhya/Physics/Upgrade/RTB/snowmass/CMSSW_11_3_0_pre4/src/step_3_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['maxEvents=-1']
config.JobType.outputFiles = ['file.root']
config.JobType.inputFiles = ["/afs/cern.ch/work/s/sandhya/Physics/Upgrade/RTB/snowmass/CMSSW_11_3_0_pre4/src/MyTools/"]
#config.JobType.inputFiles = ["MyTools"]
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset = ''
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#config.Data.totalUnits = 500
config.Data.outLFNDirBase = '/store/group/upgrade/RTB/Iter6/11_3/RAW'
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


            
    #config.Data.inputDataset = '/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v2/FEVT'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'TT_TuneCP5_14TeV-powheg-pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #
    #config.Data.inputDataset = '/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset = '/GluGluToHHTo2B2Tau_node_SM_TuneCP5_14TeV-madgraph-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext1-v3/FEVT'
    ###config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'GluGluToHHTo2B2Tau_node_SM_TuneCP5_14TeV-madgraph-pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset = '/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_castor_111X_mcRun4_realistic_T15_v1-v1/FEVT'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()


    #config.Data.inputDataset = '/MultiTau_PT15to500/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'MultiTau_PT15to500_v2'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    
    #config.Data.inputDataset =  '/MultiTau_PT15to500/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v3/FEVT'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'MultiTau_PT15to500_v3'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
   

    config.Data.inputDataset = '/QCD_Pt_20to30_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD'
    #config.Data.secondaryInputDataset = ''
    config.General.requestName = 'QCD_Pt_20to30_TuneCP5_14TeV_pythia8'
    p = Process(target=submit, args=(config,)) 
    p.start() 
    p.join()

    #config.Data.inputDataset = '/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_30to50_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset = '/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD' # do we want "withNewMB" over the other?  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_30to50_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset = '/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_50to80_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset = '/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v3/GEN-SIM-DIGI-RAW-MINIAOD' # same here  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_50to80_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  '/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_80to120_TuneCP5_14TeV_pythia'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  '/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_120to170_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  '/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_170to300_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset = '/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_300to470_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  '/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_470to600_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset = '/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD'
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = 'QCD_Pt_600oInf_TuneCP5_14TeV_pythia8'
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()

    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
    #
    #config.Data.inputDataset =  
    ##config.Data.secondaryInputDataset = ''
    #config.General.requestName = ''
    #p = Process(target=submit, args=(config,)) 
    #p.start() 
    #p.join()
