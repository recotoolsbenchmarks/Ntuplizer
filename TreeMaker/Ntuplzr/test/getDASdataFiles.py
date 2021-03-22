import argparse
import os


# Argument parser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument(
    "--getCount",
    help = "Get the number of events",
    default = False,
    action = "store_true",
)

parser.add_argument(
    "--getFiles",
    help = "Get the list of files",
    default = False,
    action = "store_true",
)

# Parse arguments
args = parser.parse_args()
d_args = vars(args)


outDir = "sourceFiles"

prefix = "root://eoscms.cern.ch//eos/cms/store/group/upgrade/RTB/Iter6/11_3/RAW/"
#prefix = "/eos/cms/store/group/upgrade/RTB/Iter6/11_3/RAW/"
toReplace = "ntuple"


l_sampleName = [
    '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/FEVT',
    # more leptons and taus, hopefully in barrel+endcap
    '/DoublePhoton_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT',
    '/DoublePhoton_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v3/FEVT',
    '/DoubleElectron_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD',
    '/DoubleMuon_gun_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT',
   
    
    # QCD flat + more QCD, although not flat in these pT ranges
    #'/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_castor_111X_mcRun4_realistic_T15_v1-v1/FEVT',
    #'/QCD_Pt_20to30_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD', # do we want "withNewMB" over the other?
    #'/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v3/GEN-SIM-DIGI-RAW-MINIAOD', # same here
    #'/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    #'/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    
    #"/QCD_Pt-15to3000_TuneCP5_Flat_14TeV-pythia8/Phase2HLTTDRWinter20DIGI-PU200_castor_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    
    #"/SingleElectron_PT2to200/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    #"/SingleElectron_PT200to500/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    #"/DoubleElectron_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD",
    
    #"/SinglePhoton_PT2to200/Phase2HLTTDRWinter20DIGI-NoPU_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/SinglePhoton_PT2to200/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    #"/SinglePhoton_PT200to500/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1/FEVT",
    #"/DoublePhoton_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT",
    
    #"/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT",
    
    #"/QCD_Pt-15to20_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-20to30_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-30to50_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-50to80_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-80to120_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-120to170_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-170to300_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-300toInf_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD",
    #"/QCD_Pt-15to3000_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD",
    
    #"/QCD_Pt-15to20_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v1/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-20to30_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-30to50_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-50to80_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-80to120_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-120to170_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-170to300_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-300toInf_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v3/GEN-SIM-DIGI-RAW",
    
    #"/QCD_Pt-15to20_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v1/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-20to30_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-30to50_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-50to80_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-80to120_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-120to170_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-170to300_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v2/GEN-SIM-DIGI-RAW",
    #"/QCD_Pt-300toInf_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRWinter20DIGI-PU200_110X_mcRun4_realistic_v3-v3/GEN-SIM-DIGI-RAW",
]


for iSample, sampleName in enumerate(l_sampleName) :
    
    print "\n"
    print "*"*50
    print "Sample %d/%d: %s" %(iSample+1, len(l_sampleName), sampleName)
    print "*"*50
    
    if (args.getCount) :
        
        command = "ls /eos/cms/store/group/upgrade/RTB/Iter6/11_3/RAW/%s " %(sampleName)
        #command = "dasgoclient -query=\"file dataset=%s | sum(file.nevents)\"" %(sampleName)
        os.system(command)
    
    if (args.getFiles) :
        
        sampleName_mod = sampleName[1:].replace("/", "_")
        
        outDir_mod = "%s/%s" %(outDir, sampleName_mod)
        
        command = "mkdir -p %s" %(outDir_mod)
        print "Command:", command
        print ""
        os.system(command)
        
        outFile = "%s/%s.txt" %(outDir_mod, sampleName_mod)
        
        #command = "dasgoclient -query=\"file dataset=%s\" > %s" %(sampleName, outFile)
        command = "ls /eos/cms/store/group/upgrade/RTB/Iter6/11_3/RAW/%s > %s"   %(sampleName_mod,outFile)
        print "Command:", command
        print ""
        os.system(command)
        
        fileContent = ""
        
        prefix2 = "%s/%s/ntuple"    %(prefix, sampleName_mod)
        print "Replacing \"%s\" with \"%s\" in file." %(toReplace, prefix2)
        print ""
        
        print "Number of lines:"
        os.system("wc -l %s" %(outFile))
        print ""
        
        with open(outFile, "r") as f :
            
            fileContent = f.read()
        
        fileContent = fileContent.replace(toReplace, prefix2)
        
        with open(outFile, "w") as f :
            
            f.write(fileContent)
