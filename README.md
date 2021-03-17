# PhaseIIAnalysis
# Repository for collecting recipes for standard physics objects and for analysis of simulated events with the CMS phase 2 detector

The recipes for accessing PhaseII objects are summarized here:

[UPG PhaseII recipes](https://twiki.cern.ch/twiki/bin/view/CMS/PhaseIIFSObjectRecipes "UPG PhaseII recipes")

Table of contents
=================

  * [Making a PR](#makingPR)
  * [Code Setup](#setup)

Making a PR
=====

```
cmsrel CMSSW_11_0_1
cd CMSSW_11_0_1/src/
cmsenv
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag
git cms-merge-topic emilbols:BTV_11_0_X

git clone -b Phase2_11_0_X_v1.02 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

mkdir new; cd new
```
If you do not attempt to contribute to this repository, simply clone it:
```
git clone -b 11_0_X https://github.com/recotoolsbenchmarks/RecoNtuplizer.git .
```

If you aim at contributing to the repository, you need to fork this repository (via the fork button) and then clone the forked repository:
```
git clone -b 11_0_X git@github.com:YOURGITUSERNAME/RecoNtuplizer.git .
cd ../; mv new/* .
cp -r new/.git .
rm -rf new
scram b -j 8
cd TreeMaker/Ntuplzr/
git remote add upstream git@github.com:recotoolsbenchmarks/RecoNtuplizer.git
```

You can then regularly update your fork via:
```
git fetch upstream && git merge upstream/master
```

If you want to submit a new feature to ```recotoolsbenchmarks/RecoNtuplizer``` you have to do it via pull-request (PR):
So, first commit and push your changes to ```YOURGITUSERNAME/RecoNtuplizer``` and then make a PR via the github interface. 


Code Setup
=====

####################################################################################
```
cmsrel CMSSW_11_0_1
cd CMSSW_11_0_1/src/
cmsenv
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git cms-addpkg RecoBTag
git cms-merge-topic emilbols:BTV_11_0_X

git clone -b Phase2_11_0_X_v1.02 --depth 1 https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements

mkdir new; cd new
git clone -b 11_0_X https://github.com/recotoolsbenchmarks/RecoNtuplizer.git .
cd ../; mv new/* .
cp -r new/.git .
rm -rf new
scram b -j 8
cd TreeMaker/Ntuplzr/
```
####################################################################################
Basic Setup : Here in plugins/Ntuplzr.cc is the EDAnalyzer that makes the ntuples, takes in information 
from python/Ntuplzr_cfi.py. These default parameters can be modified in :
test/myproduceNtuples_cfg.py

to run this file, do :  
```
cmsRun test/myproduceNtuples_cfg.py maxEvents=10 outputFile=file.root
```

###################################################################################
STEP 1 : to get the main ntuples from this Ntuplzer using test/myproduceNtuples_cfg.py, a crab setup is there. one can update 
        various parameters in test/Step1_crab/submitCrabJobs_cfgparams.py according to their choice. 
	To crab-submit, do : 

```
cd test/Step1_crab
python submitCrabJobs_cfgparams.py 
```
###################################################################################
To check on the status or resubmit the jobs, change the file Resubmit.csh accordingly and do :  
```
source Resubmit.csh
```

###################################################################################
STEP 2 : Once the main ntuples are there, C++ classes along with a proof wrapper is written to be able to run the jobs 
parallely on the system. 

Here we have multiple classes based on what needs to be done. 
SelectorClass_SIG.C/h : Mainly to plot the efficiency on a given signal sample 
SelectorClass_BKG.C/h : Mainly to plot the fake-rate  on a given bkground sample 


To get this going on the root files stored in: 
/eos/cms/store/group/upgrade/RTB/

Step1_runcreateList.sh creates the input list of files for a given sample through a script createList.sh.
And then Step2_runScript_runAll_eff.sh runs over the given class through the wrapper code of runAll.C.
To Run : 

```
cd test/Step2_PostAN
./Step1_runcreateList.sh
./Step2_runScript_runAll_eff.sh
```

This gives the root files with all the information as askd from the class that was run. 

###################################################################################

STEP 3: This is where one makes the final histograms and puts them on web page or wherever u want them to put using the proper plotter. 
for the efficiency or fake rate plots, do : 

```
cd test/Step3_MakePlots
./Step1_getEffPlots.sh
``` 	   
###################################################################################

