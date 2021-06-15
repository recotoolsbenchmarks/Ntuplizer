# PhaseIIAnalysis
# Repository for collecting recipes for standard physics objects and for analysis of simulated events with the CMS phase 2 detector

The recipes for accessing PhaseII objects are summarized here:
[UPG PhaseII recipes](https://twiki.cern.ch/twiki/bin/view/CMS/UPG#PhaseII_FS_object_recipes "UPG PhaseII recipes")


This is used for egamma, muons and jets/met parameterization in 11_3_X.

Table of contents
=================

  * [Making a PR](#makingPR)
  * [Code Setup](#setup)

Making a PR
=====

```
cmsrel CMSSW_11_3_0_pre4
cd CMSSW_11_3_0_pre4/src/
cmsenv
git cms-merge-topic SohamBhattacharya:PhaseII_forRTB_11_3_0_pre4
scram b -j10

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init
mkdir new; cd new

```
If you do not attempt to contribute to this repository, simply clone it:
```
git clone https://github.com/recotoolsbenchmarks/RecoNtuplizer.git .
```

If you aim at contributing to the repository, you need to fork this repository (via the fork button) and then clone the forked repository:
```
git clone git@github.com:YOURGITUSERNAME/RecoNtuplizer.git .
```

Continue : 
```
cd ../; mv new/* .
cp -r new/.git .
rm -rf new
scram b -j10
cd TreeMaker/Ntuplzr/
git remote add upstream git@github.com:recotoolsbenchmarks/RecoNtuplizer.git
```

You can then regularly update your fork via:
```
git fetch upstream && git merge upstream/master
```

If you want to submit a new feature to ```recotoolsbenchmarks/RecoNtuplizer``` you have to do it via pull-request (PR):
So, first commit and push your changes to ```YOURGITUSERNAME/RecoNtuplizer``` and then make a PR via the github interface. 


Condor Setup
=====

Currently in 11.3, we have re-reco setup + Basic phase-2 full-sim setup. 


##############################################################################
To submit the re-reco condor setup, please first generate your proxy and copy it in directory:  
$HOME/proxies/ 

Then in CMSSW_11_3_0_pre4/src : 

To get the files corresponding to a given dataset: 

```
mkdir sourceFiles
```

Add the corresponding dataset accrordingly in getDASdataFiles.py in the list l_samplename
and run using : 
```
python getDASdataFiles.py --getFiles
```

Then in temp.txt: add the processName, path to input files, path to putput files etc..  as the given examples:
```
./temp.txt
```
##############################################################################



To submit the ntuplizer condor setup, again please first generate your proxy and copy it in directory:  
$HOME/proxies/ 

Then in CMSSW_11_3_0_pre4/src/TreeMaker/Ntuplzr/test: 

To get the files corresponding to a given dataset: 

```
mkdir sourceFiles
```

Add the path to corresponding re-reco files in getDASdataFiles.py 
and run using : 
```
python getDASdataFiles.py --getFiles
```

Then in temp.txt: add the processName, path to input files, path to putput files etc..  as before and run : 
```
./temp.txt
```


=====