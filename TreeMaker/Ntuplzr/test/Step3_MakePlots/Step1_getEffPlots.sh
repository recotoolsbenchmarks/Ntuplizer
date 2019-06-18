#!/bin/bash 
#root -l -b -q "plotIt_Eff.C(file,outdir,pu,sig,optstat)"
root -l -b -q "plotIt_Eff.C(\"../Step2_PostAN/DYToLL_M-50_14TeV_TuneCP5_pythia8.root\", \"/eos/user/s/sandhya/www/RTB/Iter1/\",\"200\",1,0)"
root -l -b -q "plotIt_Eff.C(\"../Step2_PostAN/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8.root\", \"/eos/user/s/sandhya/www/RTB/Iter1/\",\"200\",0,0)"
#root -l -b -q "plotROC_new.C(\"../Step2_PostAN/DYToLL_M-50_14TeV_TuneCP5_pythia8.root\",\"../Step2_PostAN/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8.root\", \"/eos/user/s/sandhya/www/RTB/Iter1/\")"
#root -l -b -q "plotROC_new.C(\"DYToLL-M-50_nJ_14TeV-madgraphMLM-pythia8.root\",\"QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8.root\", \"/eos/user/s/sandhya/www/VBF/HGCal/check_OldtrackIso/\")"




