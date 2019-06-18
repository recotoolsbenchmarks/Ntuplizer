#cp SelectorClass_SIG.h SelectorClass.h 
#cp SelectorClass_SIG.C SelectorClass.C
#./script.sh DYToLL_M-50_14TeV_TuneCP5_pythia8 1


cp SelectorClass_QCD.h SelectorClass.h 
cp SelectorClass_QCD.C SelectorClass.C
./script.sh QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8 1
rm -f SelectorClass.C SelectorClass.h SelectorClass_C.d SelectorClass_C_ACLiC_dict_rdict.pcm SelectorClass_C.so


