#include "TROOT.h"
#include "TChain.h"
#include "TProof.h"
#include "TProofServ.h"
#include "TSystem.h"
#include "boost/config.hpp"
#include "boost/lexical_cast.hpp"
using namespace boost;

void runAll(string fileToOpen, string outfilename, float xs){
  
  TChain *fChain = new TChain("myana/mytree");

  ifstream file;
  file.open(fileToOpen.c_str(), ifstream::in );
  char filename[2000];
  while (true) {
    file >> filename;
    if( file.eof() ) break;
    fChain->Add(filename);   
    cout<<"Added "<<filename<<endl;
  }

  std::cout<<"Entries:" << fChain->GetEntries()<<std::endl;
  std::cout<<xs<<std::endl;
  //std::string weight= boost::lexical_cast<std::string>(xs/(float)fChain->GetEntries());
  //set it to 1 when running for efficiency
  std::string weight = "1.00";
  //std::cout<<weight<<std::endl;
  TProof *plite = TProof::Open("");
  fChain->SetProof();
  fChain->Process("SelectorClass.C+", (outfilename+" "+weight).c_str());
  //fChain->SetProof(0);

}
