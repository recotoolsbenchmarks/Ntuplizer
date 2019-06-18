#include "TLorentzVector.h"
class lepInfo{
public:
  // methods                                                                                                                         
  lepInfo();
  ~lepInfo();

  // variables
  float px;
  float pz;
  float py;
  float pt;
  float eta;
  float phi;
  int   charge;
  float iso;
  int flavour;
};

lepInfo::lepInfo()
{}
lepInfo::~lepInfo()
{}

