#ifndef _ATLASROOTSTYLE_H_
#define _ATLASROOTSTYLE_H_

#include <iostream>
#include <TROOT.h>
#include <TStyle.h>

class AtlasROOTStyle{
public:
  AtlasROOTStyle();
  virtual ~AtlasROOTStyle(){}
  void SetStyle();

protected:

  TStyle* AtlasStyle();
};
#endif
