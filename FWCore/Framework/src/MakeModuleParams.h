#ifndef FWCore_Framework_MakeModuleParams_h
#define FWCore_Framework_MakeModuleParams_h

/** ----------------------

This struct is used to communication parameters into the module factory.

---------------------- **/

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "boost/shared_ptr.hpp"

#include <string>

namespace edm {
  class ProcessConfiguration;
  class ProductRegistry;

  struct MakeModuleParams {
    MakeModuleParams() :
      pset_(nullptr), reg_(nullptr), processConfiguration_()
      {}

    MakeModuleParams(ParameterSet* pset,
                     ProductRegistry& reg,
                     boost::shared_ptr<ProcessConfiguration const> processConfiguration) :
    pset_(pset),
      reg_(&reg),
      processConfiguration_(processConfiguration) {}

    ParameterSet* pset_;
    ProductRegistry* reg_;
    boost::shared_ptr<ProcessConfiguration const> processConfiguration_;
  };
}

#endif
