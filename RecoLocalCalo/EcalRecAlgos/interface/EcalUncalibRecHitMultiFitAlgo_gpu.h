#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgoGpu_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgoGpu_HH

/** \class EcalUncalibRecHitMultiFitAlgoGpu
  * 
  *  Amplitude reconstucted by the multi-template fit on GPU!
  *
  */

// 
// #include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
// #include "FWCore/MessageLogger/interface/MessageLogger.h"
// 
// #include "CondFormats/EcalObjects/interface/EcalPedestals.h"
// #include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
// #include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS.h"
// 
// 
// #include "TMatrixDSym.h"
// #include "TVectorD.h"


#include <vector>

void EcalUncalibRecHitMultiFitAlgo_gpu_copy_run_return(std::vector<float>& h_vector_digis, std::vector<float>& d_vector_digis);

#endif
