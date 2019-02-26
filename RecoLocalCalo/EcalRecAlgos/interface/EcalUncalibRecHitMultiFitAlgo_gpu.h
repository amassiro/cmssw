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

void EcalUncalibRecHitMultiFitAlgo_gpu_copy_run_return_test(
                                                       float& h_vector_digis,      float& d_vector_digis,
                                                       float& h_vector_amplitudes, float& d_vector_amplitudes,
                                                       float &h_vector_chi2, float &d_vector_chi2
     );

void EcalUncalibRecHitMultiFitAlgo_gpu_copy_run_return(
                                                       int numRechits,
                                                       float& h_vector_vector_pulses,                       float& d_vector_vector_pulses,
                                                       float& h_vector_long_vector_correlation_matrix,      float& d_vector_long_vector_correlation_matrix,
                                                       float& h_vector_long_vector_noise_matrix,            float& d_vector_long_vector_noise_matrix,
                                                       float& h_vector_digis,      float& d_vector_digis,
                                                       float& h_vector_amplitudes, float& d_vector_amplitudes,
                                                       float &h_vector_chi2, float &d_vector_chi2
     );

#endif
