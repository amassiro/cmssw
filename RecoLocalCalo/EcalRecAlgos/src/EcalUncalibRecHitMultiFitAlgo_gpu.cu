#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo_gpu.h"

// #include "FWCore/MessageLogger/interface/MessageLogger.h"
// 
// #include "CondFormats/EcalObjects/interface/EcalPedestals.h"
// #include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
// 
// EcalUncalibRecHitMultiFitAlgo_gpu::EcalUncalibRecHitMultiFitAlgo_gpu() : 
//   _computeErrors(true),
//   _doPrefit(false),
//   _prefitMaxChiSq(1.0),
//   _dynamicPedestals(false),
//   _mitigateBadSamples(false),
//   _selectiveBadSampleCriteria(false),
//   _addPedestalUncertainty(0.),
//   _simplifiedNoiseModelForGainSwitch(true),
//   _gainSwitchUseMaxSample(false){
//     
//   _singlebx.resize(1);
//   _singlebx << 0;
//   
//   _pulsefuncSingle.disableErrorCalculation();
//   _pulsefuncSingle.setMaxIters(1);
//   _pulsefuncSingle.setMaxIterWarnings(false);
//     
// }



#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS_gpu.h"

#include <iostream>


  
__global__
void kernel_reconstruct_test( 
                        float* vector_vector_digis,
                        float* vector_vector_amplitudes,
                        float* vector_chi2
                       ) {
  //---- idx identifies the channel (the crystal!)
  
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  
  //---- do local reconstruction
  for (int iSample = 0; iSample<10; iSample++) {
    vector_vector_amplitudes[ idx*10 + iSample ] = 1000.; //---- for tests
  }
  vector_chi2[idx] = 1.0;
    
}
  
  
  
  
  
  
  

__global__
void kernel_reconstruct( 
                        float* vector_vector_pulses,
                        float* vector_long_vector_correlation_matrix,
                        float* vector_long_vector_noise_matrix,
                        float* vector_vector_digis,
                        float* vector_vector_amplitudes,
                        float* vector_chi2
                       ) {
  //---- idx identifies the channel (the crystal!)
  
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  
  //---- multifit
  ecal::multifit::PulseChiSqSNNLS multifit_fitting_NNLS;
  
  //---- prepare inputs:
  //             amplitudes_input_adc,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples
  
  SampleVector amplitudes_input_adc;
  for (int iSample = 0; iSample<10; iSample++) {
    amplitudes_input_adc[iSample] = vector_vector_digis[ idx*10 + iSample ];
  }
  
  BXVector activeBX;
  activeBX.resize(10);
  activeBX << -5,-4,-3,-2,-1,0,1,2,3,4;
  
  SampleGainVector gainsPedestal; //-1 for static pedestal
  for (int iSample=0; iSample<10; iSample++) {
    gainsPedestal[iSample] = -1;
  }
  
  SampleGainVector badSamples = SampleGainVector::Zero();
  
//   static const int TEMPLATESAMPLES = 12;
  
  //---- vector_vector_pulses ---> fullpulse
  FullSampleVector fullpulse(FullSampleVector::Zero());
  for (int iSample=0; iSample<12; iSample++) {
    fullpulse(iSample+7) = vector_vector_pulses[ idx * 12 + iSample ];
  }
  
  
  //---- vector_long_vector_correlation_matrix --> fullpulsecov
  FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
  
  //            TEMPLATESAMPLES
  for(int i=0; i<12;i++) {
    for(int j=0; j<12;j++) {
      //       fullpulsecov(i+7,j+7) = aPulseCov->covval[i][j];
      fullpulsecov(i+7,j+7) = vector_long_vector_correlation_matrix[idx * 12*12 +  12 * i + j];
    }
  }
  
  
  
  //---- vector_long_vector_noise_matrix --> noisecov
  SampleMatrix noisecov;
  int nnoise = 10;
  for (int i=0; i<nnoise; ++i) {
    for (int j=0; j<nnoise; ++j) {
      noisecov(i,j) = vector_long_vector_noise_matrix[idx * (10*10) + 10 * i + j];
    }
  }
  
  
  
  
    
    
  
  // FIXME:      multifit_fitting_NNLS.disableErrorCalculation();
  
    bool status = multifit_fitting_NNLS.DoFit(amplitudes_input_adc,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples);
//   chisq = _pulsefunc.ChiSq();
//   
//   if (!status) {
//     // TODO: Can we do better
//     printf("EcalUncalibRecHitMultiFitAlgo::makeRecHit failed fit\n");
//     //      edm::LogWarning("EcalUncalibRecHitMultiFitAlgo::makeRecHit") << "Failed Fit" << std::endl;
//   }
//   
//   unsigned int ipulseintime = 0;
//   for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
//     if (_pulsefunc.BXs().coeff(ipulse)==0) {
//       ipulseintime = ipulse;
//       break;
//     }
//   }
//   
//   amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
//   amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;
  
  
  
    //---- dump the output into something usable
    if (status) {
      for (unsigned int iPulse=0; iPulse<multifit_fitting_NNLS.BXs().rows(); iPulse++) {
        int bx = multifit_fitting_NNLS.BXs().coeff(iPulse);
        vector_vector_amplitudes[ idx*10 + bx + 5 ] = multifit_fitting_NNLS.BXs().coeff(iPulse);
      }
    }
    else {
      for (unsigned int iPulse=0; iPulse<multifit_fitting_NNLS.BXs().rows(); iPulse++) {
        vector_vector_amplitudes[ idx*10 + iPulse ] = -1.; // 0.;
      }
    }
    vector_chi2[idx] = multifit_fitting_NNLS.ChiSq();
    
  
}
  

  

void EcalUncalibRecHitMultiFitAlgo_gpu_copy_run_return(
                                                       int numRechits,
                                                       float* h_vector_vector_pulses,                       float* d_vector_vector_pulses, 
                                                       float* h_vector_long_vector_correlation_matrix,      float* d_vector_long_vector_correlation_matrix,
                                                       float* h_vector_long_vector_noise_matrix,            float* d_vector_long_vector_noise_matrix,
                                                       float* h_vector_vector_digis,                        float* d_vector_vector_digis,
                                                       float* h_vector_vector_amplitudes,                   float* d_vector_vector_amplitudes,
                                                       float* h_vector_chi2,                                float* d_vector_chi2
  ) {
  
  //---- copy on device
  
  //                                              likely 61200     EcalDataFrame::MAXSAMPLES * float
//   cudaMemcpy(&d_vector_vector_digis,      &h_vector_vector_digis,      61200 * 10 * sizeof(float), cudaMemcpyHostToDevice);
//   cudaMemcpy(&d_vector_vector_amplitudes, &h_vector_vector_amplitudes, 61200 * 10 * sizeof(float), cudaMemcpyHostToDevice);
//   cudaMemcpy(&d_vector_vector_pulses,      &h_vector_vector_pulses,      61200 * TEMPLATESAMPLES * sizeof(float), cudaMemcpyHostToDevice);


  cudaMemcpy(d_vector_vector_digis,                        h_vector_vector_digis,                      numRechits * 10   * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_vector_vector_amplitudes,                   h_vector_vector_amplitudes,                 numRechits * 10   * sizeof(float), cudaMemcpyHostToDevice);
  //                                                                                                   TEMPLATESAMPLES = 12
  cudaMemcpy(d_vector_vector_pulses,                       h_vector_vector_pulses,                     numRechits * 12   * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_vector_long_vector_correlation_matrix,      h_vector_long_vector_correlation_matrix,    numRechits * 12 *12  * sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(d_vector_long_vector_noise_matrix,            h_vector_long_vector_noise_matrix,          numRechits * 10 *10  * sizeof(float), cudaMemcpyHostToDevice);
  
  
  
  //---- then send all the multifits execution commands on the device
//   int nthreads_per_block = 256;
// //   int nblocks = (h_data.digis->size() + nthreads_per_block - 1) / nthreads_per_block;
//   int nblocks = (61200) / nthreads_per_block;
  
  
  int nthreads_per_block = 64;
  if (numRechits<64) nthreads_per_block = numRechits;
  int nblocks = (numRechits) / nthreads_per_block + 1;
  
  std::cout << " numRechits = " << numRechits << " --> nthreads_per_block (= " << nthreads_per_block << ") *  nblocks (=" << nblocks << ") = " << nthreads_per_block*nblocks << std::endl;
  
  for (int ixtal=0; ixtal<numRechits; ixtal++){
    std::cout << " ~~~ ixtal = " << ixtal << std::endl;
    for (int i=0; i<10; i++) {
      std::cout << " ~~~ ~~~ h_vector_vector_digis[ " << ixtal*10+i << " ] = " << h_vector_vector_digis[ixtal*10+i] << std::endl;
    }
  }
  
  kernel_reconstruct <<< nblocks, nthreads_per_block >>> ( 
                                                          d_vector_vector_pulses,
                                                          d_vector_long_vector_correlation_matrix,
                                                          d_vector_long_vector_noise_matrix,
                                                          d_vector_vector_digis,
                                                          d_vector_vector_amplitudes,
                                                          d_vector_chi2
                                                         );
  
  //---- wait that all the threads are done ... wait ...
  //---- in order words, synchronize
  cudaDeviceSynchronize();


  //---- now copy back from device to host the results of the multifit
  //---- and then you can access directly to the outputs of the local reconstruction, namely:
  //----    all amplitudes (in time and out-of-time, chi2, ...

  cudaMemcpy(h_vector_vector_amplitudes,      d_vector_vector_amplitudes,      61200 * 10 * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_vector_chi2,                   d_vector_chi2,                   61200 * sizeof(float),      cudaMemcpyDeviceToHost);
  
}




void EcalUncalibRecHitMultiFitAlgo_gpu_copy_run_return_test( float& h_vector_vector_digis,      float& d_vector_vector_digis,
                                                             float& h_vector_vector_amplitudes, float& d_vector_vector_amplitudes,
                                                             float &h_vector_chi2, float &d_vector_chi2
  ) {
  
  //---- copy on device
  
  //                                              likely 61200     EcalDataFrame::MAXSAMPLES * float
  cudaMemcpy(&d_vector_vector_digis,      &h_vector_vector_digis,      61200 * 10 * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(&d_vector_vector_amplitudes, &h_vector_vector_amplitudes, 61200 * 10 * sizeof(float), cudaMemcpyHostToDevice);
  
  
  //---- then send all the multifits execution commands on the device
  int nthreads_per_block = 256;
//   int nblocks = (h_data.digis->size() + nthreads_per_block - 1) / nthreads_per_block;
  int nblocks = (61200) / nthreads_per_block;

  // TEST FIXME
//   nblocks = 2;
//   nthreads_per_block = 4;
  
  kernel_reconstruct_test <<< nblocks, nthreads_per_block >>> ( 
                                                          &d_vector_vector_digis,
                                                          &d_vector_vector_amplitudes,
                                                          &d_vector_chi2
                                                         );
  
  //---- wait that all the threads are done ... wait ...
  //---- in order words, synchronize
  cudaDeviceSynchronize();


  //---- now copy back from device to host the results of the multifit
  //---- and then you can access directly to the outputs of the local reconstruction, namely:
  //----    all amplitudes (in time and out-of-time, chi2, ...

  cudaMemcpy(&h_vector_vector_amplitudes, &d_vector_vector_amplitudes, 61200 * 10 * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(&h_vector_chi2,              &d_vector_chi2,              61200 * sizeof(float),      cudaMemcpyDeviceToHost);
  
  
  
  
  
  
}






// /// compute rechits
// EcalUncalibratedRecHit EcalUncalibRecHitMultiFitAlgo_gpu::makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrixGainArray &noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX) {
// 
//   uint32_t flags = 0;
//   
//   const unsigned int nsample = EcalDataFrame::MAXSAMPLES;
//   
//   double maxamplitude = -std::numeric_limits<double>::max();
//   const unsigned int iSampleMax = 5;
//   const unsigned int iFullPulseMax = 9;
//   
//   double pedval = 0.;
//     
//   SampleVector amplitudes;
//   SampleGainVector gainsNoise;
//   SampleGainVector gainsPedestal;
//   SampleGainVector badSamples = SampleGainVector::Zero();
//   bool hasSaturation = dataFrame.isSaturated();
//   bool hasGainSwitch = hasSaturation || dataFrame.hasSwitchToGain6() || dataFrame.hasSwitchToGain1();
//   
//   //no dynamic pedestal in case of gain switch, since then the fit becomes too underconstrained
//   bool dynamicPedestal = _dynamicPedestals && !hasGainSwitch;
//   
//   for(unsigned int iSample = 0; iSample < nsample; iSample++) {
//         
//     const EcalMGPASample &sample = dataFrame.sample(iSample);
//     
//     double amplitude = 0.;
//     int gainId = sample.gainId();
//     
//     double pedestal = 0.;
//     double gainratio = 1.;
//     
//     if (gainId==0 || gainId==3) {
//       pedestal = aped->mean_x1;
//       gainratio = aGain->gain6Over1()*aGain->gain12Over6();
//       gainsNoise[iSample] = 2;
//       gainsPedestal[iSample] = dynamicPedestal ? 2 : -1;  //-1 for static pedestal
//     }
//     else if (gainId==1) {
//       pedestal = aped->mean_x12;
//       gainratio = 1.;
//       gainsNoise[iSample] = 0;
//       gainsPedestal[iSample] = dynamicPedestal ? 0 : -1; //-1 for static pedestal
//     }
//     else if (gainId==2) {
//       pedestal = aped->mean_x6;
//       gainratio = aGain->gain12Over6();
//       gainsNoise[iSample] = 1;
//       gainsPedestal[iSample] = dynamicPedestal ? 1 : -1; //-1 for static pedestals
//     }
// 
//     if (dynamicPedestal) {
//       amplitude = (double)(sample.adc())*gainratio;
//     }
//     else {
//       amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
//     }
//     
//     if (gainId == 0) {
//        edm::LogError("EcalUncalibRecHitMultiFitAlgo_gpu")<< "Saturation encountered.  Multifit is not intended to be used for saturated channels.";
//       //saturation
//       if (dynamicPedestal) {
//         amplitude = 4095.*gainratio;
//       }
//       else {
//         amplitude = (4095. - pedestal) * gainratio;
//       }
//     }
//         
//     amplitudes[iSample] = amplitude;
//     
//     if (iSample==iSampleMax) {
//       maxamplitude = amplitude;
//       pedval = pedestal;
//     }
//         
//   }
// 
//   double amplitude, amperr, chisq;
//   bool status = false;
//     
//   //special handling for gain switch, where sample before maximum is potentially affected by slew rate limitation
//   //optionally apply a stricter criteria, assuming slew rate limit is only reached in case where maximum sample has gain switched but previous sample has not
//   //option 1: use simple max-sample algorithm
//   if (hasGainSwitch && _gainSwitchUseMaxSample) {
//     double maxpulseamplitude = maxamplitude / fullpulse[iFullPulseMax];
//     EcalUncalibratedRecHit rh( dataFrame.id(), maxpulseamplitude, pedval, 0., 0., flags );
//     rh.setAmplitudeError(0.);
//     for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
//       int bx = _pulsefunc.BXs().coeff(ipulse);
//       if (bx!=0) {
//         rh.setOutOfTimeAmplitude(bx+5, 0.0);
//       }
//     }
//     return rh;
//   }
// 
//   //option2: A floating negative single-sample offset is added to the fit
//   //such that the affected sample is treated only as a lower limit for the true amplitude
//   bool mitigateBadSample = _mitigateBadSamples && hasGainSwitch && iSampleMax>0;
//   mitigateBadSample &= (!_selectiveBadSampleCriteria || (gainsNoise.coeff(iSampleMax-1)!=gainsNoise.coeff(iSampleMax)) );
//   if (mitigateBadSample) {
//     badSamples[iSampleMax-1] = 1;
//   }
//   
//   //compute noise covariance matrix, which depends on the sample gains
//   SampleMatrix noisecov;
//   if (hasGainSwitch) {
//     std::array<double,3> pedrmss = {{aped->rms_x12, aped->rms_x6, aped->rms_x1}};
//     std::array<double,3> gainratios = {{ 1., aGain->gain12Over6(), aGain->gain6Over1()*aGain->gain12Over6()}};
//     if (_simplifiedNoiseModelForGainSwitch) {
//       int gainidxmax = gainsNoise[iSampleMax];
//       noisecov = gainratios[gainidxmax]*gainratios[gainidxmax]*pedrmss[gainidxmax]*pedrmss[gainidxmax]*noisecors[gainidxmax];
//       if (!dynamicPedestal && _addPedestalUncertainty>0.) {
//         //add fully correlated component to noise covariance to inflate pedestal uncertainty
//         noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
//       }
//     }
//     else {
//       noisecov = SampleMatrix::Zero();
//       for (unsigned int gainidx=0; gainidx<noisecors.size(); ++gainidx) {
//         SampleGainVector mask = gainidx*SampleGainVector::Ones();
//         SampleVector pedestal = (gainsNoise.array()==mask.array()).cast<SampleVector::value_type>();
//         if (pedestal.maxCoeff()>0.) {
//           //select out relevant components of each correlation matrix, and assume no correlation between samples with
//           //different gain
//           noisecov += gainratios[gainidx]*gainratios[gainidx]*pedrmss[gainidx]*pedrmss[gainidx]*pedestal.asDiagonal()*noisecors[gainidx]*pedestal.asDiagonal();
//           if (!dynamicPedestal && _addPedestalUncertainty>0.) {
//             //add fully correlated component to noise covariance to inflate pedestal uncertainty
//             noisecov += gainratios[gainidx]*gainratios[gainidx]*_addPedestalUncertainty*_addPedestalUncertainty*pedestal.asDiagonal()*SampleMatrix::Ones()*pedestal.asDiagonal();
//           }
//         }
//       }
//     }
//   }
//   else {
//     noisecov = aped->rms_x12*aped->rms_x12*noisecors[0];
//     if (!dynamicPedestal && _addPedestalUncertainty>0.) {
//       //add fully correlated component to noise covariance to inflate pedestal uncertainty
//       noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
//     }
//   }
//   
//   //optimized one-pulse fit for hlt
//   bool usePrefit = false;
//   if (_doPrefit) {
//     status = _pulsefuncSingle.DoFit(amplitudes,noisecov,_singlebx,fullpulse,fullpulsecov,gainsPedestal,badSamples);
//     amplitude = status ? _pulsefuncSingle.X()[0] : 0.;
//     amperr = status ? _pulsefuncSingle.Errors()[0] : 0.;
//     chisq = _pulsefuncSingle.ChiSq();
//     
//     if (chisq < _prefitMaxChiSq) {
//       usePrefit = true;
//     }
//   }
//   
//   if (!usePrefit) {
//   
//     if(!_computeErrors) _pulsefunc.disableErrorCalculation();
//     status = _pulsefunc.DoFit(amplitudes,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples);
//     chisq = _pulsefunc.ChiSq();
//     
//     if (!status) {
//       edm::LogWarning("EcalUncalibRecHitMultiFitAlgo_gpu::makeRecHit") << "Failed Fit" << std::endl;
//     }
// 
//     unsigned int ipulseintime = 0;
//     for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
//       if (_pulsefunc.BXs().coeff(ipulse)==0) {
//         ipulseintime = ipulse;
//         break;
//       }
//     }
//     
//     amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
//     amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;
//   
//   }
//   
//   double jitter = 0.;
//   
//   EcalUncalibratedRecHit rh( dataFrame.id(), amplitude , pedval, jitter, chisq, flags );
//   rh.setAmplitudeError(amperr);
//   
//   if (!usePrefit) {
//     for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
//       int bx = _pulsefunc.BXs().coeff(ipulse);
//       if (bx!=0 && std::abs(bx)<100) {
//         rh.setOutOfTimeAmplitude(bx+5, status ? _pulsefunc.X().coeff(ipulse) : 0.);
//       }
//       else if (bx==(100+gainsPedestal[iSampleMax])) {
//         rh.setPedestal(status ? _pulsefunc.X().coeff(ipulse) : 0.);
//       }
//     }
//   }
//   
//   return rh;
// }

