//
// Builder of ECAL RecHits on GPU
//

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/DeclsForKernels.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"


namespace ecal { 
  namespace rechit {
 
    
    __global__
    void kernel_create_ecal_rehit(
      // configuration 
      int const* ChannelStatusToBeExcluded,
      uint32_t ChannelStatusToBeExcludedSize,                     
      // conditions
      float const* adc2gev,
      float const* intercalib,
      uint32_t const* status,
      float const* apdpnrefs,
      float const* alphas,
      // input for transparency corrections
      float const* p1,
      float const* p2,
      float const* p3,
      edm::TimeValue_t const* t1,
      edm::TimeValue_t const* t2,
      edm::TimeValue_t const* t3,  
      // input for linear corrections
      float const* lp1,
      float const* lp2,
      float const* lp3,
      edm::TimeValue_t const* lt1,
      edm::TimeValue_t const* lt2,
      edm::TimeValue_t const* lt3,                    
      // time, used for time dependent corrections
      edm::TimeValue_t const event_time,
      // input
      uint32_t const* did_eb,
      uint32_t const* did_ee,
      ::ecal::reco::StorageScalarType const* amplitude_eb,   // in adc counts  
      ::ecal::reco::StorageScalarType const* amplitude_ee,   // in adc counts  
      ::ecal::reco::StorageScalarType const* time_eb,   
      ::ecal::reco::StorageScalarType const* time_ee,   
      ::ecal::reco::StorageScalarType const* chi2_eb,   
      ::ecal::reco::StorageScalarType const* chi2_ee,   
      // output
      uint32_t *did,
      ::ecal::reco::StorageScalarType* energy,   // in energy [GeV]  
      ::ecal::reco::StorageScalarType* time,  
      ::ecal::reco::StorageScalarType* chi2,  
      uint32_t* flagBits,
      uint32_t* extra,
      int const nchannels,
      uint32_t const offsetForInput,
      uint32_t const offsetForHashes  
    );
      
    
    // host version, to be called by the plugin
    
    void create_ecal_rehit(
      EventInputDataGPU const& eventInputGPU,
      EventOutputDataGPU&      eventOutputGPU,
      //     eventDataForScratchGPU_,
      ConditionsProducts const& conditions, 
      ConfigurationParameters const& configParameters,
      uint32_t const offsetForInput, 
      edm::TimeValue_t const event_time,
      cuda::stream_t<>& cudaStream
    );
    
  }
  
}

