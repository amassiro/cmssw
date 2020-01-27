#ifndef RecoLocalCalo_EcalRecProducers_src_EcalChannelStatusGPU_h
#define RecoLocalCalo_EcalRecProducers_src_EcalChannelStatusGPU_h

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

#include <cuda/api_wrappers.h>

class EcalChannelStatusGPU {
public:
  struct Product {
    ~Product();
    uint16_t *status = nullptr;
  };
  
#ifndef __CUDACC__
  
  // 
  EcalChannelStatusGPU(EcalChannelStatus const&);
  
  // will call dealloation for Product thru ~Product
  ~EcalChannelStatusGPU() = default;
  
  // get device pointers
  Product const& getProduct(cuda::stream_t<>&) const;
  
  // 
  static std::string name() { return std::string{"ecalChannelStatusGPU"}; }
  
private:
  // in the future, we need to arrange so to avoid this copy on the host
  // store eb first then ee
  std::vector<uint16_t, CUDAHostAllocator<uint16_t>> status_;
    
  CUDAESProduct<Product> product_;
  
#endif
};


#endif
