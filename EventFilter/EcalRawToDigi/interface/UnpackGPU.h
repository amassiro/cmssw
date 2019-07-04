#ifndef EventFilter_EcalRawToDigi_interface_UnpackGPU_h
#define EventFilter_EcalRawToDigi_interface_UnpackGPU_h

#include "EventFilter/EcalRawToDigi/interface/DeclsForKernels.h"

#include <cuda/api_wrappers.h>

namespace ecal { namespace raw {

// FIXME: bundle up uint32_t values
void entryPoint(
        InputDataCPU const&, InputDataGPU&, 
        cuda::stream_t<>&, uint32_t const, uint32_t const);

}}

#endif // EventFilter_EcalRawToDigi_interface_UnpackGPU_h
