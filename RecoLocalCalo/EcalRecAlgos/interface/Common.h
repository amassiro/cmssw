#ifndef RecoLocalCalo_EcalRecAlgos_interface_Common_h
#define RecoLocalCalo_EcalRecAlgos_interface_Common_h

#include <cmath>

#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

// a workaround for std::abs not being a constexpr function
namespace ecal {

template<typename T>
constexpr T abs(T const& value) {
    return ::std::max(value, -value);
}

}

namespace ecal { namespace cuda {

void assert_if_error();

}}

namespace myMath {

constexpr float fast_expf(float x) { return unsafe_expf<6>(x); }
constexpr float fast_logf(float x) { return unsafe_logf<7>(x); }

}

#endif
