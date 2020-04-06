#include <iostream>
#include <limits>

#include "cuda.h"

#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/Common.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"

#include "DataFormats/Math/interface/approx_exp.h"
#include "DataFormats/Math/interface/approx_log.h"

#include "TimeComputationKernels.h"
#include "KernelHelpers.h"

//#define DEBUG

//#define ECAL_RECO_CUDA_DEBUG

namespace ecal { namespace multifit {

__device__
__forceinline__
bool use_sample(unsigned int sample_mask, unsigned int sample) {
    return sample_mask & (0x1 << (EcalDataFrame::MAXSAMPLES - (sample + 1)));
}

__global__
void kernel_time_compute_nullhypot(SampleVector::Scalar const* sample_values,
                                   SampleVector::Scalar const* sample_value_errors,
                                   bool const* useless_sample_values,
                                   SampleVector::Scalar* chi2s,
                                   SampleVector::Scalar* sum0s,
                                   SampleVector::Scalar* sumAAs,
                                   int const nchannels) {
    using ScalarType = SampleVector::Scalar;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int tx = threadIdx.x + blockDim.x*blockIdx.x;
    int ltx = threadIdx.x;
    int ch = tx / nsamples;
    int nchannels_per_block = blockDim.x / nsamples;

    // TODO: make sure that this branch plays nicely with __syncthreads inside
    // can there be a deadlock even if the thread is inactive
    if (ch < nchannels) {
        // 
        int sample = tx % nsamples;

        // shared mem inits
        extern __shared__ char sdata[];
        char* s_sum0 = sdata;
        SampleVector::Scalar* s_sum1 = reinterpret_cast<SampleVector::Scalar*>(
            s_sum0 + nchannels_per_block*nsamples);
        SampleVector::Scalar* s_sumA = s_sum1 + nchannels_per_block*nsamples;
        SampleVector::Scalar* s_sumAA = s_sumA + nchannels_per_block*nsamples;

        // TODO make sure no div by 0
        auto const inv_error = useless_sample_values[tx] 
            ? 0.0 
            : 1.0 / (sample_value_errors[tx] * sample_value_errors[tx]);
        auto const sample_value = sample_values[tx];
        s_sum0[ltx] = useless_sample_values[tx] ? 0 : 1;
        s_sum1[ltx] = inv_error;
        s_sumA[ltx] = sample_value * inv_error;
        s_sumAA[ltx] = sample_value * sample_value * inv_error;
        __syncthreads();

        // 5 threads for [0, 4] samples
        if (sample<5) {
            s_sum0[ltx] += s_sum0[ltx+5];
            s_sum1[ltx] += s_sum1[ltx+5];
            s_sumA[ltx] += s_sumA[ltx+5];
            s_sumAA[ltx] += s_sumAA[ltx+5];
        }
        __syncthreads();

        if (sample<2) {
            // note double counting of sample 3
            s_sum0[ltx] += s_sum0[ltx+2] + s_sum0[ltx+3];
            s_sum1[ltx] += s_sum1[ltx+2] + s_sum1[ltx+3];
            s_sumA[ltx] += s_sumA[ltx+2] + s_sumA[ltx+3];
            s_sumAA[ltx] += s_sumAA[ltx+2] + s_sumAA[ltx+3];
        }
        __syncthreads();

        if (sample == 0) {
            // note, subtract to remove the double counting of sample == 3
            //s_sum0[ltx] += s_sum0[ltx+1] - s_sum0[ltx+3];
            //s_sum1[ltx] += s_sum1[ltx+1] - s_sum1[ltx+3];
            //s_sumA[ltx] += s_sumA[ltx+1] - s_sumA[ltx+3];
            //s_sumAA[ltx] += s_sumAA[ltx+1] - s_sumAA[ltx+3];
            auto const sum0 = s_sum0[ltx] + s_sum0[ltx+1] - s_sum0[ltx+3];
            auto const sum1 = s_sum1[ltx] + s_sum1[ltx+1] - s_sum1[ltx+3];
            auto const sumA = s_sumA[ltx] + s_sumA[ltx+1] - s_sumA[ltx+3];
            auto const sumAA = s_sumAA[ltx] + s_sumAA[ltx+1] - s_sumAA[ltx+3];
            auto const chi2 = sum0>0 
                ? (sumAA - sumA * sumA / sum1) / sum0
                : static_cast<ScalarType>(0);
            chi2s[ch] = chi2;
            sum0s[ch] = sum0;
            sumAAs[ch] = sumAA;

#ifdef DEBUG_TC_NULLHYPOT
            if (ch == 0) {
                printf("chi2 = %f sum0 = %d sumAA = %f\n",
                    chi2, static_cast<int>(sum0), sumAA);
            }
#endif
        }
    }
}

constexpr float fast_expf(float x) { return unsafe_expf<6>(x); }
constexpr float fast_logf(float x) { return unsafe_logf<7>(x); }

//#define DEBUG_TC_MAKERATIO
//
// launch ctx parameters are 
// 45 threads per channel, X channels per block, Y blocks
// 45 comes from: 10 samples for i <- 0 to 9 and for j <- i+1 to 9
// TODO: it might be much beter to use 32 threads per channel instead of 45
// to simplify the synchronization
//
__global__
void kernel_time_compute_makeratio(SampleVector::Scalar const* sample_values,
                                   SampleVector::Scalar const* sample_value_errors,
                                   uint32_t const* dids_eb,
                                   uint32_t const* dids_ee,
                                   bool const* useless_sample_values,
                                   char const* pedestal_nums,
                                   ConfigurationParameters::type const* amplitudeFitParametersEB,
                                   ConfigurationParameters::type const* amplitudeFitParametersEE,
                                   ConfigurationParameters::type const* timeFitParametersEB,
                                   ConfigurationParameters::type const* timeFitParametersEE,
                                   SampleVector::Scalar const* sumAAsNullHypot,
                                   SampleVector::Scalar const* sum0sNullHypot,
                                   SampleVector::Scalar* tMaxAlphaBetas,
                                   SampleVector::Scalar* tMaxErrorAlphaBetas,
                                   SampleVector::Scalar* g_accTimeMax,
                                   SampleVector::Scalar* g_accTimeWgt,
                                   TimeComputationState* g_state,
                                   unsigned int const timeFitParameters_sizeEB,
                                   unsigned int const timeFitParameters_sizeEE,
                                   ConfigurationParameters::type const timeFitLimits_firstEB,
                                   ConfigurationParameters::type const timeFitLimits_firstEE,
                                   ConfigurationParameters::type const timeFitLimits_secondEB,
                                   ConfigurationParameters::type const timeFitLimits_secondEE,
                                   int const nchannels,
                                   uint32_t const offsetForInputs) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr int nthreads_per_channel = 45; // n=10, n(n-1)/2
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockDim.x*blockIdx.x;
    int const ch = gtx / nthreads_per_channel;
    int const lch = threadIdx.x / nthreads_per_channel;
    int const ltx = threadIdx.x % nthreads_per_channel;
    int const ch_start = ch*nsamples;
    int const lch_start = lch*nthreads_per_channel;
    int const nchannels_per_block = blockDim.x / nthreads_per_channel;
    auto const* dids = ch >= offsetForInputs
        ? dids_ee
        : dids_eb;
    int const inputCh = ch >= offsetForInputs
        ? ch - offsetForInputs
        : ch;
    
    // rmeove inactive threads
    // TODO: need to understand if this is 100% safe in presence of syncthreads
    if (ch >= nchannels) return;

    auto const did = DetId{dids[inputCh]};
    auto const isBarrel = did.subdetId() == EcalBarrel;
    auto const* amplitudeFitParameters = isBarrel
        ? amplitudeFitParametersEB
        : amplitudeFitParametersEE;
    auto const* timeFitParameters = isBarrel
        ? timeFitParametersEB
        : timeFitParametersEE;
    auto const timeFitParameters_size = isBarrel
        ? timeFitParameters_sizeEB
        : timeFitParameters_sizeEE;
    auto const timeFitLimits_first = isBarrel
        ? timeFitLimits_firstEB
        : timeFitLimits_firstEE;
    auto const timeFitLimits_second = isBarrel
        ? timeFitLimits_secondEB
        : timeFitLimits_secondEE;

    extern __shared__ char smem[];
    ScalarType* shr_chi2s = reinterpret_cast<ScalarType*>(smem);
    ScalarType* shr_time_wgt = shr_chi2s + blockDim.x;
    ScalarType* shr_time_max = shr_time_wgt + blockDim.x;
    ScalarType* shrTimeMax = shr_time_max + blockDim.x;
    ScalarType* shrTimeWgt = shrTimeMax + blockDim.x;

    // map tx -> (sample_i, sample_j)
    int sample_i, sample_j = 0;
    if (ltx>=0 && ltx<=8) {
        sample_i = 0;
        sample_j = 1+ltx;
    } else if (ltx<=16) {
        sample_i = 1;
        sample_j = 2+ltx-9;
    } else if (ltx<=23) {
        sample_i = 2;
        sample_j = 3 + ltx - 17;
    } else if (ltx<=29) {
        sample_i = 3;
        sample_j = 4 + ltx - 24;
    } else if (ltx<=34) {
        sample_i = 4;
        sample_j = 5 + ltx - 30;
    } else if (ltx<=38) {
        sample_i = 5;
        sample_j = 6 + ltx - 35;
    } else if (ltx<=41) {
        sample_i = 6;
        sample_j = 7 + ltx - 39;
    } else if (ltx<=43) {
        sample_i = 7;
        sample_j = 8 + ltx - 42;
    } else if (ltx <= 44) {
        sample_i = 8;
        sample_j = 9;
    } else
        assert(false);

    auto const tx_i = ch_start + sample_i;
    auto const tx_j = ch_start + sample_j;

    //
    // note, given the way we partition the block, with 45 threads per channel
    // we will end up with inactive threads which need to be dragged along
    // through the synching point
    // 
    /*
    bool const condToExit = ch >= nchannels
        ? true
        : useless_sample_values[tx_i] 
          || useless_sample_values[tx_j]
          || sample_values[tx_i]<=1 || sample_values[tx_j]<=1;
          */
    bool const condForUselessSamples = useless_sample_values[tx_i] 
        || useless_sample_values[tx_j]
        || sample_values[tx_i]<=1 || sample_values[tx_j]<=1;

    //
    // see cpu implementation for explanation
    // 
    ScalarType chi2 = std::numeric_limits<ScalarType>::max();
    ScalarType tmax = 0;
    ScalarType tmaxerr = 0;
    shrTimeMax[threadIdx.x] = 0;
    shrTimeWgt[threadIdx.x] = 0;
    bool internalCondForSkipping1 = true;
    bool internalCondForSkipping2 = true;
    if (!condForUselessSamples) {
        auto const rtmp = sample_values[tx_i] / sample_values[tx_j];
        auto const invampl_i = 1.0 / sample_values[tx_i];
        auto const relErr2_i = sample_value_errors[tx_i]*sample_value_errors[tx_i]*
            invampl_i*invampl_i;
        auto const invampl_j = 1.0 / sample_values[tx_j];
        auto const relErr2_j = sample_value_errors[tx_j]*sample_value_errors[tx_j]*
            invampl_j*invampl_j;
        auto const err1 = rtmp * rtmp * (relErr2_i + relErr2_j);
        auto err2 = sample_value_errors[tx_j]*
            (sample_values[tx_i] - sample_values[tx_j])*(invampl_j*invampl_j);
        // TODO non-divergent branch for a block if each block has 1 channel
        // otherwise non-divergent for groups of 45 threads
        // at this point, pedestal_nums[ch] can be either 0, 1 or 2
        if (pedestal_nums[ch]==2)
            err2 *= err2 * 0.5;
        auto const err3 = (0.289*0.289) * (invampl_j*invampl_j);
        auto const total_error = std::sqrt(err1 + err2 + err3);

        auto const alpha = amplitudeFitParameters[0];
        auto const beta = amplitudeFitParameters[1];
        auto const alphabeta = alpha * beta;
        auto const invalphabeta = 1.0 / alphabeta;

        // variables instead of a struct
        auto const ratio_index = sample_i;
        auto const ratio_step = sample_j - sample_i;
        auto const ratio_value = rtmp;
        auto const ratio_error = total_error;

        auto const rlim_i_j = fast_expf(
            static_cast<ScalarType>(sample_j - sample_i) / beta) - 0.001;
        internalCondForSkipping1 = !(total_error<1.0 && rtmp>0.001 && rtmp<rlim_i_j);
        if (!internalCondForSkipping1) {
            //
            // precompute.
            // in cpu version this was done conditionally
            // however easier to do it here (precompute) and then just filter out
            // if not needed
            // 
            auto const l_timeFitLimits_first = timeFitLimits_first;
            auto const l_timeFitLimits_second = timeFitLimits_second;
            if (ratio_step == 1
                && ratio_value >= l_timeFitLimits_first
                && ratio_value <= l_timeFitLimits_second) {

                auto const time_max_i = static_cast<ScalarType>(ratio_index);
                auto u = timeFitParameters[timeFitParameters_size - 1];
#pragma unroll
                for (int k=timeFitParameters_size-2; k>=0; k--)
                    u = u*ratio_value + timeFitParameters[k];

                auto du = (timeFitParameters_size - 1) *
                    (timeFitParameters[timeFitParameters_size - 1]);
                for (int k=timeFitParameters_size - 2; k>=1; k--)
                    du = du*ratio_value + k*timeFitParameters[k];

                auto const error2 = ratio_error * ratio_error * du * du;
                auto const time_max = error2 > 0
                    ? (time_max_i - u) / error2
                    : static_cast<ScalarType>(0);
                auto const time_wgt = error2 > 0
                    ? 1.0 / error2
                    : static_cast<ScalarType>(0);

                // store into shared mem
                // note, this name is essentially identical to the one used 
                // below. 
                shrTimeMax[threadIdx.x] = error2 > 0 ? time_max : 0;
                shrTimeWgt[threadIdx.x] = error2 > 0 ? time_wgt : 0;
            } else {
                shrTimeMax[threadIdx.x] = 0;
                shrTimeWgt[threadIdx.x] = 0;
            }

            // continue with ratios
            auto const stepOverBeta = static_cast<SampleVector::Scalar>(ratio_step) / beta;
            auto const offset = static_cast<SampleVector::Scalar>(ratio_index) + alphabeta;
            auto const rmin = std::max(ratio_value - ratio_error, 0.001);
            auto const rmax = std::min(ratio_value + ratio_error, 
                fast_expf(static_cast<SampleVector::Scalar>(ratio_step) / beta)
                - 0.001);
            auto const time1 = 
                offset - 
                ratio_step / 
                    (fast_expf((stepOverBeta - fast_logf(rmin)) / 
                                       alpha) - 1.0);
            auto const time2 = 
                offset - 
                ratio_step /
                    (fast_expf((stepOverBeta - fast_logf(rmax)) / 
                                       alpha) - 1.0);

            // set these guys
            tmax = 0.5 * (time1 + time2);
            tmaxerr = 0.5 * std::sqrt((time1 - time2) * (time1 - time2));
#ifdef DEBUG_TC_MAKERATIO
            if (ch == 1 || ch == 0)
                printf("ch = %d ltx = %d tmax = %f tmaxerr = %f time1 = %f time2 = %f offset = %f rmin = %f rmax = %f\n",
                    ch, ltx, tmax, tmaxerr, time1, time2, offset, rmin, rmax);
#endif

            SampleVector::Scalar sumAf = 0;
            SampleVector::Scalar sumff = 0;
            int const itmin = std::max(-1, static_cast<int>(std::floor(tmax - alphabeta)));
            auto loffset = (static_cast<ScalarType>(itmin) - tmax) * invalphabeta;
            // TODO: data dependence 
            for (int it = itmin+1; it<nsamples; it++) {
                loffset += invalphabeta;
                if (useless_sample_values[ch_start + it])
                    continue;
                auto const inverr2 = 1.0 / 
                    (sample_value_errors[ch_start + it]*sample_value_errors[ch_start + it]);
                auto const term1 = 1.0 + loffset;
                auto const f = (term1 > 1e-6)
                    ? fast_expf(alpha * (fast_logf(term1) - loffset))
                    : 0;
                sumAf += sample_values[ch_start+it] * (f * inverr2);
                sumff += f*(f*inverr2);
            }

            auto const sumAA = sumAAsNullHypot[ch];
            auto const sum0 = sum0sNullHypot[ch];
            chi2 = sumAA;
            ScalarType amp = 0;
            // TODO: sum0 can not be 0 below, need to introduce the check upfront
            if (sumff > 0) {
                chi2 = sumAA - sumAf * (sumAf / sumff);
                amp = sumAf / sumff;
            }
            chi2 /= sum0;

#ifdef DEBUG_TC_MAKERATIO
            if (ch == 1 || ch == 0)
                printf("ch = %d ltx = %d sumAf = %f sumff = %f sumAA = %f sum0 = %d tmax = %f tmaxerr = %f chi2 = %f\n",
                    ch, ltx, sumAf, sumff, sumAA, static_cast<int>(sum0), tmax, tmaxerr, chi2);
#endif

            if (chi2>0 && tmax>0 && tmaxerr>0)
                internalCondForSkipping2 = false;
            else
                chi2 = std::numeric_limits<ScalarType>::max();
        }
    }

    // store into smem
    shr_chi2s[threadIdx.x] = chi2;
    __syncthreads();

    // find min chi2 - quite crude for now
    // TODO validate/check
    char iter = nthreads_per_channel / 2 + nthreads_per_channel % 2;
    bool oddElements = nthreads_per_channel % 2;
#pragma unroll
    while (iter>=1) {
        if (ltx < iter)
            // for odd ns, the last guy will just store itself
            // exception is for ltx == 0 and iter==1
            shr_chi2s[threadIdx.x] = oddElements && (ltx==iter-1 && ltx>0)
                ? shr_chi2s[threadIdx.x] 
                : std::min(shr_chi2s[threadIdx.x], shr_chi2s[threadIdx.x+iter]);
        __syncthreads();
        oddElements = iter % 2;
        iter = iter==1 ? iter/2 : iter/2 + iter%2;
    }

    // filter out inactive or useless samples threads
    if (!condForUselessSamples && !internalCondForSkipping1 
            && !internalCondForSkipping2) {
        // min chi2, now compute weighted average of tmax measurements
        // see cpu version for more explanation
        auto const chi2min = shr_chi2s[threadIdx.x - ltx];
        auto const chi2Limit = chi2min + 1.0;
        auto const inverseSigmaSquared = 
            chi2 < chi2Limit
                ? 1.0 / (tmaxerr * tmaxerr)
                : 0.0;

#ifdef DEBUG_TC_MAKERATIO
        if (ch == 1 || ch == 0)
            printf("ch = %d ltx = %d chi2min = %f chi2Limit = %f inverseSigmaSquared = %f\n",
                ch, ltx, chi2min, chi2Limit, inverseSigmaSquared);
#endif

        // store into shared mem and run reduction
        // TODO: check if cooperative groups would be better
        // TODO: check if shuffling intrinsics are better
        shr_time_wgt[threadIdx.x] = inverseSigmaSquared;
        shr_time_max[threadIdx.x] = tmax * inverseSigmaSquared;
    } else {
        shr_time_wgt[threadIdx.x] = 0;
        shr_time_max[threadIdx.x] = 0;
    }
    __syncthreads();

    // reduce to compute time_max and time_wgt
    iter = nthreads_per_channel / 2 + nthreads_per_channel % 2;
    oddElements = nthreads_per_channel % 2;
#pragma unroll
    while (iter>=1) {
        if (ltx < iter) {
            shr_time_wgt[threadIdx.x] = oddElements && (ltx==iter-1 && ltx>0)
                ? shr_time_wgt[threadIdx.x]
                : shr_time_wgt[threadIdx.x] + shr_time_wgt[threadIdx.x+iter];
            shr_time_max[threadIdx.x] = oddElements && (ltx==iter-1 && ltx>0)
                ? shr_time_max[threadIdx.x]
                : shr_time_max[threadIdx.x] + shr_time_max[threadIdx.x+iter];
            shrTimeMax[threadIdx.x] = oddElements && (ltx==iter-1 && ltx>0)
                ? shrTimeMax[threadIdx.x]
                : shrTimeMax[threadIdx.x] + shrTimeMax[threadIdx.x+iter];
            shrTimeWgt[threadIdx.x] = oddElements && (ltx==iter-1 && ltx>0)
                ? shrTimeWgt[threadIdx.x]
                : shrTimeWgt[threadIdx.x] + shrTimeWgt[threadIdx.x+iter];
        }
        
        __syncthreads();
        oddElements = iter % 2;
        iter = iter==1 ? iter/2 : iter/2 + iter%2;
    }

    // load from shared memory the 0th guy (will contain accumulated values)
    // compute 
    // store into global mem
    if (ltx == 0) {
        auto const tmp_time_max = shr_time_max[threadIdx.x];
        auto const tmp_time_wgt = shr_time_wgt[threadIdx.x];

        // we are done if there number of time ratios is 0
        if (tmp_time_wgt==0 && tmp_time_max==0) {
            g_state[ch] = TimeComputationState::Finished;
            return ;
        }

        // no div by 0
        auto const tMaxAlphaBeta = tmp_time_max / tmp_time_wgt;
        auto const tMaxErrorAlphaBeta = 1.0 / std::sqrt(tmp_time_wgt);

        tMaxAlphaBetas[ch] = tMaxAlphaBeta;
        tMaxErrorAlphaBetas[ch] = tMaxErrorAlphaBeta;
        g_accTimeMax[ch] = shrTimeMax[threadIdx.x];
        g_accTimeWgt[ch] = shrTimeWgt[threadIdx.x];
        g_state[ch] = TimeComputationState::NotFinished;

#ifdef DEBUG_TC_MAKERATIO
            printf("ch = %d time_max = %f time_wgt = %f\n",
                ch, tmp_time_max, tmp_time_wgt);
            printf("ch = %d tMaxAlphaBeta = %f tMaxErrorAlphaBeta = %f timeMax = %f timeWgt = %f\n",
                ch, tMaxAlphaBeta, tMaxErrorAlphaBeta, 
                shrTimeMax[threadIdx.x],
                shrTimeWgt[threadIdx.x]);
#endif
    }
}

/// launch ctx parameters are 
/// 10 threads per channel, N channels per block, Y blocks
/// TODO: do we need to keep the state around or can be removed?!
//#define DEBUG_FINDAMPLCHI2_AND_FINISH
__global__
void kernel_time_compute_findamplchi2_and_finish(
        SampleVector::Scalar const* sample_values,
        SampleVector::Scalar const* sample_value_errors,
        uint32_t const* dids_eb,
        uint32_t const* dids_ee,
        bool const* useless_samples,
        SampleVector::Scalar const* g_tMaxAlphaBeta,
        SampleVector::Scalar const* g_tMaxErrorAlphaBeta,
        SampleVector::Scalar const* g_accTimeMax,
        SampleVector::Scalar const* g_accTimeWgt,
        ConfigurationParameters::type const* amplitudeFitParametersEB,
        ConfigurationParameters::type const* amplitudeFitParametersEE,
        SampleVector::Scalar const* sumAAsNullHypot,
        SampleVector::Scalar const* sum0sNullHypot,
        SampleVector::Scalar const* chi2sNullHypot,
        TimeComputationState* g_state,
        SampleVector::Scalar* g_ampMaxAlphaBeta,
        SampleVector::Scalar* g_ampMaxError,
        SampleVector::Scalar* g_timeMax,
        SampleVector::Scalar* g_timeError,
        int const nchannels,
        uint32_t const offsetForInputs) {
    using ScalarType = SampleVector::Scalar;

    // constants 
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockIdx.x*blockDim.x;
    int const ch = gtx / nsamples;
    int const sample = threadIdx.x % nsamples;
    int const ch_start = ch * nsamples;
    auto const* dids = ch >= offsetForInputs
        ? dids_ee
        : dids_eb;
    int const inputCh = ch >= offsetForInputs
        ? ch - offsetForInputs
        : ch;

    // configure shared mem
    // per block, we need #threads per block * 2 * sizeof(ScalarType)
    // we run with N channels per block
    extern __shared__ char smem[];
    ScalarType* shr_sumAf = reinterpret_cast<ScalarType*>(smem);
    ScalarType* shr_sumff = shr_sumAf + blockDim.x;

    if (ch >= nchannels) return;

    auto state = g_state[ch];
    auto const did = DetId{dids[inputCh]};
    auto const* amplitudeFitParameters = did.subdetId() == EcalBarrel
        ? amplitudeFitParametersEB
        : amplitudeFitParametersEE;


    // TODO is that better than storing into global and launching another kernel
    // for the first 10 threads
    if (state == TimeComputationState::NotFinished) {
        auto const alpha = amplitudeFitParameters[0];
        auto const beta = amplitudeFitParameters[1];
        auto const alphabeta = alpha * beta;
        auto const invalphabeta = 1.0 / alphabeta;
        auto const tMaxAlphaBeta = g_tMaxAlphaBeta[ch];
        auto const sample_value = sample_values[gtx];
        auto const sample_value_error = sample_value_errors[gtx];
        auto const inverr2 = useless_samples[gtx]
            ? static_cast<ScalarType>(0)
            : 1.0 / (sample_value_error * sample_value_error);
        auto const offset = (static_cast<ScalarType>(sample) - tMaxAlphaBeta) 
            * invalphabeta;
        auto const term1 = 1.0 + offset;
        auto const f = term1 > 1e-6 
            ? fast_expf(alpha * (fast_logf(term1) - offset))
            : static_cast<ScalarType>(0.0);
        auto const sumAf = sample_value * (f * inverr2);
        auto const sumff = f * (f * inverr2);

        // store into shared mem
        shr_sumAf[threadIdx.x] = sumAf;
        shr_sumff[threadIdx.x] = sumff;
    } else {
        shr_sumAf[threadIdx.x] = 0;
        shr_sumff[threadIdx.x] = 0;
    }
    __syncthreads();

    // reduce
    // unroll completely here (but hardcoded)
    if (sample<5) {
        shr_sumAf[threadIdx.x] += shr_sumAf[threadIdx.x+5];
        shr_sumff[threadIdx.x] += shr_sumff[threadIdx.x+5];
    }
    __syncthreads();

    if (sample<2) {
        // will need to subtract for ltx = 3, we double count here
        shr_sumAf[threadIdx.x] += shr_sumAf[threadIdx.x+2] 
            + shr_sumAf[threadIdx.x+3];
        shr_sumff[threadIdx.x] += shr_sumff[threadIdx.x+2] 
            + shr_sumff[threadIdx.x+3];
    }
    __syncthreads();

    if (sample==0) {
        // exit if the state is done
        // note, we do not exit before all __synchtreads are finished
        if (state == TimeComputationState::Finished) {
            g_timeMax[ch] = 5;
            g_timeError[ch] = -999;
            return;
        }

        // subtract to avoid double counting
        auto const sumff = shr_sumff[threadIdx.x] 
            + shr_sumff[threadIdx.x+1] 
            - shr_sumff[threadIdx.x+3];
        auto const sumAf = shr_sumAf[threadIdx.x]
            + shr_sumAf[threadIdx.x+1]
            - shr_sumAf[threadIdx.x+3];

        auto const ampMaxAlphaBeta = sumff>0 ? sumAf / sumff : 0;
        auto const sumAA = sumAAsNullHypot[ch];
        auto const sum0 = sum0sNullHypot[ch];
        auto const nullChi2 = chi2sNullHypot[ch];
        if (sumff > 0) {
            auto const chi2AlphaBeta = (sumAA - sumAf * sumAf / sumff) / sum0;
            if (chi2AlphaBeta > nullChi2) {
                // null hypothesis is better
                state = TimeComputationState::Finished;
#ifdef DEBUG_FINDAMPLCHI2_AND_FINISH
                printf("ch = %d chi2AlphaBeta = %f nullChi2 = %f sumAA = %f sumAf = %f sumff = %f sum0 = %f\n",
                    ch, chi2AlphaBeta, nullChi2, sumAA, sumAf, sumff, sum0);
#endif
            }

            // store to global
            g_ampMaxAlphaBeta[ch] = ampMaxAlphaBeta;
        } else {
#ifdef DEBUG_FINDAMPLCHI2_AND_FINISH
            printf("ch = %d sum0 = %f sumAA = %f sumff = %f sumAf = %f\n",
                ch, sum0, sumAA, sumff, sumAf);
#endif
            state = TimeComputationState::Finished;
        }

        // store the state to global and finish calcs
        g_state[ch] = state;
        if (state == TimeComputationState::Finished) {
            // store default values into global
            g_timeMax[ch] = 5;
            g_timeError[ch] = -999;
#ifdef DEBUG_FINDAMPLCHI2_AND_FINISH
            printf("ch = %d finished state\n", ch);
#endif
            return;
        }

        auto const ampMaxError = g_ampMaxError[ch];
        auto const test_ratio = ampMaxAlphaBeta / ampMaxError;
        auto const accTimeMax = g_accTimeMax[ch];
        auto const accTimeWgt = g_accTimeWgt[ch];
        auto const tMaxAlphaBeta = g_tMaxAlphaBeta[ch];
        auto const tMaxErrorAlphaBeta = g_tMaxErrorAlphaBeta[ch];
        // branch to separate large vs small pulses
        // see cpu version for more info
        if (test_ratio > 5.0 && accTimeWgt>0) {
            auto const tMaxRatio = accTimeWgt>0 
                ? accTimeMax / accTimeWgt 
                : static_cast<ScalarType>(0);
            auto const tMaxErrorRatio = accTimeWgt>0 
                ? 1.0 / std::sqrt(accTimeWgt) 
                : static_cast<ScalarType>(0);

            if (test_ratio > 10.0) {
                g_timeMax[ch] = tMaxRatio;
                g_timeError[ch] = tMaxErrorRatio;
                
#ifdef DEBUG_FINDAMPLCHI2_AND_FINISH
                    printf("ch = %d tMaxRatio = %f tMaxErrorRatio = %f\n",
                        ch, tMaxRatio, tMaxErrorRatio);
#endif
            } else {
                auto const timeMax = 
                    (tMaxAlphaBeta * (10.0 - ampMaxAlphaBeta / ampMaxError) + 
                     tMaxRatio * (ampMaxAlphaBeta / ampMaxError - 5.0)) / 5.0;
                auto const timeError = 
                    (tMaxErrorAlphaBeta * (10.0 - ampMaxAlphaBeta / ampMaxError) + 
                     tMaxErrorRatio * (ampMaxAlphaBeta / ampMaxError - 5.0)) / 5.0;
                state = TimeComputationState::Finished;
                g_state[ch] = state;
                g_timeMax[ch] = timeMax;
                g_timeError[ch] = timeError;

#ifdef DEBUG_FINDAMPLCHI2_AND_FINISH
                    printf("ch = %d timeMax = %f timeError = %f\n",
                        ch, timeMax, timeError);
#endif
            }
        }
        else {
            state = TimeComputationState::Finished;
            g_state[ch] = state;
            g_timeMax[ch] = tMaxAlphaBeta;
            g_timeError[ch] = tMaxErrorAlphaBeta;

#ifdef DEBUG_FINDAMPLCHI2_AND_FINISH
                printf("ch = %d tMaxAlphaBeta = %f tMaxErrorAlphaBeta = %f\n",
                    ch, tMaxAlphaBeta, tMaxErrorAlphaBeta);
#endif
        }
    }
}

__global__
void kernel_time_compute_fixMGPAslew(uint16_t const* digis_eb,
                                     uint16_t const* digis_ee,
                                     SampleVector::Scalar* sample_values,
                                     SampleVector::Scalar* sample_value_errors,
                                     bool* useless_sample_values,
                                     unsigned int const sample_mask,
                                     int const nchannels,
                                     uint32_t const offsetForInputs) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockIdx.x * blockDim.x;
    int const ch = gtx / nsamples;
    int const sample = threadIdx.x % nsamples;
    int const inputCh = ch >= offsetForInputs
        ? ch - offsetForInputs
        : ch;
    int const inputGtx = ch >= offsetForInputs
        ? gtx - offsetForInputs*nsamples
        : gtx;
    auto const* digis = ch >= offsetForInputs
        ? digis_ee
        : digis_eb;

    // remove thread for sample 0, oversubscribing is easier than ....
    if (ch >= nchannels || sample==0) return;

    if (!use_sample(sample_mask, sample)) return;

    auto const gainIdPrev = ecal::mgpa::gainId(digis[inputGtx-1]);
    auto const gainIdNext = ecal::mgpa::gainId(digis[inputGtx]);
    if (gainIdPrev>=1 && gainIdPrev<=3 &&
        gainIdNext>=1 && gainIdNext<=3 && gainIdPrev < gainIdNext) {
        sample_values[gtx-1] = 0;
        sample_value_errors[gtx-1] = 1e+9;
        useless_sample_values[gtx-1] = true;
    }
}

__global__
void kernel_time_compute_ampl(SampleVector::Scalar const* sample_values,
                              SampleVector::Scalar const* sample_value_errors,
                              uint32_t const* dids,
                              bool const* useless_samples,
                              SampleVector::Scalar const* g_timeMax,
                              SampleVector::Scalar const* amplitudeFitParametersEB,
                              SampleVector::Scalar const* amplitudeFitParametersEE,
                              SampleVector::Scalar *g_amplitudeMax,
                              int const nchannels) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr ScalarType corr4 = 1.;
    constexpr ScalarType corr6 = 1.;
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockIdx.x * blockDim.x;
    int const ch = gtx / nsamples;
    int const sample = threadIdx.x % nsamples;

    if (ch >= nchannels) return;

    auto const did = DetId{dids[ch]};
    auto const* amplitudeFitParameters = did.subdetId() == EcalBarrel
        ? amplitudeFitParametersEB
        : amplitudeFitParametersEE;

    // configure shared mem
    extern __shared__ char smem[];
    ScalarType* shr_sum1 = reinterpret_cast<ScalarType*>(smem);
    auto *shr_sumA = shr_sum1 + blockDim.x;
    auto *shr_sumF = shr_sumA + blockDim.x;
    auto *shr_sumAF = shr_sumF + blockDim.x;
    auto *shr_sumFF = shr_sumAF + blockDim.x;

    auto const alpha = amplitudeFitParameters[0];
    auto const beta = amplitudeFitParameters[1];
    auto const timeMax = g_timeMax[ch];
    auto const pedestalLimit = timeMax - (alpha * beta) - 1.0;
    auto const sample_value = sample_values[gtx];
    auto const sample_value_error = sample_value_errors[gtx];
    auto const inverr2 = sample_value_error > 0
        ? 1. / (sample_value_error * sample_value_error)
        : static_cast<ScalarType>(0);
    auto const termOne = 1 + (sample - timeMax) / (alpha * beta);
    auto const f = termOne > 1.e-5
        ? fast_expf(alpha * fast_logf(termOne) - 
            (sample - timeMax) / beta)
        : static_cast<ScalarType>(0.); 

    bool const cond = ((sample < pedestalLimit) ||
        (f>0.6*corr6 && sample<=timeMax) ||
        (f>0.4*corr4 && sample>=timeMax)) && !useless_samples[gtx];

    // store into shared mem
    shr_sum1[threadIdx.x] = cond ? inverr2 : static_cast<ScalarType>(0);
    shr_sumA[threadIdx.x] = cond
        ? sample_value * inverr2
        : static_cast<ScalarType>(0);
    shr_sumF[threadIdx.x] = cond 
        ? f * inverr2
        : static_cast<ScalarType>(0);
    shr_sumAF[threadIdx.x] = cond 
        ? (f*inverr2)*sample_value
        : static_cast<ScalarType>(0);
    shr_sumFF[threadIdx.x] = cond 
        ? f*(f*inverr2)
        : static_cast<ScalarType>(0);

    // reduction
    if (sample <= 4) {
        shr_sum1[threadIdx.x] += shr_sum1[threadIdx.x+5];
        shr_sumA[threadIdx.x] += shr_sumA[threadIdx.x+5];
        shr_sumF[threadIdx.x] += shr_sumF[threadIdx.x+5];
        shr_sumAF[threadIdx.x] += shr_sumAF[threadIdx.x+5];
        shr_sumFF[threadIdx.x] += shr_sumFF[threadIdx.x+5];
    }
    __syncthreads();

    if (sample < 2) {
        // note: we double count sample 3
        shr_sum1[threadIdx.x] += shr_sum1[threadIdx.x+2] + shr_sum1[threadIdx.x+3];
        shr_sumA[threadIdx.x] += shr_sumA[threadIdx.x+2] + shr_sumA[threadIdx.x+3];
        shr_sumF[threadIdx.x] += shr_sumF[threadIdx.x+2] + shr_sumF[threadIdx.x+3];
        shr_sumAF[threadIdx.x] += shr_sumAF[threadIdx.x+2] 
            + shr_sumAF[threadIdx.x+3];
        shr_sumFF[threadIdx.x] += shr_sumFF[threadIdx.x+2] 
            + shr_sumFF[threadIdx.x+3];
    }
    __syncthreads();

    if (sample == 0) {
        auto const sum1 = shr_sum1[threadIdx.x] 
            + shr_sum1[threadIdx.x+1] - shr_sum1[threadIdx.x+3];
        auto const sumA = shr_sumA[threadIdx.x] 
            + shr_sumA[threadIdx.x+1] - shr_sumA[threadIdx.x+3];
        auto const sumF = shr_sumF[threadIdx.x] 
            + shr_sumF[threadIdx.x+1] - shr_sumF[threadIdx.x+3];
        auto const sumAF = shr_sumAF[threadIdx.x] 
            + shr_sumAF[threadIdx.x+1] - shr_sumAF[threadIdx.x+3];
        auto const sumFF = shr_sumFF[threadIdx.x] 
            + shr_sumFF[threadIdx.x+1] - shr_sumFF[threadIdx.x+3];

        auto const denom = sumFF * sum1 - sumF*sumF;
        auto const condForDenom = sum1 > 0 && ecal::abs(denom)>1.e-20;
        auto const amplitudeMax = condForDenom
            ? (sumAF * sum1 - sumA * sumF) / denom
            : static_cast<ScalarType>(0.);

        // store into global mem
        g_amplitudeMax[ch] = amplitudeMax;
    }
}

//#define ECAL_RECO_CUDA_TC_INIT_DEBUG
__global__
void kernel_time_computation_init(uint16_t const* digis_eb,
                                  uint32_t const* dids_eb,
                                  uint16_t const* digis_ee,
                                  uint32_t const* dids_ee,
                                  float const* rms_x12,
                                  float const* rms_x6,
                                  float const* rms_x1,
                                  float const* mean_x12,
                                  float const* mean_x6,
                                  float const* mean_x1,
                                  float const* gain12Over6,
                                  float const* gain6Over1,
                                  SampleVector::Scalar* sample_values,
                                  SampleVector::Scalar* sample_value_errors,
                                  SampleVector::Scalar* ampMaxError,
                                  bool* useless_sample_values,
                                  char* pedestal_nums,
                                  uint32_t const offsetForHashes,
                                  uint32_t const offsetForInputs,
                                  unsigned int const sample_maskEB,
                                  unsigned int const sample_maskEE,
                                  int nchannels) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const tx = threadIdx.x + blockDim.x*blockIdx.x;
    int const ch = tx/nsamples;
    int const inputTx = ch >= offsetForInputs
        ? tx - offsetForInputs*nsamples
        : tx;
    int const inputCh = ch >= offsetForInputs
        ? ch - offsetForInputs
        : ch;
    auto const* digis = ch >= offsetForInputs
        ? digis_ee
        : digis_eb;
    auto const* dids = ch >= offsetForInputs
        ? dids_ee
        : dids_eb;

    if (ch < nchannels) {
        // indices/inits
        int const sample = tx % nsamples;
        int const ch_start = ch*nsamples;
        int const input_ch_start = inputCh*nsamples;
        SampleVector::Scalar pedestal = 0.;
        int num = 0;

        // configure shared mem
        extern __shared__ char smem[];
        ScalarType* shrSampleValues = 
            reinterpret_cast<SampleVector::Scalar*>(smem);
        ScalarType* shrSampleValueErrors = shrSampleValues + blockDim.x;

        // 0 and 1 sample values
        auto const adc0 = ecal::mgpa::adc(digis[input_ch_start]);
        auto const gainId0 = ecal::mgpa::gainId(digis[input_ch_start]);
        auto const adc1 = ecal::mgpa::adc(digis[input_ch_start+1]);
        auto const gainId1 = ecal::mgpa::gainId(digis[input_ch_start+1]);
        auto const did = DetId{dids[inputCh]};
        auto const isBarrel = did.subdetId() == EcalBarrel;
        auto const sample_mask = did.subdetId() == EcalBarrel
            ? sample_maskEB
            : sample_maskEE;
        auto const hashedId = isBarrel
            ? hashedIndexEB(did.rawId())
            : offsetForHashes + hashedIndexEE(did.rawId());

        // set pedestal
        // TODO this branch is non-divergent for a group of 10 threads
        if (gainId0 == 1 && use_sample(sample_mask, 0)) {
            pedestal = static_cast<SampleVector::Scalar>(adc0);
            num=1;

            auto const diff = adc1 - adc0;
            if (gainId1 == 1 && use_sample(sample_mask, 1)
                && std::abs(diff) < 3*rms_x12[hashedId]) {
                pedestal = 
                    (pedestal + static_cast<SampleVector::Scalar>(adc1)) / 2.0;
                num=2;
            }
        } else {
            pedestal = mean_x12[ch];
        }

        // ped subtracted and gain-renormalized samples.
        auto const gainId = ecal::mgpa::gainId(digis[inputTx]);
        auto const adc = ecal::mgpa::adc(digis[inputTx]);

        bool bad = false;
        SampleVector::Scalar sample_value, sample_value_error;
        // TODO divergent branch
        // TODO: piece below is general both for amplitudes and timing
        // potentially there is a way to reduce the amount of code...
        if (!use_sample(sample_mask, sample)) {
            bad = true;
            sample_value = 0;
            sample_value_error = 0;
        } else if (gainId == 1) {
            sample_value = static_cast<SampleVector::Scalar>(adc) - pedestal;
            sample_value_error = rms_x12[hashedId];
        } else if (gainId == 2) {
            sample_value =  (static_cast<SampleVector::Scalar>(adc) 
                - mean_x6[hashedId]) * gain12Over6[hashedId]; 
            sample_value_error = rms_x6[hashedId] * gain12Over6[hashedId];
        } else if (gainId == 3) {
            sample_value = (static_cast<SampleVector::Scalar>(adc) 
                - mean_x1[hashedId]) * gain6Over1[hashedId] * gain12Over6[hashedId];
            sample_value_error = rms_x1[hashedId] 
                * gain6Over1[hashedId] * gain12Over6[hashedId];
        } else {
            sample_value = 0;
            sample_value_error = 0;
            bad = true;
        }

        // TODO: make sure we save things correctly when sample is useless
        auto const useless_sample = (sample_value_error <= 0) | bad;
        useless_sample_values[tx] = useless_sample;
        sample_values[tx] = sample_value;
        sample_value_errors[tx] = useless_sample ? 1e+9 : sample_value_error;

        // DEBUG
#ifdef ECAL_RECO_CUDA_TC_INIT_DEBUG
        if (ch == 0) {
            printf("sample = %d sample_value = %f sample_value_error = %f useless = %c\n",
                sample, sample_value, sample_value_error, 
                useless_sample ? '1' : '0');           
        }
#endif

        // store into the shared mem
        shrSampleValues[threadIdx.x] = sample_value_error > 0
            ? sample_value
            : std::numeric_limits<ScalarType>::min();
        shrSampleValueErrors[threadIdx.x] = sample_value_error;
        __syncthreads();

        // perform the reduction with min
        if (sample < 5) {
            // note, if equal -> we keep the value with lower sample as for cpu
            shrSampleValueErrors[threadIdx.x] = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+5] 
                ? shrSampleValueErrors[threadIdx.x+5]
                : shrSampleValueErrors[threadIdx.x];
            shrSampleValues[threadIdx.x] = 
                std::max(shrSampleValues[threadIdx.x], 
                         shrSampleValues[threadIdx.x+5]);
        }
        __syncthreads();

        // a bit of an overkill, but easier than to compare across 3 values
        if (sample<3) {
            shrSampleValueErrors[threadIdx.x] = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+3]
                ? shrSampleValueErrors[threadIdx.x+3]
                : shrSampleValueErrors[threadIdx.x];
            shrSampleValues[threadIdx.x] = 
                std::max(shrSampleValues[threadIdx.x], 
                         shrSampleValues[threadIdx.x+3]);
        }
        __syncthreads();

        if (sample < 2) {
            shrSampleValueErrors[threadIdx.x] = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+2]
                ? shrSampleValueErrors[threadIdx.x+2]
                : shrSampleValueErrors[threadIdx.x];
            shrSampleValues[threadIdx.x] = 
                std::max(shrSampleValues[threadIdx.x], 
                         shrSampleValues[threadIdx.x+2]);
        }
        __syncthreads();
 
        if (sample == 0) {
            // we only needd the max error
            auto const maxSampleValueError = 
                shrSampleValues[threadIdx.x] < shrSampleValues[threadIdx.x+1]
                ? shrSampleValueErrors[threadIdx.x+1]
                : shrSampleValueErrors[threadIdx.x];

            // # pedestal samples used
            pedestal_nums[ch] = num;
            // this is used downstream
            ampMaxError[ch] = maxSampleValueError;

            // DEBUG
#ifdef ECAL_RECO_CUDA_TC_INIT_DEBUG
            if (ch == 0) {
                printf("pedestal_nums = %d ampMaxError = %f\n",
                    num, maxSampleValueError);
            }
#endif
        }
    }
}

///
/// launch context parameters: 1 thread per channel
///
//#define DEBUG_TIME_CORRECTION
__global__
void kernel_time_correction_and_finalize(
//        SampleVector::Scalar const* g_amplitude,
        ::ecal::reco::StorageScalarType const* g_amplitude,
        uint16_t const* digis_eb,
        uint32_t const* dids_eb,
        uint16_t const* digis_ee,
        uint32_t const* dids_ee,
        float const* amplitudeBinsEB,
        float const* amplitudeBinsEE,
        float const* shiftBinsEB,
        float const* shiftBinsEE,
        SampleVector::Scalar const* g_timeMax,
        SampleVector::Scalar const* g_timeError,
        float const* g_rms_x12,
        float const* timeCalibConstant,
        float *g_jitter,
        float *g_jitterError,
        uint32_t *flags,
        int const amplitudeBinsSizeEB,
        int const amplitudeBinsSizeEE,
        ConfigurationParameters::type const timeConstantTermEB,
        ConfigurationParameters::type const timeConstantTermEE,
        float const offsetTimeValueEB,
        float const offsetTimeValueEE,
        ConfigurationParameters::type const timeNconstEB,
        ConfigurationParameters::type const timeNconstEE,
        ConfigurationParameters::type const amplitudeThresholdEB,
        ConfigurationParameters::type const amplitudeThresholdEE,
        ConfigurationParameters::type const outOfTimeThreshG12pEB,
        ConfigurationParameters::type const outOfTimeThreshG12pEE,
        ConfigurationParameters::type const outOfTimeThreshG12mEB,
        ConfigurationParameters::type const outOfTimeThreshG12mEE,
        ConfigurationParameters::type const outOfTimeThreshG61pEB,
        ConfigurationParameters::type const outOfTimeThreshG61pEE,
        ConfigurationParameters::type const outOfTimeThreshG61mEB,
        ConfigurationParameters::type const outOfTimeThreshG61mEE,
        uint32_t const offsetForHashes,
        uint32_t const offsetForInputs,
        int const nchannels) {
    using ScalarType = SampleVector::Scalar;

    // constants
    constexpr int nsamples = EcalDataFrame::MAXSAMPLES;

    // indices
    int const gtx = threadIdx.x + blockIdx.x * blockDim.x;
    int const inputGtx = gtx >= offsetForInputs
        ? gtx - offsetForInputs
        : gtx;
    auto const* dids = gtx >= offsetForInputs
        ? dids_ee
        : dids_eb;
    auto const& digis = gtx >= offsetForInputs
        ? digis_ee
        : digis_eb;

    // filter out outside of range threads
    if (gtx >= nchannels) return;

    auto const did = DetId{dids[inputGtx]};
    auto const isBarrel = did.subdetId() == EcalBarrel;
    auto const hashedId = isBarrel
        ? hashedIndexEB(did.rawId())
        : offsetForHashes + hashedIndexEE(did.rawId());
    auto const* amplitudeBins = isBarrel
        ? amplitudeBinsEB
        : amplitudeBinsEE;
    auto const* shiftBins = isBarrel
        ? shiftBinsEB
        : shiftBinsEE;
    auto const amplitudeBinsSize = isBarrel
        ? amplitudeBinsSizeEB
        : amplitudeBinsSizeEE;
    auto const timeConstantTerm = isBarrel 
        ? timeConstantTermEB
        : timeConstantTermEE;
    auto const timeNconst = isBarrel 
        ? timeNconstEB
        : timeNconstEE;
    auto const offsetTimeValue = isBarrel
        ? offsetTimeValueEB
        : offsetTimeValueEE;
    auto const amplitudeThreshold = isBarrel
        ? amplitudeThresholdEB
        : amplitudeThresholdEE;
    auto const outOfTimeThreshG12p = isBarrel
        ? outOfTimeThreshG12pEB
        : outOfTimeThreshG12pEE;
    auto const outOfTimeThreshG12m = isBarrel
        ? outOfTimeThreshG12mEB
        : outOfTimeThreshG12mEE;
    auto const outOfTimeThreshG61p = isBarrel
        ? outOfTimeThreshG61pEB
        : outOfTimeThreshG61pEE;
    auto const outOfTimeThreshG61m = isBarrel
        ? outOfTimeThreshG61mEB
        : outOfTimeThreshG61mEE;
    
    // load some
    auto const amplitude = g_amplitude[gtx];
    auto const rms_x12 = g_rms_x12[hashedId];
    auto const timeCalibConst = timeCalibConstant[hashedId];

    int myBin = -1;
    for (int bin=0; bin<amplitudeBinsSize; bin++) {
        if (amplitude > amplitudeBins[bin]) 
            myBin = bin;
        else 
            break;
    }

    ScalarType correction = 0;
    if (myBin == -1) {
        correction = shiftBins[0];
    } else if (myBin == amplitudeBinsSize-1) {
        correction = shiftBins[myBin];
    } else {
        correction = shiftBins[myBin+1] - shiftBins[myBin];
        correction *= (amplitude - amplitudeBins[myBin]) / 
            (amplitudeBins[myBin+1] - amplitudeBins[myBin]);
        correction += shiftBins[myBin];
    }

    // correction * 1./25.
    correction = correction * 0.04;
    auto const timeMax = g_timeMax[gtx];
    auto const timeError = g_timeError[gtx];
    auto const jitter = timeMax - 5 + correction;
    auto const jitterError = std::sqrt(timeError*timeError + 
        timeConstantTerm*timeConstantTerm * 0.04 * 0.04); // 0.04 = 1./25.

#ifdef DEBUG_TIME_CORRECTION
//    if (gtx == 0) {
        printf("ch = %d timeMax = %f timeError = %f jitter = %f correction = %f\n",
            gtx, timeMax, timeError, jitter, correction);
//    }
#endif

    // store back to  global
    g_jitter[gtx] = jitter;
    g_jitterError[gtx] = jitterError;

    // set the flag
    // TODO: replace with something more efficient (if required), 
    // for now just to make it work
    if (amplitude > amplitudeThreshold * rms_x12) {
        auto threshP = outOfTimeThreshG12p;
        auto threshM = outOfTimeThreshG12m;
        if (amplitude > 3000.) {
            for (int isample=0; isample<nsamples; isample++) {
                int gainid = ecal::mgpa::gainId(digis[nsamples*inputGtx + isample]);
                if (gainid != 1) {
                    threshP = outOfTimeThreshG61p;
                    threshM = outOfTimeThreshG61m;
                    break;
                }
            }
        }

        auto const correctedTime = (timeMax - 5) * 25 + 
            timeCalibConst + offsetTimeValue;
        auto const nterm = timeNconst * rms_x12 / amplitude;
        auto const sigmat = std::sqrt(nterm * nterm + 
            timeConstantTerm*timeConstantTerm);
        if (correctedTime > sigmat*threshP || 
            correctedTime < -sigmat*threshM)
            flags[gtx] |= 0x1 << EcalUncalibratedRecHit::kOutOfTime;
    }
}

}}
