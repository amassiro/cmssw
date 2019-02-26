#ifndef PulseChiSqSNNLS_gpu_h
#define PulseChiSqSNNLS_gpu_h

#define EIGEN_NO_DEBUG // kill throws in eigen code
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"

#include <set>
#include <array>

namespace ecal { namespace multifit {

class PulseChiSqSNNLS {
  public:
    
    typedef BXVector::Index Index;
    
    __device__
    PulseChiSqSNNLS();
    
    
    __device__
    bool DoFit(const SampleVector &samples, const SampleMatrix &samplecov, const BXVector &bxs, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const SampleGainVector &gains = -1*SampleGainVector::Ones(), const SampleGainVector &badSamples = SampleGainVector::Zero());
    
    __device__
    const SamplePulseMatrix &pulsemat() const { return _pulsemat; }
    __device__ const SampleMatrix &invcov() const { return _invcov; }
    
    __device__
    const PulseVector &X() const { return _ampvecmin; }
    __device__
    const PulseVector &Errors() const { return _errvec; }
    __device__
    const BXVector &BXs() const { return _bxsmin; }
    
    __device__
    double ChiSq() const { return _chisq; }
    __device__
    void disableErrorCalculation() { _computeErrors = false; }
    __device__
    void setMaxIters(int n) { _maxiters = n;}
    __device__
    void setMaxIterWarnings(bool b) { _maxiterwarnings = b;}

  protected:
    
    __device__
    bool Minimize(const SampleMatrix &samplecov, const FullSampleMatrix &fullpulsecov);
    __device__
    bool NNLS();
    __device__
    void NNLSUnconstrainParameter(Index idxp);
    __device__
    void NNLSConstrainParameter(Index minratioidx);
    __device__
    bool OnePulseMinimize();
    __device__
    bool updateCov(const SampleMatrix &samplecov, const FullSampleMatrix &fullpulsecov);
    __device__
    double ComputeChiSq();
    __device__
    double ComputeApproxUncertainty(unsigned int ipulse);
    
    
    SampleVector _sampvec;
    SampleMatrix _invcov;
    SamplePulseMatrix _pulsemat;
    PulseVector _ampvec;
    PulseVector _errvec;
    PulseVector _ampvecmin;
    
    SampleDecompLLT _covdecomp;
    SampleMatrix _covdecompLinv;
    PulseMatrix _topleft_work;
    PulseDecompLDLT _pulsedecomp;

    BXVector _bxs;
    BXVector _bxsmin;
    unsigned int _npulsetot;
    unsigned int _nP;
    
    SamplePulseMatrix invcovp;
    PulseMatrix aTamat;
    PulseVector aTbvec;
    PulseVector updatework;
    
    PulseVector ampvecpermtest;
    
    double _chisq;
    bool _computeErrors;
    int _maxiters;
    bool _maxiterwarnings;
};

}}

#endif
