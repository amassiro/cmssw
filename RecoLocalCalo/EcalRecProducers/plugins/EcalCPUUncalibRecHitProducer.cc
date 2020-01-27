#include <iostream>

// framework
#include "FWCore/Framework/interface/stream/EDProducer.h"
//#include "HeterogeneousCore/Producer/interface/HeterogeneousEDProducer.h"
//#include "HeterogeneousCore/Producer/interface/HeterogeneousEvent.h"

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAScopedContext.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h" 

// algorithm specific

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"

class EcalCPUUncalibRecHitProducer
    : public edm::stream::EDProducer<edm::ExternalWork>
{
public:
    explicit EcalCPUUncalibRecHitProducer(edm::ParameterSet const& ps);
    ~EcalCPUUncalibRecHitProducer() override;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
    void acquire(edm::Event const&, 
                 edm::EventSetup const&,
                 edm::WaitingTaskWithArenaHolder) override;
    void produce(edm::Event&, edm::EventSetup const&) override;

private:
    edm::EDGetTokenT<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr>>> 
        recHitsInEBToken_, recHitsInEEToken_;
    edm::EDPutTokenT<ecal::UncalibratedRecHit<ecal::Tag::soa>>
        recHitsOutEBToken_, recHitsOutEEToken_;

    ecal::UncalibratedRecHit<ecal::Tag::soa>
        recHitsEB_, recHitsEE_;
    bool containsTimingInformation_;
};

void EcalCPUUncalibRecHitProducer::fillDescriptions(
        edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("recHitsInLabelEB", 
        edm::InputTag{"ecalUncalibRecHitProducerGPU", "EcalUncalibRecHitsEB"});
    desc.add<edm::InputTag>("recHitsInLabelEE", 
        edm::InputTag{"ecalUncalibRecHitProducerGPU", "EcalUncalibRecHitsEE"});
    desc.add<std::string>("recHitsOutLabelEB", "EcalUncalibRecHitsEB");
    desc.add<std::string>("recHitsOutLabelEE", "EcalUncalibRecHitsEE");
    desc.add<bool>("containsTimingInformation", false);

    std::string label = "ecalCPUUncalibRecHitProducer";
    confDesc.add(label, desc);
}

EcalCPUUncalibRecHitProducer::EcalCPUUncalibRecHitProducer(
        const edm::ParameterSet& ps) 
    : recHitsInEBToken_{consumes<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr>>>(
        ps.getParameter<edm::InputTag>("recHitsInLabelEB"))}
    , recHitsInEEToken_{consumes<CUDAProduct<ecal::UncalibratedRecHit<ecal::Tag::ptr>>>(
        ps.getParameter<edm::InputTag>("recHitsInLabelEE"))}
    , recHitsOutEBToken_{produces<ecal::UncalibratedRecHit<ecal::Tag::soa>>(
        ps.getParameter<std::string>("recHitsOutLabelEB"))}
    , recHitsOutEEToken_{produces<ecal::UncalibratedRecHit<ecal::Tag::soa>>(
        ps.getParameter<std::string>("recHitsOutLabelEE"))}
    , containsTimingInformation_{ps.getParameter<bool>("containsTimingInformation")}
{}

EcalCPUUncalibRecHitProducer::~EcalCPUUncalibRecHitProducer() {}

void EcalCPUUncalibRecHitProducer::acquire(
        edm::Event const& event,
        edm::EventSetup const& setup,
        edm::WaitingTaskWithArenaHolder taskHolder) 
{
    // retrieve data/ctx
    auto const& ebRecHitsProduct = event.get(recHitsInEBToken_);
    auto const& eeRecHitsProduct = event.get(recHitsInEEToken_);
    CUDAScopedContextAcquire ctx{ebRecHitsProduct, std::move(taskHolder)};
    auto const& ebRecHits = ctx.get(ebRecHitsProduct);
    auto const& eeRecHits = ctx.get(eeRecHitsProduct);

    // resize the output buffers
    recHitsEB_.resize(ebRecHits.size);
    recHitsEE_.resize(eeRecHits.size);

    // AM: FIXME : why all "uint32_t" and not "float" where needed?
    
    // enqeue transfers
    cudaCheck( cudaMemcpyAsync(recHitsEB_.did.data(),
                               ebRecHits.did,
                               recHitsEB_.did.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.did.data(),
                               eeRecHits.did,
                               recHitsEE_.did.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );

//
//     ./CUDADataFormats/EcalRecHitSoA/interface/RecoTypes.h:using ComputationScalarType = float;    
//     typename type_wrapper<reco::ComputationScalarType, L>::type amplitudesAll;
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.amplitudesAll.data(),
                               ebRecHits.amplitudesAll,
                               recHitsEB_.amplitudesAll.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.amplitudesAll.data(),
                               eeRecHits.amplitudesAll,
                               recHitsEE_.amplitudesAll.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.amplitude.data(),
                               ebRecHits.amplitude,
                               recHitsEB_.amplitude.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.amplitude.data(),
                               eeRecHits.amplitude,
                               recHitsEE_.amplitude.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.chi2.data(),
                               ebRecHits.chi2,
                               recHitsEB_.chi2.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.chi2.data(),
                               eeRecHits.chi2,
                               recHitsEE_.chi2.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.pedestal.data(),
                               ebRecHits.pedestal,
                               recHitsEB_.pedestal.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.pedestal.data(),
                               eeRecHits.pedestal,
                               recHitsEE_.pedestal.size() * sizeof(float),   // AM: FIX
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    
    cudaCheck( cudaMemcpyAsync(recHitsEB_.flags.data(),
                               ebRecHits.flags,
                               recHitsEB_.flags.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    cudaCheck( cudaMemcpyAsync(recHitsEE_.flags.data(),
                               eeRecHits.flags,
                               recHitsEE_.flags.size() * sizeof(uint32_t),
                               cudaMemcpyDeviceToHost,
                               ctx.stream().id()) );
    
    if (containsTimingInformation_) {
        cudaCheck( cudaMemcpyAsync(recHitsEB_.jitter.data(),
                                   ebRecHits.jitter,
                                   recHitsEB_.jitter.size() * sizeof(float),   // AM: FIX
                                   cudaMemcpyDeviceToHost,
                                   ctx.stream().id()) );
        cudaCheck( cudaMemcpyAsync(recHitsEE_.jitter.data(),
                                   eeRecHits.jitter,
                                   recHitsEE_.jitter.size() * sizeof(float),   // AM: FIX
                                   cudaMemcpyDeviceToHost,
                                   ctx.stream().id()) );
        
        cudaCheck( cudaMemcpyAsync(recHitsEB_.jitterError.data(),
                                   ebRecHits.jitterError,
                                   recHitsEB_.jitterError.size() * sizeof(float),   // AM: FIX
                                   cudaMemcpyDeviceToHost,
                                   ctx.stream().id()) );
        cudaCheck( cudaMemcpyAsync(recHitsEE_.jitterError.data(),
                                   eeRecHits.jitterError,
                                   recHitsEE_.jitterError.size() * sizeof(float),   // AM: FIX
                                   cudaMemcpyDeviceToHost,
                                   ctx.stream().id()) );
    }
}

void EcalCPUUncalibRecHitProducer::produce(
        edm::Event& event, 
        edm::EventSetup const& setup) 
{
    // tmp vectors
    auto recHitsOutEB = std::make_unique<ecal::UncalibratedRecHit<ecal::Tag::soa>>(
        std::move(recHitsEB_));
    auto recHitsOutEE = std::make_unique<ecal::UncalibratedRecHit<ecal::Tag::soa>>(
        std::move(recHitsEE_));

    // put into event
    event.put(recHitsOutEBToken_, std::move(recHitsOutEB));
    event.put(recHitsOutEEToken_, std::move(recHitsOutEE));
}

DEFINE_FWK_MODULE(EcalCPUUncalibRecHitProducer);
