#include "cuda.h"

#include "KernelHelpers.h"

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"

//
//
#include "EcalRecHitBuilderKernels.h"


#include "KernelHelpers.h"




namespace ecal {
  namespace rechit {


    // uncalibrecHit flags
    enum UncalibRecHitFlags {
      kGood=-1,                 // channel is good (mutually exclusive with other states)  setFlagBit(kGood) reset flags_ to zero 
      kPoorReco,                // channel has been badly reconstructed (e.g. bad shape, bad chi2 etc.)
      kSaturated,               // saturated channel
      kOutOfTime,               // channel out of time
      kLeadingEdgeRecovered,    // saturated channel: energy estimated from the leading edge before saturation
      kHasSwitchToGain6,        // at least one data frame is in G6
      kHasSwitchToGain1         // at least one data frame is in G1
    };

    
    // recHit flags
    enum RecHitFlags { 
      RecHitFlags_kGood=0,                   // channel ok, the energy and time measurement are reliable
      RecHitFlags_kPoorReco,                 // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
      RecHitFlags_kOutOfTime,                // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
      RecHitFlags_kFaultyHardware,           // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
      RecHitFlags_kNoisy,                    // the channel is very noisy
      RecHitFlags_kPoorCalib,                // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
      RecHitFlags_kSaturated,                // saturated channel (recovery not tried)
      RecHitFlags_kLeadingEdgeRecovered,     // saturated channel: energy estimated from the leading edge before saturation
      RecHitFlags_kNeighboursRecovered,      // saturated/isolated dead: energy estimated from neighbours
      RecHitFlags_kTowerRecovered,           // channel in TT with no data link, info retrieved from Trigger Primitive
      RecHitFlags_kDead,                     // channel is dead and any recovery fails
      RecHitFlags_kKilled,                   // MC only flag: the channel is killed in the real detector
      RecHitFlags_kTPSaturated,              // the channel is in a region with saturated TP
      RecHitFlags_kL1SpikeFlag,              // the channel is in a region with TP with sFGVB = 0
      RecHitFlags_kWeird,                    // the signal is believed to originate from an anomalous deposit (spike) 
      RecHitFlags_kDiWeird,                  // the signal is anomalous, and neighbors another anomalous signal  
      RecHitFlags_kHasSwitchToGain6,         // at least one data frame is in G6
      RecHitFlags_kHasSwitchToGain1,         // at least one data frame is in G1
      //
      RecHitFlags_kUnknown                   // to ease the interface with functions returning flags. 
    };
    
   
    // status code
    enum EcalChannelStatusCode_Code {
      kOk=0,
      kDAC,
      kNoLaser,
      kNoisy,
      kNNoisy,
      kNNNoisy,
      kNNNNoisy,
      kNNNNNoisy,
      kFixedG6,
      kFixedG1,
      kFixedG0,
      kNonRespondingIsolated,
      kDeadVFE,
      kDeadFE,
      kNoDataNoTP      
    };
    
    
    
    
    
    __global__
    void kernel_create_ecal_rehit(
                     // configuration 
                     int const* ChannelStatusToBeExcluded,
                     uint32_t ChannelStatusToBeExcludedSize,   
                     bool const killDeadChannels,
                     float const EBLaserMIN,
                     float const EELaserMIN,
                     float const EBLaserMAX,
                     float const EELaserMAX,
                     // for flags setting
                     int const* expanded_v_DB_reco_flags,    // FIXME AM: to be checked
                     uint32_t const* expanded_Sizes_v_DB_reco_flags,
                     uint32_t const* expanded_flagbit_v_DB_reco_flags,
                     uint32_t expanded_v_DB_reco_flagsSize,
                     uint32_t flagmask,
                     // conditions
                     float const* adc2gev,
                     float const* intercalib,
                     uint16_t const* status,
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
                     uint32_t const* flags_eb,
                     uint32_t const* flags_ee,
                     // output
                     uint32_t *did,
                     ::ecal::reco::StorageScalarType* energy,   // in energy [GeV]  
                     ::ecal::reco::StorageScalarType* time,  
                     ::ecal::reco::StorageScalarType* chi2,  
                     uint32_t* flagBits,
                     uint32_t* extra,
                     // other
                     int const nchannels,
                     uint32_t const offsetForInput,
                     uint32_t const offsetForHashes                     
    ) {
      
      
//       
//    NB: energy   "type_wrapper<reco::StorageScalarType, L>::type" most likely std::vector<float>
//       
      
      int ch = threadIdx.x + blockDim.x*blockIdx.x;

      if (ch < nchannels) {
        
        int const inputCh = ch >= offsetForInput
                          ? ch - offsetForInput
                          : ch;
                            
        uint32_t const * didCh = ch >= offsetForInput
                          ? did_ee
                          : did_eb;
        
        // only two values, EB or EE
        // AM : FIXME : why not using "isBarrel" ?    isBarrel ? adc2gev[0] : adc2gev[1]
        float adc2gev_to_use = ch >= offsetForInput
                          ? adc2gev[1]  // ee
                          : adc2gev[0]; // eb
                          
  
        // first EB and then EE
                          
        ::ecal::reco::StorageScalarType const* amplitude = ch >= offsetForInput
                          ? amplitude_ee
                          : amplitude_eb;
     
        ::ecal::reco::StorageScalarType const* time_in = ch >= offsetForInput
                          ? time_ee
                          : time_eb;
  
        ::ecal::reco::StorageScalarType const* chi2_in = ch >= offsetForInput
                          ? chi2_ee
                          : chi2_eb;
        
        uint32_t const* flags_in = ch >= offsetForInput
                          ? flags_ee
                          : flags_eb;
  
        // simple copy
        did[ch] = didCh[inputCh];
        
        auto const did_to_use = DetId{didCh[inputCh]};
        
        //
        // AM: move hashedIndexEB outside ecal::multifit   ---> change into reconstruction   [done]
        //
        
        auto const isBarrel = did_to_use.subdetId() == EcalBarrel;
        auto const hashedId = isBarrel
                            ? ecal::reconstruction::hashedIndexEB(did_to_use.rawId())
                            : offsetForHashes + ecal::reconstruction::hashedIndexEE(did_to_use.rawId());

        float const intercalib_to_use = intercalib[hashedId];
       
        
        // get laser coefficient
        float lasercalib = 1.;

        //
        // AM: ideas
        //
        //    One possibility is to create the map of laser corrections once on CPU
        //    for all crystals and push them on GPU.
        //    Then only if the LS is different, update the laser correction
        //    The variation within a LS is not worth pursuing (<< 0.1% !!)
        //    and below the precision we can claim on the laser corrections (right?).
        //    This will save quite some time (also for the CPU version?)    
        //
        
        //       
        //       https://github.com/cms-sw/cmssw/blob/5d27b6509171a20b6e6a4bbbaf29ca471d612913/CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
        //       ./CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
        //       
        
        int iLM = 1;
        
        if (isBarrel) {
          iLM = ecal::reconstruction::laser_monitoring_region_EB (did_to_use.rawId());
        }
        else {
          iLM = ecal::reconstruction::laser_monitoring_region_EE (did_to_use.rawId());
        }
        
        
        long long t_i = 0, t_f = 0;
        float p_i = 0, p_f = 0;
        long long lt_i = 0, lt_f = 0;
        float lp_i = 0, lp_f = 0;
        
        // laser
        if (event_time >= t1[iLM - 1] && event_time < t2[iLM - 1]) {
          t_i = t1[iLM - 1];
          t_f = t2[iLM - 1];
          p_i = p1[hashedId];
          p_f = p2[hashedId];
        } else if (event_time >= t2[iLM - 1] && event_time <= t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[hashedId];
          p_f = p3[hashedId];
        } else if (event_time < t1[iLM - 1]) {
          t_i = t1[iLM - 1];
          t_f = t2[iLM - 1];
          p_i = p1[hashedId];
          p_f = p2[hashedId];
          
        } else if (event_time > t3[iLM - 1]) {
          t_i = t2[iLM - 1];
          t_f = t3[iLM - 1];
          p_i = p2[hashedId];
          p_f = p3[hashedId];
        }

        
        // linear corrections
        if (event_time >= lt1[iLM - 1] && event_time < lt2[iLM - 1]) {
          lt_i = lt1[iLM - 1];
          lt_f = lt2[iLM - 1];
          lp_i = lp1[hashedId];
          lp_f = lp2[hashedId];
        } else if (event_time >= lt2[iLM - 1] && event_time <= lt3[iLM - 1]) {
          lt_i = lt2[iLM - 1];
          lt_f = lt3[iLM - 1];
          lp_i = lp2[hashedId];
          lp_f = lp3[hashedId];
        } else if (event_time < lt1[iLM - 1]) {
          lt_i = lt1[iLM - 1];
          lt_f = lt2[iLM - 1];
          lp_i = lp1[hashedId];
          lp_f = lp2[hashedId];
          
        } else if (event_time > lt3[iLM - 1]) {
          lt_i = lt2[iLM - 1];
          lt_f = lt3[iLM - 1];
          lp_i = lp2[hashedId];
          lp_f = lp3[hashedId];
        }
        
        
        // apdpnref and alpha 
        float apdpnref = apdpnrefs[hashedId];
        float alpha = alphas[hashedId];
        
        // now calculate transparency correction
        if (apdpnref != 0 && (t_i - t_f) != 0 && (lt_i - lt_f) != 0) {
          long long tt = event_time;  // never subtract two unsigned!
          float interpolatedLaserResponse =   p_i / apdpnref + float(tt - t_i)  * (p_f - p_i)   / (apdpnref * float(t_f - t_i));
          //           float interpolatedLaserResponse =   iLM+10; p_f*=2; p_i*=2;
//           float interpolatedLaserResponse =   iLM+float(t_f - t_i)+float(tt - t_i); p_f*=2; p_i*=2;
//           float interpolatedLaserResponse =   iLM+float(t_f - t_i)+float(tt - t_i)+apdpnref+alpha+100; p_f*=2; p_i*=2;
// float interpolatedLaserResponse =   iLM+float(t_f - t_i)+float(tt - t_i)+apdpnref+alpha+100+p_i; p_f*=2; p_i*=2;
// float interpolatedLaserResponse =   p_i+apdpnref; p_f*=2; p_i*=2;
// float interpolatedLaserResponse =   apdpnref; p_f*=2; p_i*=2; // good
// float interpolatedLaserResponse =   p_i; p_f*=2; p_i*=2; // bad
// float interpolatedLaserResponse =   iLM; p_f*=2; p_i*=2; // bad
// float interpolatedLaserResponse =   did_to_use.rawId()/1000000000.; p_f*=2; p_i*=2; // good
// float interpolatedLaserResponse =   iLM; p_f*=2; p_i*=2; // bad ??  ---> now good??????

// float interpolatedLaserResponse =   tt/6589055619675208704.; p_f*=2; p_i*=2; // good ... but not 100% sure
// float interpolatedLaserResponse =   1.+float(tt-t_i)*(p_f-p_i)/(apdpnref*float(t_f-t_i)); p_f*=2; p_i*=2; // bad ?


          float interpolatedLinearResponse = lp_i / apdpnref + float(tt - lt_i) * (lp_f - lp_i) / (apdpnref * float(lt_f - lt_i));  // FIXED BY FC
          
          if (interpolatedLinearResponse > 2.f || interpolatedLinearResponse < 0.1f) {
            interpolatedLinearResponse = 1.f;
          }
          if (interpolatedLaserResponse <= 0.) {
            // AM :  how the heck is it possible?
//             interpolatedLaserResponse = 0.0001;
            lasercalib = 1.;
            
          }
          else {
            
            float interpolatedTransparencyResponse = interpolatedLaserResponse / interpolatedLinearResponse;
            
            // ... and now this:
            lasercalib = 1.f / ( std::pow(interpolatedTransparencyResponse, alpha) * interpolatedLinearResponse);
            
//             alpha *= 2;
//             
//             lasercalib = interpolatedTransparencyResponse;
//             
//             lasercalib = 10.;
//             
//             lasercalib = interpolatedLaserResponse;

//             lasercalib = interpolatedLinearResponse; // this is ok
            
          }
        }
        
        
             
        //
        // Check for channels to be excluded from reconstruction
        //        
        //
        // Default energy? Not to be updated if "ChannelStatusToBeExcluded"
        // Exploited later by the module "EcalRecHitConvertGPU2CPUFormat"
        //
        energy[ch] = -1; //---- -1 or -100 FIXME AM
        
        //
        // From EcalRecHitWorkerSimple.cc
        //
        //             std::vector<int> v_chstatus_;
        //
        //     EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
        //     EcalChannelStatusCode::Code  dbstatus = chit->getStatusCode();
        //     
        //     // check for channels to be excluded from reconstruction
        //     if ( !v_chstatus_.empty()) {
        //       
        //       std::vector<int>::const_iterator res = 
        //       std::find( v_chstatus_.begin(), v_chstatus_.end(), dbstatus );
        //       if ( res != v_chstatus_.end() ) return false;
        //       
        //     }
        //
        //         /// return decoded status
        //             Code  getStatusCode() const { return Code(status_&chStatusMask); }
        //         
        //              static const int chStatusMask      = 0x1F;
        //
        //
        static const int chStatusMask      = 0x1F;
        // ChannelStatusToBeExcluded is a "int" then I put "dbstatus" to be the same
        int dbstatus = EcalChannelStatusCode_Code( (status[hashedId]) & chStatusMask );
//         printf (" ChannelStatusToBeExcludedSize = %d \n", ChannelStatusToBeExcludedSize) ; 
        if (ChannelStatusToBeExcludedSize != 0) {
          for (int ich_to_check = 0; ich_to_check<ChannelStatusToBeExcludedSize; ich_to_check++) {
//             printf (" ChannelStatusToBeExcluded[%d] = %d  --> dbstatus = %d \n", ich_to_check, ChannelStatusToBeExcluded[ich_to_check], dbstatus) ; 
            if ( ChannelStatusToBeExcluded[ich_to_check] == dbstatus ) {
//               printf (" killed @ GPU %d -->  ChannelStatusToBeExcluded[%d] = %d  --> dbstatus = %d \n", hashedId, ich_to_check, ChannelStatusToBeExcluded[ich_to_check], dbstatus) ; 
              return; 
            }
          }
        }
        
 
        // Take our association map of dbstatuses-> recHit flagbits and return the apporpriate flagbit word
        
        //
        // AM: get the smaller "flagbit_counter" with match
        //
        
        uint32_t temporary_flagBits = 0;
        
        int iterator_flags = 0;
        bool need_to_exit = false;
        int flagbit_counter = 0;
        while (!need_to_exit) {
          iterator_flags = 0;
          for (unsigned int i = 0; i != expanded_v_DB_reco_flagsSize; ++i) { 
            // check the correct "flagbit"
            if (expanded_flagbit_v_DB_reco_flags[i] == flagbit_counter) {
              
              for (unsigned int j = 0; j < expanded_Sizes_v_DB_reco_flags[i]; j++) {
                
                if ( expanded_v_DB_reco_flags[iterator_flags] == dbstatus ) {
                  temporary_flagBits =  0x1 << expanded_flagbit_v_DB_reco_flags[i];      
                  need_to_exit = true;
                  break; // also from the big loop!!!
                  
                }
                iterator_flags++;
              }
            }
            else {
              // if not, got to the next bunch directly
              iterator_flags += expanded_Sizes_v_DB_reco_flags[i];
            }
            
            if (need_to_exit) {
              break;
            }
            
          }
          flagbit_counter+=1;
        }
        
        
        //         flagBits[ch] =  dbstatus; // AM: FIXME : just a test
        
        
        //         for (unsigned int i = 0; i != v_DB_reco_flags_.size(); ++i) {
        //           if (std::find(v_DB_reco_flags_[i].begin(), v_DB_reco_flags_[i].end(), dbstatus) != v_DB_reco_flags_[i].end()) {
        //             flagBits[ch] =  0x1 << i;
        //             break;
        //           }
        //         }
        
        
        
        
//         class EcalRecHit {
//         public:
//           typedef DetId key_type;
//           
//           // recHit flags
//           enum Flags { 
//             kGood=0,                   // channel ok, the energy and time measurement are reliable
//             kPoorReco,                 // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
//             kOutOfTime,                // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
//             kFaultyHardware,           // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
//             kNoisy,                    // the channel is very noisy
//             kPoorCalib,                // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
//             kSaturated,                // saturated channel (recovery not tried)
//             kLeadingEdgeRecovered,     // saturated channel: energy estimated from the leading edge before saturation
//             kNeighboursRecovered,      // saturated/isolated dead: energy estimated from neighbours
//             kTowerRecovered,           // channel in TT with no data link, info retrieved from Trigger Primitive
//             kDead,                     // channel is dead and any recovery fails
//             kKilled,                   // MC only flag: the channel is killed in the real detector
//             kTPSaturated,              // the channel is in a region with saturated TP
//             kL1SpikeFlag,              // the channel is in a region with TP with sFGVB = 0
//             kWeird,                    // the signal is believed to originate from an anomalous deposit (spike) 
//             kDiWeird,                  // the signal is anomalous, and neighbors another anomalous signal  
//             kHasSwitchToGain6,         // at least one data frame is in G6
//             kHasSwitchToGain1,         // at least one data frame is in G1
//             //
//             kUnknown                   // to ease the interface with functions returning flags. 
//           };
//           
//           
//           
//           // Associate reco flagbit ( outer vector) to many db status flags (inner vector)
//           std::vector<std::vector<uint32_t> > v_DB_reco_flags_;
//           
//           // Traslate string representation of flagsMapDBReco into enum values 
//           const edm::ParameterSet & p=ps.getParameter< edm::ParameterSet >("flagsMapDBReco");
//           std::vector<std::string> recoflagbitsStrings = p.getParameterNames();
//           v_DB_reco_flags_.resize(32); 
//         
//         for (unsigned int i=0;i!=recoflagbitsStrings.size();++i){
//           EcalRecHit::Flags recoflagbit = (EcalRecHit::Flags)  StringToEnumValue<EcalRecHit::Flags>(recoflagbitsStrings[i]);
//           std::vector<std::string> dbstatus_s =  p.getParameter<std::vector<std::string> >(recoflagbitsStrings[i]);
//           std::vector<uint32_t> dbstatuses;
//           for (unsigned int j=0; j!= dbstatus_s.size(); ++j){
//             EcalChannelStatusCode::Code  dbstatus  = (EcalChannelStatusCode::Code)
//             StringToEnumValue<EcalChannelStatusCode::Code>(dbstatus_s[j]);
//             dbstatuses.push_back(dbstatus);
//           }
//           
//           v_DB_reco_flags_[recoflagbit]=dbstatuses;
//         }
//         
//         uint32_t flagBits = setFlagBits(v_DB_reco_flags_, dbstatus);
//         
//         // Take our association map of dbstatuses-> recHit flagbits and return the apporpriate flagbit word
//         uint32_t EcalRecHitWorkerSimple::setFlagBits(const std::vector<std::vector<uint32_t> >& map, 
//                                                      const uint32_t& status  ){
//           
//           for (unsigned int i = 0; i!=map.size(); ++i){
//             if (std::find(map[i].begin(), map[i].end(),status)!= map[i].end()) 
//               return 0x1 << i;
//           }
//           
//           return 0;
//                                                      }
                                                     
                                                     
                                                     
                                                     
        
        
        // AM: FIXME
        // why are these channels killed in a different way?
        // and not via configuration, as the other bits?
        // is it like this from the past development????
        // ----> make it uniform!
        //
        //         flagmask_ = 0;
        //         flagmask_ |= 0x1 << EcalRecHit::kNeighboursRecovered;
        //         flagmask_ |= 0x1 << EcalRecHit::kTowerRecovered;
        //         flagmask_ |= 0x1 << EcalRecHit::kDead;
        //         flagmask_ |= 0x1 << EcalRecHit::kKilled;
        //         flagmask_ |= 0x1 << EcalRecHit::kTPSaturated;
        //         flagmask_ |= 0x1 << EcalRecHit::kL1SpikeFlag;
        //       
        //         if (!(flagmask_ & flagBits) || !killDeadChannels_) {
        //
        //         ----> from EcalRecHitWorkerSimple.cc
        //
        //         -----> I have filled "flagmask" in the EcalRecHitProducerGPU.cc
        //
        //         bool killDeadChannels = true;
        
        //
        // defined here: https://github.com/amassiro/cmssw/blob/amassiro_patatrack_ECALrechit_on_GPU/RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitProducerGPU.cc#L265-L273
        // 
//          
//     In EcalRecHitWorkerSimple.cc as:
// 
//         if (! (flagmask_ & flagBits ) || !killDeadChannels_) {
//           
        
        if ( (flagmask & temporary_flagBits) && killDeadChannels ) {
                    return;
        }
        
        
        //
        // AM: Seriously, was it useless to set it before???
//         flagBits[ch] = 0;
        // FIXME
        // AM : no no ... now we store in flagBits[ch] what we had before in temporary_flagBits 
        flagBits[ch] = temporary_flagBits;
        // we don't leave it to 0 !!!
        
        
        
        
        
        
        //
        // multiply the adc counts with factors to get GeV
        //

//         energy[ch] = amplitude[inputCh] * adc2gev_to_use * intercalib_to_use ;
        energy[ch] = amplitude[inputCh] * adc2gev_to_use * intercalib_to_use * lasercalib;
        // AM: lasercalib *was* the wrong one!  ---> now fixed
        
        // Time is not saved so far, FIXME
//         time[ch] = time_in[inputCh];
 
        
        //      setChi2
        //         https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L126-L133
        //                 
        //         
        //         // bound the max value of the chi2
        //         if (chi2 > 64) chi2 = 64;
        //         
        // -> introduce bound also here in order to get cpu = gpu
        //
//         
//         chi2[ch] = chi2_in[inputCh];
//         
        if (chi2_in[inputCh] > 64) chi2[ch] = 64;
        else chi2[ch] = chi2_in[inputCh];
        
        
        
        
        // FIXME: calculate the "flagBits extra"  --> not really "flags", but actually an encoded version of energy uncertainty, time unc., ...
        extra[ch] = 0;  // -1 or 0 as default? ---> it is an unsgned int, it cannot be -1 !
        
        //
        // extra packing ...
        //

//      setChi2
//         https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L126-L133
//                 
//         
//         // bound the max value of the chi2
//         if (chi2 > 64) chi2 = 64;
//         
//         // use 7 bits
//         uint32_t rawChi2 = lround(chi2 / 64. * ((1<<7)-1));
//         extra_ = setMasked(extra_, rawChi2, 0, 7);
//         
//         

//         
//         static inline uint32_t setMasked(uint32_t value, uint32_t x, uint32_t offset, uint32_t width) {
//           const uint32_t mask = ((1 << width) - 1) << offset;
//           value &= ~mask;
//           value |= (x & ((1U << width) - 1)) << offset;
//           return value;
//         }
//         
        
        uint32_t offset;
        uint32_t width;
        uint32_t value;
        
        float chi2_temp = chi2[ch];
        if (chi2_temp > 64) chi2_temp = 64;
        // use 7 bits
        uint32_t rawChi2 = lround(chi2_temp / 64. * ((1<<7)-1));
        
//         extra_ = setMasked(extra_, rawChi2, 0, 7);
        
        offset = 0;
        width = 7;
        value = 0; // default: https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L65
        
        uint32_t mask = ((1 << width) - 1) << offset;
        value &= ~mask;
        value |= (rawChi2 & ((1U << width) - 1)) << offset;
        
//         extra[ch] = value;
//         
     
        
//      setEnergyError
//         rh.setEnergyError( uncalibRH.amplitudeError()*adcToGeVConstant_*intercalibConstant);
//         https://github.com/cms-sw/cmssw/blob/266e21cfc9eb409b093e4cf064f4c0a24c6ac293/DataFormats/EcalRecHit/interface/EcalRecHit.h#L141-L163
//         
        
        // rawEnergy is actually "error" !!!
        uint32_t rawEnergy = 0;
        

        // AM: FIXME: this is not propagated currently to the uncalibrecHit collection SOA 
        //            if you want to store this in "extra", we need first to add it to the uncalibrecHit results
        //            then it will be something like the following
//         amplitudeError[inputCh] * adc2gev_to_use * intercalib_to_use * lasercalib
//         
//         
        
        float amplitudeError_ch = 0. ; // amplitudeError[ch];
        
        if (amplitudeError_ch > 0.001) {
//           uint16_t exponent = getPower10(amplitudeError_ch);

          static constexpr float p10[] = {1.e-2f,1.e-1f,1.f,1.e1f,1.e2f,1.e3f,1.e4f,1.e5f,1.e6f};
          int b = amplitudeError_ch<p10[4] ? 0 : 5;
          for (;b<9;++b) if (amplitudeError_ch<p10[b]) break;
          
          uint16_t exponent = b;
          
          static constexpr float ip10[] = {1.e5f,1.e4f,1.e3f,1.e2f,1.e1f,1.e0f,1.e-1f,1.e-2f,1.e-3f,1.e-4};
          uint16_t significand = lround( amplitudeError_ch * ip10[exponent]);
          // use 13 bits (3 exponent, 10 significand)
          rawEnergy = exponent << 10 | significand;
          /* here for doc and regression test
           *             u int16_t exponent_old = lround(floor(log10(*energy))) + 3;  
           *             uint16_t significand_old = lround(energy/pow(10, exponent - 5));
           *             std::cout << energy << ' ' << exponent << ' ' << significand 
           *             << ' ' << exponent_old <<     ' ' << significand_old << std::endl;
           *             assert(exponent==exponent_old);
           *             assert(significand==significand_old);
           */
        }
        
//         extra_ = setMasked(extra_, rawEnergy, 8, 13);
        
        offset = 8;
        width = 13;
        // value from last change, ok
        
        mask = ((1 << width) - 1) << offset;
        value &= ~mask;
        value |= (rawEnergy & ((1U << width) - 1)) << offset;


//         // set the energy uncertainty
//         // (only energy >= 0 will be stored)
//         void setEnergyError(float energy) {
//           uint32_t rawEnergy = 0;
//           if (energy > 0.001) {
//             uint16_t exponent = getPower10(energy);
//             static constexpr float ip10[] = {1.e5f,1.e4f,1.e3f,1.e2f,1.e1f,1.e0f,1.e-1f,1.e-2f,1.e-3f,1.e-4};
//             uint16_t significand = lround(energy*ip10[exponent]);
//             // use 13 bits (3 exponent, 10 significand)
//             rawEnergy = exponent << 10 | significand;
//           }
//           
//           extra_ = setMasked(extra_, rawEnergy, 8, 13);
//         }
        

//         
//         void setTimeError(uint8_t timeErrBits) {
//           extra_ = setMasked(extra_, timeErrBits & 0xFF, 24, 8);
//         }
//         
//         uint8_t uncalibRH.jitterErrorBits() ---> set to 0 if not available? FIXME: check default value
//        
//         EcalRecHit rh( uncalibRH.id(), energy, time );
//         rh.setChi2( uncalibRH.chi2() );
//         rh.setEnergyError( uncalibRH.amplitudeError()*adcToGeVConstant_*intercalibConstant);
//         rh.setTimeError(uncalibRH.jitterErrorBits());
//         
//         
        
        uint32_t jitterErrorBits = 0;
        jitterErrorBits = jitterErrorBits & 0xFF;
        
        
        offset = 24;
        width = 8;
        // value from last change, ok
        
        mask = ((1 << width) - 1) << offset;
        value &= ~mask;
        value |= (jitterErrorBits & ((1U << width) - 1)) << offset;
        
        //
        // now finally set "extra[ch]"
        //
        extra[ch] = value;
            
          
        
        //
        // additional flags setting
        //
        // using correctly the flags as calculated at the UncalibRecHit stage
        //
        // Now fill flags
        
        
        
        
        bool good = true;

//
//         from: EcalUncalibratedRecHit.cc
//         
//         
//         bool EcalUncalibratedRecHit::checkFlag(EcalUncalibratedRecHit::Flags flag) const {
//           if(flag == kGood){ if ( ! flags_ ) return true;else return false;} // if all flags are unset, then hit is good
//           return  flags_ & ( 0x1<<flag);
//         }
//         
//         
//         
//         void EcalUncalibratedRecHit::setFlagBit(EcalUncalibratedRecHit::Flags flag){
//           if  (flag == kGood) {
//             //then set all bits to zero;
//             flags_  = 0;
//             return;
//           }
//           // else set the flagbit
//           flags_|= 0x1 <<  flag;  
//         }
//         
//         
//  
//         enum Flags {
//           kGood=-1,                 // channel is good (mutually exclusive with other states)  setFlagBit(kGood) reset flags_ to zero 
//           kPoorReco,                // channel has been badly reconstructed (e.g. bad shape, bad chi2 etc.)
//           kSaturated,               // saturated channel
//           kOutOfTime,               // channel out of time
//           kLeadingEdgeRecovered,    // saturated channel: energy estimated from the leading edge before saturation
//           kHasSwitchToGain6,        // at least one data frame is in G6
//           kHasSwitchToGain1         // at least one data frame is in G1
//           
//         };
//         
        
      
        if ( flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kLeadingEdgeRecovered) ) ) {
          flagBits[ch]  |=  (0x1 <<  (RecHitFlags::RecHitFlags_kLeadingEdgeRecovered));  
          good = false;          
        }
        
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kSaturated) ) ) {
          // leading edge recovery failed - still keep the information
          // about the saturation and do not flag as dead
          flagBits[ch]  |=  (0x1 <<  (RecHitFlags::RecHitFlags_kSaturated));  
          good = false;
        }

        //
        // AM: why do we have two tests one after the other checking almost the same thing??? 
        // Please clean up the code, ... also the original one!
        //        
        // uncalibRH.isSaturated() ---> 
        //         
        //                                   bool EcalUncalibratedRecHit::isSaturated() const {
        //                                     return EcalUncalibratedRecHit::checkFlag(kSaturated);
        //                                   }
        //
        //
        
//         if (( 0x1 << (UncalibRecHitFlags::kSaturated) ) ) {
        if ( flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kSaturated) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kSaturated));  
          good = false;
        }
        
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kOutOfTime) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kOutOfTime));
          good = false;
        }
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kPoorReco) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kPoorReco));
          good = false;
        }
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kHasSwitchToGain6) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kHasSwitchToGain6));
        }
        if (flags_in[inputCh] & ( 0x1 << (UncalibRecHitFlags::kHasSwitchToGain1) ) ) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kHasSwitchToGain1));
        }
        
        
        if (good) {
          flagBits[ch] |= (0x1 << (RecHitFlags::RecHitFlags_kGood));
        }
        
        
          
//         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kLeadingEdgeRecovered)) {
//           rh.setFlag(EcalRecHit::kLeadingEdgeRecovered);
//           good = false;
//         }
//         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kSaturated)) {
//           // leading edge recovery failed - still keep the information
//           // about the saturation and do not flag as dead
//           rh.setFlag(EcalRecHit::kSaturated);
//           good = false;
//         }
//         if (uncalibRH.isSaturated()) {
//           rh.setFlag(EcalRecHit::kSaturated);
//           good = false;
//         }
//         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kOutOfTime)) {
//           rh.setFlag(EcalRecHit::kOutOfTime);
//           good = false;
//         }
//         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kPoorReco)) {
//           rh.setFlag(EcalRecHit::kPoorReco);
//           good = false;
//         }
//         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kHasSwitchToGain6)) {
//           rh.setFlag(EcalRecHit::kHasSwitchToGain6);
//         }
//         if (uncalibRH.checkFlag(EcalUncalibratedRecHit::kHasSwitchToGain1)) {
//           rh.setFlag(EcalRecHit::kHasSwitchToGain1);
//         }
//         
//         if (good)
//           rh.setFlag(EcalRecHit::kGood);
        
        
//        
//        @ ecalrechit
//        void setFlag(int flag) {flagBits_|= (0x1 << flag);}
//        
        
        
        
        
        
        
        
        //         if (detid.subdetId() == EcalBarrel && (lasercalib < EBLaserMIN_ || lasercalib > EBLaserMAX_)) 
        //           myrechit.setFlag(EcalRecHit::kPoorCalib);
        //         if (detid.subdetId() == EcalEndcap && (lasercalib < EELaserMIN_ || lasercalib > EELaserMAX_)) 
        //           myrechit.setFlag(EcalRecHit::kPoorCalib);
        //         
        //         void setFlag(int flag) {flagBits_|= (0x1 << flag);}
        
         
        if (isBarrel  && (lasercalib < EBLaserMIN || lasercalib > EBLaserMAX)) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kPoorCalib));
          
        }
        if (!isBarrel && (lasercalib < EELaserMIN || lasercalib > EELaserMAX)) {
          flagBits[ch]  |= (0x1 <<  (RecHitFlags::RecHitFlags_kPoorCalib));
        }
        
        
        
        
        //-- AM TEST just to test
//         flagBits[ch] = 42;
//         flagBits[ch] = 500;
        //--
         
        
        
      }
      
    }
    
    
    
    // host version, to be called by the plugin
    void create_ecal_rehit(
                  EventInputDataGPU const& eventInputGPU,
                  EventOutputDataGPU&      eventOutputGPU,
                  //     eventDataForScratchGPU_,
                  ConditionsProducts const& conditions, 
                  ConfigurationParameters const& configParameters,
                  uint32_t const  offsetForInput,
                  edm::TimeValue_t const event_time,
                  cuda::stream_t<>& cudaStream
             ){
    
      int nchannels = eventInputGPU.ebUncalibRecHits.size + eventInputGPU.eeUncalibRecHits.size ;
            
      unsigned int nchannels_per_block = 32;
      unsigned int threads_1d = nchannels_per_block;
      unsigned int blocks_1d = (nchannels + threads_1d) / threads_1d; // TEST : to be optimized (AM)
      
      
      // 
      // kernel create rechit
      //
      kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
// configuration 
        configParameters.ChannelStatusToBeExcluded,
        configParameters.ChannelStatusToBeExcludedSize,
        configParameters.killDeadChannels,
        configParameters.EBLaserMIN,
        configParameters.EELaserMIN,
        configParameters.EBLaserMAX,
        configParameters.EELaserMAX,
// for flags setting
        configParameters.expanded_v_DB_reco_flags,
        configParameters.expanded_Sizes_v_DB_reco_flags,
        configParameters.expanded_flagbit_v_DB_reco_flags,
        configParameters.expanded_v_DB_reco_flagsSize,
        configParameters.flagmask,
        // conditions
        conditions.ADCToGeV.adc2gev,
        conditions.Intercalib.values,  
        conditions.ChannelStatus.status,  
        conditions.LaserAPDPNRatiosRef.values,  
        conditions.LaserAlphas.values,  
// input for transparency corrections
        conditions.LaserAPDPNRatios.p1,
        conditions.LaserAPDPNRatios.p2,
        conditions.LaserAPDPNRatios.p3,
        conditions.LaserAPDPNRatios.t1,
        conditions.LaserAPDPNRatios.t2,
        conditions.LaserAPDPNRatios.t3,
// input for linear corrections
        conditions.LinearCorrections.p1,
        conditions.LinearCorrections.p2,
        conditions.LinearCorrections.p3,
        conditions.LinearCorrections.t1,
        conditions.LinearCorrections.t2,
        conditions.LinearCorrections.t3,
// time, used for time dependent corrections
        event_time,
// input
        eventInputGPU.ebUncalibRecHits.did,
        eventInputGPU.eeUncalibRecHits.did,
        eventInputGPU.ebUncalibRecHits.amplitude, 
        eventInputGPU.eeUncalibRecHits.amplitude, 
        eventInputGPU.ebUncalibRecHits.jitter, 
        eventInputGPU.eeUncalibRecHits.jitter, 
        eventInputGPU.ebUncalibRecHits.chi2, 
        eventInputGPU.eeUncalibRecHits.chi2, 
        eventInputGPU.ebUncalibRecHits.flags, 
        eventInputGPU.eeUncalibRecHits.flags, 
// output
        eventOutputGPU.did,
        eventOutputGPU.energy,
        eventOutputGPU.time,
        eventOutputGPU.chi2,
        eventOutputGPU.flagBits,
        eventOutputGPU.extra,
// other
        nchannels,
        offsetForInput,
        conditions.offsetForHashes
      );
      
      
      
    }
    
    

    
  }
  
}

