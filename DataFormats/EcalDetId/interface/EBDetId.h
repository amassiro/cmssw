#ifndef ECALDETID_EBDETID_H
#define ECALDETID_EBDETID_H

#include <iosfwd>
#include <cmath>
#include <cstdlib>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "FWCore/Utilities/interface/Exception.h"


/** \class EBDetId
 *  Crystal identifier class for the ECAL barrel
 *
 *
 */



class EBDetId : public DetId {
 public:
  enum { Subdet=EcalBarrel};
  /** Constructor of a null id */
  constexpr EBDetId() {}
  /** Constructor from a raw value */
  constexpr EBDetId(uint32_t rawid) : DetId(rawid) {}
  /** Constructor from crystal ieta and iphi 
      or from SM# and crystal# */
  // fast
  constexpr EBDetId(  int crystal_ieta, int crystal_iphi) : DetId(Ecal,EcalBarrel) {
    id_|=((crystal_ieta>0)?(0x10000|(crystal_ieta<<9)):((-crystal_ieta)<<9))|(crystal_iphi&0x1FF);
  }
  // slow
  constexpr EBDetId(int index1, int index2, int mode)
    : DetId(Ecal,EcalBarrel)
    {
      int crystal_ieta = 0;
      int crystal_iphi = 0;
      if (mode == ETAPHIMODE) {
        crystal_ieta = index1;
        crystal_iphi = index2;  
      } else if (mode == SMCRYSTALMODE) {
        int SM = index1;
        int crystal = index2;
        int i = (int)  floor((crystal-1) / kCrystalsInPhi);
        int j = ((crystal-1) - (kCrystalsInPhi*i));
        if (SM <= 18) {
          crystal_ieta = i + 1;
          crystal_iphi = ((SM-1) * kCrystalsInPhi) + (kCrystalsInPhi-j);
        } else {
          crystal_ieta = -(i+1);
          crystal_iphi = ((SM-19) * kCrystalsInPhi) + j+1;
        }
      } else {
        throw cms::Exception("InvalidDetId") << "EBDetId:  Cannot create object.  Unknown mode for (int, int) constructor."; 
      }
      
      if ( !validDetId(crystal_ieta, crystal_iphi) ) {
        //    std::cout << "crystal_eta " << crystal_ieta << "crystal_phi " << crystal_iphi << std::endl;
        throw cms::Exception("InvalidDetId") << "EBDetId:  Cannot create object.  Indexes out of bounds \n" 
        << "eta = " << crystal_ieta << " phi = " << crystal_iphi;
      }
      id_|=((crystal_ieta>0)?(0x10000|(crystal_ieta<<9)):((-crystal_ieta)<<9))|(crystal_iphi&0x1FF);
  }
  /** Constructor from a generic cell id */
  constexpr EBDetId(const DetId& id) : DetId(id){}
  /** Assignment operator from cell id */
  constexpr EBDetId& operator=(const DetId& id) {
    id_=id.rawId();
  return *this;
  }

  /// get the subdetector .i.e EcalBarrel (what else?)
  constexpr EcalSubdetector subdet() const { return (EcalSubdetector)(subdetId()); }
//   constexpr EcalSubdetector subdet() { return EcalBarrel;}

  /// get the z-side of the crystal (1/-1)
  constexpr int zside() const { return (id_&0x10000)?(1):(-1); }
  /// get the absolute value of the crystal ieta
  constexpr int ietaAbs() const { return (id_>>9)&0x7F; }
  /// get the crystal ieta
  constexpr int ieta() const { return zside()*ietaAbs(); }
  /// get the crystal iphi
  constexpr int iphi() const { return id_&0x1FF; }
  /// get the HCAL/trigger ieta of this crystal
  constexpr int tower_ieta() const { return ((ietaAbs()-1)/5+1)*zside(); }
  /// get the HCAL/trigger iphi of this crystal
  constexpr int tower_iphi() const {
    int iphi_simple=((iphi()-1)/5)+1; 
    iphi_simple-=2;
    return ((iphi_simple<=0)?(iphi_simple+72):(iphi_simple));
  }
  /// get the HCAL/trigger iphi of this crystal
//   constexpr
  EcalTrigTowerDetId tower() const { return EcalTrigTowerDetId(zside(),EcalBarrel,abs(tower_ieta()),tower_iphi()); }
  /// get the ECAL/SM id
  constexpr int ism() const  {
    int id = ( iphi() - 1 ) / kCrystalsInPhi + 1;
    return positiveZ() ? id : id+18;
  }
  /// get the number of module inside the SM (1-4)
  constexpr int im() const {
    int ii = (ietaAbs()-26);
    return ii<0 ? 1 : (ii/20 +2);
  }
  /// get ECAL/crystal number inside SM
  // Following TB 2004  numbering scheme 
  constexpr int ic() const {
    int ie = ietaAbs() -1;
    return  (ie * kCrystalsInPhi)
    +  ( positiveZ() ?
    ( kCrystalsInPhi - ( (iphi() -1 ) % kCrystalsInPhi ) )
    : ( ( iphi() -1 ) % kCrystalsInPhi  + 1)  
    );
  }
  
  /// get the crystal ieta in the SM convention (1-85)
  constexpr int ietaSM() const { return ietaAbs(); }
  /// get the crystal iphi (1-20)
  constexpr int iphiSM() const { return (( ic() -1 ) % kCrystalsInPhi ) + 1; }
  
  // is z positive?
  constexpr bool positiveZ() const { return id_&0x10000;}
  // crystal number in eta-phi grid
  constexpr int numberByEtaPhi() const { 
    return (MAX_IETA + (positiveZ() ? ietaAbs()-1 : -ietaAbs()) )*MAX_IPHI+ iphi()-1;
  }
  // index numbering crystal by SM
  // Maintains SM crystals in bunch of 1700 indices
  constexpr int numberBySM() const {
    return (ism()-1) * kCrystalsPerSM + ic() -1;   
  }
  /// get a compact index for arrays
  constexpr int hashedIndex() const { return numberByEtaPhi(); }

  constexpr uint32_t denseIndex() const { return hashedIndex() ; }

  /** returns a new EBDetId offset by nrStepsEta and nrStepsPhi (can be negative), 
    * returns EBDetId(0) if invalid */
  constexpr EBDetId offsetBy( int nrStepsEta, int nrStepsPhi ) const {
    int newEta = ieta()+nrStepsEta;
    if( newEta*ieta() <= 0 ) {
      if( ieta() < 0 ) {
        newEta++;
      } else if ( ieta() > 0 ) {
        newEta--;
      }
    }
    int newPhi = iphi() + nrStepsPhi;
    while ( newPhi>360 ) newPhi -= 360;
    while ( newPhi<=0  ) newPhi += 360;
    
    if( validDetId( newEta, newPhi ) ) {
      return EBDetId( newEta, newPhi);
    } else {
      return EBDetId(0);
    }
  }

  /** returns a new EBDetId on the other zside of barrel (ie iEta*-1), 
    * returns EBDetId(0) if invalid (shouldnt happen) */
  constexpr EBDetId switchZSide() const {
    int newEta = ieta()*-1;
    if( validDetId( newEta, iphi() ) ) {
      return EBDetId( newEta, iphi() );
    } else {
      return EBDetId(0);
    }
  }
 
  /** following are static member functions of the above two functions
    * which take and return a DetId, returns DetId(0) if invalid 
    */
  constexpr static DetId offsetBy( const DetId startId, int nrStepsEta, int nrStepsPhi ) {
    if( startId.det() == DetId::Ecal && startId.subdetId() == EcalBarrel ) {
      EBDetId ebStartId(startId);
      return ebStartId.offsetBy( nrStepsEta, nrStepsPhi ).rawId();
    } else {
      return DetId(0);
    }
  }
  constexpr static DetId switchZSide( const DetId startId ) {
    if( startId.det() == DetId::Ecal && startId.subdetId() == EcalBarrel ) {
      EBDetId ebStartId(startId);
      return ebStartId.switchZSide().rawId();
    } else {
      return DetId(0);
    }
  }

  /** return an approximate values of eta (~0.15% precise)
   */
//   constexpr float approxEta() const { return ieta() * (crystalUnitToEta); }
  constexpr float approxEta() const { return ieta() * (0.017453292519943295); }

  constexpr static float approxEta( const DetId id ) {
    if( id.subdetId() == EcalBarrel ) {
      EBDetId ebId( id );
      return ebId.approxEta();
    } else {
      return 0;
    }
  }

  constexpr static bool validDenseIndex( uint32_t din ) { return ( din < kSizeForDenseIndexing ) ; }

  constexpr static EBDetId detIdFromDenseIndex( uint32_t di ) { return unhashIndex( di ) ; }

  /// get a DetId from a compact index for arrays
  constexpr static EBDetId unhashIndex( int hi ) {
    const int pseudo_eta = hi/MAX_IPHI - MAX_IETA;
    return ( validHashIndex( hi ) ?
	     EBDetId(pseudo_eta<0 ? pseudo_eta :  pseudo_eta+1, hi%MAX_IPHI+1) :
	     EBDetId() ) ;
  }

  constexpr static bool validHashIndex(int i) { return !(i<MIN_HASH || i>MAX_HASH); }

  /// check if a valid index combination
  constexpr static bool validDetId(int i, int j) {
    return i!=0 && (std::abs(i) <= MAX_IETA)
      && (j>=MIN_IPHI) && (j <= MAX_IPHI); 
  }

  constexpr static bool isNextToBoundary(EBDetId id) {
    return isNextToEtaBoundary( id ) || isNextToPhiBoundary( id );
  }

  constexpr static bool isNextToEtaBoundary(EBDetId id) {
    int ieta = id.ietaSM();
//     return ieta == 1 || (kModuleBoundaries + 4)!=std::find( kModuleBoundaries, kModuleBoundaries + 4, ieta );
    return ieta == 1 || (ieta == 25) || (ieta == 45) || (ieta == 65) || (ieta == 85);
  }

  constexpr static bool isNextToPhiBoundary(EBDetId id) {
    int iphi = id.iphiSM();
    return iphi == 1 || iphi == 20;
  }

  //return the distance in eta units between two EBDetId
  constexpr static int distanceEta(const EBDetId& a,const EBDetId& b) {
    if (a.ieta() * b.ieta() > 0)
      return abs(a.ieta()-b.ieta());
    else
      return abs(a.ieta()-b.ieta())-1;
  }
  //return the distance in phi units between two EBDetId
  constexpr static int distancePhi(const EBDetId& a,const EBDetId& b) {
    int PI = 180;
    int  result = a.iphi() - b.iphi();
    
    while  (result > PI)    result -= 2*PI;
    while  (result <= -PI)  result += 2*PI;
    return abs(result);
  }

  /// range constants
  static const int MIN_IETA = 1;
  static const int MIN_IPHI = 1;
  static const int MAX_IETA = 85;
  static const int MAX_IPHI = 360;
  static const int kChannelsPerCard = 5;
  static const int kTowersInPhi = 4;  // per SM
  static const int kModulesPerSM = 4;
  static const int kModuleBoundaries[4] ;
  static const int kCrystalsInPhi = 20; // per SM
  static const int kCrystalsInEta = 85; // per SM
  static const int kCrystalsPerSM = 1700;
  static const int MIN_SM = 1;
  static const int MAX_SM = 36;
  static const int MIN_C = 1;
  static const int MAX_C = kCrystalsPerSM;
  static const int MIN_HASH =  0; // always 0 ...
  static const int MAX_HASH =  2*MAX_IPHI*MAX_IETA-1;

  // eta coverage of one crystal (approximate)
  static const float crystalUnitToEta;
//   static const float crystalUnitToEta = 0.017453292519943295;
  
  enum { kSizeForDenseIndexing = MAX_HASH + 1 } ;
  

  // function modes for (int, int) constructor
  static const int ETAPHIMODE = 0;
  static const int SMCRYSTALMODE = 1;
};

std::ostream& operator<<(std::ostream& s,const EBDetId& id);

// const int EBDetId::kModuleBoundaries[4] = { 25, 45, 65, 85 };

// pi / 180.
// const float EBDetId::crystalUnitToEta = 0.017453292519943295;



#endif
