// -*- Mode: C++; c-basic-offset: 2; indent-tabs-mode: t; tab-width: 8; -*-
//
// \author Philippe Gras (CEA/Saclay). Code adapted from EEDetId.
//
#ifndef EcalDetId_EcalScDetId_h
#define EcalDetId_EcalScDetId_h

#include <iosfwd>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "FWCore/Utilities/interface/thread_safety_macros.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <ostream>
#include <cassert>
#include <mutex>


static std::once_flag initializedFlag;



/** \class EcalScDetId
 *  Supercrystal identifier class for the ECAL endcap.
 *  <P>Note: internal representation of ScDetId:
 *  \verbatim
 *  31              .               15              .              0
 *  |-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-|-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-| 
 *  |  det  |sudet|         0       |1|z|     ix      |     iy      |
 *  +-------+-----+-----------------+-+-+-------------+-------------+
 *  \endverbatim
 */

class EcalScDetId : public DetId {
 public:

  /** Constructor of a null id
   */
  constexpr EcalScDetId() {
  }
  
  /** Constructor from a raw value
   * @param rawid det ID number of the supecrystal, as defined in this class
   * description.
   */
  constexpr EcalScDetId(uint32_t rawid) {
  }

  /** Constructor from supercrystal ix,iy,iz (iz=+1/-1)
   * ix x-index runs from 1 to 20 along x-axis of standard CMS coordinates
   * iy y-index runs from 1 to 20 along y-axis of standard CMS coordinates
   * iz z-index (also called "z-side") is -1 for EE- and +1 for EE+
   * @param ix x-index
   * @param iy y-index
   * @param iz z-side /z-index: -1 for EE-, +1 for EE+
   */
  constexpr EcalScDetId(int ix, int iy, int iz)
  {
    if(!validDetId(ix,iy,iz))
    {
      throw cms::Exception("InvalidDetId") << "EcalScDetId:  Cannot create object.  Indexes out of bounds \n" 
      << "x = " << ix << " y = " << iy << " z = " << iz;
    }
    const int scBit = 1<<15; //bit set to 1 to distinguish from crystal id (EEDetId)
    //                         and for a reasonale behaviour of DetId ccomparison operators.
    id_|=(iy&0x7f)|((ix&0x7f)<<7)|((iz>0)?(1<<14):(0))|scBit;
  }
  
  /** Constructor from a raw value
   * @param id det ID number
   */
  constexpr EcalScDetId(const DetId& gen)
  {
    if (!gen.null() && (gen.det()!=Ecal || gen.subdetId()!=EcalEndcap)) {
      throw cms::Exception("InvalidDetId"); 
    }
    id_=gen.rawId();
  }

  /** Assignment operator
   * @param id source det id
   */ 
  constexpr EcalScDetId& operator=(const DetId& gen)
  {
    if (!gen.null() && ( gen.det()!=Ecal || gen.subdetId()!=EcalEndcap )) {
      throw cms::Exception("InvalidDetId"); 
    }
    id_=gen.rawId();
    return *this;
  }

  /** Gets the subdetector
   * @return subdetectot ID, that is EcalEndcap
   */
  constexpr EcalSubdetector subdet() const { return EcalSubdetector(subdetId()); }
  
  /** Gets the z-side of the crystal (1/-1)
   * @return -1 for EE-, +1 for EE+
   */
  constexpr int zside() const { return (id_&0x4000)?(1):(-1); }
  
  /** Gets the crystal x-index.
   * @see EcalDetId(int, int, int) for x-index definition
   * @return x-index
   */
  constexpr int ix() const { return (id_>>7)&0x7F; }
  
  /** Get the crystal y-index
   * @see EcalDetId(int, int, int) for y-index definition.
   * @return y-index
   */
  constexpr int iy() const { return id_&0x7F; }
  
  /** Gets the quadrant of the DetId.
   *
   * Quadrant number definition for EE+, x and y in std CMS coordinates:
   *
   * \verbatim
   *                 A y
   *                 |
   *          Q2     |    Q1
   *                 |
   *       ----------o---------> x
   *                 |
   *          Q3     |    Q4
   *                 |
   * \endverbatim
   * This method will return the same quadrant number independently of
   * z: that is two supercrystals which are face to face will be considered
   * will have the same quadrant number. It is not clear it is the correct
   * or usual definition.
   * @see EEDetId::iquadrant()
   * @return quadrant number, from 1 to 4.
   * @deprecated This method might be withdraw in a future release
   */
  constexpr int iquadrant() const 
  {
    const int xMiddle = IX_MAX/2; //y = 0 between xMiddle and xMiddle+1
    const int yMiddle = IY_MAX/2; //x = 0 between yMiddle and yMiddle+1
    if (iy()>yMiddle){// y>0
      if(ix()>xMiddle)   //             A y            
        return 1;        //             |              
        else               //      Q2     |    Q1        
          return 2;        //             |              
    } else{// y<0        //   ----------o---------> x   
      if(ix()>xMiddle)   //             |               
        return 4;        //      Q3     |    Q4       
        else               //             |               
          return 3;
    }
    //Should never be reached
    return -1;
  }  
  
  
  /** Gets a compact index for arrays. Index runs from 0 to 623.
   * They are ordered by increasing z (EE- then EE+), then for
   * same z by increasing y. then for same z and y by increasing x
   */
  int hashedIndex() const{
    checkHashedIndexMap();
    if(!validDetId(ix(),iy(),zside())) return -1;
    return xyz2HashedIndex[ix()-IX_MIN][iy()-IY_MIN][zside()>0?1:0];
  }

  /** Gets EcalScDetId from hasedIndex as defined by hashedIndex method
   * @param hi hashed index
   * @return the EcalScDetId. If hi is invalid return a null EcalScDetId.
   */
  static EcalScDetId unhashIndex(int hi){
    checkHashedIndexMap();
    if(hi < 0 || hi >= kSizeForDenseIndexing) return EcalScDetId();
    return hashedIndex2DetId[hi];
  }
  
  /** Same as hashed index.
   * @return the dense/hashed index
   */
  uint32_t denseIndex() const { return hashedIndex() ; }

  /** Validates a hashed index.
   * @param din hashed index to validate
   * @return true if the index is valid, false if it is invalid.
   */
  constexpr static bool validDenseIndex(uint32_t din) { return din < kSizeForDenseIndexing; }

  
  /** Validates a hashed index.
   * @param hi hashed index to validate
   * @return true if the index is valid, false if it is invalid.
   */
  constexpr static bool validHashIndex(int hi) { return validDenseIndex(hi) ; }

  /** Number of supercrystals per endcap
   */
  constexpr static const int SC_PER_EE_CNT = 312;
  
  /** Lower bound of EE supercrystal x-index
   */
  static const int IX_MIN=1;

  /** Lower bound of EE supercrystal y-index
   */
  static const int IY_MIN=1;

  /** Upper bound of EE crystal y-index
   */
  static const int IX_MAX=20;

  /** Upper bound of EE crystal y-index
   */
  static const int IY_MAX=20;

  /** Lower bound for hashed/dense index
   */
  static const int IHASHED_MIN = 0;

  /** Upper bound for hashed/dense index
   */
  static const int IHASHED_MAX = SC_PER_EE_CNT*2 - 1;
  
  /** Checks validity of a crystal (x,y.z) index triplet.
   * @param ix supercrystal x-index
   * @param iy supercrystal y-index
   * @param iz supercrystal z-index (aka z-side)
   * @see EEDetId(int, int, int) for index definition
   * @return true if valid, false otherwise
   */
  static bool validDetId(int iX, int iY, int iZ) 
  {
    const char endcapMap[401] = {
      "       XXXXXX       "
      "    XXXXXXXXXXXX    "
      "   XXXXXXXXXXXXXX   "
      "  XXXXXXXXXXXXXXXX  "
      " XXXXXXXXXXXXXXXXXX "
      " XXXXXXXXXXXXXXXXXX "             //    Z
      " XXXXXXXXXXXXXXXXXX "             //     x-----> X
      "XXXXXXXXXXXXXXXXXXXX"             //     |
      "XXXXXXXXX  XXXXXXXXX"             //     |
      "XXXXXXXX    XXXXXXXX"//_          //     |
      "XXXXXXXX    XXXXXXXX"             //     V Y
      "XXXXXXXXX  XXXXXXXXX"
      "XXXXXXXXXXXXXXXXXXXX"
      " XXXXXXXXXXXXXXXXXX "
      " XXXXXXXXXXXXXXXXXX "
      " XXXXXXXXXXXXXXXXXX "
      "  XXXXXXXXXXXXXXXX  "
      "   XXXXXXXXXXXXXX   "
      "    XXXXXXXXXXXX    "
      "       XXXXXX       "};
      
      return std::abs(iZ)==1 && endcapMap[iX-1+(iY-1)*20]!=' ';
  }

private:
  /** Initializes x,y,z <-> hashed index map if not yet done.
   */
  
  static void checkHashedIndexMap()
  {
    std::call_once(initializedFlag, []() 
    {
      int hashedIndex = -1;
      for(int iZ = -1; iZ <= +1; iZ+=2){
        for(int iY = IY_MIN; iY <= IY_MAX; ++iY){
          for(int iX = IX_MIN; iX <= IX_MAX; ++iX){
            if(validDetId(iX,iY,iZ)){
              xyz2HashedIndex[iX-IX_MIN][iY-IY_MIN][iZ>0?1:0] = ++hashedIndex;
              assert((unsigned)hashedIndex < sizeof(hashedIndex2DetId)/sizeof(hashedIndex2DetId[0]));
              hashedIndex2DetId[hashedIndex] = EcalScDetId(iX, iY, iZ);
            }
          }
        }
      }
    });
  }
  

  //fields
public:
  enum {
    /** Number of dense supercrystal indices.
     */
    kSizeForDenseIndexing = SC_PER_EE_CNT * 2
  };

private:
  static const int nEndcaps = 2;
  
  /** Map of z,x,y index to hashed index. See hashedIndex/
   */
  CMS_THREAD_SAFE static short xyz2HashedIndex[IX_MAX][IY_MAX][nEndcaps];
//   static const short xyz2HashedIndex[IX_MAX][IY_MAX][nEndcaps];
  
  /** Map of hased index to x,y,z. See hashedIndex/
   */
  CMS_THREAD_SAFE static EcalScDetId hashedIndex2DetId[kSizeForDenseIndexing];
//   static const EcalScDetId hashedIndex2DetId[kSizeForDenseIndexing];
  
  /*The two arrays are thread safe since they are filled safely using std::call_once and
    then only read and never modified.
   */
};


std::ostream& operator<<(std::ostream& s,const EcalScDetId& id);


#endif //EcalDetId_EcalScDetId_h not defined
