#ifndef ECALDETID_EEDETID_H
#define ECALDETID_EEDETID_H

#include <iosfwd>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalScDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <algorithm>

/** \class EEDetId
 *  Crystal/cell identifier class for the ECAL endcap
 *
 *
 */
class EEDetId : public DetId {
public:
  enum {
    /** Sudetector type. Here it is ECAL endcap.
     */
    Subdet=EcalEndcap
  };
  
  /** Constructor of a null id
   */
  constexpr EEDetId() {}
  
  /** Constructor from a raw value
   * @param rawid det ID number
   */
  constexpr EEDetId(uint32_t rawid) : DetId(rawid) {}
  
  /** Constructor from crystal ix,iy,iz (iz=+1/-1) (mode = XYMODE)
   * or from sc,cr,iz (mode = SCCRYSTALMODE).
   * <p>ix runs from 1 to 100 along x-axis of standard CMS coordinates<br>
   * iy runs from 1 to 100 along y-axis of standard CMS coordinates<br>
   * iz is -1 for EE- and +1 for EE+<br>
   * <p>For isc see isc(), for ic see ic()
   * @see isc(), ic()
   * @param i ix or isc index
   * @param j iy or isc index
   * @param iz iz/zside index: -1 for EE-, +1 for EE+
   * @param mode pass XYMODE if i j refer to ix, iy, SCCRYSTALMODE if thery refer to isc, ic
   */
  // fast  
  constexpr EEDetId(int crystal_ix, int crystal_iy, int iz) : 
    DetId( Ecal, EcalEndcap ) {
    id_|=(crystal_iy&0x7f)|((crystal_ix&0x7f)<<7)|((iz>0)?(0x4000):(0));
  }
  // slow
  constexpr EEDetId(int index1, int index2, int iz, int mode) {
    int crystal_ix=0;
    int crystal_iy=0;
    if (mode == XYMODE) 
    {
      crystal_ix = index1;
      crystal_iy = index2;  
    } 
    else if (mode == SCCRYSTALMODE) 
    {
      int SC = index1;
      int crystal = index2;
      //      std::cout << "iz " << iz << " SC " << index1 << "crystal " << index2  << std::endl;
      
      crystal_ix=iz*ix(SC,crystal);
      if (crystal_ix<0)
        crystal_ix++;
      crystal_ix+=50;
      crystal_iy=iy(SC,crystal);
      if (crystal_iy<0)
        crystal_iy++;
      crystal_iy+=50;
      
    } 
    else 
    {
      throw cms::Exception("InvalidDetId") << "EEDetId:  Cannot create object.  Unknown mode for (int, int, int) constructor.";
    }
    
    if (!validDetId(crystal_ix,crystal_iy,iz))
    {
      throw cms::Exception("InvalidDetId") << "EEDetId:  Cannot create object.  Indexes out of bounds \n"
      << "x = " << crystal_ix << " y = " << crystal_iy << " z = " << iz;
    }
    
    id_|=(crystal_iy&0x7f)|((crystal_ix&0x7f)<<7)|((iz>0)?(0x4000):(0));    
  }
  
  /** Constructor from a generic cell id
   * @param id source detid
   */
  constexpr EEDetId(const DetId& id) : DetId(id){}
  
  /** Assignment operator
   * @param id source det id
   */ 
  constexpr EEDetId& operator=(const DetId& id) {id_ = id.rawId(); return *this;}
  
  /** Gets the subdetector
   * @return subdetectot ID, that is EcalEndcap
   */
  constexpr EcalSubdetector subdet() { return EcalEndcap;}
  
  /** Gets the z-side of the crystal (1/-1)
   * @return -1 for EE-, +1 for EE+
   */
  constexpr int zside() const { return (id_&0x4000)?(1):(-1); }
  
  /** Gets the crystal x-index.
   * @see EEDetId(int, int, int, int) for x-index definition
   * @return x-index
   */
  constexpr int ix() const { return (id_>>7)&0x7F; }
  
  /** Get the crystal y-index
   * @see EEDetId(int, int, int, int) for y-index definition.
   * @return y-index
   */
  constexpr int iy() const { return id_&0x7F; }
  
  /** Gets the DetId of the supercrystal the crystal belong to.
   * @return the supercrystal det id
   * @throw cms::Exception if the crystal det id is invalid 
   */ 
  constexpr EcalScDetId sc() const {
    const int scEdge = 5;
    return EcalScDetId(1+(ix()-1)/scEdge, 1+(iy()-1)/scEdge, zside());
  }
  
  /** Gets the SuperCrystal number within the endcap. This number runs from 1 to 316,
   * numbers 70 149 228 307 are not used.
   *
   * BEWARE: This number is not consistent with indices used in constructor:  see details below.
   *
   * Numbering in quadrant 1 of EE+ is the following
   * \verbatim 
   *  08 17 27        
   *  07 16 26 36 45 54     
   *  06 15 25 35 44 53 62    
   *  05 14 24 34 43 52 61 69   
   *  04 13 23 33 42 51 60 68 76  
   *  03 12 22 32 41 50 59 67 75  
   *  02 11 21 31 40 49 58 66 74  
   *  01 10 20 30 39 48 57 65 73 79 
   *     09 19 29 38 47 56 64 72 78 
   *        18 28 37 46 55 63 71 77
   *  
   *        == THERE IS NO INDEX 70! ==
   * \endverbatim
   *
   * Quadrant 2 indices are deduced by a symetry about y-axis and by adding an offset
   * of 79.<br>
   * Quadrant 3 and 4 indices are deduced from quadrant 1 and 2 by a symetry
   * about x-axis and adding an offset. Quadrant N starts with index 1 + (N-1)*79.
   *
   * <p>EE- indices are deduced from EE+ by a symetry about (x,y)-plane (mirrored view). <b>It is
   * inconsistent with indices used in constructor EEDetId(int, int,int) in
   * SCCRYSTALMODE</b>. Indices of constructor uses a symetry along y-axis: in principal it
   * considers the isc as a local index. The discrepancy is most probably due to a bug in the
   * implementation of this isc() method.
   */
  int isc() const {
    return isc( 1 + ( ix() - 1 )/nCrys,
                1 + ( iy() - 1 )/nCrys ) ; 
  }
  
  /** Gets crystal number inside SuperCrystal.
   * Crystal numbering withing a supercrystal in each quadrant:
   * \verbatim
   *                       A y
   *  (Q2)                 |                    (Q1)
   *       25 20 15 10 5   |     5 10 15 20 25
   *       24 19 14  9 4   |     4  9 14 19 24
   *       23 18 13  8 3   |     3  8 13 18 23
   *       22 17 12  7 2   |     2  7 12 17 22
   *       21 16 11  6 1   |     1  6 11 16 21
   *                       |
   * ----------------------o---------------------------> x
   *                       |
   *       21 16 11  6 1   |     1  6 11 16 21
   *       22 17 12  7 2   |     2  7 12 17 22
   *       23 18 13  8 3   |     3  8 13 18 23
   *       24 19 14  9 4   |     4  9 14 19 24
   *       25 20 15 10 5   |     5 10 15 20 25
   *  (Q3)                                       (Q4)
   * \endverbatim
   *
   * @return crystal number from 1 to 25
   */
  constexpr int ic() const {
    /*
     *  Return crystal number from (x,y) coordinates.
     *
     *  Author    : B W Kennedy
     *  Version   : 1.00
     *  Created   : 5 May 2006
     *  Last Mod  :
     *
     *  Input     : ix, iy - (x,y) position of crystal
     */
    
    /*  Useful constants . */
    int iQuadrant = iquadrant();
    int icrCol=-1;
    int icrRow=-1;
    
    if (iQuadrant == 1 || iQuadrant == 3)
    {
      icrCol=(ixQuadrantOne()-1) % nCrys;
      icrRow=(iyQuadrantOne()-1) % nCrys;
    }
    
    else if (iQuadrant == 2 || iQuadrant == 4)
    {
      icrRow=(ixQuadrantOne()-1) % nCrys;
      icrCol=(iyQuadrantOne()-1) % nCrys;
    } 
    
    int icrys = 5*icrCol + icrRow + 1;
    
    return icrys;
  }
  
  /** Gets the quadrant of the DetId.
   * Quadrant number definition, x and y in std CMS coordinates, for EE+:
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
   *
   * @return quadrant number
   */
  constexpr int iquadrant() const {
    if (ix()>50)
    {
      if(iy()>50)
        return 1;
      else
        return 4;
    }
    else
    {
      if(iy()>50)
        return 2;
      else
        return 3;
    }
    //Should never be reached
    return -1;
  }
  
  /** Checks if crystal is in EE+
   * @return true for EE+, false for EE-
   */
  constexpr bool positiveZ() const { return id_&0x4000;}
  
  constexpr int iPhiOuterRing() const {
    int returnValue ( 0 ) ;
    if( isOuterRing() )
    {
      const int ax ( std::abs( ix() - IX_MAX/2 ) ) ;
      const int ay ( std::abs( iy() - IY_MAX/2 ) ) ;
      returnValue = ax + 50 - ay ;
      if( ay <= 47 ) --returnValue ;
      if( ay <= 45 ) --returnValue ;
      if( ay <= 42 ) --returnValue ;
      if( ay <= 37 ) --returnValue ;
      if( ay <= 35 ) --returnValue ;
      if( ay <= 30 ) --returnValue ;
      if( ay <= 25 ) --returnValue ;
      if( ay <= 15 ) --returnValue ;
      if( ay <= 10 ) --returnValue ;
      const int iq ( iquadrant() ) ;
      if( 1==iq )
      {
        returnValue = 91 - returnValue ;
      }
      else
      {
        if( 2==iq )
        {
          returnValue += 90 ;
        }
        else
        {
          if( 3==iq )
          {
            returnValue = 271 - returnValue ;
          }
          else
          {
            returnValue += 270 ;
          }
        }
      }
      returnValue = 1 + ( 360 + returnValue - 10 -1 )%360 ;
    }
    //   if( positiveZ() ) returnValue += 360 ;
    return returnValue ;
  } // 1-360 else==0 if not on outer ring!
  
  constexpr static EEDetId idOuterRing( int iPhi , int zEnd ) {
    iPhi -= 10 ; // phi=1 in barrel is at -10deg
    while( iPhi <   1 ) iPhi+=360 ;
    while( iPhi > 360 ) iPhi-=360 ;
    
    const int index1 ( iPhi - 1 ) ;
    const int quad   ( index1/90 ) ;
    int       indexq (  index1 - quad*90 + 1 ) ;
    if( 0==quad || 2==quad ) indexq = 91 - indexq ;
    const int indexh ( indexq > 45 ? 91 - indexq : indexq ) ;
    const int axh    ( indexh<=10 ? indexh :
    ( indexh<=12 ? 10 :
    ( indexh<=17 ? indexh - 2 :
    ( indexh<=18 ? 15 :
    ( indexh<=28 ? indexh - 3 :
    ( indexh<=30 ? 25 :
    ( indexh<=35 ? indexh - 5 :
    ( indexh<=39 ? 30 :
    ( indexh<=44 ? indexh - 9 : 35 ))))))))) ;
    const int ayh    ( indexh<=10 ? 50 :
    ( indexh<=12 ? 60 - indexh :
    ( indexh<=17 ? 47 :
    ( indexh<=18 ? 64 - indexh : 
    ( indexh<=28 ? 45 :
    ( indexh<=30 ? 73 - indexh :
    ( indexh<=35 ? 42 :
    ( indexh<=39 ? 77 - indexh :
    ( indexh<=44 ? 37 : 36 ))))))))) ;
    const int bxh ( indexq>45 ? ayh : axh ) ;
    const int byh ( indexq>45 ? axh : ayh ) ;
    const int cx  ( ( quad==0 || quad==3 ? bxh : -bxh+1 ) + IX_MAX/2 ) ;
    const int cy  ( ( quad==0 || quad==1 ? byh : -byh+1 ) + IY_MAX/2 ) ;
    
    return EEDetId( cx, cy, ( zEnd > 0 ? 1 : -1 ) ) ;
  }
  
  /** Gets a compact index for arrays
   * @return compact index from 0 to kSizeForDenseIndexing-1
   */
  int hashedIndex() const 
  {
    const uint32_t jx ( ix() ) ;
    const uint32_t jd ( 2*( iy() - 1 ) + ( jx - 1 )/50 ) ;
    return (  ( positiveZ() ? kEEhalf : 0) + kdi[jd] + jx - kxf[jd] ) ;
  }
  
  /** Same as hashedIndex()
   * @return compact index from 0 to kSizeForDenseIndexing-1
   */
  uint32_t denseIndex() const { return hashedIndex() ; }
  
  /** returns a new EEDetId offset by nrStepsX and nrStepsY (can be negative),
   * returns EEDetId(0) if invalid */
  constexpr EEDetId offsetBy( int nrStepsX, int nrStepsY ) const {
    int newX = ix() + nrStepsX;
    int newY = iy() + nrStepsY;
    
    if( validDetId( newX, newY, zside() ) ) {
      return EEDetId( newX, newY, zside() );
    } else {
      return EEDetId(0);
    }
  }
  
      /** returns a new EEDetId swapped (same iX, iY) to the other endcap, 
       * returns EEDetId(0) if invalid (shouldnt happen) */
  constexpr EEDetId switchZSide() const {
    int newZSide = -1 * zside();
    if( validDetId(ix(), iy(), newZSide ) ) {
      return EEDetId( ix(), iy(), newZSide );
    } else {
      return EEDetId(0);
    }
  }
  
  /** following are static member functions of the above two functions
   *  which take and return a DetId, returns DetId(0) if invalid 
   */
  constexpr static DetId offsetBy( const DetId startId, int nrStepsX, int nrStepsY ) {
    if( startId.det() == DetId::Ecal && startId.subdetId() == EcalEndcap ) {
      EEDetId eeStartId( startId );
      return eeStartId.offsetBy( nrStepsX, nrStepsY ).rawId();
    } else {
      return DetId(0);
    }
  }
  constexpr static DetId switchZSide( const DetId startId ) {
    if( startId.det() == DetId::Ecal && startId.subdetId() == EcalEndcap ) {
      EEDetId eeStartId(startId);
      return eeStartId.switchZSide().rawId();
    } else {
      return DetId(0);
    }
  }
  
  /** Checks validity of a dense/hashed index
   * @param din dense/hashed index as returned by hashedIndex() or denseIndex()
   * method
   * @return true if index is valid, false otherwise
   */
  constexpr static bool validDenseIndex( uint32_t din ) { return validHashIndex( din ) ; }
  
  /** Converts a hashed/dense index as defined in hashedIndex() and denseIndex()
   * methods to a det id.
   * @param din hashed/dense index
   * @return det id
   */
  static EEDetId detIdFromDenseIndex( uint32_t din ) { return unhashIndex( din ) ; }
  
  constexpr static bool isNextToBoundary(     EEDetId id ) {
    return isNextToDBoundary( id ) || isNextToRingBoundary( id ) ;
  }
  
  constexpr static bool isNextToDBoundary(    EEDetId id ) {
    // hardcoded values for D boundary
    return id.ix() == 50 || id.ix() == 51 ;
  }
  
  constexpr static bool isNextToRingBoundary( EEDetId id ) {
    for (int i = -1; i <= 1; ++i) {
      for (int j = -1; j <= 1; ++j) {
        if ( ! validDetId( id.ix() + i, id.iy() + j, id.zside() ) ) {
          return true;
        }
      }
    }
    return false;
  }
  
  /** Gets a DetId from a compact index for arrays. Converse of hashedIndex() method.
   * @param hi dense/hashed index
   * @return det id
   */
  static EEDetId unhashIndex( int hi ) {
    if( validHashIndex( hi ) )
    {
      const int iz ( hi<kEEhalf ? -1 : 1 ) ;
      const uint32_t di ( hi%kEEhalf ) ;
      const int ii ( ( std::upper_bound( kdi, kdi+(2*IY_MAX), di ) - kdi ) - 1 ) ;
      const int iy ( 1 + ii/2 ) ;
      const int ix ( kxf[ii] + di - kdi[ii] ) ;
      return EEDetId( ix, iy, iz ) ;
    }
    else
    {
      return EEDetId() ;
    }
  }
  
  /** Checks if a hashed/dense index is valid
   * @see hashedIndex(), denseIndex()
   * @param i hashed/dense index
   * @return true if the index is valid, false otherwise
   */
  constexpr static bool validHashIndex( int i ) { return ( i < kSizeForDenseIndexing ) ; }
  
  /** Checks validity of a crystal (x,y.z) index triplet.
   * @param crystal_ix crystal x-index
   * @param crystal_iy crystal y-index
   * @param iz crystal z-index
   * @see EEDetId(int, int, int, int) for index definition
   * @return true if valid, false otherwise
   */
  constexpr static bool validDetId(int crystal_ix, int crystal_iy, int iz) {
    return 
      crystal_ix >= IX_MIN && crystal_ix <= IX_MAX &&
      crystal_iy >= IY_MIN && crystal_iy <= IY_MAX &&  
      std::abs(iz)==1 && 
      ( fastValidDetId(crystal_ix,crystal_iy) ||
	slowValidDetId(crystal_ix,crystal_iy) );
  }
  constexpr static bool slowValidDetId(int crystal_ix, int crystal_iy) {
    return // negative logic!
    !( 
    (crystal_ix >= 1 && crystal_ix <= 3 && (crystal_iy <= 40 || crystal_iy > 60) ) ||
    (crystal_ix >= 4 && crystal_ix <= 5 && (crystal_iy <= 35 || crystal_iy > 65) ) || 
    (crystal_ix >= 6 && crystal_ix <= 8 && (crystal_iy <= 25 || crystal_iy > 75) ) || 
    (crystal_ix >= 9 && crystal_ix <= 13 && (crystal_iy <= 20 || crystal_iy > 80) ) || 
    (crystal_ix >= 14 && crystal_ix <= 15 && (crystal_iy <= 15 || crystal_iy > 85) ) || 
    (crystal_ix >= 16 && crystal_ix <= 20 && (crystal_iy <= 13 || crystal_iy > 87) ) || 
    (crystal_ix >= 21 && crystal_ix <= 25 && (crystal_iy <= 8 || crystal_iy > 92) ) || 
    (crystal_ix >= 26 && crystal_ix <= 35 && (crystal_iy <= 5 || crystal_iy > 95) ) || 
    (crystal_ix >= 36 && crystal_ix <= 39 && (crystal_iy <= 3 || crystal_iy > 97) ) || 
    (crystal_ix >= 98 && crystal_ix <= 100 && (crystal_iy <= 40 || crystal_iy > 60) ) ||
    (crystal_ix >= 96 && crystal_ix <= 97 && (crystal_iy <= 35 || crystal_iy > 65) ) || 
    (crystal_ix >= 93 && crystal_ix <= 95 && (crystal_iy <= 25 || crystal_iy > 75) ) || 
    (crystal_ix >= 88 && crystal_ix <= 92 && (crystal_iy <= 20 || crystal_iy > 80) ) || 
    (crystal_ix >= 86 && crystal_ix <= 87 && (crystal_iy <= 15 || crystal_iy > 85) ) || 
    (crystal_ix >= 81 && crystal_ix <= 85 && (crystal_iy <= 13 || crystal_iy > 87) ) || 
    (crystal_ix >= 76 && crystal_ix <= 80 && (crystal_iy <= 8 || crystal_iy > 92) ) || 
    (crystal_ix >= 66 && crystal_ix <= 75 && (crystal_iy <= 5 || crystal_iy > 95) ) || 
    (crystal_ix >= 62 && crystal_ix <= 65 && (crystal_iy <= 3 || crystal_iy > 97) ) ||
    ( (crystal_ix == 40 || crystal_ix == 61) && ( (crystal_iy >= 46 && crystal_iy <= 55 ) || crystal_iy <= 3 || crystal_iy > 97 )) ||
    ( (crystal_ix == 41 || crystal_ix == 60) && crystal_iy >= 44 && crystal_iy <= 57 ) ||
    ( (crystal_ix == 42 || crystal_ix == 59) && crystal_iy >= 43 && crystal_iy <= 58 ) ||
    ( (crystal_ix == 43 || crystal_ix == 58) && crystal_iy >= 42 && crystal_iy <= 59 ) ||
    ( (crystal_ix == 44 || crystal_ix == 45 || crystal_ix == 57 || crystal_ix == 56) && crystal_iy >= 41 && crystal_iy <= 60 ) ||
    ( crystal_ix >= 46 && crystal_ix <= 55 && crystal_iy >= 40 && crystal_iy <= 61 ) 
    );
  }

  /**  check if ix and iy is in a "ring" inscribed in EE
   *   if is inside is valid for sure
   *   if not the slow version shall be called
   */
  constexpr static bool fastValidDetId(int crystal_ix, int crystal_iy) {
    float x =  crystal_ix; float y =  crystal_iy;
    float r = (x - 50.5f) * (x - 50.5f) + (y - 50.5f) * (y - 50.5f);
    return r > 12.f * 12.f && r < 48.f * 48.f;
  }

  /** Returns the distance along x-axis in crystal units between two EEDetId
   * @param a det id of first crystal
   * @param b det id of second crystal
   * @return distance
   */
  constexpr static int distanceX(const EEDetId& a,const EEDetId& b){
//     return std::abs(a.ix()-b.ix());
    return ((a.ix() >b.ix()) ? (a.ix()-b.ix()) : (b.ix()-a.ix()));
  }
  
  /** Returns the distance along y-axis in crystal units between two EEDetId
   * @param a det id of first crystal
   * @param b det id of second crystal
   * @return distance
   */
  constexpr static int distanceY(const EEDetId& a,const EEDetId& b) {
//     return std::abs(a.iy() - b.iy()); 
    return ((a.iy() >b.iy()) ? (a.iy()-b.iy()) : (b.iy()-a.iy()));
  }
  
  
  /** Gives supercrystal index from endcap *supercrystal* x and y indexes.
   * @see isc() for the index definition
   * @param iscCol supercrystal column number: supecrystal x-index for EE+
   * @param iscRow: supecrystal y-index
   * @return supercystal index
   */
  static int isc( int iscCol,   // output is 1-316
		  int iscRow ) 
  {
    if( 0  < iscCol &&
      21 > iscCol &&
      0  < iscRow &&
      21 > iscRow    )
    {
      const int iquad (  ( 10<iscCol && 10<iscRow ? 1 :
      ( 11>iscCol && 10<iscRow ? 2 :
      ( 11>iscCol && 11>iscRow ? 3 : 4 ) ) ) ) ;
      
      const int iCol = ( 1 == iquad || 4 == iquad ? iscCol - 10 : 11 - iscCol ) ;
      const int iRow = ( 1 == iquad || 2 == iquad ? iscRow - 10 : 11 - iscRow ) ;
      
      const int nSCinQuadrant = ISC_MAX/4;
      
      const int yOff ( iYoffset[iCol] ) ;
      
      const int qOff ( nSCinQuadrant*( iquad - 1 ) ) ;
      
      const int iscOne ( QuadColLimits[iCol-1] + iRow - yOff ) ;
      
      return ( yOff                >= iRow   ? -1 : 
      ( QuadColLimits[iCol] <  iscOne ? -2 :
      iscOne + qOff ) ) ;
    }
    else
    {
      return -3 ; // bad inputs
    }
  }
  
  /** Lower bound of EE crystal x-index
   */
  static const int IX_MIN =1;
  
  /** Lower bound of EE crystal y-index
   */
  static const int IY_MIN =1;
  
  /** Upper bound of EE crystal y-index
   */
  static const int IX_MAX =100;
  
  /** Upper bound of EE crystal y-index
   */
  static const int IY_MAX =100;
  
  /** Lower bound of supercystal index as defined in isc()
   */
  static const int ISC_MIN=1;
  
  /** Lower bound of crystal index within a supercrystal
   */
  static const int ICR_MIN=1;
  
  /** Upper bound of supercystal index defined in isc()
   * <p>Beware it differs from the number of supercrystals in one endcap,
   * which is 312, because the numbering is not dense.
   */
  static const int ISC_MAX=316;
  
  /** Upper bound of crystal index within a supercrystal
   */
  static const int ICR_MAX=25;
  
  enum {
    /** Number of crystals per Dee
     */
    kEEhalf = 7324 ,
    /** Number of dense crystal indices, that is number of
     * crystals per endcap.
     */
    kSizeForDenseIndexing = 2*kEEhalf
  };
  
  /*@{*/
  /** function modes for EEDetId(int, int, int, int) constructor
   */
  static const int XYMODE        = 0;
  static const int SCCRYSTALMODE = 1;
  /*@}*/
  
private:
  
  constexpr bool        isOuterRing() const 
  {
    const int kx ( ix() ) ;
    const int ky ( iy() ) ;
    const int ax ( kx>IX_MAX/2 ? kx-IX_MAX/2 : IX_MAX/2 + 1 - kx ) ;
    const int ay ( ky>IY_MAX/2 ? ky-IY_MAX/2 : IY_MAX/2 + 1 - ky ) ;
    return ( isOuterRingXY( ax, ay ) ||
    isOuterRingXY( ay, ax )    ) ;
  }
  
  constexpr static bool isOuterRingXY( int ax, int ay ) 
  {
    return ( ( ax<=10 &&           ay==50 ) ||
    ( ax==10 &&           ay>=48 ) ||
    ( ax<=15 && ax>=11 && ay==47 ) ||
    ( ax==15 &&           ay==46 ) ||
    ( ax<=25 && ax>=16 && ay==45 ) ||
    ( ax==25 &&           ay<=44 && ay>=43 ) ||
    ( ax<=30 && ax>=26 && ay==42 ) ||
    ( ax==30 &&           ay<=41 && ay>=38 ) ||
    ( ax<=35 && ax>=31 && ay==37 ) ||
    ( ax==35 &&           ay==36 )              ) ;
  }
  
  
  //Functions from B. Kennedy to retrieve ix and iy from SC and Crystal number
  
  static const int nCols = 10;
  static const int nCrys = 5; /* Number of crystals per row in SC */
  static const int QuadColLimits[nCols+1];
  static const int iYoffset[nCols+1];
  
  static const unsigned short kxf[2*IY_MAX] ;
  static const unsigned short kdi[2*IY_MAX] ;
  
  int ix( int iSC, int iCrys ) const
  {
    /*
     *  ix() return individual crystal x-coordinate
     *
     *  Author    : B W Kennedy
     *  Version   : 1.00
     *  Created   : 21 December 2005
     *  Last Mod  : 31 January 2006
     *
     *  Input     : iSC, iCrys - Supercrystal and crystal ids
     */
    
    
    int nSCinQuadrant = QuadColLimits[nCols];
    
    if (iSC > 4*nSCinQuadrant || iSC < 1) 
    {
      throw new std::exception();
    }
    
    //  Map SC number into (x>0,y>0) quadrant.
    int iSCmap = 0;
    int iqx = 0;
    int iq = 0;
    if (iSC > 3*nSCinQuadrant) 
    {
      iSCmap = iSC - 3*nSCinQuadrant;
      iqx =  1;
      iq=4;
    } 
    else if (iSC > 2*nSCinQuadrant) 
    {
      iSCmap = iSC - 2*nSCinQuadrant;
      iqx = -1;
      iq=3;
    } 
    else if (iSC > nSCinQuadrant) 
    {
      iSCmap = iSC - nSCinQuadrant;
      iqx = -1;
      iq=2;
    } 
    else 
    {
      iSCmap = iSC;
      iqx = 1;
      iq=1;
    }
    
    // Decide which column the SC is in
    int iCol = 0 ;
    while (iSCmap > QuadColLimits[iCol++]) ;
    iCol-- ;
    
    int ixCrys=-1;
    if (iq == 1 || iq == 3) 
      ixCrys = iqx*(5*(iCol-1) + (int)(iCrys+4)/5);
    else   if (iq == 2 || iq == 4) 
      ixCrys = iqx*(5*(iCol-1) + (iCrys-1)%5 + 1);
    
    // returning a value from 1 to 100  
    
    return ixCrys;
  }
  
  int iy( int iSC, int iCrys ) const
  {
    /*
     *  iy() return individual crystal y-coordinate
     *
     *  Author    : B W Kennedy
     *  Version   : 1.00
     *  Created   : 21 December 2005
     *  Last Mod  : 31 January 2006
     *
     *  Input     : iSC, iCrys - Supercrystal and crystal ids
     */
    
    int nSCinQuadrant = QuadColLimits[nCols];
    if (iSC > 4*nSCinQuadrant || iSC < 1) 
    {
      throw new std::exception();
    }
    
    //  Map SC number into (x>0,y>0) quadrant
    int iSCmap, iqy,iq;
    if (iSC > 3*nSCinQuadrant) 
    {
      iSCmap = iSC - 3*nSCinQuadrant;
      iqy = -1;
      iq=4;
    } 
    else if (iSC > 2*nSCinQuadrant) 
    {
      iSCmap = iSC - 2*nSCinQuadrant;
      iqy = -1;
      iq=3;
    } 
    else if (iSC > nSCinQuadrant) 
    {
      iSCmap = iSC - nSCinQuadrant;
      iqy = 1;
      iq=2;
    } else 
    {
      iSCmap = iSC;
      iqy = 1;
      iq=1;
    }
    
    // Decide which column the SC is in
    int iCol = 0;
    while (iSCmap > QuadColLimits[iCol++]) ;
    iCol--;
    
    int iSCy = iSCmap - QuadColLimits[iCol-1] + iYoffset[iCol];
    
    int iyCrys=-1;
    if (iq == 1 || iq == 3)
      iyCrys = iqy*(5*(iSCy-1) + (iCrys-1)%5 + 1);
    else if (iq == 2 || iq == 4)
      iyCrys = iqy*(5*(iSCy-1) + (int)(iCrys+4)/5 );
    return iyCrys;
  }
  
  constexpr int ixQuadrantOne() const
  { 
    int iQuadrant = iquadrant();
    if ( iQuadrant == 1 || iQuadrant == 4)
      return (ix() - 50);
    else if ( iQuadrant == 2 || iQuadrant == 3)
      return (51 - ix());
    //Should never be reached
    return -1;
  }
  
  constexpr int iyQuadrantOne() const
  { 
    int iQuadrant = iquadrant();
    if ( iQuadrant == 1 || iQuadrant == 2)
      return (iy() - 50);
    else if ( iQuadrant == 3 || iQuadrant == 4)
      return 51 - iy();
    //Should never be reached
    return -1;
  }

  
};


std::ostream& operator<<(std::ostream& s,const EEDetId& id);

#endif
