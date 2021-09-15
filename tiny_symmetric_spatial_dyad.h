/*
 * Copyright 2020 Google LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TINY_SYMMETRIC_SPATIAL_DYAD_H
#define TINY_SYMMETRIC_SPATIAL_DYAD_H

#include "tiny_spatial_transform.h"

template <typename TinyScalar, typename TinyConstants>
struct TinySymmetricSpatialDyad {
  typedef ::TinyMatrix3x3<TinyScalar, TinyConstants> TinyMatrix3x3;
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinySpatialMotionVector<TinyScalar, TinyConstants>
      TinySpatialMotionVector;
  typedef ::TinySpatialTransform<TinyScalar, TinyConstants>
      TinySpatialTransform;

  TinyMatrix3x3 m_topLeftMat, m_topRightMat, m_bottomLeftMat, m_bottomRightMat;

  TinyVector3 m_center_of_mass;
  //
  TinySymmetricSpatialDyad() { setIdentity(); }
  TinySymmetricSpatialDyad(const TinyMatrix3x3& topLeftMat,
                           const TinyMatrix3x3& topRightMat,
                           const TinyMatrix3x3& bottomLeftMat,
                           const TinyMatrix3x3& bottomRightMat) {
    m_topLeftMat = topLeftMat;
    m_topRightMat = topRightMat;
    m_bottomLeftMat = bottomLeftMat;
    m_bottomRightMat = bottomRightMat;
  }
  //
  void setIdentity() {
    m_topLeftMat.set_identity();
    m_topRightMat.set_zero();
    m_bottomLeftMat.set_zero();
    m_bottomRightMat.set_identity();
    m_center_of_mass.set_zero();
  }
  void set_zero() {
    m_topLeftMat.set_zero();
    m_topRightMat.set_zero();
    m_bottomLeftMat.set_zero();
    m_bottomRightMat.set_zero();
    m_center_of_mass.set_zero();
  }
  // void setZero() {
  //   m_topLeftMat.set_zero();
  //   m_topRightMat.set_zero();
  //   m_bottomLeftMat.set_zero();
  //   m_bottomRightMat.set_zero();
  //   m_center_of_mass.set_zero();
  // }

  //

  TinySymmetricSpatialDyad operator-() const {
    TinySymmetricSpatialDyad m;
    m.m_topLeftMat = -m_topLeftMat;
    m.m_topRightMat = -m_topRightMat;
    m.m_bottomLeftMat = -m_bottomLeftMat;
    m.m_bottomRightMat = -m_bottomRightMat;
    return m;
  }

  TinySymmetricSpatialDyad operator+(
      const TinySymmetricSpatialDyad& m1) {
    TinySymmetricSpatialDyad m;
    m.m_topLeftMat = m1.m_topLeftMat + m_topLeftMat;
    m.m_topRightMat = m1.m_topRightMat + m_topRightMat;
    m.m_bottomLeftMat = m1.m_bottomLeftMat + m_bottomLeftMat;
    m.m_bottomRightMat = m1.m_bottomRightMat + m_bottomRightMat;
    return m;
  }

  TinySymmetricSpatialDyad operator-(
      const TinySymmetricSpatialDyad& m1) {
    TinySymmetricSpatialDyad m;
    m.m_topLeftMat = m_topLeftMat - m1.m_topLeftMat;
    m.m_topRightMat = m_topRightMat - m1.m_topRightMat;
    m.m_bottomLeftMat = m_bottomLeftMat - m1.m_bottomLeftMat;
    m.m_bottomRightMat = m_bottomRightMat - m1.m_bottomRightMat;
    return m;
  }

  TinySymmetricSpatialDyad& operator-=(const TinySymmetricSpatialDyad& mat) {
    m_topLeftMat -= mat.m_topLeftMat;
    m_topRightMat -= mat.m_topRightMat;
    m_bottomLeftMat -= mat.m_bottomLeftMat;
    m_bottomRightMat -= mat.m_bottomRightMat;
    return *this;
  }

  TinySymmetricSpatialDyad& operator+=(const TinySymmetricSpatialDyad& mat) {
    m_topLeftMat += mat.m_topLeftMat;
    m_topRightMat += mat.m_topRightMat;
    m_bottomLeftMat += mat.m_bottomLeftMat;
    m_bottomRightMat += mat.m_bottomRightMat;
    return *this;
  }

  static TinySymmetricSpatialDyad computeInertiaDyad(
      TinyScalar mass, const TinyVector3& com, const TinyMatrix3x3& inertia_C, TinyScalar armature) {
    TinySymmetricSpatialDyad result;
    TinyVector3 h = com * mass;
    TinyScalar o = TinyConstants::zero();
    TinyMatrix3x3 I = inertia_C + TinyVectorCrossMatrix(com) *
                                      TinyVectorCrossMatrix(com).transpose() *
                                      mass;
    TinySymmetricSpatialDyad& mat = result;
    for (int r = 0; r < 3; r++) {
      for (int c = 0; c < 3; c++) {
        mat(r, c) = I(r, c);
      }
    }

    {
      mat(3, 0) = o;
      mat(3, 1) = h[2];
      mat(3, 2) = -h[1];
      mat(4, 0) = -h[2];
      mat(4, 1) = o;
      mat(4, 2) = h[0];
      mat(5, 0) = h[1];
      mat(5, 1) = -h[0];
      mat(5, 2) = o;
    }
    {
      mat(0, 3) = o;
      mat(0, 4) = -h[2];
      mat(0, 5) = h[1];
      mat(1, 3) = h[2];
      mat(1, 4) = o;
      mat(1, 5) = -h[0];
      mat(2, 3) = -h[1];
      mat(2, 4) = h[0];
      mat(2, 5) = o;
    }
    {
      mat(3, 3) = mass;
      mat(3, 4) = o;
      mat(3, 5) = o;
      mat(4, 3) = o;
      mat(4, 4) = mass;
      mat(4, 5) = o;
      mat(5, 3) = o;
      mat(5, 4) = o;
      mat(5, 5) = mass;
    }
    for (int i = 0; i < 6; ++i)
      mat(i, i) = mat(i, i) + armature;


    result.m_center_of_mass = com;

    return result;
  }
  TinyScalar operator()(int r, int c) const {
    TinyConstants::FullAssert(r >= 0);
    TinyConstants::FullAssert(c >= 0);
    TinyConstants::FullAssert(r < 6);
    TinyConstants::FullAssert(c < 6);

    if (r < 3) {
      if (c < 3) {
        return m_topLeftMat(r, c);
      } else {
        return m_topRightMat(r, c - 3);
      }
    } else {
      if (c < 3) {
        return m_bottomLeftMat(r - 3, c);
      } else {
        return m_bottomRightMat(r - 3, c - 3);
      }
    }

    return TinyConstants::zero();
  }
  TinyScalar& operator()(int r, int c) {
    TinyConstants::FullAssert(r >= 0);
    TinyConstants::FullAssert(c >= 0);
    TinyConstants::FullAssert(r < 6);
    TinyConstants::FullAssert(c < 6);

    if (r < 3) {
      if (c < 3) {
        return m_topLeftMat(r, c);
      } else {
        return m_topRightMat(r, c - 3);
      }
    } else {
      if (c < 3) {
        return m_bottomLeftMat(r - 3, c);
      } else {
        return m_bottomRightMat(r - 3, c - 3);
      }
    }

    return m_bottomRightMat(0, 0);
  }

  void print(const char* txt) const {
    printf("%s\n", txt);
    for (int r = 0; r < 6; r++) {
      for (int c = 0; c < 6; c++) {
        TinyScalar val = (*this)(r, c);
        double v = TinyConstants::getDouble(val);
        printf("%f, ", v);
      }
      printf("\n");
    }
  }

  TinySpatialMotionVector mul_org(const TinySpatialMotionVector& vec) const {
    return TinySpatialMotionVector(
        m_bottomLeftMat * vec.m_topVec +
            m_topLeftMat.transpose() * vec.m_bottomVec,
        m_topLeftMat * vec.m_topVec + m_topRightMat * vec.m_bottomVec);
  }
  
  void adj_sd_mul_org(const TinySpatialMotionVector& R, const TinySpatialMotionVector& a,
                      TinySymmetricSpatialDyad& adj_m, TinySpatialMotionVector& adj_a) const {
    const TinyMatrix3x3& m11 = m_topLeftMat;
    const TinyMatrix3x3& m12 = m_topRightMat;
    const TinyMatrix3x3& m21 = m_bottomLeftMat;
    const TinyMatrix3x3& m22 = m_bottomRightMat;

    const TinyVector3& w = a.m_topVec;
    const TinyVector3& v = a.m_bottomVec;

    const TinyVector3& w0 = R.m_topVec;
    const TinyVector3& v0 = R.m_bottomVec;

    adj_m.m_topLeftMat += TinyMatrix3x3::vvt(v0, w)+TinyMatrix3x3::vvt(v, w0);
    adj_m.m_topRightMat += TinyMatrix3x3::vvt(w0, v);
    adj_m.m_bottomLeftMat += TinyMatrix3x3::vvt(w0, w);

    adj_a.m_topVec += m21.transpose()*w0+m11.transpose()*v0;
    adj_a.m_bottomVec += m11*w0+m12.transpose()*v0;

  }
  
  TinySpatialMotionVector mul_inv(const TinySpatialMotionVector& vec) const {
    return TinySpatialMotionVector(
        m_topLeftMat * vec.m_topVec + m_topRightMat * vec.m_bottomVec,
        m_bottomLeftMat * vec.m_topVec + m_bottomRightMat * vec.m_bottomVec);
  }
  
  void adj_sd_mul_inv(const TinySpatialMotionVector& R, const TinySpatialMotionVector& a,
                      TinySymmetricSpatialDyad& adj_m, TinySpatialMotionVector& adj_a) const {
    const TinyMatrix3x3& m11 = m_topLeftMat;
    const TinyMatrix3x3& m12 = m_topRightMat;
    const TinyMatrix3x3& m21 = m_bottomLeftMat;
    const TinyMatrix3x3& m22 = m_bottomRightMat;

    const TinyVector3& w = a.m_topVec;
    const TinyVector3& v = a.m_bottomVec;

    const TinyVector3& w0 = R.m_topVec;
    const TinyVector3& v0 = R.m_bottomVec;

    adj_m.m_topLeftMat += TinyMatrix3x3::vvt(w0, w);
    adj_m.m_topRightMat += TinyMatrix3x3::vvt(w0, v);
    adj_m.m_bottomLeftMat += TinyMatrix3x3::vvt(v0, w);
    adj_m.m_bottomRightMat += TinyMatrix3x3::vvt(v0, v);
    adj_a.m_topVec += m11.transpose()*w0+m21.transpose()*v0;
    adj_a.m_bottomVec += m12.transpose()*w0+m22.transpose()*v0;

  }


  static TinySymmetricSpatialDyad vTimesvTranspose(
      const TinySpatialMotionVector& vecA,
      const TinySpatialMotionVector& vecB) {
    TinySymmetricSpatialDyad diad;
    for (int i = 0; i < 3; i++) {
      diad.m_topLeftMat[i] = vecA.m_topVec[i] * vecB.m_topVec;
      diad.m_bottomLeftMat[i] = vecA.m_bottomVec[i] * vecB.m_topVec;
      diad.m_topRightMat[i] = vecA.m_topVec[i] * vecB.m_bottomVec;
      diad.m_bottomRightMat[i] = vecA.m_bottomVec[i] * vecB.m_bottomVec;
    }

    return diad;
  }

  static void adj_sd_vTimesvTranspose(
      const TinySymmetricSpatialDyad& m,
      const TinySpatialMotionVector& vecA,
      const TinySpatialMotionVector& vecB,
      TinySpatialMotionVector& adj_A,
      TinySpatialMotionVector& adj_B) {
    adj_A.m_topVec += m.m_topLeftMat*vecB.m_topVec 
                    + m.m_topRightMat*vecB.m_bottomVec;
    adj_A.m_bottomVec += m.m_bottomLeftMat*vecB.m_topVec 
                    + m.m_bottomRightMat*vecB.m_bottomVec;
    adj_B.m_topVec += m.m_topLeftMat.transpose() * vecA.m_topVec 
                    + m.m_bottomLeftMat.transpose() * vecA.m_bottomVec;
    adj_B.m_bottomVec += m.m_topRightMat.transpose() * vecA.m_topVec 
                    + m.m_bottomRightMat.transpose() * vecA.m_bottomVec;
  }

  TinySymmetricSpatialDyad transposed() const {
    TinySymmetricSpatialDyad mT;
    mT.m_topLeftMat = this->m_topLeftMat.transpose();
    mT.m_bottomRightMat = this->m_bottomRightMat.transpose();
    mT.m_topRightMat = this->m_bottomLeftMat.transpose();
    mT.m_bottomLeftMat = this->m_topRightMat.transpose();
    return mT;
  }

  static TinySymmetricSpatialDyad mul(const TinySymmetricSpatialDyad& a,
                                      const TinySymmetricSpatialDyad& b) {
    TinySymmetricSpatialDyad res;
    res.m_topLeftMat = (a.m_topLeftMat * b.m_topLeftMat) +
                       (a.m_topRightMat * b.m_bottomLeftMat);
    res.m_topRightMat = (a.m_topLeftMat * b.m_topRightMat) +
                        (a.m_topRightMat * b.m_bottomRightMat);
    res.m_bottomLeftMat = (a.m_bottomLeftMat * b.m_topLeftMat) +
                          (a.m_bottomRightMat * b.m_bottomLeftMat);
    res.m_bottomRightMat = (a.m_bottomLeftMat * b.m_topRightMat) +
                           (a.m_bottomRightMat * b.m_bottomRightMat);
    return res;
  }

  static TinySymmetricSpatialDyad shift(const TinySymmetricSpatialDyad& ia,
                                        const TinySpatialTransform& trans) {
    // toMatrix
    TinyMatrix3x3 rx = TinyVectorCrossMatrix(trans.m_translation);
    // rx.print("rx");
    // TinyMatrix3x3 rxT = rx.transpose();
    // trans.m_rotation.print("E:");
    TinyMatrix3x3 Erx = trans.m_rotation.transpose() * rx;
    // Erx.print("Erx");
    TinySymmetricSpatialDyad m;
    m.m_topLeftMat = trans.m_rotation.transpose();
    m.m_topRightMat.set_zero();
    m.m_bottomLeftMat = -Erx;
    m.m_bottomRightMat = trans.m_rotation.transpose();
    // m.print("m");
    // trans
    TinySymmetricSpatialDyad mT = m.transposed();
    // mT.print("mT");

    TinySymmetricSpatialDyad mTia = TinySymmetricSpatialDyad::mul(mT, ia);
    // mTia.print("mTia");
    TinySymmetricSpatialDyad mTiam = TinySymmetricSpatialDyad::mul(mTia, m);
    // mTiam.print("mTiam");
    return mTiam;
  }


  static void adj_sd_shift(const TinySymmetricSpatialDyad& R, const TinySymmetricSpatialDyad& ia, 
              const TinySpatialTransform& trans, TinySymmetricSpatialDyad& adj_ia, TinySpatialTransform& adj_t)  {
    


    TinyMatrix3x3 rx = TinyVectorCrossMatrix(trans.m_translation);
    TinyMatrix3x3 Erx = trans.m_rotation.transpose() * rx;
    TinySymmetricSpatialDyad m;
    m.m_topLeftMat = trans.m_rotation.transpose();
    m.m_topRightMat.set_zero();
    m.m_bottomLeftMat = -Erx;
    m.m_bottomRightMat = trans.m_rotation.transpose();


    TinySymmetricSpatialDyad adj_a = TinySymmetricSpatialDyad::mul(
        TinySymmetricSpatialDyad::mul(ia.transposed(), m), R);
    adj_a += TinySymmetricSpatialDyad::mul(
        TinySymmetricSpatialDyad::mul(ia, m), R.transposed());    
    TinySymmetricSpatialDyad adj_b = TinySymmetricSpatialDyad::mul(
                TinySymmetricSpatialDyad::mul(m, R), m.transposed());

    adj_ia += adj_b;
    adj_t.m_rotation += adj_a.m_topLeftMat.transpose() 
                  + adj_a.m_bottomRightMat.transpose() 
                  - rx * adj_a.m_bottomLeftMat.transpose();

    adj_t.m_translation += - adj_TinyVectorCrossMatrix(trans.m_rotation 
                          * adj_a.m_bottomLeftMat);
  }


  TinySymmetricSpatialDyad inverseOld() const {
    TinyMatrix3x3 Binv =
        m_bottomRightMat.inverse() * TinyConstants::fraction(-1, 1);
    TinyMatrix3x3 tmp = m_topRightMat * Binv;
    TinyMatrix3x3 temp = tmp * m_bottomLeftMat + m_topLeftMat;
    TinyMatrix3x3 invI_upper_right = temp.inverse();
    tmp = invI_upper_right * m_topRightMat;
    TinyMatrix3x3 invI_upper_left = (tmp * Binv);
    TinyMatrix3x3 invI_lower_right = (invI_upper_left).transpose();
    tmp = m_bottomLeftMat * invI_upper_left;
    tmp[0][0] -= TinyConstants::one();
    tmp[1][1] -= TinyConstants::one();
    tmp[2][2] -= TinyConstants::one();
    TinyMatrix3x3 invI_lower_left = (Binv * tmp);

    return TinySymmetricSpatialDyad(invI_upper_right, invI_upper_left,
                                    invI_lower_right, invI_lower_left);
  }

  // Inverse of a symmetric block matrix
  // according to (4.1) in
  // http://msvlab.hre.ntou.edu.tw/grades/now/inte/Inverse%20&%20Border/border-LuTT.pdf
  TinySymmetricSpatialDyad inverse() const {
    const TinyMatrix3x3& A = m_topLeftMat;
    const TinyMatrix3x3& B = m_topRightMat;
    const TinyMatrix3x3& C = m_bottomLeftMat;  // denoted as B* in paper
    const TinyMatrix3x3& D = m_bottomRightMat;
    TinyMatrix3x3 Ainv = A.inverse();
    TinyMatrix3x3 Dinv = D.inverse();
    TinyMatrix3x3 DCAB = (D - C * Ainv * B).inverse();
    TinySymmetricSpatialDyad result;
    result.m_topLeftMat = Ainv + Ainv * B * DCAB * C * Ainv;
    result.m_topRightMat = -Ainv * B * DCAB;
    result.m_bottomLeftMat = -DCAB * C * Ainv;
    result.m_bottomRightMat = DCAB;
    return result;
  }

  
  // Inverse of a symmetric block matrix
  // according to (4.1) in
  // http://msvlab.hre.ntou.edu.tw/grades/now/inte/Inverse%20&%20Border/border-LuTT.pdf
  void adj_sd_inverse(const TinySymmetricSpatialDyad& R,
                      TinySymmetricSpatialDyad &adj_m ) const {

    const TinyMatrix3x3& A = m_topLeftMat;
    const TinyMatrix3x3& B = m_topRightMat;
    const TinyMatrix3x3& C = m_bottomLeftMat;  // denoted as B* in paper
    const TinyMatrix3x3& D = m_bottomRightMat;
    const TinyMatrix3x3& m11 = R.m_topLeftMat;
    const TinyMatrix3x3& m12 = R.m_topRightMat;
    const TinyMatrix3x3& m21 = R.m_bottomLeftMat;
    const TinyMatrix3x3& m22 = R.m_bottomRightMat;
    TinyMatrix3x3 Ainv = A.inverse();
    TinyMatrix3x3 AinvB = Ainv * B;
    TinyMatrix3x3 CinvA = C * Ainv;
    TinyMatrix3x3 DCAB = (D - CinvA * B).inverse();

    TinyMatrix3x3 RAinv;
    TinyMatrix3x3 RinvAB;
    TinyMatrix3x3 RCinvA;
    TinyMatrix3x3 RDCAB;
    TinyMatrix3x3 oriRDCAB;

    RAinv += m11;
    RinvAB += m11*CinvA.transpose()*DCAB.transpose() - m12*DCAB.transpose();
    RCinvA += DCAB.transpose()*AinvB.transpose()*m11-DCAB.transpose()*m21;
    RDCAB += AinvB.transpose()*m11*CinvA.transpose() - AinvB.transpose()*m12
             - m21*CinvA.transpose() + m22;
    oriRDCAB += -DCAB.transpose() * RDCAB * DCAB.transpose();
    
    RCinvA += -oriRDCAB*B.transpose();

    RAinv += C.transpose()*RCinvA 
            + RinvAB * B.transpose();


    adj_m.m_topLeftMat += -Ainv.transpose() * RAinv * Ainv.transpose();

    adj_m.m_topRightMat += -CinvA.transpose() * oriRDCAB;
    adj_m.m_topRightMat += Ainv.transpose() * RinvAB;


    adj_m.m_bottomLeftMat += RCinvA * Ainv.transpose();
    adj_m.m_bottomRightMat += oriRDCAB;
  }

  /*
  // Inverse of a symmetric block matrix
  // according to (4.1) in
  // http://msvlab.hre.ntou.edu.tw/grades/now/inte/Inverse%20&%20Border/border-LuTT.pdf
  void adj_sd_inverse(const TinySymmetricSpatialDyad& R,
                      TinySymmetricSpatialDyad &adj_m ) const {

    const TinyMatrix3x3& A = m_topLeftMat;
    const TinyMatrix3x3& B = m_topRightMat;
    const TinyMatrix3x3& C = m_bottomLeftMat;  // denoted as B* in paper
    const TinyMatrix3x3& D = m_bottomRightMat;
    const TinyMatrix3x3& m11 = R.m_topLeftMat;
    const TinyMatrix3x3& m12 = R.m_topRightMat;
    const TinyMatrix3x3& m21 = R.m_bottomLeftMat;
    const TinyMatrix3x3& m22 = R.m_bottomRightMat;
    TinyMatrix3x3 Ainv = A.inverse();
    TinyMatrix3x3 AinvB = Ainv * B;
    TinyMatrix3x3 CinvA = C * Ainv;
    TinyMatrix3x3 DCAB = (D - CinvA * B).inverse();

    TinyMatrix3x3 RAinv; //
    TinyMatrix3x3 RinvAB;//
    TinyMatrix3x3 RCinvA;//
    TinyMatrix3x3 RDCAB;
    TinyMatrix3x3 oriRDCAB;//


    RDCAB += AinvB.transpose()*m11*CinvA.transpose() - AinvB.transpose()*m12
             - m21*CinvA.transpose() + m22;
    oriRDCAB += DCAB.transpose() * RDCAB * DCAB.transpose();
    
    RCinvA += DCAB.transpose()*AinvB.transpose()*m11-DCAB.transpose()*m21
              -oriRDCAB*B.transpose();

    RinvAB += m11*CinvA.transpose()*DCAB.transpose() 
            - m12*DCAB.transpose();

    RAinv += m11 + C.transpose()*RCinvA 
            + RinvAB * B.transpose();

  

    // A
    adj_m.m_topLeftMat += Ainv.transpose() * RAinv * Ainv.transpose();
    // B
    adj_m.m_topRightMat += Ainv.transpose() * RinvAB
                        -CinvA.transpose() * oriRDCAB;
    // C
    adj_m.m_bottomLeftMat += RCinvA * Ainv.transpose();
    // D
    adj_m.m_bottomRightMat += oriRDCAB;
  }
  */
};

#endif  // TINY_SYMMETRIC_SPATIAL_DYAD_H
