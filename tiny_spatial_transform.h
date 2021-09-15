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

#ifndef TINY_SPATIAL_TRANSFORM_H
#define TINY_SPATIAL_TRANSFORM_H

#include "tiny_matrix3x3.h"
#include "tiny_spatial_motion_vector.h"
#include "tiny_vector3.h"

/**
 * We follow the spatial algebra from Featherstone, but use right-handed
 * transformation matrix applications so that the rotations need to be
 * transposed and the transformation matrices are multiplied from right to left.
 */
template <typename TinyScalar, typename TinyConstants>
class TinySpatialTransform {
 public:
  enum eOutputOperation { None = 0, Add = 1, Subtract = 2 };

  typedef ::TinyMatrix3x3<TinyScalar, TinyConstants> TinyMatrix3x3;
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinySpatialMotionVector<TinyScalar, TinyConstants>
      TinySpatialMotionVector;

  TinyVector3 m_translation;
  TinyMatrix3x3 m_rotation;

  TinySpatialTransform() { set_identity(); }

  void set_identity() {
    m_translation.set_zero();
    m_rotation.set_identity();
  }
  void set_zero() {
    m_translation.set_zero();
    m_rotation.set_zero();
  }
  /**
   * X1*X2 = plx(E1*E2, r2 + E2T*r1)
   */
  TinySpatialTransform operator*(const TinySpatialTransform& t) const {
    TinySpatialTransform tr = *this;
    tr.m_translation += m_rotation * t.m_translation;
    tr.m_rotation *= t.m_rotation;
    return tr;
  }


  TinySpatialTransform& operator+=(const TinySpatialTransform& t) {
    m_rotation += t.m_rotation;
    m_translation += t.m_translation;
    return *this;
  }

  // template <typename SpatialVectorType>
  // void adj_st_multiply(const SpatialVectorType& R, const TinySpatialTransform& a,
  //                                 TinySpatialTransform& adj_t1, SpatialVectorType& adj_t2) const {
  //   const TinyMatrix3x3& E1 = m_rotation;
  //   const TinyVector3& r1 = m_translation;
  //   const TinyMatrix3x3& E2 = a.m_rotation;
  //   const TinyVector3& r2 = a.m_translation;

  //   const TinyMatrix3x3& E0 = R.m_rotation;
  //   const TinyVector3& r0 = R.m_translation;

  //   adj_t1.m_rotation += E0*E2.transpose();
  //   adj_t1.m_translation += E2.transpose()*r0;
  //   adj_t2.m_rotation += E1.transpose()*E0+TinyMatrix3x3::vvt(r0, r1);
  //   adj_t2.m_translation += r0;
  // }

  template <typename SpatialVectorType>
  void adj_st_multiply(const SpatialVectorType& R, const TinySpatialTransform& a,
                                  TinySpatialTransform& adj_t1, SpatialVectorType& adj_t2) const {
    const TinyMatrix3x3& E1 = m_rotation;
    const TinyVector3& r1 = m_translation;
    const TinyMatrix3x3& E2 = a.m_rotation;
    const TinyVector3& r2 = a.m_translation;

    const TinyMatrix3x3& E0 = R.m_rotation;
    const TinyVector3& r0 = R.m_translation;

    adj_t1.m_rotation += E0*E2.transpose()+TinyMatrix3x3::vvt(r0, r2);
    adj_t1.m_translation += r0;
    adj_t2.m_rotation += E1.transpose()*E0;
    adj_t2.m_translation += E1.transpose()*r0;
  }

  void print(const char* name) const {
    printf("%s\n", name);
    double x = TinyConstants::getDouble(m_translation.getX());
    double y = TinyConstants::getDouble(m_translation.getY());
    double z = TinyConstants::getDouble(m_translation.getZ());
    printf("translation: %f,%f,%f\n", x, y, z);
    printf("rotation:\n");
    for (int r = 0; r < 3; r++) {
      for (int c = 0; c < 3; c++) {
        double v = TinyConstants::getDouble(m_rotation[r][c]);
        printf("%f, ", v);
      }
      printf("\n");
    }
  }

  template <typename SpatialVectorType>
  SpatialVectorType transformRotateOnly(const SpatialVectorType& vecIn) const {
    SpatialVectorType outVec;
    outVec.m_topVec = m_rotation * vecIn.m_topVec;
    outVec.m_bottomVec = m_rotation * vecIn.m_bottomVec;
    return outVec;
  }

  TinyVector3 apply(const TinyVector3& point) const {
    return m_rotation * point + m_translation;
  }
  void adj_apply(const TinyVector3& R, const TinyVector3& point,
    TinySpatialTransform& Rst, TinyVector3& Rpoint) const {
    Rpoint += m_rotation.transpose() * R;
    Rst.m_rotation += TinyMatrix3x3::vvt(R, point);
    Rst.m_translation += R;
  }

  TinyVector3 apply_inverse(const TinyVector3& point) const {
    return m_rotation.transpose() * (point - m_translation);
  }

  TinySpatialTransform get_inverse() const {
    TinySpatialTransform inv;
    inv.m_rotation = m_rotation.transpose();
    inv.m_translation = inv.m_rotation * -m_translation;
    return inv;
  }

  /**
   * V = mv(w, v)
   * X*V = mv(E*w, E*(v - r x w))
   */
  template <typename SpatialVectorType>
  SpatialVectorType apply(const SpatialVectorType& inVec) const {
    SpatialVectorType outVec;

    TinyVector3 rxw = inVec.m_topVec.cross(m_translation);
    TinyVector3 v_rxw = inVec.m_bottomVec + rxw;

    TinyVector3 tmp3 = m_rotation.transpose() * v_rxw;
    TinyVector3 tmp4 = m_rotation.transpose() * inVec.m_topVec;

    outVec.m_topVec = tmp4;
    outVec.m_bottomVec = tmp3;

    return outVec;
  }

  /**
   * V = mv(w, v)
   * X*V = mv(E*w, E*(v - r x w))
   */
  template <typename SpatialVectorType>
  void adj_st_apply(const SpatialVectorType& R, const SpatialVectorType& a,
                                  TinySpatialTransform& adj_t, SpatialVectorType& adj_a) const {
    const TinyVector3& a1 = a.m_topVec;
    const TinyVector3& a2 = a.m_bottomVec;
    const TinyVector3& b1 = R.m_topVec;
    const TinyVector3& b2 = R.m_bottomVec;
    const TinyMatrix3x3& E = m_rotation;
    const TinyVector3& r = m_translation;


    adj_t.m_rotation += TinyMatrix3x3::vvt(a1, b1)+TinyMatrix3x3::vvt(a2-r.cross(a1), b2);
    adj_t.m_translation += (E*b2).cross(a1);
    adj_a.m_topVec += E*b1-(E*b2).cross(r);
    adj_a.m_bottomVec += E*b2;
  }

  /**
   * V = mv(w, v)
   * inv(X)*V = mv(ET*w, ET*v + r x (ET*w))
   */
  template <typename SpatialVectorType>
  SpatialVectorType apply_inverse(const SpatialVectorType& inVec) const {
    SpatialVectorType outVec;
    outVec.m_topVec = m_rotation * inVec.m_topVec;
    outVec.m_bottomVec =
        m_rotation * inVec.m_bottomVec + m_translation.cross(outVec.m_topVec);
    return outVec;
  }

  template <typename SpatialVectorType>
  void adj_st_apply_inverse(const SpatialVectorType& R, const SpatialVectorType& a,
                                  TinySpatialTransform& adj_t, SpatialVectorType& adj_a) const {
    const TinyVector3& a1 = a.m_topVec;
    const TinyVector3& a2 = a.m_bottomVec;
    const TinyVector3& b1 = R.m_topVec;
    const TinyVector3& b2 = R.m_bottomVec;
    const TinyMatrix3x3& E = m_rotation;
    const TinyVector3& r = m_translation;


    adj_t.m_rotation += TinyMatrix3x3::vvt(b1+b2.cross(r), a1)
                        +TinyMatrix3x3::vvt(b2, a2);
    adj_t.m_translation += -(b2).cross(E*a1);
    adj_a.m_topVec += E.transpose()*(b1+r.cross(b2));
    adj_a.m_bottomVec += E.transpose()*b2;
  }



  /**
   * F = fv(n, f)
   * XT*F = fv(ETn + rxETf, ETf)
   */
  template <typename SpatialVectorType>
  SpatialVectorType apply_transpose(const SpatialVectorType& inVec) const {
    SpatialVectorType outVec;
    outVec.m_bottomVec = m_rotation * inVec.m_bottomVec;
    outVec.m_topVec = m_rotation * inVec.m_topVec;
    outVec.m_topVec +=
        TinyVectorCrossMatrix(m_translation) * outVec.m_bottomVec;

    return outVec;
  }

  /**
   * F = fv(n, f)
   * XT*F = fv(ETn + rxETf, ETf)
   */
  template <typename SpatialVectorType>
  void adj_st_apply_trans(const SpatialVectorType& R, const SpatialVectorType& a,
                                  TinySpatialTransform& adj_t, SpatialVectorType& adj_a) const {
    const TinyVector3& a1 = a.m_topVec;
    const TinyVector3& a2 = a.m_bottomVec;
    const TinyVector3& b1 = R.m_topVec;
    const TinyVector3& b2 = R.m_bottomVec;
    const TinyMatrix3x3& E = m_rotation;
    const TinyVector3& r = m_translation;


    adj_t.m_rotation += TinyMatrix3x3::vvt(b1, a1)+TinyMatrix3x3::vvt(b1.cross(r)+b2, a2);
    adj_t.m_translation += -b1.cross(E*a2);
    adj_a.m_topVec += E.transpose()*b1;
    adj_a.m_bottomVec += E.transpose()*(b1.cross(r)+b2);

  }

  /**
   * F = fv(n, f)
   * X^* F = fv(E(n - rxf), Ef)
   */
  template <typename SpatialVectorType>
  SpatialVectorType apply_inverse_transpose(
      const SpatialVectorType& inVec) const {
    const TinyVector3& n = inVec.m_topVec;
    const TinyVector3& f = inVec.m_bottomVec;
    SpatialVectorType outVec;
    outVec.m_topVec = m_rotation.transpose() * (n - m_translation.cross(f));
    outVec.m_bottomVec = m_rotation.transpose() * f;

    return outVec;
  }

  /**
   * F = fv(n, f)
   * X^* F = fv(E(n - rxf), Ef)
   */
  // template <typename SpatialVectorType>
  // SpatialVectorType adj_st_apply_inverse_transpose(
  //   const SpatialVectorType& R, const SpatialVectorType& a,
  //   TinySpatialTransform& adj_t, SpatialVectorType& adj_a) const {

  //   const TinyVector3& n = inVec.m_topVec;
  //   const TinyVector3& f = inVec.m_bottomVec;
  //   SpatialVectorType outVec;
  //   outVec.m_topVec = m_rotation.transpose() * (n - m_translation.cross(f));
  //   outVec.m_bottomVec = m_rotation.transpose() * f;

  //   return outVec;
  // }
};

#endif  // TINY_SPATIAL_TRANSFORM_H
