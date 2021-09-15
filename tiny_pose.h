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

#ifndef TINY_POSE_H
#define TINY_POSE_H

#include "tiny_quaternion.h"
#include "tiny_vector3.h"

/// TinyPose is a coordinate frame specified as position and orientation
/// quaternion.
template <typename TinyScalar, typename TinyConstants>
struct TinyPose {
  TinyVector3<TinyScalar, TinyConstants> m_position, adjm_Rposition;
  TinyQuaternion<TinyScalar, TinyConstants> m_orientation, adjm_Rorientation;

  TinyPose() {
    adjm_Rposition.set_zero();
    adjm_Rorientation.set_zero();
  }

  TinyPose(const TinyVector3<TinyScalar, TinyConstants>& position,
           const TinyQuaternion<TinyScalar, TinyConstants>& orientation)
      : m_position(position), m_orientation(orientation) {
    adjm_Rposition.set_zero();
    adjm_Rorientation.set_zero();
  }

  TinyVector3<TinyScalar, TinyConstants> transform(
      const TinyVector3<TinyScalar, TinyConstants>& point) const {
    return m_orientation.rotate(point) + m_position;
  }

  TinyVector3<TinyScalar, TinyConstants> inverse_transform(
      const TinyVector3<TinyScalar, TinyConstants>& point) const {
    TinyVector3<TinyScalar, TinyConstants> point_out;
    point_out = point - m_position;
    return m_orientation.inversed().rotate(point_out);
  }

  void adj_pose_mul(TinyPose& res, TinyPose& a, TinyPose& b) {
    a.m_orientation.adj_mul_q(res.adjm_Rorientation, b.m_orientation, a.adjm_Rorientation, b.adjm_Rorientation);
    a.adjm_Rposition += res.adjm_Rposition;
    a.m_orientation.adj_rotate(res.adjm_Rposition, b.m_position, a.adjm_Rorientation, b.adjm_Rposition);
    res.adjm_Rposition.set_zero();
    res.adjm_Rorientation.set_zero();
  }

  TinyPose operator*(const TinyPose& b) const {
    const TinyPose& a = *this;
    TinyPose res;
    res.m_position = a.m_position + a.m_orientation.rotate(b.m_position);
    res.m_orientation = a.m_orientation * b.m_orientation;
    return res;
  }

  void set_identity() {
    m_position.set_zero();
    m_orientation.set_identity();
  }

  void inverse() {
    m_orientation = m_orientation.inversed();
    m_position = m_orientation.rotate(-m_position);
  }
};

#endif  // TINY_POSE_H
