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

#ifndef TINY_GEOMETRY_H
#define TINY_GEOMETRY_H

#include <vector>

#include "tiny_pose.h"

enum TinyJointType {
  JOINT_FIXED = -1,
  JOINT_PRISMATIC_X = 0,
  JOINT_PRISMATIC_Y,
  JOINT_PRISMATIC_Z,
  JOINT_PRISMATIC_AXIS,
  JOINT_REVOLUTE_X,
  JOINT_REVOLUTE_Y,
  JOINT_REVOLUTE_Z,
  JOINT_REVOLUTE_AXIS,
  JOINT_INVALID,

};

enum TinyGeometryTypes {
  TINY_SPHERE_TYPE = 0,
  TINY_PLANE_TYPE,
  TINY_CAPSULE_TYPE,
  TINY_MESH_TYPE,      // only for visual shapes at the moment
  TINY_BOX_TYPE,       // only for visual shapes at the moment
  TINY_CYLINDER_TYPE,  // unsupported
  TINY_MAX_GEOM_TYPE,
};

template <typename TinyScalar, typename TinyConstants>
class TinyGeometry {
  int m_type;

 public:
  explicit TinyGeometry(int type) : m_type(type) {}
  virtual ~TinyGeometry() {}
  int get_type() const { return m_type; }
};

template <typename TinyScalar, typename TinyConstants>
class TinyBox : public TinyGeometry<TinyScalar, TinyConstants> {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  TinyVector3 m_extents;

 public:
  explicit TinyBox(TinyVector3 extents)
      : TinyGeometry<TinyScalar, TinyConstants>(TINY_BOX_TYPE),
        m_extents(extents) {}

  TinyVector3 get_extents() const { return m_extents; }

  // TinyVector3 compute_local_inertia(
  //     TinyScalar mass) const {
  //   TinyScalar elem =
  //       TinyConstants::fraction(4, 10) * mass * m_radius * m_radius;
  //   return TinyVector3<TinyScalar, TinyConstants>(elem, elem, elem);
  // }
};

template <typename TinyScalar, typename TinyConstants>
class TinySphere : public TinyGeometry<TinyScalar, TinyConstants> {
  TinyScalar m_radius;

 public:
  explicit TinySphere(TinyScalar radius)
      : TinyGeometry<TinyScalar, TinyConstants>(TINY_SPHERE_TYPE),
        m_radius(radius) {}

  TinyScalar get_radius() const { return m_radius; }

  TinyVector3<TinyScalar, TinyConstants> compute_local_inertia(
      TinyScalar mass) const {
    TinyScalar elem =
        TinyConstants::fraction(4, 10) * mass * m_radius * m_radius;
    return TinyVector3<TinyScalar, TinyConstants>(elem, elem, elem);
  }
};

// capsule aligned with the Z axis
template <typename TinyScalar, typename TinyConstants>
class TinyCapsule : public TinyGeometry<TinyScalar, TinyConstants> {
  TinyScalar m_radius;
  TinyScalar m_length;

 public:
  explicit TinyCapsule(TinyScalar radius, TinyScalar length)
      : TinyGeometry<TinyScalar, TinyConstants>(TINY_CAPSULE_TYPE),
        m_radius(radius),
        m_length(length) {}

  TinyScalar get_radius() const { return m_radius; }
  TinyScalar get_length() const { return m_length; }

  TinyVector3<TinyScalar, TinyConstants> compute_local_inertia(
      TinyScalar mass) const {
    TinyScalar lx = TinyConstants::fraction(2, 1) * (m_radius);
    TinyScalar ly = TinyConstants::fraction(2, 1) * (m_radius);
    TinyScalar lz = m_length + TinyConstants::fraction(2, 1) * (m_radius);
    TinyScalar x2 = lx * lx;
    TinyScalar y2 = ly * ly;
    TinyScalar z2 = lz * lz;
    TinyScalar scaledmass = mass * TinyConstants::fraction(1, 12);

    TinyVector3<TinyScalar, TinyConstants> inertia;
    inertia[0] = scaledmass * (y2 + z2);
    inertia[1] = scaledmass * (x2 + z2);
    inertia[2] = scaledmass * (x2 + y2);
    return inertia;
  }
};

template <typename TinyScalar, typename TinyConstants>
class TinyPlane : public TinyGeometry<TinyScalar, TinyConstants> {
  TinyVector3<TinyScalar, TinyConstants> m_normal;
  TinyScalar m_constant;

 public:
  TinyPlane()
      : TinyGeometry<TinyScalar, TinyConstants>(TINY_PLANE_TYPE),
        m_normal(TinyConstants::zero(), TinyConstants::zero(),
                 TinyConstants::one()),
        m_constant(TinyConstants::zero()) {}

  const TinyVector3<TinyScalar, TinyConstants>& get_normal() const {
    return m_normal;
  }
  TinyScalar get_constant() const { return m_constant; }
};

template <typename TinyScalar, typename TinyConstants>
struct TinyContactPoint {
  TinyVector3<TinyScalar, TinyConstants> m_world_normal_on_b, 
                                          adjm_Rworld_normal_on_b;
  TinyVector3<TinyScalar, TinyConstants> m_world_point_on_a,
                                          adjm_Rworld_point_on_a;
  TinyVector3<TinyScalar, TinyConstants> m_world_point_on_b,
                                          adjm_Rworld_point_on_b;
  TinyScalar m_distance, adjm_Rdistance;
};

template <typename TinyScalar, typename TinyConstants>
int contactBoxSphere(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    const TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    const TinyPose<TinyScalar, TinyConstants>& poseB,
    std::vector<TinyContactPoint<TinyScalar, TinyConstants> >& contactsOut,
    bool is_tape,
    std::vector<TinyScalar>& ax) {
  TinyScalar CONTACT_EPSILON = TinyConstants::fraction(1, 100);
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  typedef ::TinyBox<TinyScalar, TinyConstants> TinyBox;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  assert(geomA->get_type() == TINY_BOX_TYPE);
  assert(geomB->get_type() == TINY_SPHERE_TYPE);
 
  TinyBox* boxA = (TinyBox*)geomA;
  TinySphere* sphereB = (TinySphere*)geomB;

  TinyQuaternion<TinyScalar, TinyConstants> box_orientation = poseA.m_orientation;
  auto tmp_posb = poseB.m_position;

  TinyVector3 sphereRelPos = box_orientation.inversed().rotate(
      tmp_posb - poseA.m_position);
  TinyVector3 closestPoint = sphereRelPos;
  TinyVector3 extents = boxA->get_extents();

  closestPoint.setX(std::min(extents.getX(), closestPoint.getX()));
  closestPoint.setX(std::max(-extents.getX(), closestPoint.getX()));
  closestPoint.setY(std::min(extents.getY(), closestPoint.getY()));
  closestPoint.setY(std::max(-extents.getY(), closestPoint.getY()));
  closestPoint.setZ(std::min(extents.getZ(), closestPoint.getZ()));
  closestPoint.setZ(std::max(-extents.getZ(), closestPoint.getZ()));

  TinyVector3 diff = sphereRelPos - closestPoint;

  TinyScalar length = diff.length();
  TinyScalar distance = length - sphereB->get_radius();

  if (distance < CONTACT_EPSILON) {
    TinyVector3 normal_on_a = TinyConstants::one() / length * diff;
    auto tmp_normal_on_a = normal_on_a;

    normal_on_a = box_orientation.rotate(tmp_normal_on_a);

    TinyVector3 normal_on_b = - normal_on_a;
    TinyVector3 point_a_world = box_orientation.rotate(closestPoint) + poseA.m_position;
    


    TinyVector3 point_b_world = point_a_world - distance * normal_on_b;
    TinyContactPoint pt;
    pt.m_world_normal_on_b = normal_on_b;
    pt.m_world_point_on_a = point_a_world;
    pt.m_world_point_on_b = point_b_world;
    pt.m_distance = distance;
    contactsOut.push_back(pt);
    return 1;

  }
  return 0;
}

template <typename TinyScalar, typename TinyConstants>
void adj_contactBoxSphere(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    TinyPose<TinyScalar, TinyConstants>& poseB,
    TinyContactPoint<TinyScalar, TinyConstants>& pt) {
  TinyScalar CONTACT_EPSILON = TinyConstants::fraction(1, 100);
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinyQuaternion<TinyScalar, TinyConstants> TinyQuaternion;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  typedef ::TinyBox<TinyScalar, TinyConstants> TinyBox;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  assert(geomA->get_type() == TINY_BOX_TYPE);
  assert(geomB->get_type() == TINY_SPHERE_TYPE);

  poseA.adjm_Rposition.set_zero();
  poseA.adjm_Rorientation.set_zero();
  poseB.adjm_Rorientation.set_zero();
  poseB.adjm_Rorientation.set_zero();
  
  TinyBox* boxA = (TinyBox*)geomA;
  TinySphere* sphereB = (TinySphere*)geomB;

  TinyQuaternion box_orientation = poseA.m_orientation;
  TinyVector3 sphereRelPos = box_orientation.inversed().rotate(
      poseB.m_position - poseA.m_position);
  TinyQuaternion inv_box_orientation = box_orientation.inversed();
  TinyVector3 closestPoint = sphereRelPos;
  TinyVector3 extents = boxA->get_extents();

  closestPoint.setX(std::min(extents.getX(), closestPoint.getX()));
  closestPoint.setX(std::max(-extents.getX(), closestPoint.getX()));
  closestPoint.setY(std::min(extents.getY(), closestPoint.getY()));
  closestPoint.setY(std::max(-extents.getY(), closestPoint.getY()));
  closestPoint.setZ(std::min(extents.getZ(), closestPoint.getZ()));
  closestPoint.setZ(std::max(-extents.getZ(), closestPoint.getZ()));

  TinyVector3 diff = sphereRelPos - closestPoint;
  TinyScalar length = diff.length();
  TinyScalar distance = length - sphereB->get_radius();

  TinyVector3 normal_on_a = TinyConstants::one() / length * diff;
  auto normal_on_a2 = normal_on_a;
  normal_on_a = box_orientation.rotate(normal_on_a2);

  // adjoint
  TinyScalar Rdistance = TinyConstants::zero();
  TinyScalar Rlength = TinyConstants::zero();
  TinyVector3 Rdiff, RclosetPoint;
  TinyVector3 Rpoint_b_world, Rpoint_a_world, Rnormal_on_b;
  TinyVector3 Rnormal_on_a, Rnormal_on_a2, RsphereRelPos;
  TinyQuaternion Rbox_orientation;

  RsphereRelPos.set_zero();
  Rdiff.set_zero(); Rbox_orientation.set_zero(); RclosetPoint.set_zero();
  Rnormal_on_b.set_zero(); Rpoint_b_world.set_zero(); 
  Rnormal_on_a.set_zero(); Rnormal_on_a2.set_zero(); 
  Rpoint_a_world.set_zero();

  Rdistance += pt.adjm_Rdistance;
  Rpoint_b_world += pt.adjm_Rworld_point_on_b;
  Rpoint_a_world += pt.adjm_Rworld_point_on_a;
  Rnormal_on_b += pt.adjm_Rworld_normal_on_b;


  Rpoint_a_world += Rpoint_b_world;
  Rdistance += -Rpoint_b_world.dot(pt.m_world_normal_on_b);
  Rnormal_on_b += -pt.m_distance * Rpoint_b_world;

  poseA.adjm_Rposition += Rpoint_a_world;
  box_orientation.adj_rotate(Rpoint_a_world, closestPoint, Rbox_orientation, RclosetPoint);
  Rnormal_on_a = -Rnormal_on_b;

  box_orientation.adj_rotate(Rnormal_on_a, 
    normal_on_a2, Rbox_orientation, Rnormal_on_a2);
  
  Rlength += -TinyConstants::one() / length / length * diff.dot(Rnormal_on_a2);
  
  Rdiff += TinyConstants::one() / length * Rnormal_on_a2;
  Rlength += Rdistance;
  Rdiff += diff.adj_length(Rlength);

  RsphereRelPos += Rdiff;
  RclosetPoint += -Rdiff;

  if (sphereRelPos.getX() <= extents.getX() && sphereRelPos.getX() >= -extents.getX()) 
    RsphereRelPos.m_x += RclosetPoint.m_x;
  if (sphereRelPos.getY() <= extents.getY() && sphereRelPos.getY() >= -extents.getY()) 
    RsphereRelPos.m_y += RclosetPoint.m_y;
  if (sphereRelPos.getZ() <= extents.getZ() && sphereRelPos.getZ() >= -extents.getZ()) 
    RsphereRelPos.m_z += RclosetPoint.m_z;

  TinyQuaternion Rinv_box_orientation;
  TinyVector3 RB_A, B_A;
  B_A = poseB.m_position - poseA.m_position;
  Rinv_box_orientation.set_zero(); RB_A.set_zero();
  inv_box_orientation.adj_rotate(RsphereRelPos, B_A, Rinv_box_orientation, RB_A);
  box_orientation.adj_inversed(Rinv_box_orientation, Rbox_orientation);
  poseB.adjm_Rposition += RB_A;
  poseA.adjm_Rposition -= RB_A;
  poseA.adjm_Rorientation += Rbox_orientation;
}

template <typename TinyScalar, typename TinyConstants>
int contactSphereSphere(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    const TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    const TinyPose<TinyScalar, TinyConstants>& poseB,
    std::vector<TinyContactPoint<TinyScalar, TinyConstants> >& contactsOut,
    bool is_tape,
    std::vector<TinyScalar>& ax) {
  TinyScalar CONTACT_EPSILON = TinyConstants::fraction(1, 100);
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  assert(geomA->get_type() == TINY_SPHERE_TYPE);
  assert(geomB->get_type() == TINY_SPHERE_TYPE);
  TinySphere* sphereA = (TinySphere*)geomA;
  TinySphere* sphereB = (TinySphere*)geomB;

  TinyVector3 diff = poseA.m_position - poseB.m_position;
  TinyScalar length = diff.length();
  TinyScalar distance =
      length - (sphereA->get_radius() + sphereB->get_radius());
  TinyVector3 normal_on_b;
  normal_on_b.setValue(TinyConstants::one(), TinyConstants::zero(),
                       TinyConstants::zero());
  if (TinyConstants::getBool(
    distance < CONTACT_EPSILON)) {
    TinyVector3 normal_on_b = TinyConstants::one() / length * diff;
    TinyVector3 point_a_world =
        poseA.m_position - sphereA->get_radius() * normal_on_b;
    TinyVector3 point_b_world = point_a_world - distance * normal_on_b;
    TinyContactPoint pt;
    pt.m_world_normal_on_b = normal_on_b;
    pt.m_world_point_on_a = point_a_world;
    pt.m_world_point_on_b = point_b_world;
    pt.m_distance = distance;
    contactsOut.push_back(pt);
    return 1;
  }
  return 0;
}


template <typename TinyScalar, typename TinyConstants>
void adj_contactSphereSphere(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    TinyPose<TinyScalar, TinyConstants>& poseB,
    TinyContactPoint<TinyScalar, TinyConstants>& pt) {
  TinyScalar CONTACT_EPSILON = TinyConstants::fraction(1, 100);
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  assert(geomA->get_type() == TINY_SPHERE_TYPE);
  assert(geomB->get_type() == TINY_SPHERE_TYPE);
  TinySphere* sphereA = (TinySphere*)geomA;
  TinySphere* sphereB = (TinySphere*)geomB;

  TinyVector3 diff = poseA.m_position - poseB.m_position;
  TinyScalar length = diff.length();
  TinyScalar distance =
      length - (sphereA->get_radius() + sphereB->get_radius());

  TinyScalar Rdistance = TinyConstants::zero();
  TinyScalar Rlength = TinyConstants::zero();
  TinyVector3 Rdiff;
  TinyVector3 Rpoint_b_world, Rpoint_a_world, Rnormal_on_b;
  Rnormal_on_b.set_zero();
  Rpoint_b_world.set_zero();
  Rpoint_a_world.set_zero();
  Rdiff.set_zero();


  Rdistance += pt.adjm_Rdistance;
  Rpoint_b_world += pt.adjm_Rworld_point_on_b;
  Rpoint_a_world += pt.adjm_Rworld_point_on_a;
  Rnormal_on_b += pt.adjm_Rworld_normal_on_b;

  Rpoint_a_world += Rpoint_b_world;
  Rdistance += -Rpoint_b_world.dot(pt.m_world_normal_on_b);
  Rnormal_on_b += -pt.m_distance * Rpoint_b_world;

  poseA.adjm_Rposition += Rpoint_a_world;
  Rnormal_on_b += -sphereA->get_radius() * Rpoint_a_world;
  Rlength += -TinyConstants::one() / length / length * diff.dot(Rnormal_on_b);
  Rdiff += TinyConstants::one() / length * Rnormal_on_b;

  Rlength += Rdistance;
  Rdiff += diff.adj_length(Rlength);
  poseA.adjm_Rposition += Rdiff;
  poseB.adjm_Rposition += - Rdiff;
   
}

template <typename TinyScalar, typename TinyConstants>
void adj_contactPlaneSphere(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    TinyPose<TinyScalar, TinyConstants>& poseB,
    TinyContactPoint<TinyScalar, TinyConstants>& pt) {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  typedef ::TinyPlane<TinyScalar, TinyConstants> TinyPlane;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  assert(geomA->get_type() == TINY_PLANE_TYPE);
  assert(geomB->get_type() == TINY_SPHERE_TYPE);
  TinyPlane* planeA = (TinyPlane*)geomA;
  TinySphere* sphereB = (TinySphere*)geomB;


  TinyScalar Rdistance = TinyConstants::zero();
  TinyScalar Rlength = TinyConstants::zero();
  TinyVector3 Rdiff;
  Rdiff.set_zero();
  TinyVector3 RpointBWorld, RpointAWorld, Rnormal_on_b;
  Rnormal_on_b.set_zero();
  RpointBWorld.set_zero();
  RpointAWorld.set_zero();

  Rdistance += pt.adjm_Rdistance;
  RpointAWorld += pt.adjm_Rworld_point_on_a;
  RpointBWorld += pt.adjm_Rworld_point_on_b;
  poseB.adjm_Rposition += RpointBWorld;
  poseB.adjm_Rposition += Rdistance * planeA->get_normal();

}

template <typename TinyScalar, typename TinyConstants>
int contactPlaneSphere(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    const TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    const TinyPose<TinyScalar, TinyConstants>& poseB,
    std::vector<TinyContactPoint<TinyScalar, TinyConstants> >& contactsOut,
    bool is_tape,
    std::vector<TinyScalar>& ax) {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  typedef ::TinyPlane<TinyScalar, TinyConstants> TinyPlane;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  assert(geomA->get_type() == TINY_PLANE_TYPE);
  assert(geomB->get_type() == TINY_SPHERE_TYPE);
  TinyPlane* planeA = (TinyPlane*)geomA;
  TinySphere* sphereB = (TinySphere*)geomB;

  TinyScalar t =
      -(poseB.m_position.dot(-planeA->get_normal()) + planeA->get_constant());
  TinyVector3 pointAWorld = poseB.m_position + t * -planeA->get_normal();
  TinyScalar distance = t - sphereB->get_radius();
  TinyVector3 pointBWorld =
      poseB.m_position - sphereB->get_radius() * planeA->get_normal();
  TinyContactPoint pt;
  pt.m_world_normal_on_b = -planeA->get_normal();
  pt.m_world_point_on_a = pointAWorld;
  pt.m_world_point_on_b = pointBWorld;
  pt.m_distance = distance;
  TinyScalar CONTACT_EPSILON = TinyConstants::fraction(1, 100);
  if (TinyConstants::getBool(distance > CONTACT_EPSILON))
    return 0;
  contactsOut.push_back(pt);
  return 1;
}

template <typename TinyScalar, typename TinyConstants>
int contactPlaneCapsule(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    const TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    const TinyPose<TinyScalar, TinyConstants>& poseB,
    std::vector<TinyContactPoint<TinyScalar, TinyConstants> >& contactsOut,
    bool is_tape,
    std::vector<TinyScalar>& ax) {
  typedef ::TinyPose<TinyScalar, TinyConstants> TinyPose;
  typedef ::TinyPlane<TinyScalar, TinyConstants> TinyPlane;
  typedef ::TinyCapsule<TinyScalar, TinyConstants> TinyCapsule;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  assert(geomA->get_type() == TINY_PLANE_TYPE);
  assert(geomB->get_type() == TINY_CAPSULE_TYPE);
  TinyCapsule* capsule = (TinyCapsule*)geomB;

  // create twice a plane-sphere contact
  TinySphere sphere(capsule->get_radius());
  // shift the sphere to each end-point
  TinyPose offset;
  offset.m_orientation.set_identity();
  offset.m_position.setValue(
      TinyConstants::zero(), TinyConstants::zero(),
      TinyConstants::fraction(1, 2) * capsule->get_length());
  TinyPose poseEndSphere = poseB * offset;
  int old = contactsOut.size();
  contactPlaneSphere<TinyScalar, TinyConstants>(geomA, poseA, &sphere,
                                                poseEndSphere, contactsOut, is_tape, ax);
  offset.m_position.setValue(
      TinyConstants::zero(), TinyConstants::zero(),
      TinyConstants::fraction(-1, 2) * capsule->get_length());
  poseEndSphere = poseB * offset;
  contactPlaneSphere<TinyScalar, TinyConstants>(geomA, poseA, &sphere,
                                                poseEndSphere, contactsOut, is_tape, ax);

  return contactsOut.size() - old;
}


template <typename TinyScalar, typename TinyConstants>
void adj_contactPlaneCapsule(
    const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    TinyPose<TinyScalar, TinyConstants>& poseA,
    const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    TinyPose<TinyScalar, TinyConstants>& poseB,
    TinyContactPoint<TinyScalar, TinyConstants>& pt) {
  typedef ::TinyPose<TinyScalar, TinyConstants> TinyPose;
  typedef ::TinyPlane<TinyScalar, TinyConstants> TinyPlane;
  typedef ::TinyCapsule<TinyScalar, TinyConstants> TinyCapsule;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  assert(geomA->get_type() == TINY_PLANE_TYPE);
  assert(geomB->get_type() == TINY_CAPSULE_TYPE);
  TinyCapsule* capsule = (TinyCapsule*)geomB;

  // create twice a plane-sphere contact
  TinySphere sphere(capsule->get_radius());
  // shift the sphere to each end-point
  TinyPose offset;
  offset.m_orientation.set_identity();
  offset.m_position.setValue(
      TinyConstants::zero(), TinyConstants::zero(),
      TinyConstants::fraction(1, 2) * capsule->get_length());
  TinyPose poseEndSphere = poseB * offset;

  poseEndSphere.adjm_Rposition.set_zero();
  poseEndSphere.adjm_Rorientation.set_zero();

  std::vector<TinyContactPoint > contactsOut1;
  
  bool is_tape = false;
  std::vector<TinyScalar> ax;
  contactPlaneSphere<TinyScalar, TinyConstants>(geomA, poseA, &sphere,
                                                poseEndSphere, contactsOut1, is_tape, ax);
  if (contactsOut1.size() > 0 && 
    (contactsOut1[0].m_world_point_on_b - pt.m_world_point_on_b).length() < TinyConstants::fraction(1,100000) ) {
    adj_contactPlaneSphere(geomA, poseA, &sphere, poseEndSphere, pt);
    poseB.adj_pose_mul(poseEndSphere, poseB, offset);
    return;
  } 


  offset.m_position.setValue(
      TinyConstants::zero(), TinyConstants::zero(),
      TinyConstants::fraction(-1, 2) * capsule->get_length());
  poseEndSphere = poseB * offset;
  poseEndSphere.adjm_Rposition.set_zero();
  poseEndSphere.adjm_Rorientation.set_zero();
  std::vector<TinyContactPoint > contactsOut2;
  contactPlaneSphere<TinyScalar, TinyConstants>(geomA, poseA, &sphere,
                                                poseEndSphere, contactsOut2, is_tape, ax);
  if (contactsOut2.size() > 0 && 
    (contactsOut2[0].m_world_point_on_b - pt.m_world_point_on_b).length() < TinyConstants::fraction(1,100000) )  {
    adj_contactPlaneSphere(geomA, poseA, &sphere, poseEndSphere, pt);
    poseB.adj_pose_mul(poseEndSphere, poseB, offset);
    return;
  }
}

    // const TinyGeometry<TinyScalar, TinyConstants>* geomA,
    // TinyPose<TinyScalar, TinyConstants>& poseA,
    // const TinyGeometry<TinyScalar, TinyConstants>* geomB,
    // TinyPose<TinyScalar, TinyConstants>& poseB,
    // TinyContactPoint<TinyScalar, TinyConstants>& pt

template <typename TinyScalar, typename TinyConstants>
struct TinyCollisionDispatcher {
  typedef ::TinyGeometry<TinyScalar, TinyConstants> TinyGeometry;
  typedef ::TinyPose<TinyScalar, TinyConstants> TinyPose;
  typedef ::TinyContactPoint<TinyScalar, TinyConstants> TinyContactPoint;

  typedef int (*contact_func)(const TinyGeometry* geomA, const TinyPose& poseA,
                              const TinyGeometry* geomB, const TinyPose& poseB,
                              std::vector<TinyContactPoint>& contactsOut,
    bool is_tape,
    std::vector<TinyScalar>& ax);

  typedef void (*adj_contact_func)(const TinyGeometry* geomA, TinyPose& poseA,
                              const TinyGeometry* geomB, TinyPose& poseB,
                              TinyContactPoint& pt);

  contact_func m_contactFuncs[TINY_MAX_GEOM_TYPE][TINY_MAX_GEOM_TYPE];
  adj_contact_func adjm_contactFuncs[TINY_MAX_GEOM_TYPE][TINY_MAX_GEOM_TYPE];

  TinyCollisionDispatcher() {
    for (int i = 0; i < TINY_MAX_GEOM_TYPE; i++) {
      for (int j = 0; j < TINY_MAX_GEOM_TYPE; j++) {
        m_contactFuncs[i][j] = 0;
        adjm_contactFuncs[i][j] = 0;
      }
    }
    m_contactFuncs[TINY_SPHERE_TYPE][TINY_SPHERE_TYPE] = contactSphereSphere;
    m_contactFuncs[TINY_PLANE_TYPE][TINY_SPHERE_TYPE] = contactPlaneSphere;
    m_contactFuncs[TINY_PLANE_TYPE][TINY_CAPSULE_TYPE] = contactPlaneCapsule;
    m_contactFuncs[TINY_BOX_TYPE][TINY_SPHERE_TYPE] = contactBoxSphere;

    adjm_contactFuncs[TINY_SPHERE_TYPE][TINY_SPHERE_TYPE] = adj_contactSphereSphere;
    adjm_contactFuncs[TINY_PLANE_TYPE][TINY_SPHERE_TYPE] = adj_contactPlaneSphere;
    adjm_contactFuncs[TINY_PLANE_TYPE][TINY_CAPSULE_TYPE] = adj_contactPlaneCapsule;
    adjm_contactFuncs[TINY_BOX_TYPE][TINY_SPHERE_TYPE] = adj_contactBoxSphere;

  }

  int computeContacts(const TinyGeometry* geomA, const TinyPose& poseA,
                      const TinyGeometry* geomB, const TinyPose& poseB,
                      std::vector<TinyContactPoint>& contactsOut) {
    contact_func f = m_contactFuncs[geomA->get_type()][geomB->get_type()];
    
    std::vector<TinyScalar> ax;
    if (f) {
      return f(geomA, poseA, geomB, poseB, contactsOut, false, ax);
    }
    return 0;
  }

  void adj_computeContacts(const TinyGeometry* geomA, TinyPose& poseA,
                      const TinyGeometry* geomB, TinyPose& poseB,
                      TinyContactPoint& contact_point) {
    adj_contact_func adjf = 
              adjm_contactFuncs[geomA->get_type()][geomB->get_type()];
    if (adjf) {
      adjf(geomA, poseA, geomB, poseB, contact_point);
    }
    
  }

  std::vector<TinyContactPoint> compute_contacts(const TinyGeometry* geomA,
                                                 const TinyPose& poseA,
                                                 const TinyGeometry* geomB,
                                                 const TinyPose& poseB) {
    std::vector<TinyContactPoint> pts;
    // int num = computeContacts(geomA, poseA, geomB, poseB, pts);
    int num = 0; // is_tape
    return pts;
  }
};

#endif  // TINY_GEOMETRY_H