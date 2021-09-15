// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
     
#include <pybind11/operators.h>  
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdio.h> 
    
#include "fix64_scalar.h" 
#include "tiny_double_utils.h"
#include "tiny_matrix3x3.h"
#include "tiny_matrix_x.h"
#include "tiny_mb_constraint_solver.h"
#include "tiny_multi_body.h"
#include "tiny_pose.h"
#include "tiny_quaternion.h"
#include "tiny_urdf_structures.h"
#include "tiny_raycast.h"  
#include "tiny_rigid_body.h"
#include "tiny_urdf_to_multi_body.h"
#include "tiny_vector3.h"
#include "tiny_world.h"
#include "examples/motion_import.h"
#include "examples/tiny_urdf_parser.h"
#include "tiny_metrics.h"
// #include "examples/pybullet_urdf_import.h"
// #include "examples/pybullet_visualizer_api.h"

  
  
// typedef PyBulletVisualizerAPI VisualizerAPI;
typedef double Scalar; 
typedef DoubleUtils Utils;
 
std::string to_string(const Scalar &x) { 
  return std::to_string(Utils::getDouble(x));
} 
    
template <typename TinyScalar, typename TinyConstants>
struct UrdfToMultiBody2 { 
  typedef ::TinyUrdfStructures<TinyScalar, TinyConstants> TinyUrdfStructures;

  void convert(TinyUrdfStructures *urdf_structures,
               TinyWorld<TinyScalar, TinyConstants> *world,
               TinyMultiBody<TinyScalar, TinyConstants> *mb) {
 
    char search_path[4096];
    sprintf(search_path, "./data/");
    world->meshcat_viz.m_path_prefix = search_path;
    std::string texture_path = "checker_purple.jpg";
    world->meshcat_viz.convert_visuals(*urdf_structures, texture_path);
 
                   
    TinyUrdfToMultiBody<TinyScalar, TinyConstants>::convert_to_multi_body(
        *urdf_structures, *world, *mb);
    mb->initialize();
  }
}; 

using std::vector;
template <typename TinyScalar, typename TinyConstants>
void do_mix(TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &state,
  vector<double> &q, vector<double> &qd) {
  int curr_q = 0, curr_qd = 0;
  for (auto *mb : world.m_multi_bodies) {
    for (int i = 0; i < mb->dof(); ++i)
      state.push_back(TinyConstants::scalar_from_double(q[curr_q++]));
    for (int i = 0; i < mb->dof_qd(); ++i)
      state.push_back(TinyConstants::scalar_from_double(qd[curr_qd++]));
  }
} 
 
   
template <typename TinyScalar, typename TinyConstants>
vector<TinyScalar> get_joints(TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &q, vector<TinyScalar> &tans) {
  int curr_q = 0;
  for (auto *mb : world.m_multi_bodies) {
    vector<TinyScalar> tmp_v = vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof());
    mb->forward_kinematics(tmp_v);

    curr_q += mb->dof();
    auto root = mb->get_world_transform(-1).m_translation;
    tans.push_back(root.x());
    tans.push_back(root.y());
    tans.push_back(root.z());
    for (int i = 0; i < mb->m_links.size(); ++i) {
      auto tr = mb->get_world_transform(i);
      auto pos = tr.m_translation;
      tans.push_back(pos.x());
      tans.push_back(pos.y());
      tans.push_back(pos.z());
    }
  }
  return tans;
}

template <typename TinyScalar, typename TinyConstants>
void adj_get_joints(
  TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &Rtans,
  vector<TinyScalar> &q, vector<double> &dldq) {
  int curr_q = 0; 
  int i_tan = 0;
  for (auto *mb : world.m_multi_bodies) {
    vector<TinyScalar> tmp_v = vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof());
    mb->adj_constraint_set_zero();
    mb->adj_f_kinematics(tmp_v);

    TinyVector3<TinyScalar, TinyConstants> Rroot; 
    Rroot.set_zero();
    Rroot.m_x = Rtans[i_tan++];
    Rroot.m_y = Rtans[i_tan++];
    Rroot.m_z = Rtans[i_tan++];
    mb->adj_get_world_transform(-1).m_translation += Rroot;

    for (int i = 0; i < mb->m_links.size(); ++i) {
      TinyVector3<TinyScalar, TinyConstants> Rpos;
      Rpos.set_zero();
      Rpos.m_x = Rtans[i_tan++];
      Rpos.m_y = Rtans[i_tan++];
      Rpos.m_z = Rtans[i_tan++];
      mb->adj_get_world_transform(i).m_translation += Rpos;
    }

    TinyVectorX<TinyScalar, TinyConstants> Rq(tmp_v.size()), Rqd(0), Rqdd(0);
    Rq.set_zero();
    mb->adj_fk(Rq, Rqd, Rqdd, tmp_v);

    for (int i = 0; i < Rq.m_size; i++)
      dldq.push_back(TinyConstants::getDouble(Rq[i]));
    mb->adj_constraint_set_zero();
 
    curr_q += mb->dof();
  }
}

//input: dq: [sum_dof]
//output: tans: [(num_links + num_mb)*3]
template <typename TinyScalar, typename TinyConstants>
vector<double> forward_get_joints(vector<double> &dq, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq)
    q.push_back(TinyConstants::scalar_from_double(ele));
  get_joints(world, q, tans);
  vector<double> ans;
  for (auto ele : tans)
    ans.push_back(TinyConstants::getDouble(ele));
  return ans;
}
   
//input: dq: [sum_dof], dldans: [(num_links + num_mb)*3]
//output: dldq: [sum_dof]
template <typename TinyScalar, typename TinyConstants>
vector<double> backward_get_joints(vector<double> &dq, vector<double> &dldans, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq)
    q.push_back(TinyConstants::scalar_from_double(ele));

  vector<double> dldq;
  vector<TinyScalar> Rtans;
  for (auto ele : dldans)
    Rtans.push_back(TinyConstants::scalar_from_double(ele));
  adj_get_joints(world, Rtans, q, dldq);
  return dldq;  
}   
 

template <typename TinyScalar, typename TinyConstants>
vector<TinyScalar> get_com(TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &q, vector<TinyScalar> &tans) {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  int curr_q = 0;
  for (auto *mb : world.m_multi_bodies) {
    TinyVector3 ans(TinyConstants::zero(),TinyConstants::zero(),TinyConstants::zero());
    TinyScalar totm = TinyConstants::zero();
    mb->forward_kinematics(vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof()));
    curr_q += mb->dof();
    for (int i = 0; i < mb->m_links.size(); ++i) {
      TinyVector3 tmp = mb->get_world_com(i);
      TinyScalar mass = mb->m_links[i].m_I(3,3);
      ans = ans + tmp * mass;
      totm = totm + mass;
    }
    ans = ans * (1./ totm);
    tans.push_back(ans[0]);
    tans.push_back(ans[1]);
    tans.push_back(ans[2]);
  }
  return tans;
}
 
//input: dq: [sum_dof]
//output: tans: [3*num_mb]
template <typename TinyScalar, typename TinyConstants>
vector<double> forward_get_com(vector<double> &dq, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq)
    q.push_back(TinyConstants::scalar_from_double(ele));
  get_com(world, q, tans);
  vector<double> ans;
  for (auto ele : tans)
    ans.push_back(TinyConstants::getDouble(ele));
  return ans;
}

 

template <typename TinyScalar, typename TinyConstants>
void adj_get_com(
  TinyWorld<TinyScalar, TinyConstants> &world, vector<TinyScalar> &Rtans,
  vector<TinyScalar> &q, vector<double> &dldq) {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;

  for (int ii = 0; ii < Rtans.size(); ii++) {
  }
  int curr_q = 0;
  int curr_mb_i = 0;
  for (auto *mb : world.m_multi_bodies) { 
 
    TinyVector3 Rtans_vec(Rtans[3*curr_mb_i], Rtans[3*curr_mb_i+1], Rtans[3*curr_mb_i+2]);
    curr_mb_i++; 

    vector<TinyScalar> tmp_v = vector<TinyScalar>(q.begin()+curr_q, q.begin()+curr_q+mb->dof());
    mb->adj_constraint_set_zero();
    mb->adj_f_kinematics(tmp_v);
    TinyScalar totm = TinyConstants::zero();
    for (int i = 0; i < mb->m_links.size(); ++i) {
      TinyScalar mass = mb->m_links[i].m_I(3,3);
      totm = totm + mass;
    }

    for (int i = mb->m_links.size()-1; i >= 0; --i) {
      TinyScalar mass = mb->m_links[i].m_I(3,3);
      TinyVector3 Rv = Rtans_vec * mass * (1./totm);
      mb->adj_get_world_com(i, Rv); 
    }

    TinyVectorX<TinyScalar, TinyConstants> Rq(tmp_v.size()), Rqd(0), Rqdd(0);
    Rq.set_zero();
    mb->adj_fk(Rq, Rqd, Rqdd, tmp_v);
 
    for (int i = 0; i < Rq.m_size; i++)
      dldq.push_back(TinyConstants::getDouble(Rq[i]));
    mb->adj_constraint_set_zero();

    curr_q += mb->dof();
  }

} 

//input: dq: [sum_dof], dldans: [3*num_mb]
//output: dldq: [sum_dof]
template <typename TinyScalar, typename TinyConstants>
vector<double> backward_get_com(vector<double> &dq, vector<double> &dldans, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> q, tans;
  for (auto ele : dq) {
    q.push_back(TinyConstants::scalar_from_double(ele));
  }

  vector<double> dldq;
  vector<TinyScalar> Rtans;
  for (auto ele : dldans) {
    Rtans.push_back(TinyConstants::scalar_from_double(ele));
  }
  adj_get_com(world, Rtans, q, dldq);

  return dldq;
}

template <typename TinyScalar, typename TinyConstants>
vector<vector<double> > forward_step(vector<double> &q, vector<double> &qd, 
  vector<double> &tau, TinyWorld<TinyScalar, TinyConstants> &world) {
  vector<TinyScalar> state, taut;
  do_mix(world, state, q, qd);
  for (auto ele : tau)
    taut.push_back(TinyConstants::scalar_from_double(ele));
  TinyVectorX<TinyScalar, TinyConstants> taux(taut);
  world.adj_step_torch(state, taux);  
  vector<double> ansq, ansqd;
  for (auto *mb : world.m_multi_bodies) {
    auto state = mb->adjm_this_state;
    int dof1 = mb->dof(), dof2 = mb->dof_qd(); 
    for (int i = 0; i < dof1; ++i)
      ansq.push_back(TinyConstants::getDouble(state[i]));
    for (int i = 0; i < dof2; ++i)
      ansqd.push_back(TinyConstants::getDouble(state[i+dof1])); 
  } 
  return {ansq, ansqd};
}      
  
template <typename TinyScalar, typename TinyConstants>
vector<vector<double> > backward_step(vector<double> &q, vector<double> &qd, vector<double> &tau, 
  vector<double> &dldq, vector<double> &dldqd,  
  TinyWorld<TinyScalar, TinyConstants> &world) { 
  static int n_step = 0;  

  vector<TinyScalar> Rt, taut, state;
  do_mix(world, Rt, dldq, dldqd);  
  do_mix(world, state, q, qd);
  for (auto ele : tau)
    taut.push_back(TinyConstants::scalar_from_double(ele)); 
  TinyVectorX<TinyScalar, TinyConstants> R(Rt), taux(taut);     
                  
  std::vector<TinyVectorX<TinyScalar, TinyConstants>> ans, our_ans;
  TinyVectorX<TinyScalar, TinyConstants> tR(Rt), ttaux(taut), tstate(state.size());   
            
  ans = world.adj_back_step_man(R, taux, state); 
 
  auto ans_state = ans.back();
  vector<double> ansq, ansqd, anstau;  
  int curr_state = 0, curr_body = 0;
  for (auto *mb : world.m_multi_bodies) {
    int dof1 = mb->dof(), dof2 = mb->dof_qd(), dof3 = mb->dof_u();;
    for (int i = 0; i < dof1; ++i)
      ansq.push_back(TinyConstants::getDouble(ans_state[curr_state++]));
    for (int i = 0; i < dof2; ++i)
      ansqd.push_back(TinyConstants::getDouble(ans_state[curr_state++]));
    for (int i = 0; i < dof3; ++i)
      anstau.push_back(TinyConstants::getDouble(ans[curr_body][i]));
    curr_body++;
  }
   
  // if (n_step==0)
  //   exit(0);
  n_step++;
      

  return {ansq, ansqd, anstau};
}

namespace py = pybind11;
  
PYBIND11_MODULE(pydiffarti, m) {
  m.doc() = R"pbdoc(
        tiny differentiable physics python plugin
        -----------------------

        .. currentmodule:: pydiffarti

        .. autosummary::
           :toctree: _generate

    )pbdoc";

  m.def("forward_step", ::forward_step<Scalar, Utils>, py::return_value_policy::copy)
    .def("backward_step", ::backward_step<Scalar, Utils>, py::return_value_policy::copy)
    .def("forward_get_joints", ::forward_get_joints<Scalar, Utils>, py::return_value_policy::copy)
    .def("backward_get_joints", ::backward_get_joints<Scalar, Utils>, py::return_value_policy::copy)
    .def("forward_get_com", ::forward_get_com<Scalar, Utils>, py::return_value_policy::copy)
    .def("backward_get_com", ::backward_get_com<Scalar, Utils>, py::return_value_policy::copy);
  // py::class_<VisualizerAPI>(m, "VisualizerAPI")
  //   .def(py::init())
  //   .def("setAdditionalSearchPath", &VisualizerAPI::setAdditionalSearchPath)
  //   .def("connect", &VisualizerAPI::connect)
  //   .def("disconnect", &VisualizerAPI::disconnect) 
  //   .def("canSubmitCommand", &VisualizerAPI::canSubmitCommand);
    // .def("resetSimulation", &VisualizerAPI::resetSimulation);
 
  // py::class_<PyBulletUrdfImport<Scalar, Utils>>(m, "PyBulletUrdfImport") 
  //   // .def_static("convert_visuals", &PyBulletUrdfImport<Scalar, Utils>::convert_visuals1)
  //   .def_static("sync_graphics_transforms", &PyBulletUrdfImport<Scalar, Utils>::sync_graphics_transforms);
  py::class_<Scalar>(m, "Scalar");
  py::class_<Utils>(m, "Utils")
    .def_static("getDouble", &Utils::getDouble)
    .def_static("scalar_from_double", &Utils::scalar_from_double)
    .def_static("zero", &Utils::zero)
    .def_static("fraction", &Utils::fraction);
  
  py::class_<TinyVector3<Scalar, Utils>>(m, "TinyVector3")
      .def(py::init<Scalar, Scalar, Scalar>())
      .def("set_zero", &TinyVector3<Scalar, Utils>::set_zero)
      .def_readwrite("x", &TinyVector3<Scalar, Utils>::m_x)
      .def_readwrite("y", &TinyVector3<Scalar, Utils>::m_y)
      .def_readwrite("z", &TinyVector3<Scalar, Utils>::m_z)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self += py::self)
      .def(py::self -= py::self)
      .def(-py::self)
      .def("__repr__",
           [](const TinyVector3<Scalar, Utils> &a) {
             return "[" + to_string(a.m_x) + " " + to_string(a.m_y) +
                    " " + to_string(a.m_z) + "]";
           })
      .def("__getitem__", [](const TinyVector3<Scalar, Utils> &a,
                             int i) { return a[i]; })
      .def("__setitem__", [](TinyVector3<Scalar, Utils> &a, int i,
                             Scalar v) { a[i] = v; });

  py::class_<TinyGeometry<Scalar, Utils>,
             std::unique_ptr<TinyGeometry<Scalar, Utils>>>
      geom(m, "TinyGeometry");

  geom.def(py::init<int>())
      .def("get_type", &TinyGeometry<Scalar, Utils>::get_type);

  py::class_<TinySphere<Scalar, Utils>,
             std::unique_ptr<TinySphere<Scalar, Utils>>>(m, "TinySphere",
                                                               geom)
      .def(py::init<Scalar>())
      .def("get_radius", &TinySphere<Scalar, Utils>::get_radius);

  py::class_<TinyPlane<Scalar, Utils>,
             std::unique_ptr<TinyPlane<Scalar, Utils>>>(m, "TinyPlane",
                                                              geom)
      .def(py::init<>())
      .def("get_normal", &TinyPlane<Scalar, Utils>::get_normal);

  py::class_<TinyRigidBody<Scalar, Utils>,
             std::unique_ptr<TinyRigidBody<Scalar, Utils>>>(
      m, "TinyRigidBody")
      .def(py::init<Scalar, TinyGeometry<Scalar, Utils> *>(),
           py::return_value_policy::reference_internal)
      .def_readwrite("world_pose",
                     &TinyRigidBody<Scalar, Utils>::m_world_pose)
      .def_readwrite("collision_geometry",
                     &TinyRigidBody<Scalar, Utils>::m_geometry)
      .def_readwrite("linear_velocity",
                     &TinyRigidBody<Scalar, Utils>::m_linear_velocity)
      .def_readwrite("angular_velocity",
                     &TinyRigidBody<Scalar, Utils>::m_angular_velocity)
      .def_readwrite("local_inertia",
                     &TinyRigidBody<Scalar, Utils>::m_local_inertia)
      .def_readwrite("total_force",
                     &TinyRigidBody<Scalar, Utils>::m_total_force)
      .def_readwrite("total_torque",
                     &TinyRigidBody<Scalar, Utils>::m_total_torque)
      .def_readwrite("user_index",
                     &TinyRigidBody<Scalar, Utils>::m_user_index)

      .def("apply_gravity", &TinyRigidBody<Scalar, Utils>::apply_gravity)
      .def("apply_force_impulse",
           &TinyRigidBody<Scalar, Utils>::apply_force_impulse)

      .def("apply_central_force",
           &TinyRigidBody<Scalar, Utils>::apply_central_force)
      .def("apply_impulse", &TinyRigidBody<Scalar, Utils>::apply_impulse)
      .def("clear_forces", &TinyRigidBody<Scalar, Utils>::clear_forces)
      .def("integrate", &TinyRigidBody<Scalar, Utils>::integrate);

  py::class_<TinyMatrixXxX<Scalar, Utils>>(m, "TinyMatrixXxX")
      .def(py::init<int, int>())
      .def("inversed", &TinyMatrixXxX<Scalar, Utils>::inversed)
      .def("set_zero", &TinyMatrixXxX<Scalar, Utils>::set_zero)
      .def("print", &TinyMatrixXxX<Scalar, Utils>::print)
      .def("__getitem__", [](const TinyMatrixXxX<Scalar, Utils>& a,
          int row, int col) { return a(row,col); })
      .def_readonly("num_rows", &TinyMatrixXxX<Scalar, Utils>::m_rows)
      .def_readonly("num_columns", &TinyMatrixXxX<Scalar, Utils>::m_cols)
      
      ;

  py::class_<TinyMatrix3x3<Scalar, Utils>>(m, "TinyMatrix3x3")
      .def(py::init<>())
      .def(py::init<TinyQuaternion<Scalar, Utils>>())
      .def("get_row", &TinyMatrix3x3<Scalar, Utils>::getRow)
      .def("set_identity", &TinyMatrix3x3<Scalar, Utils>::set_identity);

  py::class_<TinyQuaternion<Scalar, Utils>>(m, "TinyQuaternion")
      .def(py::init<Scalar, Scalar, Scalar, Scalar>())
      .def("set_identity", &TinyQuaternion<Scalar, Utils>::set_identity)
      // .def("get_euler_rpy", &TinyQuaternion<Scalar, Utils>::get_euler_rpy)
      .def("get_euler_rpy2",
           &TinyQuaternion<Scalar, Utils>::get_euler_rpy2)
      .def("set_euler_rpy", &TinyQuaternion<Scalar, Utils>::set_euler_rpy)

      .def_readwrite("x", &TinyQuaternion<Scalar, Utils>::m_x)
      .def_readwrite("y", &TinyQuaternion<Scalar, Utils>::m_y)
      .def_readwrite("z", &TinyQuaternion<Scalar, Utils>::m_z)
      .def_readwrite("w", &TinyQuaternion<Scalar, Utils>::m_w)
      .def("__repr__",
           [](const TinyQuaternion<Scalar, Utils> &q) {
             return "[" + to_string(q.m_x) + " " + to_string(q.m_y) +
                    " " + to_string(q.m_z) + " " + to_string(q.m_w) +
                    "]";
           })
      .def("__getitem__", [](const TinyQuaternion<Scalar, Utils> &a,
                             int i) { return a[i]; })
      .def("__setitem__", [](TinyQuaternion<Scalar, Utils> &a, int i,
                             Scalar v) { a[i] = v; });

  py::class_<TinyPose<Scalar, Utils>>(m, "TinyPose")
      .def(py::init<TinyVector3<Scalar, Utils>,
                    TinyQuaternion<Scalar, Utils>>())
      .def_readwrite("position", &TinyPose<Scalar, Utils>::m_position)
      .def_readwrite("orientation",
                     &TinyPose<Scalar, Utils>::m_orientation)
      .def("inverse_transform",
           &TinyPose<Scalar, Utils>::inverse_transform);

  py::class_<TinySpatialTransform<Scalar, Utils>>(m,
                                                        "TinySpatialTransform")
      .def(py::init<>())
      .def("set_identity",
           &TinySpatialTransform<Scalar, Utils>::set_identity)
      .def_readwrite("translation",
                     &TinySpatialTransform<Scalar, Utils>::m_translation)
      .def_readwrite("rotation",
                     &TinySpatialTransform<Scalar, Utils>::m_rotation)
      .def(py::self * py::self)
      .def("get_inverse",
           &TinySpatialTransform<Scalar, Utils>::get_inverse);

  py::class_<TinySpatialMotionVector<Scalar, Utils>>(
      m, "TinySpatialMotionVector")
      .def(py::init<int>())
      .def_readwrite("topVec",
                     &TinySpatialMotionVector<Scalar, Utils>::m_topVec)
      .def_readwrite(
          "bottomVec",
          &TinySpatialMotionVector<Scalar, Utils>::m_bottomVec);

  py::class_<TinySymmetricSpatialDyad<Scalar, Utils>>(
      m, "TinySymmetricSpatialDyad")
      .def(py::init<>())
      .def("set_identity",
           &TinySymmetricSpatialDyad<Scalar, Utils>::setIdentity)
      .def("compute_inertia_dyad",
           &TinySymmetricSpatialDyad<Scalar, Utils>::computeInertiaDyad)
      .def("mul", &TinySymmetricSpatialDyad<Scalar, Utils>::mul)
      .def("shift", &TinySymmetricSpatialDyad<Scalar, Utils>::shift)
      .def("inverse", &TinySymmetricSpatialDyad<Scalar, Utils>::inverse)
      .def_readwrite(
          "topLeftMat",
          &TinySymmetricSpatialDyad<Scalar, Utils>::m_topLeftMat)
      .def_readwrite(
          "topRightMat",
          &TinySymmetricSpatialDyad<Scalar, Utils>::m_topRightMat)
      .def_readwrite(
          "bottomLeftMat",
          &TinySymmetricSpatialDyad<Scalar, Utils>::m_bottomLeftMat)
      .def_readwrite(
          "bottomRightMat",
          &TinySymmetricSpatialDyad<Scalar, Utils>::m_bottomRightMat)
      .def_readwrite(
          "center_of_mass",
          &TinySymmetricSpatialDyad<Scalar, Utils>::m_center_of_mass)
      .def(py::self -= py::self);

  py::enum_<TinyJointType>(m, "TinyJointType")
      .value("JOINT_FIXED", JOINT_FIXED, "JOINT_FIXED")
      .value("JOINT_PRISMATIC_X", JOINT_PRISMATIC_X, "JOINT_PRISMATIC_X")
      .value("JOINT_PRISMATIC_Y", JOINT_PRISMATIC_Y, "JOINT_PRISMATIC_Y")
      .value("JOINT_PRISMATIC_Z", JOINT_PRISMATIC_Z, "JOINT_PRISMATIC_Z")
      .value("JOINT_PRISMATIC_AXIS", JOINT_PRISMATIC_AXIS,
             "JOINT_PRISMATIC_AXIS")
      .value("JOINT_REVOLUTE_X", JOINT_REVOLUTE_X, "JOINT_REVOLUTE_X")
      .value("JOINT_REVOLUTE_Y", JOINT_REVOLUTE_Y, "JOINT_REVOLUTE_Y")
      .value("JOINT_REVOLUTE_Z", JOINT_REVOLUTE_Z, "JOINT_REVOLUTE_Z")
      .value("JOINT_REVOLUTE_AXIS", JOINT_REVOLUTE_AXIS, "JOINT_REVOLUTE_AXIS")
      .value("JOINT_INVALID", JOINT_INVALID, "JOINT_INVALID")
      .export_values();
          
  py::class_<TinyLink<Scalar, Utils>,
             std::unique_ptr<TinyLink<Scalar, Utils>>>(m, "TinyLink")
      .def(py::init<TinyJointType, TinySpatialTransform<Scalar, Utils> &,
                    const TinySymmetricSpatialDyad<Scalar, Utils> &>())
      .def("jcalc", &TinyLink<Scalar, Utils>::jcalc1)
      .def("set_joint_type", &TinyLink<Scalar, Utils>::set_joint_type)
      .def_readwrite("stiffness", &TinyLink<Scalar, Utils>::m_stiffness)
      .def_readwrite("joint_type", &TinyLink<Scalar, Utils>::m_joint_type)
      .def_readwrite("damping", &TinyLink<Scalar, Utils>::m_damping);
                                        
  py::class_<TinyMultiBody<Scalar, Utils>,
             std::unique_ptr<TinyMultiBody<Scalar, Utils>>>(
      m, "TinyMultiBody")
      .def(py::init<bool>(), py::return_value_policy::reference)
      .def("initialize", &TinyMultiBody<Scalar, Utils>::initialize)
      .def("set_position",
           &TinyMultiBody<Scalar, Utils>::set_position)
      .def("get_world_transform",
           &TinyMultiBody<Scalar, Utils>::get_world_transform)
      .def("mass_matrix", &TinyMultiBody<Scalar, Utils>::mass_matrix1)
      .def("attach_link", &TinyMultiBody<Scalar, Utils>::attach_link)
      .def("forward_kinematics",
           &TinyMultiBody<Scalar, Utils>::forward_kinematics1)
      .def("forward_dynamics", 
           py::overload_cast<const TinyVector3<Scalar, Utils> &>(
               &TinyMultiBody<Scalar, Utils>::forward_dynamics))
      .def("integrate", py::overload_cast<Scalar>(
                            &TinyMultiBody<Scalar, Utils>::integrate))
      .def("integrate_q", &TinyMultiBody<Scalar, Utils>::integrate_q)
      .def("dof", &TinyMultiBody<Scalar, Utils>::dof)
      .def("dof_qd", &TinyMultiBody<Scalar, Utils>::dof_qd)
      .def("dof_u", &TinyMultiBody<Scalar, Utils>::dof_u)
      .def("dof_state", &TinyMultiBody<Scalar, Utils>::dof_state)
      .def("body_to_world", &TinyMultiBody<Scalar, Utils>::body_to_world)
      .def("world_to_body", &TinyMultiBody<Scalar, Utils>::world_to_body)
      .def("point_jacobian",
           &TinyMultiBody<Scalar, Utils>::point_jacobian1)
      // .def("bias_forces", &TinyMultiBody<Scalar, Utils>::bias_forces)
      .def_property_readonly("num_dofs", &TinyMultiBody<Scalar, Utils>::dof)
      .def_property_readonly("num_dofs_qd", &TinyMultiBody<Scalar, Utils>::dof_qd)
      .def_readwrite("q", &TinyMultiBody<Scalar, Utils>::m_q)
      .def_readwrite("links", &TinyMultiBody<Scalar, Utils>::m_links)
      .def_readwrite("qd", &TinyMultiBody<Scalar, Utils>::m_qd)
      .def_readwrite("qdd", &TinyMultiBody<Scalar, Utils>::m_qdd)
      .def_readwrite("tau", &TinyMultiBody<Scalar, Utils>::m_tau)
      .def_readwrite("isFloating", &TinyMultiBody<Scalar, Utils>::m_isFloating)
      .def_readwrite("adjm_this_state", &TinyMultiBody<Scalar, Utils>::adjm_this_state);

  py::class_<TinyCollisionDispatcher<Scalar, Utils>>(
      m, "TinyCollisionDispatcher")
      .def(py::init<>())
      .def("compute_contacts", 
           &TinyCollisionDispatcher<Scalar, Utils>::compute_contacts);

  py::class_<TinyContactPoint<Scalar, Utils>> contact(m,
                                                            "TinyContactPoint");
  contact.def(py::init<>())
      .def_readwrite(
          "world_normal_on_b",
          &TinyContactPoint<Scalar, Utils>::m_world_normal_on_b)
      .def_readwrite("world_point_on_a",
                     &TinyContactPoint<Scalar, Utils>::m_world_point_on_a)
      .def_readwrite("world_point_on_b",
                     &TinyContactPoint<Scalar, Utils>::m_world_point_on_b)
      .def_readwrite("distance",
                     &TinyContactPoint<Scalar, Utils>::m_distance);

  py::class_<TinyContactPointRigidBody<Scalar, Utils>>(
      m, "TinyContactPointRigidBody", contact)
      .def(py::init<>())
      .def_readwrite(
          "rigid_body_a",
          &TinyContactPointRigidBody<Scalar, Utils>::m_rigid_body_a)
      .def_readwrite(
          "rigid_body_b",
          &TinyContactPointRigidBody<Scalar, Utils>::m_rigid_body_b)
      .def_readwrite(
          "restitution",
          &TinyContactPointRigidBody<Scalar, Utils>::m_restitution)
      .def_readwrite(
          "friction",
          &TinyContactPointRigidBody<Scalar, Utils>::m_friction);

  py::class_<TinyContactPointMultiBody<Scalar, Utils>>(
      m, "TinyContactPointMultiBody", contact)
      .def(py::init<>())
      .def_readwrite(
          "multi_body_a",
          &TinyContactPointMultiBody<Scalar, Utils>::m_multi_body_a)
      .def_readwrite(
          "multi_body_b",
          &TinyContactPointMultiBody<Scalar, Utils>::m_multi_body_b)
      .def_readwrite(
          "restitution",
          &TinyContactPointMultiBody<Scalar, Utils>::m_restitution)
      .def_readwrite(
          "friction",
          &TinyContactPointMultiBody<Scalar, Utils>::m_friction)
      .def_readwrite("link_a",
                     &TinyContactPointMultiBody<Scalar, Utils>::m_link_a)
      .def_readwrite("link_b",
                     &TinyContactPointMultiBody<Scalar, Utils>::m_link_b);

  py::class_<TinyConstraintSolver<Scalar, Utils>>(m,
                                                        "TinyConstraintSolver")
      .def(py::init<>())
      .def("resolve_collision",
           &TinyConstraintSolver<Scalar, Utils>::resolveCollision);

  py::class_<TinyMultiBodyConstraintSolver<Scalar, Utils>>(
      m, "TinyMultiBodyConstraintSolver")
      .def(py::init<>())
      .def("resolve_collision",
           &TinyMultiBodyConstraintSolver<Scalar,
                                          Utils>::resolveCollision);
  py::class_<TinyUrdfParser<Scalar, Utils>>(m, "TinyUrdfParser")
      .def(py::init<>())
      .def("load_urdf", &TinyUrdfParser<Scalar, Utils>::load_urdf);

  // py::enum_<TinyVelocitySmoothingMethod>(m, "TinyVelocitySmoothingMethod",
  //                                        py::arithmetic())
  //     .value("SMOOTH_VEL_NONE", SMOOTH_VEL_NONE)
  //     .value("SMOOTH_VEL_SIGMOID", SMOOTH_VEL_SIGMOID)
  //     .value("SMOOTH_VEL_TANH", SMOOTH_VEL_TANH)
  //     .value("SMOOTH_VEL_ABS", SMOOTH_VEL_ABS)
  //     .export_values();

  // typedef TinyMultiBodyConstraintSolverSpring<Scalar, Utils> TMBCSS;
  // py::class_<TMBCSS>(m, "TinyMultiBodyConstraintSolverSpring")
  //     .def(py::init<>())
  //     .def("resolve_collision", &TMBCSS::resolveCollision)
  //     .def_readwrite("spring_k", &TMBCSS::spring_k)
  //     .def_readwrite("damper_d", &TMBCSS::damper_d)
  //     .def_readwrite("hard_contact_condition", &TMBCSS::hard_contact_condition)
  //     .def_readwrite("exponent_n", &TMBCSS::exponent_n)
  //     .def_readwrite("exponent_n_air", &TMBCSS::exponent_n_air)
  //     .def_readwrite("exponent_vel_air", &TMBCSS::exponent_vel_air)
  //     .def_readwrite("smoothing_method", &TMBCSS::smoothing_method)
  //     .def_readwrite("smooth_alpha_vel", &TMBCSS::smooth_alpha_vel)
  //     .def_readwrite("smooth_alpha_normal", &TMBCSS::smooth_alpha_normal)
  //     .def_readwrite("mu_static", &TMBCSS::mu_static)
  //     .def_readwrite("andersson_vs", &TMBCSS::andersson_vs)
  //     .def_readwrite("andersson_p", &TMBCSS::andersson_p)
  //     .def_readwrite("andersson_ktanh", &TMBCSS::andersson_ktanh)
  //     .def_readwrite("v_transition", &TMBCSS::v_transition)
  //     .def_readwrite("friction_model", &TMBCSS::friction_model)
  //     .def("compute_contact_force", &TMBCSS::compute_contact_force)
  //     .def("compute_friction_force", &TMBCSS::compute_friction_force);
 
  py::class_<TinyWorld<Scalar, Utils>>(m, "TinyWorld")
      .def(py::init<Scalar, bool>(), py::return_value_policy::reference, py::arg("gravity_z")=Utils::fraction(-98, 10),
        py::arg("do_vis")=true)
      .def("step", &TinyWorld<Scalar, Utils>::step)
      .def_property("gravity", &TinyWorld<Scalar, Utils>::get_gravity,
                    &TinyWorld<Scalar, Utils>::set_gravity)
      .def("compute_contacts_rigid_body",
           &TinyWorld<Scalar, Utils>::compute_contacts_rigid_body)
      .def("compute_contacts_multi_body",
           &TinyWorld<Scalar, Utils>::compute_contacts_multi_body)
      .def("create_multi_body", 
           &TinyWorld<Scalar, Utils>::create_multi_body)
      .def("get_collision_dispatcher",
           &TinyWorld<Scalar, Utils>::get_collision_dispatcher)
      .def("adj_initialize",
           &TinyWorld<Scalar, Utils>::adj_initialize)
      .def("adj_step",
           &TinyWorld<Scalar, Utils>::adj_step)
      .def("sync_visual_meshcat",  
           &TinyWorld<Scalar, Utils>::sync_visual_meshcat)
      .def("do_convert_visuals",
           &TinyWorld<Scalar, Utils>::do_convert_visuals)
      .def("sync",
           &TinyWorld<Scalar, Utils>::sync)
      .def_readwrite("R_friction", 
                     &TinyWorld<Scalar, Utils>::adjm_Rfriction)
      .def_readwrite("R_restitution", 
                     &TinyWorld<Scalar, Utils>::adjm_Rrestituition)
      .def_readwrite("friction", 
                     &TinyWorld<Scalar, Utils>::default_friction)
      .def_readwrite("restitution",
                     &TinyWorld<Scalar, Utils>::default_restitution)
      .def_readwrite("dt",
                     &TinyWorld<Scalar, Utils>::dt)
      .def_readwrite("m_multi_bodies",
                     &TinyWorld<Scalar, Utils>::m_multi_bodies, py::return_value_policy::reference);

  py::class_<TinyRaycastResult<Scalar, Utils>>(m, "TinyRaycastResult")
      .def(py::init<>())
      .def_readwrite("hit_fraction",
                     &TinyRaycastResult<Scalar, Utils>::m_hit_fraction)
      .def_readwrite("collider_index",
                     &TinyRaycastResult<Scalar, Utils>::m_collider_index);

  py::class_<TinyRaycast<Scalar, Utils>>(m, "TinyRaycast")
      .def(py::init<>())
      .def("cast_rays", &TinyRaycast<Scalar, Utils>::cast_rays)
      .def("volume", &TinyRaycast<Scalar, Utils>::volume)
      .def("intersection_volume",
           &TinyRaycast<Scalar, Utils>::intersection_volume);

  ///////////////////////////////////////////////////////////////////////////////////////////

  py::class_<Fix64Scalar>(m, "Fix64Scalar")
      .def(py::init<>())
      .def("zero", &Fix64Scalar::zero)
      .def("fraction", &Fix64Scalar::fraction_internal)
      .def("getScalar", &Fix64Scalar::getScalar)
      .def("sqrt", &Fix64Scalar::sqrt1)
      .def("sin", &Fix64Scalar::sin1)
      .def("cos", &Fix64Scalar::cos1)
      .def("abs", &Fix64Scalar::Fix64Abs)
      .def("max", &Fix64Scalar::Fix64Max)
      .def("half", &Fix64Scalar::half)
      .def("pi", &Fix64Scalar::pi)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * py::self)
      .def(py::self / py::self)
      .def(py::self < py::self)
      .def(py::self == py::self)
      .def("__repr__",
           [](const Fix64Scalar &a) { return to_string(a.getScalar()); });
  py::class_<TinyVector3<Fix64Scalar, Fix64Scalar>>(m, "Fix64Vector3")
      .def(py::init<Fix64Scalar, Fix64Scalar, Fix64Scalar>())
      .def("set_zero", &TinyVector3<Fix64Scalar, Fix64Scalar>::set_zero)
      .def_readwrite("x", &TinyVector3<Fix64Scalar, Fix64Scalar>::m_x)
      .def_readwrite("y", &TinyVector3<Fix64Scalar, Fix64Scalar>::m_y)
      .def_readwrite("z", &TinyVector3<Fix64Scalar, Fix64Scalar>::m_z)
      .def(py::self + py::self)
      .def(py::self - py::self)

      .def("__repr__", [](const TinyVector3<Fix64Scalar, Fix64Scalar> &a) {
        Scalar x = Fix64Scalar::getDouble(a.m_x);
        Scalar y = Fix64Scalar::getDouble(a.m_y);
        Scalar z = Fix64Scalar::getDouble(a.m_z);

        return "[" + to_string(x) + "," + to_string(y) + "," +
               to_string(z) + "]";
      });

  py::class_<TinyUrdfInertial<Scalar, Utils>>(m, "TinyUrdfInertial")
      .def(py::init<>())
      .def_readwrite("mass", &TinyUrdfInertial<Scalar, Utils>::mass)
      .def_readwrite("inertia_xxyyzz",
                     &TinyUrdfInertial<Scalar, Utils>::inertia_xxyyzz)
      .def_readwrite("origin_xyz",
                     &TinyUrdfInertial<Scalar, Utils>::origin_xyz)
      .def_readwrite("origin_rpy",
                     &TinyUrdfInertial<Scalar, Utils>::origin_rpy);

  py::class_<TinyUrdfCollisionSphere<Scalar, Utils>>(
      m, "TinyUrdfCollisionSphere")
      .def(py::init<>())
      .def_readwrite("radius",
                     &TinyUrdfCollisionSphere<Scalar, Utils>::m_radius);

  py::class_<TinyUrdfCollisionBox<Scalar, Utils>>(m,
                                                        "TinyUrdfCollisionBox")
      .def(py::init<>())
      .def_readwrite("extents",
                     &TinyUrdfCollisionBox<Scalar, Utils>::m_extents);

  py::class_<TinyUrdfCollisionCapsule<Scalar, Utils>>(
      m, "TinyUrdfCollisionCapsule")
      .def(py::init<>())
      .def_readwrite("radius",
                     &TinyUrdfCollisionCapsule<Scalar, Utils>::m_radius)
      .def_readwrite("length",
                     &TinyUrdfCollisionCapsule<Scalar, Utils>::m_length);

  py::class_<TinyUrdfCollisionPlane<Scalar, Utils>>(
      m, "TinyUrdfCollisionPlane")
      .def(py::init<>())
      .def_readwrite("constant",
                     &TinyUrdfCollisionPlane<Scalar, Utils>::m_constant)
      .def_readwrite("normal",
                     &TinyUrdfCollisionPlane<Scalar, Utils>::m_normal);

  py::class_<TinyUrdfCollisionMesh<Scalar, Utils>>(
      m, "TinyUrdfCollisionMesh")
      .def(py::init<>())
      .def_readwrite("file_name",
                     &TinyUrdfCollisionMesh<Scalar, Utils>::m_file_name)
      .def_readwrite("scale",
                     &TinyUrdfCollisionMesh<Scalar, Utils>::m_scale);

  py::class_<TinyUrdfGeometry<Scalar, Utils>>(m, "TinyUrdfGeometry")
      .def(py::init<>())
      .def_readwrite("geom_type",
                     &TinyUrdfGeometry<Scalar, Utils>::geom_type)
      .def_readwrite("sphere", &TinyUrdfGeometry<Scalar, Utils>::m_sphere)
      .def_readwrite("box", &TinyUrdfGeometry<Scalar, Utils>::m_box)
      .def_readwrite("plane", &TinyUrdfGeometry<Scalar, Utils>::m_plane)
      .def_readwrite("capsule",
                     &TinyUrdfGeometry<Scalar, Utils>::m_capsule)
      .def_readwrite("mesh", &TinyUrdfGeometry<Scalar, Utils>::m_mesh);

  py::class_<TinyUrdfCollision<Scalar, Utils>>(m, "TinyUrdfCollision")
      .def(py::init<>())

      .def_readwrite("origin_xyz",
                     &TinyUrdfCollision<Scalar, Utils>::origin_xyz)
      .def_readwrite("origin_rpy",
                     &TinyUrdfCollision<Scalar, Utils>::origin_rpy)
      .def_readwrite("geometry",
                     &TinyUrdfCollision<Scalar, Utils>::geometry);

  py::class_<TinyUrdfVisual<Scalar, Utils>>(m, "TinyUrdfVisual")
      .def(py::init<>())
      .def_readwrite("origin_xyz",
                     &TinyUrdfVisual<Scalar, Utils>::origin_xyz)
      .def_readwrite("origin_rpy",
                     &TinyUrdfVisual<Scalar, Utils>::origin_rpy)
      .def_readwrite("geometry", &TinyUrdfVisual<Scalar, Utils>::geometry)
      .def_readwrite(
          "sync_visual_body_uid1",
          &TinyUrdfVisual<Scalar, Utils>::sync_visual_body_uid1)
      .def_readwrite(
          "sync_visual_body_uid2",
          &TinyUrdfVisual<Scalar, Utils>::sync_visual_body_uid2);

  py::class_<TinyUrdfJoint<Scalar, Utils>>(m, "TinyUrdfJoint")
      .def(py::init<>())

      .def_readwrite("joint_name",
                     &TinyUrdfJoint<Scalar, Utils>::joint_name)
      .def_readwrite("joint_type",
                     &TinyUrdfJoint<Scalar, Utils>::joint_type)
      .def_readwrite("joint_lower_limit",
                     &TinyUrdfJoint<Scalar, Utils>::joint_lower_limit)
      .def_readwrite("joint_upper_limit",
                     &TinyUrdfJoint<Scalar, Utils>::joint_upper_limit)
      .def_readwrite("parent_name",
                     &TinyUrdfJoint<Scalar, Utils>::parent_name)
      .def_readwrite("child_name",
                     &TinyUrdfJoint<Scalar, Utils>::child_name)
      .def_readwrite("joint_origin_xyz",
                     &TinyUrdfJoint<Scalar, Utils>::joint_origin_xyz)
      .def_readwrite("joint_origin_rpy",
                     &TinyUrdfJoint<Scalar, Utils>::joint_origin_rpy)
      .def_readwrite("joint_axis_xyz",
                     &TinyUrdfJoint<Scalar, Utils>::joint_axis_xyz);

  py::class_<TinyUrdfLink<Scalar, Utils>>(m, "TinyUrdfLink")
      .def(py::init<>())
      .def_readwrite("link_name", &TinyUrdfLink<Scalar, Utils>::link_name)
      .def_readwrite("urdf_inertial",
                     &TinyUrdfLink<Scalar, Utils>::urdf_inertial)
      .def_readwrite("urdf_visual_shapes",
                     &TinyUrdfLink<Scalar, Utils>::urdf_visual_shapes)
      .def_readwrite("urdf_collision_shapes",
                     &TinyUrdfLink<Scalar, Utils>::urdf_collision_shapes)
      .def_readwrite("parent_index",
                     &TinyUrdfLink<Scalar, Utils>::m_parent_index);

  py::class_<TinyUrdfStructures<Scalar, Utils>>(m, "TinyUrdfStructures")
      .def(py::init<>(), py::return_value_policy::reference)
      .def_readwrite("robot_name",
                     &TinyUrdfStructures<Scalar, Utils>::m_robot_name)
      .def_readwrite("base_links",
                     &TinyUrdfStructures<Scalar, Utils>::m_base_links)
      .def_readwrite("links", &TinyUrdfStructures<Scalar, Utils>::m_links)
      .def_readwrite("joints",
                     &TinyUrdfStructures<Scalar, Utils>::m_joints)
      .def_readwrite(
          "name_to_link_index",
          &TinyUrdfStructures<Scalar, Utils>::m_name_to_link_index);

  py::class_<UrdfToMultiBody2<Scalar, Utils>>(m, "UrdfToMultiBody2")
      .def(py::init<>(), py::return_value_policy::reference)
      .def("convert2", &UrdfToMultiBody2<Scalar, Utils>::convert);

  py::class_<Motion>(m, "Motion")
      .def(py::init([](const std::string &filename) {
        Motion motion;
        Motion::load_from_file(filename, &motion);
        return motion;
      }))
      .def_readwrite("frames", &Motion::frames)
      .def_readonly("frame_duration", &Motion::frame_duration)
      .def_property_readonly("total_duration", &Motion::total_duration)
      .def("calculate_frame", &Motion::calculate_frame);

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev"; 
#endif
  m.attr("SPHERE_TYPE") = py::int_(int(TINY_SPHERE_TYPE));
  m.attr("BOX_TYPE") = py::int_(int(TINY_BOX_TYPE));
  m.attr("PLANE_TYPE") = py::int_(int(TINY_PLANE_TYPE));
  m.attr("CAPSULE_TYPE") = py::int_(int(TINY_CAPSULE_TYPE));
  m.attr("MESH_TYPE") = py::int_(int(TINY_MESH_TYPE));

  m.attr("JOINT_FIXED") = py::int_(int(JOINT_FIXED));
  m.attr("JOINT_PRISMATIC_X") = py::int_(int(JOINT_PRISMATIC_X));
  m.attr("JOINT_PRISMATIC_Y") = py::int_(int(JOINT_PRISMATIC_Y));
  m.attr("JOINT_PRISMATIC_Z") = py::int_(int(JOINT_PRISMATIC_Z));
  m.attr("JOINT_PRISMATIC_AXIS") = py::int_(int(JOINT_PRISMATIC_AXIS));
  m.attr("JOINT_REVOLUTE_X") = py::int_(int(JOINT_REVOLUTE_X));
  m.attr("JOINT_REVOLUTE_Y") = py::int_(int(JOINT_REVOLUTE_Y));
  m.attr("JOINT_REVOLUTE_Z") = py::int_(int(JOINT_REVOLUTE_Z));
  m.attr("JOINT_REVOLUTE_AXIS") = py::int_(int(JOINT_REVOLUTE_AXIS));
}
