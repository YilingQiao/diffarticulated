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

#ifndef TINY_WORLD_H
#define TINY_WORLD_H

#include <string>
#include <vector>

#include "tiny_constraint_solver.h"
#include "tiny_geometry.h"
#include "tiny_mb_constraint_solver.h"
#include "tiny_multi_body.h"
#include "tiny_rigid_body.h"
#include "tiny_metrics.h"
// #include "examples/pybullet_urdf_import.h"
// #include "examples/pybullet_visualizer_api.h"


#include "examples/meshcat_urdf_visualizer.h"
// typedef PyBulletVisualizerAPI VisualizerAPI;

template <typename TinyScalar, typename TinyConstants>
class TinyWorld {
  typedef ::TinyRigidBody<TinyScalar, TinyConstants> TinyRigidBody;
  typedef ::TinyMultiBody<TinyScalar, TinyConstants> TinyMultiBody;
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinyGeometry<TinyScalar, TinyConstants> TinyGeometry;
  typedef ::TinySpatialTransform<TinyScalar, TinyConstants>
      TinySpatialTransform;
  typedef ::TinySpatialMotionVector<TinyScalar, TinyConstants>
      TinySpatialMotionVector;

  typedef ::TinyCapsule<TinyScalar, TinyConstants> TinyCapsule;
  typedef ::TinySphere<TinyScalar, TinyConstants> TinySphere;
  typedef ::TinyPlane<TinyScalar, TinyConstants> TinyPlane;
  typedef ::TinyMatrixXxX<TinyScalar, TinyConstants> TinyMatrixXxX;
  typedef ::TinyVectorX<TinyScalar, TinyConstants> TinyVectorX;
  typedef ::TinyMetric<TinyScalar, TinyConstants> TinyMetric;
  typedef ::TinyUrdfStructures<TinyScalar, TinyConstants> TinyUrdfStructures;
  typedef ::TinyBox<TinyScalar, TinyConstants> TinyBox;


  

  std::vector<TinyRigidBody*> m_bodies;

  TinyVector3 m_gravity_acceleration;

  std::vector<TinyGeometry*> m_geoms;

  TinyCollisionDispatcher<TinyScalar, TinyConstants> m_dispatcher;

 public:
  MeshcatUrdfVisualizer<TinyScalar, TinyConstants> meshcat_viz;
  // VisualizerAPI *vis;
  TinyUrdfStructures cur_stu;
  std::vector<TinyMultiBody*> m_multi_bodies;
  TinySubmitProfileTiming m_profileTimingFunc{nullptr};
  TinyConstraintSolver<TinyScalar, TinyConstants>* m_constraint_solver{nullptr};
  TinyMultiBodyConstraintSolver<TinyScalar, TinyConstants>*
      m_mb_constraint_solver{nullptr};

  int m_num_solver_iterations{50};

  // quantities related to adjoint method
  int adjm_step_count = 0;
  int adjm_total_num_state = 0;
  std::vector<std::vector<TinyScalar> > adjm_dphi_du;
  TinyScalar dt;

  std::vector<TinyScalar> adjm_global_ax, adjm_global_ay; 

  // contact settings
  TinyScalar default_friction{TinyConstants::fraction(2, 10)};
  TinyScalar default_restitution{TinyConstants::zero()};
  TinyScalar adjm_Rfriction{TinyConstants::zero()};
  TinyScalar adjm_Rrestituition{TinyConstants::zero()};

  void sync() {
    // for (auto *mb : m_multi_bodies)
    //   PyBulletUrdfImport<TinyScalar, TinyConstants>::sync_graphics_transforms(mb, *vis);
  }

  void do_convert_visuals(TinyUrdfStructures stu) {
    // cur_stu = stu;
    // PyBulletUrdfImport<TinyScalar, TinyConstants>::convert_visuals(cur_stu, cur_stu.m_base_links[0], *vis);
    // for (int i = 0; i < cur_stu.m_links.size(); ++i)
    //   PyBulletUrdfImport<TinyScalar, TinyConstants>::convert_visuals(cur_stu, cur_stu.m_links[i], *vis);
  }

  explicit TinyWorld(TinyScalar gravity_z = TinyConstants::fraction(-98, 10), bool do_vis=true)
      : m_gravity_acceleration(TinyConstants::zero(), TinyConstants::zero(),
                               gravity_z),
        m_constraint_solver(
            new TinyConstraintSolver<TinyScalar, TinyConstants>),
        m_mb_constraint_solver(
            new TinyMultiBodyConstraintSolver<TinyScalar, TinyConstants>) {
          if (do_vis) {
            // vis = new VisualizerAPI;
            // vis->connect(1); //1:gui; 2:direct; 3:shared
            // vis->setAdditionalSearchPath("/scratch1/diffarti/data");
            // vis->resetSimulation();

            meshcat_viz.delete_all();
            // MeshcatUrdfVisualizer<double, DoubleUtils> meshcat_viz;
          }
        }

  inline void submitProfileTiming(const std::string& name) {
    if (m_profileTimingFunc) {
      m_profileTimingFunc(name);
    }
  }
  virtual ~TinyWorld() { 
    clear(); 
  }

  void clear() {
    for (int i = 0; i < m_geoms.size(); i++) {
      delete m_geoms[i];
    }
    m_geoms.clear();

    for (int i = 0; i < m_bodies.size(); i++) {
      delete m_bodies[i];
    }
    m_bodies.clear();

    for (int i = 0; i < m_multi_bodies.size(); i++) {
      delete m_multi_bodies[i];
    }
    m_multi_bodies.clear();

    if (m_constraint_solver) {
      delete m_constraint_solver;
      m_constraint_solver = nullptr;
    }
  }

  const TinyVector3& get_gravity() const { return m_gravity_acceleration; }

  void set_gravity(const TinyVector3& gravity) {
    m_gravity_acceleration = gravity;
  }

  TinyConstraintSolver<TinyScalar, TinyConstants>* get_constraint_solver() {
    return m_constraint_solver;
  }

  TinyCapsule* create_capsule(TinyScalar radius, TinyScalar length) {
    TinyCapsule* capsule = new TinyCapsule(radius, length);
    m_geoms.push_back(capsule);
    return capsule;
  }

  TinyPlane* create_plane() {
    TinyPlane* plane = new TinyPlane();
    m_geoms.push_back(plane);
    return plane;
  }

  TinySphere* create_sphere(TinyScalar radius) {
    TinySphere* sphere = new TinySphere(radius);
    m_geoms.push_back(sphere);
    return sphere;
  }

  TinyBox* create_box(TinyVector3 extents) {
    TinyBox* box = new TinyBox(extents);
    m_geoms.push_back(box);
    return box;
  }

  TinyCollisionDispatcher<TinyScalar, TinyConstants>
  get_collision_dispatcher() {
    return m_dispatcher;
  }

  TinyRigidBody* create_rigid_body(TinyScalar mass, const TinyGeometry* geom) {
    TinyRigidBody* body = new TinyRigidBody(mass, geom);
    this->m_bodies.push_back(body);
    return body;
  }

  TinyMultiBody* create_multi_body() {
    TinyMultiBody* body = new TinyMultiBody();
    this->m_multi_bodies.push_back(body);
    return body;
  }

  std::vector<TinyContactPointRigidBody<TinyScalar, TinyConstants>>
      m_allContacts;
  std::vector<std::vector<TinyContactPointMultiBody<TinyScalar, TinyConstants>>>
      m_allMultiBodyContacts;

  std::vector<const TinyGeometry*> vec_geom_a, vec_geom_b;
  std::vector<TinyPose<TinyScalar, TinyConstants>> vec_pose_a, vec_pose_b;
  std::vector<int> vec_ii, vec_iii, vec_jj, vec_jjj;
  std::vector<TinyMultiBody*> vec_mb_b, vec_mb_a;

  std::vector<TinyContactPoint<TinyScalar, TinyConstants>> m_contacts;

  static void compute_contacts_rigid_body_internal(
      std::vector<TinyRigidBody*> bodies,
      TinyCollisionDispatcher<TinyScalar, TinyConstants>* dispatcher,
      std::vector<TinyContactPointRigidBody<TinyScalar, TinyConstants>>&
          contactsOut,
      const TinyScalar& restitution, const TinyScalar& friction) {
    std::vector<TinyContactPoint<TinyScalar, TinyConstants>> contacts;
    {
      for (int i = 0; i < bodies.size(); i++) {
        for (int j = i + 1; j < bodies.size(); j++) {
          contacts.reserve(1);
          contacts.resize(0);

          std::vector<TinyScalar> ax(0);
          int numContacts = dispatcher->computeContacts(
              bodies[i]->m_geometry, bodies[i]->m_world_pose,
              bodies[j]->m_geometry, bodies[j]->m_world_pose, contacts);
          for (int c = 0; c < numContacts; c++) {
            TinyContactPointRigidBody<TinyScalar, TinyConstants> rb_pt;
            TinyContactPoint<TinyScalar, TinyConstants>& pt = rb_pt;
            pt = contacts[c];
            rb_pt.m_rigid_body_a = bodies[i];
            rb_pt.m_rigid_body_b = bodies[j];
            // TODO(erwincoumans): combine friction and restitution based on
            // material properties of the two touching bodies
            rb_pt.m_restitution = restitution;
            rb_pt.m_friction = friction;
            contactsOut.push_back(rb_pt);
          }
        }
      }
    }
  }

  std::vector<TinyContactPointRigidBody<TinyScalar, TinyConstants>>
  compute_contacts_rigid_body(
      std::vector<TinyRigidBody*> bodies,
      TinyCollisionDispatcher<TinyScalar, TinyConstants>* dispatcher) {
    std::vector<TinyContactPointRigidBody<TinyScalar, TinyConstants>>
        contactsOut;
    compute_contacts_rigid_body_internal(bodies, dispatcher, contactsOut,
                                         default_restitution, default_friction);
    return contactsOut;
  }

  void compute_contacts_multi_body_internal(
      std::vector<TinyMultiBody*> multi_bodies,
      TinyCollisionDispatcher<TinyScalar, TinyConstants>* dispatcher,
      std::vector<
          std::vector<TinyContactPointMultiBody<TinyScalar, TinyConstants>>>&
          contacts_out,
      const TinyScalar& restitution, const TinyScalar& friction) {
    int num_multi_bodies = multi_bodies.size();

    vec_geom_a.clear();
    vec_geom_b.clear();
    vec_pose_a.clear();
    vec_pose_b.clear();
    vec_mb_b.clear();
    vec_mb_a.clear();
    vec_ii.clear();
    vec_iii.clear();
    vec_jj.clear();
    vec_jjj.clear();


    for (int i = 0; i < num_multi_bodies; i++) {
      TinyMultiBody* mb_a = multi_bodies[i];
      int num_links_a = mb_a->m_links.size();
      for (int j = i + 1; j < multi_bodies.size(); j++) {
        std::vector<TinyContactPoint<TinyScalar, TinyConstants>> contacts;
        TinyMultiBody* mb_b = multi_bodies[j];
        int num_links_b = mb_b->m_links.size();
        std::vector<TinyContactPointMultiBody<TinyScalar, TinyConstants>>
            contacts_ab;

        for (int ii = -1; ii < num_links_a; ii++) {
          const TinySpatialTransform& world_transform_a =
              mb_a->get_world_transform(ii);
          int num_geoms_a = mb_a->get_collision_geometries(ii).size();
          for (int iii = 0; iii < num_geoms_a; iii++) {
            const TinyGeometry* geom_a =
                mb_a->get_collision_geometries(ii)[iii];
            TinyPose<TinyScalar, TinyConstants> pose_a;
            const TinySpatialTransform& local_a =
                mb_a->get_collision_transforms(ii)[iii];
            TinySpatialTransform tr_a = world_transform_a * local_a;
            pose_a.m_position = tr_a.m_translation;
            tr_a.m_rotation.getRotation(pose_a.m_orientation);

            for (int jj = -1; jj < num_links_b; jj++) {
              const TinySpatialTransform& world_transform_b =
                  mb_b->get_world_transform(jj);
              int num_geoms_b = mb_b->get_collision_geometries(jj).size();
              for (int jjj = 0; jjj < num_geoms_b; jjj++) {
                const TinyGeometry* geom_b =
                    mb_b->get_collision_geometries(jj)[jjj];
                TinyPose<TinyScalar, TinyConstants> pose_b;
                const TinySpatialTransform& local_b =
                    mb_b->get_collision_transforms(jj)[jjj];

                TinySpatialTransform tr_b = world_transform_b * local_b;
                pose_b.m_position = tr_b.m_translation;
                tr_b.m_rotation.getRotation(pose_b.m_orientation);

                contacts.reserve(1);
                contacts.resize(0);

                bool is_tape = false;    

                int numContacts = dispatcher->computeContacts(
                    geom_a, pose_a, geom_b, pose_b, contacts);
                for (int c = 0; c < numContacts; c++) {
                  TinyContactPointMultiBody<TinyScalar, TinyConstants> mb_pt;
                  TinyContactPoint<TinyScalar, TinyConstants>& pt = mb_pt;
                  pt = contacts[c];
                  mb_pt.m_multi_body_a = multi_bodies[i];
                  mb_pt.m_multi_body_b = multi_bodies[j];
                  mb_pt.m_link_a = ii;
                  mb_pt.m_link_b = jj;
                  mb_pt.m_restitution = restitution;
                  mb_pt.m_friction = friction;
                  
                  contacts_ab.push_back(mb_pt);
                  vec_geom_a.push_back(geom_a);
                  vec_geom_b.push_back(geom_b);
                  vec_pose_a.push_back(pose_a);
                  vec_pose_b.push_back(pose_b);
                  vec_mb_b.push_back(mb_b);
                  vec_mb_a.push_back(mb_a);
                  vec_ii.push_back(ii);
                  vec_iii.push_back(iii);
                  vec_jj.push_back(jj);
                  vec_jjj.push_back(jjj);
                }
              }
            }
          }
        }

        contacts_out.push_back(contacts_ab);
      }
    }
  }


  void adj_compute_contacts_multi_body_internal(
      std::vector<TinyMultiBody*> multi_bodies,
      TinyCollisionDispatcher<TinyScalar, TinyConstants>* dispatcher,
      std::vector<
          std::vector<TinyContactPointMultiBody<TinyScalar, TinyConstants>>>&
          contacts_out,
      const TinyScalar& restitution, const TinyScalar& friction) {

    int ic = 0; 

    for (int i1 = 0; i1 < m_allMultiBodyContacts.size(); ++i1) {
        std::vector<TinyContactPointMultiBody<TinyScalar, TinyConstants>>
            contacts_ab = contacts_out[i1];

      // for (int i2 = contacts_ab.size() - 1; i2 >= 0; --i2)  {
      for (int i2 = 0; i2  < contacts_ab.size(); i2++)  {

        int ii = vec_ii[ic];
        int iii = vec_iii[ic];
        int jj = vec_jj[ic];
        int jjj = vec_jjj[ic];


        TinyContactPointMultiBody<TinyScalar, TinyConstants> mb_pt =
            contacts_ab[i2];
        adjm_Rfriction += mb_pt.adjm_Rfriction;
        adjm_Rrestituition += mb_pt.adjm_Rrestituition;

        dispatcher->adj_computeContacts(vec_geom_a[ic], vec_pose_a[ic],
                    vec_geom_b[ic], vec_pose_b[ic], mb_pt);
        TinySpatialTransform Rtr_a, Rtr_b, Rlocal_b, Rlocal_a,
          Rworld_transform_b, Rworld_transform_a;
        Rtr_a.set_zero();
        Rtr_b.set_zero();
        Rworld_transform_b.set_zero();
        Rworld_transform_a.set_zero();
        Rlocal_a.set_zero();
        Rlocal_b.set_zero();

        TinyMultiBody* mb_b = vec_mb_b[ic];
        TinyMultiBody* mb_a = vec_mb_a[ic];


        const TinySpatialTransform& world_transform_a =
                      mb_a->get_world_transform(ii);
        const TinySpatialTransform& local_a =
                mb_a->get_collision_transforms(ii)[iii];

        const TinySpatialTransform& world_transform_b =
                  mb_b->get_world_transform(jj);
        const TinySpatialTransform& local_b =
                  mb_b->get_collision_transforms(jj)[jjj];

        TinySpatialTransform tr_b = world_transform_b * local_b;
        TinySpatialTransform tr_a = world_transform_a * local_a;
        
        tr_b.m_rotation.adj_getRotation(vec_pose_b[ic].adjm_Rorientation, Rtr_b.m_rotation);
        Rtr_b.m_translation += vec_pose_b[ic].adjm_Rposition;
        world_transform_b.adj_st_multiply(Rtr_b, local_b,
              Rworld_transform_b, Rlocal_b);
        
        tr_a.m_rotation.adj_getRotation(vec_pose_a[ic].adjm_Rorientation, Rtr_a.m_rotation);
        Rtr_a.m_translation += vec_pose_a[ic].adjm_Rposition;

        world_transform_a.adj_st_multiply(Rtr_a, local_a,
              Rworld_transform_a, Rlocal_a);
        
        mb_b->adj_get_world_transform(jj) += Rworld_transform_b;
        // mb_b->adj_get_collision_transforms(jj)[jjj] += Rlocal_b;
        mb_a->adj_get_world_transform(ii) += Rworld_transform_a;
        // mb_a->adj_get_collision_transforms(ii)[iii] += Rlocal_a;
        ic++;
      }
    }

  }

  std::vector<std::vector<TinyContactPointMultiBody<TinyScalar, TinyConstants>>>
  compute_contacts_multi_body(
      std::vector<TinyMultiBody*> bodies,
      TinyCollisionDispatcher<TinyScalar, TinyConstants>* dispatcher) {
    std::vector<
        std::vector<TinyContactPointMultiBody<TinyScalar, TinyConstants>>>
        contactsOut;
    compute_contacts_multi_body_internal(bodies, dispatcher, contactsOut,
                                         default_restitution, default_friction);
    return contactsOut;
  }


  void adj_initialize(const TinyVector3& gravity,
                      int n_step, int dof_u) {
    TinySpatialMotionVector spatial_gravity(
        TinyVector3(TinyConstants::zero(), TinyConstants::zero(),
                    TinyConstants::zero()),
        gravity);

    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      // adjm_total_num_state += mb->dof_state();
      int this_dof_u = mb->dof_state() > 0 ? mb->m_dof : 0;
      mb->adj_initialize(gravity, n_step, this_dof_u);

      mb->adjm_this_state.resize(mb->dof_state());
      mb->adjm_next_state.resize(mb->dof_state());

      for (int cc = 0; cc < mb->dof(); cc++) 
          mb->adjm_this_state[cc] = mb->m_q[cc];
      for (int cc = 0; cc < mb->dof_qd(); cc++) 
          mb->adjm_this_state[cc+mb->dof()] = mb->m_qd[cc];
    }

    for (auto *mb : m_multi_bodies)
      adjm_dphi_du.push_back(std::vector<TinyScalar>(mb->dof_u()*n_step));
  }

  // void sync_visual_meshcat() {
  //   for (auto *mb : m_multi_bodies) {
  //     meshcat_viz.sync_visual_transforms(mb);
  //   }
  // }

  void sync_visual_meshcat(int step) {
    for (auto *mb : m_multi_bodies) {
      meshcat_viz.sync_visual_transforms(mb, step);
    }
  }

  void adj_man_constraints(const TinyVectorX& R,
                          TinyScalar dt){
    // forward
    {
      for (int i = 0; i < m_multi_bodies.size(); i++) {
        TinyMultiBody* mb = m_multi_bodies[i]; 
        for (int cc = 0; cc < mb->dof(); cc++) 
            mb->m_q[cc] = mb->adjm_next_state[cc];
        for (int cc = 0; cc < mb->dof_qd(); cc++) 
            mb->m_qd[cc] = mb->adjm_next_state[cc+mb->dof()];

      }
    }

    {
      m_allMultiBodyContacts.reserve(1024);
      m_allMultiBodyContacts.resize(0);
    }
    
    {
      submitProfileTiming("compute multi body contacts");
      compute_contacts_multi_body_internal(
          m_multi_bodies, &m_dispatcher, m_allMultiBodyContacts,
          default_restitution, default_friction);
      submitProfileTiming("");
    }
    // adjoint begin

    int curr_state = 0;
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      if (mb->dof()==0 || mb->dof_qd()==0)
        continue;
      mb->adj_constraint_set_zero();

      for (int cc = 0; cc < mb->dof(); cc++) 
          mb->adjm_Rq[cc] += R[curr_state+cc];
      for (int cc = 0; cc < mb->dof_qd(); cc++) 
          mb->adjm_Rqd[cc] += R[curr_state + cc + mb->dof()];


      mb->adj_b_integrate(dt);

      curr_state += mb->dof_state();

    }

    {
      submitProfileTiming("adj solve constraints");
      // use outer loop in case the multi-body constraint solver requires it
      // (e.g. sequential impulse method)
        for (int c = m_allMultiBodyContacts.size() - 1; c >= 0; c--) {
          bool is_tape = false;
          m_mb_constraint_solver->adj_resolveCollision(
                        m_allMultiBodyContacts[c], dt);
        }
      
      submitProfileTiming("");
    }
    

    {
      submitProfileTiming("adj compute multi body contacts");
      adj_compute_contacts_multi_body_internal(
          m_multi_bodies, &m_dispatcher, m_allMultiBodyContacts,
          default_restitution, default_friction);
      submitProfileTiming("");
    }

  }

  void adj_step_torch(std::vector<TinyScalar> state, TinyVectorX &u) {
    int curr_u = 0, curr_state = 0;
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      
      auto bodystate = std::vector<TinyScalar>(state.begin()+curr_state, state.begin()+curr_state+mb->dof_state());
      
      mb->adj_forward_dynamics(u.segment(curr_u, mb->dof_u()), 
                               bodystate,
                               mb->adjm_next_state, dt);

      curr_u += mb->dof_u();
      curr_state += mb->dof_state();

    }
    // change state
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];

      for (int cc = 0; cc < mb->dof_state(); cc++) 
        mb->adjm_this_state[cc] = mb->adjm_next_state[cc];
    }

    adj_forward_constraint(dt);

    // change state
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      for (int cc = 0; cc < mb->dof_state(); cc++) 
        mb->adjm_this_state[cc] = mb->adjm_next_state[cc];
    }

  }



  std::vector<TinyVectorX> adj_back_step_man(TinyVectorX& R, TinyVectorX &u, std::vector<TinyScalar> &state) {
    int curr_u = 0, curr_state = 0;
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      auto bodystate = std::vector<TinyScalar>(state.begin()+curr_state, state.begin()+curr_state+mb->dof_state());

      mb->adj_forward_dynamics(u.segment(curr_u, mb->dof_u()), 
                               bodystate,
                               mb->adjm_next_state, dt);
      curr_u += mb->dof_u();
      curr_state += mb->dof_state();
    }
    adj_man_constraints(R, dt);

    TinyVectorX pf_px_R(R.m_size);
    int ns = 0;
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      for (int cc = 0; cc < mb->dof(); cc++) 
          pf_px_R[ns++] = mb->adjm_Rq[cc];
      for (int cc = 0; cc < mb->dof_qd(); cc++) 
          pf_px_R[ns++] = mb->adjm_Rqd[cc];
    }

    std::vector<TinyVectorX> ans;

    curr_state = 0;
    curr_u = 0;
    TinyVectorX new_R(R.m_size);
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];

      if (mb->dof_state() == 0) {
        TinyVectorX dldu(mb->dof_u());
        ans.push_back(dldu);
        continue;
      }
      auto bodystate = std::vector<TinyScalar>(state.begin()+curr_state, state.begin()+curr_state+mb->dof_state());
      auto bodyu = u.segment(curr_u, mb->dof_u());

      TinyVectorX mb_pf_px_R(mb->dof_state()), mb_R_pf_pu(mb->dof_u()), 
        mb_dphi_du(mb->dof_u()), mb_R;

      mb_pf_px_R = pf_px_R.segment(curr_state, mb->dof_state());
      mb_R = R.segment(curr_state, mb->dof_state());
      mb->adj_bb_integrate_grad(mb_pf_px_R, mb_R);
      mb_pf_px_R.set_zero();
      mb->adj_bb_dynamics(bodyu, 
                          bodystate, 
                          mb_R, mb_pf_px_R, mb_R_pf_pu, dt);
      ans.push_back(mb_R_pf_pu);

      for (int cc = 0; cc < mb->dof_state(); cc++) {
        new_R[curr_state + cc] = mb_pf_px_R[cc];
      }
   
      curr_u += mb->dof_u();
      curr_state += mb->dof_state();

    }

    ans.push_back(new_R);
    return ans;
  }


  void adj_step(TinyScalar dt, const TinyVectorX& u) {
    // forward dynamics
    int curr_u = 0;
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      // save record
      mb->adjm_ckpt_state[adjm_step_count] = mb->adjm_this_state;
      mb->adjm_ckpt_dt[adjm_step_count] = dt;
      mb->adjm_ckpt_u[adjm_step_count] = u.segment(curr_u, mb->dof_u());
      
      // forward dynamics
      mb->adj_forward_dynamics(mb->adjm_ckpt_u[adjm_step_count], 
                               mb->adjm_this_state, 
                               mb->adjm_next_state, dt);
      curr_u += mb->dof_u();
    }

    // change state
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i];
      for (int cc = 0; cc < mb->dof_state(); cc++) 
        mb->adjm_this_state[cc] = mb->adjm_next_state[cc];
    }

    // forward constraint
    for (int i = 0; i < m_multi_bodies.size(); i++) {
      TinyMultiBody* mb = m_multi_bodies[i]; 
      mb->adjm_ckpt_state_cons[adjm_step_count] = mb->adjm_this_state;
    }
    adj_forward_constraint(dt);
    adjm_step_count++;
  }

  void adj_forward_constraint(TinyScalar dt) {
    {
      for (int i = 0; i < m_multi_bodies.size(); i++) {
        TinyMultiBody* mb = m_multi_bodies[i]; 
        for (int cc = 0; cc < mb->dof(); cc++) {
          mb->m_q[cc] = mb->adjm_this_state[cc];
        }
        for (int cc = 0; cc < mb->dof_qd(); cc++) 
          mb->m_qd[cc] = mb->adjm_this_state[cc+mb->dof()];
      }
    }

    {
      m_allMultiBodyContacts.reserve(1024);
      m_allMultiBodyContacts.resize(0);
    }
    {
      submitProfileTiming("compute multi body contacts");
      compute_contacts_multi_body_internal(
          m_multi_bodies, &m_dispatcher, m_allMultiBodyContacts,
          default_restitution, default_friction);
      submitProfileTiming("");
    }
    

    {
      submitProfileTiming("solve constraints");
      // use outer loop in case the multi-body constraint solver requires it
      // (e.g. sequential impulse method)

      for (int c = 0; c < m_allMultiBodyContacts.size(); c++) {
        m_mb_constraint_solver->resolveCollision(m_allMultiBodyContacts[c], dt);
      }
      
      submitProfileTiming("");
    }

    {
      submitProfileTiming("integrate");
      for (int i = 0; i < m_multi_bodies.size(); i++) {
        TinyMultiBody* mb = m_multi_bodies[i];
        
        mb->adj_f_integrate(mb->m_q, mb->m_qd, mb->m_qdd, dt);
      }
      submitProfileTiming("");
    }
    // update state
    {
      for (int i = 0; i < m_multi_bodies.size(); i++) {
        TinyMultiBody* mb = m_multi_bodies[i];  
        for (int cc = 0; cc < mb->dof(); cc++) 
            mb->adjm_next_state[cc] = mb->m_q[cc];
        for (int cc = 0; cc < mb->dof_qd(); cc++) 
            mb->adjm_next_state[cc+mb->dof()] = mb->m_qd[cc];
      }
    }
  }

  void step(TinyScalar dt) {
    {
      m_allContacts.reserve(1024);
      m_allContacts.resize(0);
      m_allMultiBodyContacts.reserve(1024);
      m_allMultiBodyContacts.resize(0);
      submitProfileTiming("apply forces");
      for (int i = 0; i < m_bodies.size(); i++) {
        TinyRigidBody* b = m_bodies[i];
        b->apply_gravity(m_gravity_acceleration);
        b->apply_force_impulse(dt);
        b->clear_forces();
      }
      submitProfileTiming("");
    }
    {
      submitProfileTiming("compute contacts");
      compute_contacts_rigid_body_internal(m_bodies, &m_dispatcher,
                                           m_allContacts, default_restitution,
                                           default_friction);
      submitProfileTiming("");
    }
    {
      submitProfileTiming("compute multi body contacts");
      compute_contacts_multi_body_internal(
          m_multi_bodies, &m_dispatcher, m_allMultiBodyContacts,
          default_restitution, default_friction);
      submitProfileTiming("");
    }

    {
      submitProfileTiming("solve constraints");
      for (int i = 0; i < m_num_solver_iterations; i++) {
        for (int c = 0; c < m_allContacts.size(); c++) {
          m_constraint_solver->resolveCollision(m_allContacts[c], dt);
        }
      }
      // use outer loop in case the multi-body constraint solver requires it
      // (e.g. sequential impulse method)
      int mb_solver_iters;
      if (!m_mb_constraint_solver->needs_outer_iterations) {
        mb_solver_iters = 1;
      } else {
        mb_solver_iters = m_num_solver_iterations;
      }
      for (int i = 0; i < mb_solver_iters; i++) {
        for (int c = 0; c < m_allMultiBodyContacts.size(); c++) { 
          m_mb_constraint_solver->resolveCollision(m_allMultiBodyContacts[c],
                                                   dt);
        }
      }
      submitProfileTiming("");
    }

    {
      submitProfileTiming("integrate");
      for (int i = 0; i < m_bodies.size(); i++) {
        TinyRigidBody* b = m_bodies[i];
        b->integrate(dt);
      }
      submitProfileTiming("");
    }
  }
};

#endif  // TINY_WORLD_H
