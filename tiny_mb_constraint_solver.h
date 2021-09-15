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

#ifndef TINY_MB_CONSTRAINT_SOLVER_H
#define TINY_MB_CONSTRAINT_SOLVER_H

#include "tiny_constraint_solver.h"
#include "tiny_multi_body.h"

template <typename TinyScalar, typename TinyConstants>
struct TinyContactPointMultiBody
    : public TinyContactPoint<TinyScalar, TinyConstants> {
  typedef ::TinyMultiBody<TinyScalar, TinyConstants> TinyMultiBody;
  TinyMultiBody* m_multi_body_a;
  TinyMultiBody* m_multi_body_b;
  TinyScalar m_restitution, m_friction;
  TinyScalar adjm_Rrestituition{TinyConstants::zero()}; 
  TinyScalar adjm_Rfriction{TinyConstants::zero()};
  int m_link_a;
  int m_link_b;
};

template <typename TinyScalar, typename TinyConstants>
struct TinyMultiBodyConstraintSolver {
  typedef ::TinyVector3<TinyScalar, TinyConstants> TinyVector3;
  typedef ::TinyVectorX<TinyScalar, TinyConstants> TinyVectorX;
  typedef ::TinyMatrix3xX<TinyScalar, TinyConstants> TinyMatrix3xX;
  typedef ::TinyMatrixXxX<TinyScalar, TinyConstants> TinyMatrixXxX;
  typedef ::TinyMultiBody<TinyScalar, TinyConstants> TinyMultiBody;
  typedef ::TinyContactPointMultiBody<TinyScalar, TinyConstants>
      TinyContactPoint;

  /**
   * Whether this method requires an outer iteration loop (such as the
   * sequential impulse method).
   */

  TinyScalar m_hi_friction = TinyConstants::fraction(1000, 1);
  TinyScalar adjm_Rhi_friction;

  bool needs_outer_iterations{false};

  int m_pgs_iterations{50};
  double m_least_squares_residual_threshold{0};
  std::vector<int> m_limit_dependency;

  /**
   * Projected Gauss-Seidel solver for a MLCP defined by coefficient matrix A
   * and vector b.
   *
   *   Ax + b >= 0
   *   s.t. x >= lo and x <= hi and xT(Ax + b) = 0
   *
   * where lo and hi are the respective lower and upper bounds on x.
   *
   * Reference: Jakub Stępień, PhD Thesis, 2013, p. 91.
   * Code reference: Bullet physics engine, Erwin Coumans.
   * https://github.com/bulletphysics/bullet3/blob/master/src/BulletDynamics/MLCPSolvers/btSolveProjectedGaussSeidel.h#L35
   */
  void solve_pgs(const TinyMatrixXxX& A, const TinyVectorX& b, TinyVectorX& x,
                 int num_iterations, double least_squares_residual_threshold,
                 const TinyVectorX* lo = nullptr,
                 const TinyVectorX* hi = nullptr) const {
    assert(A.m_rows == b.m_size);
    assert(A.m_cols == x.m_size);

    TinyScalar delta;
    double least_squares_residual;
    for (int k = 0; k < num_iterations; ++k) {
      least_squares_residual = 0;
      for (int i = 0; i < A.m_rows; ++i) {
        delta = TinyConstants::zero();
        for (int j = 0; j < i; j++) delta = delta + A(i, j) * x[j];
        for (int j = i + 1; j < A.m_rows; j++) delta = delta + A(i, j) * x[j];

        TinyScalar a_diag = A(i, i);
        TinyScalar x_old = x[i];
        x[i] = (b[i] - delta) / a_diag;
        TinyScalar s = TinyConstants::one();
        if (!m_limit_dependency.empty() && m_limit_dependency[i] >= 0) {
          s = x[m_limit_dependency[i]];
          if (TinyConstants::getBool(
            s <= TinyConstants::zero())) {
            s = TinyConstants::one();
          }
        }

        if (lo && TinyConstants::getBool(
            x[i] < (*lo)[i] * s)) {
          x[i] = (*lo)[i] * s;
        }
        if (hi && TinyConstants::getBool(
          x[i] > (*hi)[i] * s)) {
          x[i] = (*hi)[i] * s;
        }
        TinyScalar diff = x[i] - x_old;
      }
    }
  }

  void adj_solve_pgs(TinyVectorX& Rx, const int num_iterations,
                    const TinyMatrixXxX& A, 
                    const TinyVectorX& b, 
                    TinyVectorX& x, 
                    TinyMatrixXxX& RA,
                    TinyVectorX& Rb,
                    TinyVectorX& Rcon_lo,
                    TinyVectorX& Rcon_hi,
                    const TinyVectorX* lo = nullptr,
                    const TinyVectorX* hi = nullptr) const {

    std::vector<TinyVectorX> vec_x(num_iterations);
    std::vector<TinyVectorX> vec_x2(num_iterations);
    std::vector<TinyVectorX> vec_x3(num_iterations);
    std::vector<TinyVectorX> vec_s(num_iterations);    
    std::vector<TinyVectorX> vec_s0(num_iterations); 
    std::vector<TinyVectorX> vec_delta(num_iterations);

    TinyScalar delta;
    for (int k = 0; k < num_iterations; ++k) {
      vec_s[k] = TinyVectorX(A.m_rows);
      vec_s0[k] = TinyVectorX(A.m_rows);
      vec_delta[k] = TinyVectorX(A.m_rows);
      vec_x[k] = x;
      vec_x2[k] = x;
      vec_x3[k] = x;

      for (int i = 0; i < A.m_rows; ++i) {
        delta = TinyConstants::zero();
        for (int j = 0; j < i; j++) delta = delta + A(i, j) * x[j];
        for (int j = i + 1; j < A.m_rows; j++) delta = delta + A(i, j) * x[j];
        vec_delta[k][i] = delta;

        TinyScalar a_diag = A(i, i);
        x[i] = (b[i] - delta) / a_diag;
        vec_x2[k][i] = x[i];
        TinyScalar s = TinyConstants::one();
        if (!m_limit_dependency.empty() && m_limit_dependency[i] >= 0) {
          s = x[m_limit_dependency[i]];
          vec_s0[k][i] = s;
          if (TinyConstants::getBool(
            s <= TinyConstants::zero())) {
            s = TinyConstants::one();
          }
        }
        vec_s[k][i] = s;

        if (lo && TinyConstants::getBool(
            x[i] < (*lo)[i] * s)) {
          x[i] = (*lo)[i] * s;
        }
        if (hi && TinyConstants::getBool(
          x[i] > (*hi)[i] * s)) {
          x[i] = (*hi)[i] * s;
        }
        vec_x3[k][i] = x[i];
      }
    }

    for (int k = num_iterations - 1; k >= 0 ; --k) {
      for (int i = A.m_rows - 1; i >= 0 ; --i) {

        TinyScalar Rs = TinyConstants::zero();

        if (lo && TinyConstants::getBool(
            vec_x2[k][i] < (*lo)[i] * vec_s[k][i])) {
          Rcon_lo[i] += Rx[i] * vec_s[k][i];
          Rs += Rx[i] * (*lo)[i];
          Rx[i] = TinyConstants::zero();
        }

        if (hi && TinyConstants::getBool(
            vec_x2[k][i] > (*hi)[i] * vec_s[k][i])) {
          Rcon_hi[i] += Rx[i] * vec_s[k][i];
          Rs += Rx[i] * (*hi)[i];
          Rx[i] = TinyConstants::zero();
        }


        if (!m_limit_dependency.empty() && m_limit_dependency[i] >= 0) {
          if (TinyConstants::getBool(
            vec_s0[k][i] > TinyConstants::zero())) {
            Rx[m_limit_dependency[i]] += Rs;
          }
        }

        TinyScalar a_diag = A(i, i);
        TinyScalar Ra_diag, Rdelta;

        Rb[i] += Rx[i] / a_diag;
        Rdelta = -Rx[i] / a_diag;
        Ra_diag = - 1 / a_diag / a_diag * Rx[i] * (b[i] - vec_delta[k][i]);
        Rx[i] = TinyConstants::zero();

        RA(i, i) += Ra_diag;

        for (int j = i + 1; j < A.m_rows; j++) {
          RA(i, j) += Rdelta * vec_x[k][j];
          Rx[j] += Rdelta * A(i, j);
        }

        for (int j = 0; j < i; j++) {
          RA(i, j) += Rdelta * vec_x3[k][j];
          Rx[j] += Rdelta * A(i, j);
        }

      }
    }

  }

  // Solve impulse-based collision response using MLCP formulation from
  // Jakub Stępień, PhD Thesis, 2013.
  //
  // Args:
  // cps: contact points with distances < 0
  // dt : delta time (in seconds)
  virtual void resolveCollision(std::vector<TinyContactPoint>& cps,
                                TinyScalar dt) {
    if (cps.empty()) return;
    const int n_c = static_cast<int>(cps.size());

    const TinyContactPoint& cp0 = cps[0];

    TinyMultiBody* mb_a = cp0.m_multi_body_a;
    TinyMultiBody* mb_b = cp0.m_multi_body_b;


    const int n_a = mb_a->dof_qd();
    const int n_b = mb_b->dof_qd();
    const int n_ab = n_a + n_b;
    if (n_ab == 0) return;


    TinyMatrixXxX mass_matrix_a(n_a, n_a);
    mb_a->mass_matrix(&mass_matrix_a);
    bool is_positive_definite_a = true;
    bool is_positive_definite_b = true;

    TinyMatrixXxX mass_matrix_a_inv(mass_matrix_a.m_rows, mass_matrix_a.m_cols);
    if (mass_matrix_a.m_cols * mass_matrix_a.m_rows > 0) {
      mb_a->submitProfileTiming("inverse_mass_matrix_a");
      is_positive_definite_a = mass_matrix_a.inversed(mass_matrix_a_inv);
      mb_a->submitProfileTiming("");
    }

    TinyMatrixXxX mass_matrix_b(n_b, n_b);

    mb_b->mass_matrix(&mass_matrix_b);



    TinyMatrixXxX mass_matrix_b_inv(mass_matrix_b.m_rows, mass_matrix_b.m_cols);
    if (mass_matrix_b.m_cols * mass_matrix_b.m_rows > 0) {
      mb_b->submitProfileTiming("inverse_mass_matrix_b");
      is_positive_definite_b = mass_matrix_b.inversed(mass_matrix_b_inv);
      mb_b->submitProfileTiming("");
    }
    if (!is_positive_definite_a) {
      printf("LCP: mass matrix a is not positive definite");
    }
    if (!is_positive_definite_b) {
      printf("LCP: mass matrix b is not positive definite");
    }
    assert(is_positive_definite_a);
    assert(is_positive_definite_b);


    TinyMatrixXxX mass_matrix_inv(n_ab, n_ab);
    mass_matrix_inv.set_zero();
    mass_matrix_inv.assign_matrix(0, 0, mass_matrix_a_inv);
    mass_matrix_inv.assign_matrix(n_a, n_a, mass_matrix_b_inv);
    // Assemble constraint Jacobian J_C for a and b
    // The convention for constructing the constraint Jacobian is as follows:
    // For each contact point i the rows are as follows:
    //  i    is for the contact normal
    //  c+i  is for the friction direction towards the lateral velocity

    int num_friction_dir = 1;
    int dof_per_contact = 1 + num_friction_dir;

    TinyMatrixXxX jac_con(dof_per_contact * n_c, n_ab);
    jac_con.set_zero();
    TinyVectorX lcp_b(dof_per_contact * n_c);
    lcp_b.set_zero();

    TinyScalar erp =
        TinyConstants::fraction(1, 100);  // BAUMGARTE_ERROR_REDUCTION_PARAMETER
    TinyScalar cfm =
        TinyConstants::fraction(1, 100000);  // Constraint Force Mixing

    for (int i = 0; i < n_c; ++i) {
      TinyContactPoint& cp = cps[i];
      // all contact points are already assumed to have distance < 0
      if (TinyConstants::getBool(
        cp.m_distance > TinyConstants::zero())) continue;
      TinyVector3& world_point_a = cp.m_world_point_on_a;

      TinyMatrix3xX jac_a = 
        cp.m_multi_body_a->point_jacobian(cp.m_link_a, world_point_a);
      TinyVectorX jac_a_i = jac_a.mul_transpose(cp.m_world_normal_on_b);
      jac_con.assign_vector_horizontal(i, 0, jac_a_i);

      TinyVector3& world_point_b = cp.m_world_point_on_b;
      
      TinyMatrix3xX jac_b =
          cp.m_multi_body_b->point_jacobian(cp.m_link_b, world_point_b);
      TinyVectorX jac_b_i = jac_b.mul_transpose(cp.m_world_normal_on_b);
      std::vector<TinyScalar> qd_empty;
      int szb = cp.m_multi_body_b->m_qd.size();
      qd_empty.resize(szb, TinyConstants::zero());
      std::vector<TinyScalar> tau_jac;
      tau_jac.resize(szb);
      for (int i = 0; i < szb; i++) {
        tau_jac[i] = -jac_b_i[i];
      }

      // compare with unit impulse method
      // std::vector<TinyScalar> qdd_delta_unit_impulse;
      // qdd_delta_unit_impulse.resize(szb);
      // cp.m_multi_body_b->forward_dynamics(
      //    cp.m_multi_body_b->m_q, qd_empty, tau_jac, qdd_delta_unit_impulse,
      //    TinyConstants::fraction(100, 10));  // TinyConstants::zero());

      jac_con.assign_vector_horizontal(i, n_a, jac_b_i);

      TinyVectorX qd_a(cp.m_multi_body_a->m_qd);
      TinyVectorX qd_b(cp.m_multi_body_b->m_qd);
      TinyVector3 vel_a = jac_a * qd_a;
      TinyVector3 vel_b = jac_b * qd_b;
      TinyVector3 rel_vel = vel_a - vel_b;
      TinyScalar normal_rel_vel = cp.m_world_normal_on_b.dot(rel_vel);
      // Baumgarte stabilization

      TinyScalar baumgarte_rel_vel = erp * cp.m_distance / dt;
      lcp_b[i] = -(TinyConstants::one() + cp.m_restitution) * normal_rel_vel -
                 baumgarte_rel_vel;
      // friction direction
      TinyVector3 lateral_rel_vel =
          rel_vel - normal_rel_vel * cp.m_world_normal_on_b;
      const TinyScalar lateral = lateral_rel_vel.length();
      TinyVector3 fr_direction1, fr_direction2;
      
      if (TinyConstants::getBool(
        lateral < TinyConstants::fraction(1, 10000))) {
        // use the plane space of the contact normal as friction directions
        cp.m_world_normal_on_b.plane_space(fr_direction1, fr_direction2);
      } else {
        // use the negative lateral velocity and its orthogonal as friction
        // directions
        fr_direction1 = lateral_rel_vel * (TinyConstants::one() / lateral);
        fr_direction2 = fr_direction1.cross(cp.m_world_normal_on_b);
      }
      TinyScalar l1 = fr_direction1.dot(rel_vel);
      lcp_b[n_c + i] = -l1;
      if (num_friction_dir > 1) {
        TinyScalar l2 = fr_direction2.dot(rel_vel);
        lcp_b[2 * n_c + i] = -l2;
      }
      TinyVectorX jac_a_i_fr1 = jac_a.mul_transpose(fr_direction1);
      jac_con.assign_vector_horizontal(n_c + i, 0, jac_a_i_fr1);
      TinyVectorX jac_b_i_fr1 = jac_b.mul_transpose(fr_direction1);
      jac_con.assign_vector_horizontal(n_c + i, n_a, jac_b_i_fr1);
      if (num_friction_dir > 1) {
        TinyVectorX jac_a_i_fr2 = jac_a.mul_transpose(fr_direction2);
        jac_con.assign_vector_horizontal(2 * n_c + i, 0, jac_a_i_fr2);
        TinyVectorX jac_b_i_fr2 = jac_b.mul_transpose(fr_direction2);
        jac_con.assign_vector_horizontal(2 * n_c + i, n_a, jac_b_i_fr2);
      }
    }

    TinyMatrixXxX jac_con_t = jac_con.transpose();
    TinyMatrixXxX lcp_A;
    {
      mb_b->submitProfileTiming("lcpA");
      lcp_A = jac_con * mass_matrix_inv * jac_con_t;
      mb_b->submitProfileTiming("");
    }
    // apply CFM
    // This is _necessary_ for fixed base systems where the degrees of freedom
    // would otherwise not allow for the lateral friction directions, leading
    // to zero constraint Jacobian rows and eventually zeros on the diagonal
    // of the A matrix.
    {
      mb_b->submitProfileTiming("cfm");
      for (int i = 0; i < dof_per_contact * n_c; ++i) {
        lcp_A(i, i) = lcp_A(i, i) + cfm;
      }
      mb_b->submitProfileTiming("");
    }
    TinyVectorX lcp_p(dof_per_contact * n_c);
    lcp_p.set_zero();
    TinyVectorX con_lo(dof_per_contact * n_c);
    TinyVectorX con_hi(dof_per_contact * n_c);
    con_lo.set_zero();
    m_limit_dependency.reserve(dof_per_contact * n_c);
    m_limit_dependency.resize(dof_per_contact * n_c);
    for (int i = 0; i < n_c; ++i) {
      m_limit_dependency[i] = -1;
      con_hi[i] = TinyConstants::fraction(1000, 1);
      // ||friction impulse|| <= mu * ||normal impulse||
      con_hi[n_c + i] = cps[i].m_friction;
      con_lo[n_c + i] = -cps[i].m_friction;
      m_limit_dependency[n_c + i] = i;
      if (num_friction_dir > 1) {
        con_hi[2 * n_c + i] = cps[i].m_friction;
        con_lo[2 * n_c + i] = -cps[i].m_friction;
        m_limit_dependency[2 * n_c + i] = i;
      }
    }
    
    {
      //      fflush(stdout);
      mb_b->submitProfileTiming("solve_pgs");
      solve_pgs(lcp_A, lcp_b, lcp_p, m_pgs_iterations,
                m_least_squares_residual_threshold, &con_lo, &con_hi);
      mb_b->submitProfileTiming("");
    }
    if (n_a > 0) {
      TinyVectorX p_a = lcp_p.segment(0, n_c);
      TinyMatrixXxX jac_con_a = jac_con.block(0, 0, n_c, n_a);
      TinyVectorX delta_qd_a = mass_matrix_a_inv * jac_con_a.mul_transpose(p_a);
      // add friction impulse
      TinyVectorX p_a_fr = lcp_p.segment(0, n_c);
      TinyMatrixXxX jac_con_a_fr = jac_con.block(n_c, 0, n_c, n_a);
      delta_qd_a += mass_matrix_a_inv * jac_con_a_fr.mul_transpose(p_a_fr);
      for (int i = 0; i < n_a; ++i) {
        mb_a->m_qd[i] = mb_a->m_qd[i] + delta_qd_a[i];
      }
    }
    
    if (n_b > 0) {
      TinyVectorX p_b = lcp_p.segment(0, n_c);
      TinyMatrixXxX jac_con_b = jac_con.block(0, n_a, n_c, n_b);
      TinyVectorX delta_qd_b = mass_matrix_b_inv * jac_con_b.mul_transpose(p_b);
      // add friction impulse
      if (1) {
        TinyVectorX p_b_fr = lcp_p.segment(n_c, n_c);
        TinyMatrixXxX jac_con_b_fr = jac_con.block(n_c, n_a, n_c, n_b);
        TinyVectorX tmp_mul_trans = jac_con_b_fr.mul_transpose(p_b_fr);
        TinyVectorX fr_qd =
            mass_matrix_b_inv * tmp_mul_trans;
        delta_qd_b += fr_qd;
      }
      if (num_friction_dir > 1) {
        TinyVectorX p_b_fr = lcp_p.segment(2 * n_c, n_c);
        TinyMatrixXxX jac_con_b_fr = jac_con.block(2 * n_c, n_a, n_c, n_b);
        TinyVectorX fr_qd =
            mass_matrix_b_inv * jac_con_b_fr.mul_transpose(p_b_fr);
        delta_qd_b += fr_qd;
      }

      for (int i = 0; i < n_b; ++i) {
        mb_b->m_qd[i] = mb_b->m_qd[i] - delta_qd_b[i];
      }
    }
  }


  // Solve impulse-based collision response using MLCP formulation from
  // Jakub Stępień, PhD Thesis, 2013.
  //
  // Args:
  // cps: contact points with distances < 0
  // dt : delta time (in seconds)
  void adj_resolveCollision(std::vector<TinyContactPoint>& cps,
                                TinyScalar dt) {
    if (cps.empty()) return;
    const int n_c = static_cast<int>(cps.size());
    const TinyContactPoint& cp0 = cps[0];
    // start forward
    // -------------------------------------------------
    // if (cps.empty()) return;
    TinyMultiBody* mb_a = cp0.m_multi_body_a;
    TinyMultiBody* mb_b = cp0.m_multi_body_b;

    const int n_a = mb_a->dof_qd();
    const int n_b = mb_b->dof_qd();
    const int n_ab = n_a + n_b;
    if (n_ab == 0) return;
    bool is_positive_definite_a = true;
    bool is_positive_definite_b = true;

    TinyMatrixXxX mass_matrix_a(n_a, n_a), Rmass_matrix_a(n_a, n_a);
    mb_a->mass_matrix(&mass_matrix_a);
    TinyMatrixXxX mass_matrix_a_inv(mass_matrix_a.m_rows, mass_matrix_a.m_cols),
                Rmass_matrix_a_inv(mass_matrix_a.m_rows, mass_matrix_a.m_cols);
    if (mass_matrix_a.m_cols * mass_matrix_a.m_rows > 0) {
      mb_a->submitProfileTiming("inverse_mass_matrix_a");
      is_positive_definite_a = mass_matrix_a.inversed(mass_matrix_a_inv);
      mb_a->submitProfileTiming("");
    }

    TinyMatrixXxX mass_matrix_b(n_b, n_b), Rmass_matrix_b(n_b, n_b);
    mb_b->mass_matrix(&mass_matrix_b);
    TinyMatrixXxX mass_matrix_b_inv(mass_matrix_b.m_rows, mass_matrix_b.m_cols),
                Rmass_matrix_b_inv(mass_matrix_b.m_rows, mass_matrix_b.m_cols);
    if (mass_matrix_b.m_cols * mass_matrix_b.m_rows > 0) {
      mb_b->submitProfileTiming("inverse_mass_matrix_b");
      is_positive_definite_b = mass_matrix_b.inversed(mass_matrix_b_inv);
      mb_b->submitProfileTiming("");
    }
    if (!is_positive_definite_a) {
      printf("LCP: mass matrix a is not positive definite");
    }
    if (!is_positive_definite_b) {
      printf("LCP: mass matrix b is not positive definite");
    }
    assert(is_positive_definite_a);
    assert(is_positive_definite_b);

    TinyMatrixXxX mass_matrix_inv(n_ab, n_ab), Rmass_matrix_inv(n_ab, n_ab);
    mass_matrix_inv.assign_matrix(0, 0, mass_matrix_a_inv);
    mass_matrix_inv.assign_matrix(n_a, n_a, mass_matrix_b_inv);

    int num_friction_dir = 1;
    int dof_per_contact = 1 + num_friction_dir;

    TinyMatrixXxX jac_con(dof_per_contact * n_c, n_ab),
                Rjac_con(dof_per_contact * n_c, n_ab);
    TinyVectorX lcp_b(dof_per_contact * n_c), Rlcp_b(dof_per_contact * n_c);

    TinyScalar erp =
        TinyConstants::fraction(1, 100);  // BAUMGARTE_ERROR_REDUCTION_PARAMETER
    TinyScalar cfm =
        TinyConstants::fraction(1, 100000);  // Constraint Force Mixing

    std::vector<TinyMatrix3xX> vec_jac_a(n_c);
    std::vector<TinyVectorX> vec_jac_a_i(n_c);
    std::vector<TinyMatrix3xX> vec_jac_b(n_c);
    std::vector<TinyVectorX> vec_jac_b_i(n_c);
    std::vector<TinyVector3> vec_vel_a(n_c);
    std::vector<TinyVector3> vec_vel_b(n_c);
    std::vector<TinyVector3> vec_rel_vel(n_c);
    std::vector<TinyScalar> vec_normal_rel_vel(n_c);
    std::vector<TinyScalar> vec_baumgarte_rel_vel(n_c);
    std::vector<TinyVector3> vec_lateral_rel_vel(n_c);
    std::vector<TinyScalar> vec_lateral(n_c);
    std::vector<TinyVector3> vec_fr_direction1(n_c);
    std::vector<TinyVector3> vec_fr_direction2(n_c);
    std::vector<TinyVectorX> vec_jac_a_i_fr1(n_c);
    std::vector<TinyVectorX> vec_jac_b_i_fr1(n_c);
    std::vector<TinyVectorX> vec_jac_a_i_fr2(n_c);
    std::vector<TinyVectorX> vec_jac_b_i_fr2(n_c);

    for (int i = 0; i < n_c; ++i) {
      TinyContactPoint& cp = cps[i];
      // all contact points are already assumed to have distance < 0
      if (TinyConstants::getBool(
        cp.m_distance > TinyConstants::zero())) continue;

      vec_jac_a[i] =
          cp.m_multi_body_a->point_jacobian(cp.m_link_a, cp.m_world_point_on_a);
      vec_jac_a_i[i] = vec_jac_a[i].mul_transpose(cp.m_world_normal_on_b);
      
      jac_con.assign_vector_horizontal(i, 0, vec_jac_a_i[i]);

      vec_jac_b[i] =
          cp.m_multi_body_b->point_jacobian(cp.m_link_b, cp.m_world_point_on_b);
      vec_jac_b_i[i] = vec_jac_b[i].mul_transpose(cp.m_world_normal_on_b);
 

      jac_con.assign_vector_horizontal(i, n_a, vec_jac_b_i[i]);

      TinyVectorX qd_a(cp.m_multi_body_a->m_qd);
      TinyVectorX qd_b(cp.m_multi_body_b->m_qd);
      vec_vel_a[i] = vec_jac_a[i] * qd_a;
      vec_vel_b[i] = vec_jac_b[i] * qd_b;
      vec_rel_vel[i] = vec_vel_a[i] - vec_vel_b[i];
      vec_normal_rel_vel[i] = cp.m_world_normal_on_b.dot(vec_rel_vel[i]);
      // Baumgarte stabilization
      vec_baumgarte_rel_vel[i] = erp * cp.m_distance / dt;

      lcp_b[i] = -(TinyConstants::one() + cp.m_restitution) * vec_normal_rel_vel[i] -
                 vec_baumgarte_rel_vel[i];
      vec_lateral_rel_vel[i] =
          vec_rel_vel[i] - vec_normal_rel_vel[i] * cp.m_world_normal_on_b;
      vec_lateral[i] = vec_lateral_rel_vel[i].length();
      vec_fr_direction1[i], vec_fr_direction2[i];

      if (TinyConstants::getBool(
        vec_lateral[i] < TinyConstants::fraction(1, 10000))) {
        cp.m_world_normal_on_b.plane_space(vec_fr_direction1[i], vec_fr_direction2[i]);
      } else {
        vec_fr_direction1[i] = vec_lateral_rel_vel[i] * (TinyConstants::one() / vec_lateral[i]);
        vec_fr_direction2[i] = vec_fr_direction1[i].cross(cp.m_world_normal_on_b);
      }

      lcp_b[n_c + i] = -vec_fr_direction1[i].dot(vec_rel_vel[i]);
      if (num_friction_dir > 1) {
        lcp_b[2 * n_c + i] = -vec_fr_direction2[i].dot(vec_rel_vel[i]);
      }

      vec_jac_a_i_fr1[i] = vec_jac_a[i].mul_transpose(vec_fr_direction1[i]);
      jac_con.assign_vector_horizontal(n_c + i, 0, vec_jac_a_i_fr1[i]);
      vec_jac_b_i_fr1[i] = vec_jac_b[i].mul_transpose(vec_fr_direction1[i]);
      jac_con.assign_vector_horizontal(n_c + i, n_a, vec_jac_b_i_fr1[i]);
      if (num_friction_dir > 1) {
        vec_jac_a_i_fr2[i] = vec_jac_a[i].mul_transpose(vec_fr_direction2[i]);
        jac_con.assign_vector_horizontal(2 * n_c + i, 0, vec_jac_a_i_fr2[i]);
        vec_jac_b_i_fr2[i] = vec_jac_b[i].mul_transpose(vec_fr_direction2[i]);
        jac_con.assign_vector_horizontal(2 * n_c + i, n_a, vec_jac_b_i_fr2[i]);
      }
    }
    TinyMatrixXxX jac_con_t = jac_con.transpose();
    TinyMatrixXxX lcp_A;

    {
      mb_b->submitProfileTiming("lcpA");
      lcp_A = jac_con * mass_matrix_inv * jac_con_t;
      mb_b->submitProfileTiming("");
    }

    // apply CFM
    // This is _necessary_ for fixed base systems where the degrees of freedom
    // would otherwise not allow for the lateral friction directions, leading
    // to zero constraint Jacobian rows and eventually zeros on the diagonal
    // of the A matrix.
    {
      mb_b->submitProfileTiming("cfm");
      for (int i = 0; i < dof_per_contact * n_c; ++i) {
        lcp_A(i, i) = lcp_A(i, i) + cfm;
      }
      mb_b->submitProfileTiming("");
    }
    TinyVectorX lcp_p(dof_per_contact * n_c);
    lcp_p.set_zero();
    TinyVectorX con_lo(dof_per_contact * n_c);
    TinyVectorX con_hi(dof_per_contact * n_c);
    con_lo.set_zero();
    m_limit_dependency.reserve(dof_per_contact * n_c);
    m_limit_dependency.resize(dof_per_contact * n_c);
    for (int i = 0; i < n_c; ++i) {
      m_limit_dependency[i] = -1;
      con_hi[i] = TinyConstants::fraction(1000, 1);
      con_hi[n_c + i] = cps[i].m_friction;
      con_lo[n_c + i] = -cps[i].m_friction;
      m_limit_dependency[n_c + i] = i;
      if (num_friction_dir > 1) {
        con_hi[2 * n_c + i] = cps[i].m_friction;
        con_lo[2 * n_c + i] = -cps[i].m_friction;
        m_limit_dependency[2 * n_c + i] = i;
      }
    }

    {
      mb_b->submitProfileTiming("solve_pgs");
      solve_pgs(lcp_A, lcp_b, lcp_p, m_pgs_iterations,
                m_least_squares_residual_threshold, &con_lo, &con_hi);
      mb_b->submitProfileTiming("");
    }

    // start adjoint
    // -------------------------------------------------
    Rjac_con.set_zero();
    Rmass_matrix_a_inv.set_zero();
    Rmass_matrix_b_inv.set_zero();
    TinyVectorX Rlcp_p(lcp_p.m_size);

    if (n_b > 0) {
      TinyVectorX Rdelta_qd_b(n_b);
      for (int i = 0; i < n_b; ++i) {
        Rdelta_qd_b[i] += -mb_b->adjm_Rqd[i];
      }
      TinyVectorX tmp_mul_trans;
      TinyVectorX Rtmp_mul_trans;
      TinyMatrixXxX jac_con_b_fr;
      TinyMatrixXxX Rjac_con_b_fr;
      TinyVectorX Rp_b_fr(n_c);
      TinyVectorX p_b_fr;

      if (num_friction_dir > 1) {
        jac_con_b_fr = jac_con.block(2 * n_c, n_a, n_c, n_b);
        Rjac_con_b_fr = TinyMatrixXxX(jac_con_b_fr.m_rows, jac_con_b_fr.m_cols);
      
        p_b_fr = lcp_p.segment(2 * n_c, n_c);

        tmp_mul_trans = jac_con_b_fr.mul_transpose(p_b_fr);
        Rmass_matrix_b_inv += TinyMatrixXxX::vvt(Rdelta_qd_b, tmp_mul_trans);
        Rtmp_mul_trans = mass_matrix_b_inv.transpose() * Rdelta_qd_b;
        jac_con_b_fr.adj_mx_mul_transpose(Rtmp_mul_trans, p_b_fr, 
                                    Rjac_con_b_fr, Rp_b_fr);
        Rjac_con.assign_matrix_add(2 * n_c, n_a, Rjac_con_b_fr);
        Rlcp_p.assign_vector_add(2 * n_c, Rp_b_fr);
      }
      p_b_fr = lcp_p.segment(n_c, n_c);
      jac_con_b_fr = jac_con.block(n_c, n_a, n_c, n_b);
      Rp_b_fr.set_zero();
      Rjac_con_b_fr = TinyMatrixXxX(jac_con_b_fr.m_rows, jac_con_b_fr.m_cols);
      tmp_mul_trans = jac_con_b_fr.mul_transpose(p_b_fr);
      
      Rmass_matrix_b_inv += TinyMatrixXxX::vvt(Rdelta_qd_b, tmp_mul_trans);
      Rtmp_mul_trans = mass_matrix_b_inv.transpose()*Rdelta_qd_b;
      jac_con_b_fr.adj_mx_mul_transpose(Rtmp_mul_trans, p_b_fr, 
                                  Rjac_con_b_fr, Rp_b_fr);
      Rjac_con.assign_matrix_add(n_c, n_a, Rjac_con_b_fr);
      Rlcp_p.assign_vector_add(n_c, Rp_b_fr);

      TinyVectorX p_b = lcp_p.segment(0, n_c);
      TinyVectorX Rp_b(n_c);
      TinyMatrixXxX jac_con_b = jac_con.block(0, n_a, n_c, n_b);
      TinyMatrixXxX Rjac_con_b(jac_con_b);
      Rjac_con_b.set_zero();

      tmp_mul_trans = jac_con_b.mul_transpose(p_b);
      Rmass_matrix_b_inv += TinyMatrixXxX::vvt(Rdelta_qd_b, tmp_mul_trans);
      Rtmp_mul_trans = mass_matrix_b_inv.transpose()*Rdelta_qd_b;
      jac_con_b.adj_mx_mul_transpose(Rtmp_mul_trans, p_b, 
                                  Rjac_con_b, Rp_b);

      Rjac_con.assign_matrix_add(0, n_a, Rjac_con_b);
      Rlcp_p.assign_vector_add(0, Rp_b);
    }
    if (n_a > 0) {
      TinyVectorX Rdelta_qd_a(n_a);

      for (int i = 0; i < n_a; ++i) {
        Rdelta_qd_a[i] += mb_a->adjm_Rqd[i];
      }

      // ++++++++++++
      TinyMatrixXxX jac_con_a_fr = jac_con.block(n_c, 0, n_c, n_a);
      TinyMatrixXxX Rjac_con_a_fr(jac_con_a_fr.m_rows, jac_con_a_fr.m_cols);
      TinyVectorX p_a_fr = lcp_p.segment(0, n_c);
      TinyVectorX Rp_a_fr(n_c);

      TinyVectorX tmp_mul_trans = jac_con_a_fr.mul_transpose(p_a_fr);
      TinyVectorX Rtmp_mul_trans(n_a);
      Rmass_matrix_a_inv += TinyMatrixXxX::vvt(Rdelta_qd_a, tmp_mul_trans);
      Rtmp_mul_trans = mass_matrix_a_inv.transpose() * Rdelta_qd_a;
      jac_con_a_fr.adj_mx_mul_transpose(Rtmp_mul_trans, p_a_fr, 
                                  Rjac_con_a_fr, Rp_a_fr);
      Rjac_con.assign_matrix_add(n_c, 0, Rjac_con_a_fr);
      Rlcp_p.assign_vector_add(0, Rp_a_fr);

      // ++++++++++++
      TinyMatrixXxX jac_con_a = jac_con.block(0, 0, n_c, n_a);
      TinyMatrixXxX Rjac_con_a(jac_con_a.m_rows, jac_con_a.m_cols);
      TinyVectorX p_a = lcp_p.segment(0, n_c);
      TinyVectorX Rp_a(n_c);

      tmp_mul_trans = jac_con_a.mul_transpose(p_a);
      Rtmp_mul_trans.set_zero();
      Rmass_matrix_a_inv += TinyMatrixXxX::vvt(Rdelta_qd_a, tmp_mul_trans);
      Rtmp_mul_trans = mass_matrix_a_inv.transpose() * Rdelta_qd_a;
      jac_con_a.adj_mx_mul_transpose(Rtmp_mul_trans, p_a, 
                                  Rjac_con_a, Rp_a);
      Rjac_con.assign_matrix_add(0, 0, Rjac_con_a);
      Rlcp_p.assign_vector_add(0, Rp_a);
    }


    TinyMatrixXxX Rlcp_A(lcp_A);
    Rlcp_A.set_zero();

    TinyVectorX Rcon_lo(dof_per_contact * n_c);
    TinyVectorX Rcon_hi(dof_per_contact * n_c);
    TinyScalar Rhi_fraction;

    {
      mb_b->submitProfileTiming("adj_solve_pgs");
      lcp_p.set_zero();

      adj_solve_pgs(Rlcp_p, m_pgs_iterations,
                    lcp_A, lcp_b, lcp_p,
                    Rlcp_A, Rlcp_b,
                    Rcon_lo, Rcon_hi, &con_lo, &con_hi); // TODO , Rcon_lo, Rcon_hi
      mb_b->submitProfileTiming("");
    }
    for (int i = n_c -1 ; i >= 0; --i) {
      adjm_Rhi_friction += Rcon_hi[i];
      cps[i].adjm_Rfriction += Rcon_hi[n_c + i];
      cps[i].adjm_Rfriction += -Rcon_lo[n_c + i];
    }

    {
      mb_b->submitProfileTiming("adj cfm");
      for (int i = 0; i < dof_per_contact * n_c; ++i) {
        lcp_A(i, i) = lcp_A(i, i) - cfm;
      }
      mb_b->submitProfileTiming("");
    }

    {
      mb_b->submitProfileTiming("adj lcpA");
      Rmass_matrix_inv += jac_con.transpose() * Rlcp_A * jac_con;
      Rjac_con += Rlcp_A * jac_con *mass_matrix_inv.transpose() 
              + Rlcp_A.transpose() * jac_con * mass_matrix_inv;
      mb_b->submitProfileTiming("");
    }

    TinyVector3 Rfr_direction2, Rfr_direction1;
    TinyVector3 Rrel_vel,Rlateral_rel_vel;
    TinyScalar Rnormal_rel_vel, Rbaumgarte_rel_vel;
    TinyVector3 Rvel_a, Rvel_b;

    for (int i = n_c - 1; i >= 0; i--) {
      TinyContactPoint& cp = cps[i];
      // all contact points are already assumed to have distance < 0
      if (TinyConstants::getBool(
        cp.m_distance > TinyConstants::zero())) continue;

      TinyVectorX Rjac_b_i_fr2(vec_jac_b_i_fr2[i].m_size), 
                  Rjac_a_i_fr2(vec_jac_a_i_fr2[i].m_size);
      TinyVectorX Rjac_b_i_fr1(vec_jac_b_i_fr1[i].m_size), 
                  Rjac_a_i_fr1(vec_jac_a_i_fr1[i].m_size);
      TinyVectorX Rjac_b_i(vec_jac_b_i[i].m_size),
                  Rjac_a_i(vec_jac_a_i[i].m_size);
      TinyMatrix3xX Rjac_b(vec_jac_b[i].m_rows, vec_jac_b[i].m_cols), 
                    Rjac_a(vec_jac_a[i].m_rows, vec_jac_a[i].m_cols);
      Rfr_direction2.set_zero();
      Rfr_direction1.set_zero();
      Rrel_vel.set_zero();
      Rvel_a.set_zero(); 
      Rvel_b.set_zero();
      Rlateral_rel_vel.set_zero();
      Rnormal_rel_vel = TinyConstants::zero();
      Rbaumgarte_rel_vel = TinyConstants::zero();

      TinyVectorX Rqd_a(cp.m_multi_body_a->m_qd.size());
      TinyVectorX Rqd_b(cp.m_multi_body_b->m_qd.size());
      if (num_friction_dir > 1) {
        Rjac_b_i_fr2 += Rjac_con.get_vector_horizontal(2 * n_c + i, n_a, 
                                          Rjac_b_i_fr2.size());
        Rjac_a_i_fr2 += Rjac_con.get_vector_horizontal(2 * n_c + i, 0, 
                                        Rjac_a_i_fr2.size());
        vec_jac_b[i].adj_mx_mul_transpose(Rjac_b_i_fr2, vec_fr_direction2[i],
                                    Rjac_b, Rfr_direction2);
        vec_jac_a[i].adj_mx_mul_transpose(Rjac_a_i_fr2, vec_fr_direction2[i],
                                    Rjac_a, Rfr_direction2);
      }

      Rjac_b_i_fr1 += Rjac_con.get_vector_horizontal(n_c + i, n_a, 
                                        Rjac_b_i_fr1.size());
      Rjac_a_i_fr1 += Rjac_con.get_vector_horizontal(n_c + i, 0, 
                                        Rjac_a_i_fr1.size());
      vec_jac_b[i].adj_mx_mul_transpose(Rjac_b_i_fr1, vec_fr_direction1[i],
                                    Rjac_b, Rfr_direction1);
      vec_jac_a[i].adj_mx_mul_transpose(Rjac_a_i_fr1, vec_fr_direction1[i],
                                    Rjac_a, Rfr_direction1);
      if (num_friction_dir > 1) {
  
        Rfr_direction2 += -vec_rel_vel[i] * Rlcp_b[2 * n_c + i];
        Rrel_vel += -vec_fr_direction2[i] * Rlcp_b[2 * n_c + i];
      }

      Rfr_direction1 += -vec_rel_vel[i] * Rlcp_b[n_c + i];
      Rrel_vel += -vec_fr_direction1[i] * Rlcp_b[n_c + i];
      
      TinyScalar Rlateral = TinyConstants::zero();
      if (TinyConstants::getBool(
        vec_lateral[i] < TinyConstants::fraction(1, 10000))) {
        cp.m_world_normal_on_b.adj_plane_space(Rfr_direction1, Rfr_direction2,
                                        vec_fr_direction1[i], vec_fr_direction2[i],
                                       cp.adjm_Rworld_normal_on_b);
      } else {
        Rfr_direction1 += -Rfr_direction2.cross(cp.m_world_normal_on_b);
        cp.adjm_Rworld_normal_on_b += Rfr_direction2.cross(vec_fr_direction1[i]);
        Rfr_direction2.set_zero();
        Rlateral_rel_vel += Rfr_direction1 * 
                            (TinyConstants::one() / vec_lateral[i]);
        Rlateral += -(TinyConstants::one() / vec_lateral[i] /vec_lateral[i]) * 
                                  Rfr_direction1.dot(vec_lateral_rel_vel[i]);
        Rfr_direction2.set_zero();
        Rlateral_rel_vel += Rlateral / vec_lateral[i] * vec_lateral_rel_vel[i];
      }
      Rrel_vel += Rlateral_rel_vel;
      Rnormal_rel_vel += -cp.m_world_normal_on_b.dot(Rlateral_rel_vel);
      cp.adjm_Rworld_normal_on_b += -vec_normal_rel_vel[i] * Rlateral_rel_vel;
      Rbaumgarte_rel_vel += -Rlcp_b[i];
      Rnormal_rel_vel += -(TinyConstants::one() + cp.m_restitution) * Rlcp_b[i];
      cp.adjm_Rdistance += Rbaumgarte_rel_vel * erp / dt;

      cp.adjm_Rworld_normal_on_b += vec_rel_vel[i] * Rnormal_rel_vel;
      Rrel_vel += Rnormal_rel_vel * cp.m_world_normal_on_b;

      Rvel_a += Rrel_vel;
      Rvel_b += -Rrel_vel;

      TinyVectorX qd_a(cp.m_multi_body_a->m_qd);
      TinyVectorX qd_b(cp.m_multi_body_b->m_qd);
      Rjac_b += TinyMatrixXxX::vvt(Rvel_b, qd_b);
      Rqd_b += vec_jac_b[i].transpose() * Rvel_b;


      Rjac_a += TinyMatrixXxX::vvt(Rvel_a, qd_a);
      Rqd_a += vec_jac_a[i].transpose() * Rvel_a;

      cp.m_multi_body_a->adjm_Rqd += Rqd_a;
      cp.m_multi_body_b->adjm_Rqd += Rqd_b;

      Rqd_a.set_zero();
      Rqd_b.set_zero();

      Rjac_b_i += Rjac_con.get_vector_horizontal(i, n_a, Rjac_b_i.m_size);
      vec_jac_b[i].adj_mx_mul_transpose(Rjac_b_i, cp.m_world_normal_on_b,
                                Rjac_b, cp.adjm_Rworld_normal_on_b);
      cp.m_multi_body_b->adj_point_jacobian(Rjac_b, cp.m_link_b, 
                        cp.m_world_point_on_b, cp.adjm_Rworld_point_on_b);
      Rjac_a_i += Rjac_con.get_vector_horizontal(i, 0, Rjac_a_i.m_size);
      vec_jac_a[i].adj_mx_mul_transpose(Rjac_a_i, cp.m_world_normal_on_b,
                                Rjac_a, cp.adjm_Rworld_normal_on_b);
      cp.m_multi_body_a->adj_point_jacobian(Rjac_a, cp.m_link_a, 
                        cp.m_world_point_on_a, cp.adjm_Rworld_point_on_a);
    }
    
    Rmass_matrix_a_inv += Rmass_matrix_inv.block(0, 0, 
                Rmass_matrix_a_inv.m_rows, Rmass_matrix_a_inv.m_cols);
    Rmass_matrix_b_inv += Rmass_matrix_inv.block(n_a, n_a, 
                Rmass_matrix_b_inv.m_rows, Rmass_matrix_b_inv.m_cols);
    if (mass_matrix_b.m_cols * mass_matrix_b.m_rows > 0) {
      mb_b->submitProfileTiming("adj inverse_mass_matrix_b");
      mass_matrix_b.adj_inversed(Rmass_matrix_b_inv, mass_matrix_b_inv, Rmass_matrix_b);
      mb_b->submitProfileTiming("");
    }

    mb_b->adj_mass_matrix(Rmass_matrix_b, mass_matrix_b);

    if (mass_matrix_a.m_cols * mass_matrix_a.m_rows > 0) {
      mass_matrix_a.adj_inversed(Rmass_matrix_a_inv, mass_matrix_a_inv, Rmass_matrix_a);
      mb_a->submitProfileTiming("");
    }
    mb_a->adj_mass_matrix(Rmass_matrix_a, mass_matrix_a);
  }
};

#endif  // TINY_MB_CONSTRAINT_SOLVER_H
