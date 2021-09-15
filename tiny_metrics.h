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

#ifndef TINY_METRICS_H
#define TINY_METRICS_H

#include <string>
#include <vector>
#include <stdio.h>

#include "tiny_matrix_x.h"
#include "tiny_vector_x.h"

template <typename TinyScalar, typename TinyConstants>
class TinyMetric {
  typedef ::TinyMatrixXxX<TinyScalar, TinyConstants> TinyMatrixXxX;
  typedef ::TinyVectorX<TinyScalar, TinyConstants> TinyVectorX;
 public:

  TinyVectorX m_x, m_u;
  TinyVectorX m_px, m_pu;
  
  explicit TinyMetric() {}

  virtual ~TinyMetric() { clear(); }

  void clear() { }

  TinyScalar forward(TinyVectorX x, TinyVectorX u) {
    m_x = x;
    m_u = u;

    TinyScalar ans = TinyConstants::zero();
    for (int i = 0; i < m_x.size(); i++)
      ans = ans + m_x[i];
    return ans;
  }

  TinyScalar forward(TinyVectorX x, std::vector<TinyScalar> u) {
    TinyVectorX input_u(u.size());
    for (int i = 0; i < u.size(); i++)
      input_u[i] = u[i];

    m_x = x;
    m_u = input_u;

    TinyScalar ans = TinyConstants::zero();
    for (int i = 0; i < m_x.size(); i++)
      ans = ans + m_x[i];
    return ans;
  }

  void backward() {
    TinyVectorX px(m_x.size());
    TinyVectorX pu(m_u.size());
    m_px = px;
    m_pu = pu;
    for (int i = 0; i < m_px.size(); i++)
      m_px[i] = TinyConstants::one();
    for (int i = 0; i < m_pu.size(); i++)
      m_pu[i] = TinyConstants::zero();
  }

  TinyVectorX get_px() {
    return m_px;
  }

  TinyVectorX get_pu() {
    return m_pu;
  }


  TinyVectorX get_px(int start, int length) const {
    assert(start >= 0);
    assert(start + length <= m_px.size());
    return m_px.segment(start, length);
  }

  TinyVectorX get_pu(int start, int length) const {
    assert(start >= 0);
    assert(start + length <= m_pu.size());
    return m_pu.segment(start, length);
  }

  
};

#endif  // TINY_METRICS_H
