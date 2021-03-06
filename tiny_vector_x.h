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

#ifndef TINY_VECTOR_X_H
#define TINY_VECTOR_X_H

#include <stdio.h>

#include <vector>

/**
 * Represents a vector with arbitrary number of dimensions.
 */
template <typename TinyScalar, typename TinyConstants>
class TinyVectorX {
  TinyScalar* m_data{nullptr};

 public:
  int m_size{0};

  TinyVectorX() = default;

  explicit TinyVectorX(int size) : m_size(size) {
    m_data = new TinyScalar[m_size];
    set_zero();
  }
  TinyVectorX(int size, TinyScalar* data) : m_size(size), m_data(data) {}

  TinyVectorX(const TinyVectorX& v) : m_size(v.m_size) {
    m_data = new TinyScalar[m_size];
    for (int i = 0; i < m_size; ++i) m_data[i] = v.m_data[i];
  }

  TinyVectorX(const std::vector<TinyScalar>& v) : m_size(v.size()) {
    m_data = new TinyScalar[m_size];
    for (int i = 0; i < m_size; ++i) {
      m_data[i] = v.at(i);
    }
  }

  inline TinyScalar* data() { return m_data; }
  inline const TinyScalar* data() const { return m_data; }

  TinyVectorX& operator=(const TinyVectorX& v) {
    delete[] m_data;
    m_size = v.m_size;
    m_data = new TinyScalar[m_size];
    for (int i = 0; i < m_size; ++i) m_data[i] = v.m_data[i];
    return *this;
  }

  virtual ~TinyVectorX() { delete[] m_data; }

  void set_zero() {
    for (int i = 0; i < m_size; ++i) m_data[i] = TinyConstants::zero();
  }
  void resize(int n) {
    if (m_size != n) {
      m_size = n;
      delete[] m_data;
      m_data = new TinyScalar[m_size];
    }
    set_zero();
  }

  std::vector<TinyScalar> std() const {
    std::vector<TinyScalar> v(m_data, m_data + m_size);
    return v;
  }

  template <template <typename, typename> typename VectorType>
  inline TinyScalar dot(
      const VectorType<TinyScalar, TinyConstants>& other) const {
    assert(m_size == other.m_size);
    TinyScalar res = TinyConstants::zero();
    for (int i = 0; i < m_size; ++i) res = res + m_data[i] * other[i];
    return res;
  }

  inline TinyScalar sqnorm() const {
    TinyScalar res = TinyConstants::zero();
    for (int i = 0; i < m_size; ++i) res += m_data[i] * m_data[i];
    return res;
  }

  inline TinyScalar length() const {
    TinyScalar res = sqnorm();
    res = TinyConstants::sqrt1(res);
    return res;
  }

  inline int size() const {
    return m_size;
  }

  template <template <typename, typename> typename VectorType>
  inline TinyVectorX& operator+=(
      const VectorType<TinyScalar, TinyConstants>& v) {
    assert(m_size == v.m_size);
    for (int i = 0; i < m_size; ++i) m_data[i] = m_data[i] + v[i];
    return *this;
  }
  template <template <typename, typename> typename VectorType>
  inline TinyVectorX& operator-=(
      const VectorType<TinyScalar, TinyConstants>& v) {
    assert(m_size == v.m_size);
    for (int i = 0; i < m_size; ++i) m_data[i] -= v[i];
    return *this;
  }

  inline TinyVectorX operator-() const {
    TinyVectorX v(m_size);
    for (int i = 0; i < m_size; ++i) v.m_data[i] = -m_data[i];
    return v;
  }

  inline TinyScalar& operator[](int i) {
    TinyConstants::FullAssert(0 <= i && i < m_size);
    return m_data[i];
  }
  inline const TinyScalar& operator[](int i) const {
    TinyConstants::FullAssert(0 <= i && i < m_size);
    return m_data[i];
  }

  inline friend TinyVectorX operator+(const TinyVectorX& a, TinyScalar s) {
    TinyVectorX v(a.m_size);
    for (int i = 0; i < a.m_size; ++i) v.m_data[i] = a.m_data[i] + s;
    return v;
  }
  inline friend TinyVectorX operator*(const TinyVectorX& a, TinyScalar s) {
    TinyVectorX v(a.m_size);
    for (int i = 0; i < a.m_size; ++i) v.m_data[i] = a.m_data[i] * s;
    return v;
  }
  inline friend TinyVectorX operator/(const TinyVectorX& a, TinyScalar s) {
    TinyVectorX v(a.m_size);
    for (int i = 0; i < a.m_size; ++i) v.m_data[i] = a.m_data[i] / s;
    return v;
  }
  inline friend TinyVectorX operator+(TinyScalar s, const TinyVectorX& a) {
    return a + s;
  }
  inline friend TinyVectorX operator*(TinyScalar s, const TinyVectorX& a) {
    return a * s;
  }
  inline friend TinyVectorX operator/(TinyScalar s, const TinyVectorX& a) {
    return a / s;
  }

  template <template <typename, typename> typename VectorType>
  inline friend TinyVectorX operator+(
      const TinyVectorX& a, const VectorType<TinyScalar, TinyConstants>& b) {
    assert(a.m_size == b.m_size);
    TinyVectorX v(a.m_size);
    for (int i = 0; i < a.m_size; ++i) v.m_data[i] = a.m_data[i] + b[i];
    return v;
  }
  template <template <typename, typename> typename VectorType>
  inline friend TinyVectorX operator-(
      const TinyVectorX& a, const VectorType<TinyScalar, TinyConstants>& b) {
    assert(a.m_size == b.m_size);
    TinyVectorX v(a.m_size);
    for (int i = 0; i < a.m_size; ++i) v.m_data[i] = a.m_data[i] - b[i];
    return v;
  }
  template <template <typename, typename> typename VectorType>
  inline friend TinyVectorX operator*(
      const TinyVectorX& a, const VectorType<TinyScalar, TinyConstants>& b) {
    assert(a.m_size == b.m_size);
    TinyVectorX v(a.m_size);
    for (int i = 0; i < a.m_size; ++i) v.m_data[i] = a.m_data[i] * b[i];
    return v;
  }

  void print(const char* txt) const {
    printf("%s\n", txt);
    for (int i = 0; i < m_size; ++i) {
      double v = TinyConstants::getDouble(m_data[i]);
      printf("%f, ", v);
    }
    printf("\n");
  }

  TinyVectorX segment(int start, int length) const {
    assert(start >= 0);
    assert(start + length <= m_size);
    TinyVectorX v(length);
    for (int i = 0; i < length; ++i) {
      v[i] = m_data[i + start];
    }
    return v;
  }


  void assign_vector_add(int start, const TinyVectorX& v) {
    TinyConstants::FullAssert(0 <= start);
    TinyConstants::FullAssert(start + v.m_size <= m_size);

    for (int j = 0; j < v.m_size; ++j) {
      m_data[start + j] += v[j];
    }
  }
};

#endif  // TINY_VECTOR_X_H
