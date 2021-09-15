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

#ifndef DOUBLE_UTILS_H
#define DOUBLE_UTILS_H

#define _USE_MATH_DEFINES 1
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "math.h"

struct DoubleUtils {
  static double zero() { return 0.; }
  static double one() { return 1.; }

  static double two() { return 2.; }
  static double half() { return 0.5; }
  static double pi() { return M_PI; }
  static double half_pi() { return M_PI / 2.; }

  static double cos1(double v) { return ::cos(v); }
  static double sin1(double v) { return ::sin(v); }
  static double atan2(double dy, double dx) { return ::atan2(dy, dx); }
  static double asin(double v) { return ::asin(v); }
  static double copysign(double x, double y) { return ::copysign(x, y); }
  static double abs(double v) { return ::fabs(v); }
  static double pow(double a, double b) { return ::pow(a, b); }
  static double exp(double v) { return ::exp(v); }
  static double log(double v) { return ::log(v); }
  static double tanh(double v) { return ::tanh(v); }
  static double min1(double a, double b) { return a < b ? a : b; }
  static double max1(double a, double b) { return a > b ? a : b; }
  static double max(double a, double b) { return a > b ? a : b; }
  static double min(double a, double b) { return a < b ? a : b; }

  template <class T>
  static T sqrt1(T v) {
    // std::cout << "sqrt1 " << v << " " << sqrt(v) << "\n";
    return sqrt(v);
  }


  static bool getBool(bool v) {  return v; }

  // template <class T>
  // static double getDouble(T v) {
  //   return (double)v;
  // }

  static double getDouble(double v) {
    // std::cout << "getDouble " << v << " " << v << "\n";
    return v;
  }

  template <class T>
  static double convert(T) = delete;  // C++11

  static double convert(int value) { 
    // std::cout << "convert " << value << " " << double(value) << "\n";
    return double(value); 
  }

  // template <class T>
  // static double fraction(T, T) = delete;  // C++11

  static double fraction(int num, int denom) {
    // std::cout << "fraction " << num << " " << denom << " " << double(num) / double(denom) << "\n";
    return double(num) / double(denom);
  }

  static double scalar_from_string(const std::string& txt) {
    double result = atof(txt.c_str());
    // std::cout << "scalar_from_string " << txt << " " << result << "\n";
    return result;
  }

  static double scalar_from_double(double value) { return value; }

  static void FullAssert(bool a) {
    if (!a) {
      printf("!");
      assert(0);
      exit(0);
    }
  }
};

#endif  // DOUBLE_UTILS_H
