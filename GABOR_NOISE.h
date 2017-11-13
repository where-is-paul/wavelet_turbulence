#ifndef GABOR_NOISE_H
#define GABOR_NOISE_H

/* Copyright (c) 2009 by Ares Lagae, Sylvain Lefebvre,
 * George Drettakis and Philip Dutre
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// -----------------------------------------------------------------------------

#include <climits>
#include <cmath>
#include <ctime>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

#include <VEC3.h>
#include <MERSENNETWISTER.h>

using namespace BasicVector;
using namespace std;

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

class pseudo_random_number_generator
{
public:
  void seed(unsigned s) { x_ = s; }
  unsigned operator()() { x_ *= 3039177861u; return x_; }
  float uniform_0_1() { return float(operator()()) / float(UINT_MAX); }
  float uniform(float min, float max)
    { return min + (uniform_0_1() * (max - min)); }
  unsigned poisson(float mean)
  {
    float g_ = std::exp(-mean);
    unsigned em = 0;
    double t = uniform_0_1();
    while (t > g_) {
      ++em;
      t *= uniform_0_1();
    }
    return em;
  }
private:
  unsigned x_;
};

static inline float gabor(float K, float a, float F_0, float omega_0, float omega_1, float x, float y, float z)
{
  float gaussian_envelop = K * std::exp(-M_PI * (a * a) * ((x * x) + (y * y) + (z * z)));
  float sinusoidal_carrier = std::cos(2.0 * M_PI * F_0 * 
    ((x * std::cos(omega_0) * std::sin(omega_1)) + 
     (y * std::sin(omega_0) * std::sin(omega_1)) + 
     (z * std::cos(omega_1))));
  return gaussian_envelop * sinusoidal_carrier;
}

static inline float gabor_derivative(float K, float a, float F_0, float omega_0, float omega_1, float x, float y, float z, int dir)
{
  float gaussian_envelop = K * std::exp(-M_PI * (a * a) * ((x * x) + (y * y) + (z * z)));
  float c_arg = 2.0 * M_PI * F_0 * ((x * std::cos(omega_0) * std::sin(omega_1)) + (y * std::sin(omega_0) * std::sin(omega_1)) + (z * std::cos(omega_1)));
  float sinusoidal_carrier = std::cos(c_arg);
  float derivative;
  if (dir == 0) {
    derivative = -2 * M_PI * (a * a * x + F_0 * cos(omega_0) * sin(omega_1) * tan(c_arg));
  } else if (dir == 1) {
    derivative = -2 * M_PI * (a * a * y + F_0 * sin(omega_0) * sin(omega_1) * tan(c_arg));
  } else {
    derivative = -2 * M_PI * (a * a * z + F_0 * cos(omega_1) * tan(c_arg));
  }
  return gaussian_envelop * sinusoidal_carrier * derivative;
}

static double g_t = 0;
class noise3d
{
public:
  noise3d(float K, float a, float F_0, float omega_0, float omega_1, float number_of_impulses_per_kernel, unsigned period, unsigned random_offset)
    :  K_(K), a_(a), F_0_(F_0), omega_0_(omega_0), omega_1_(omega_1), period_(period), random_offset_(random_offset)
  {
    kernel_radius_ = std::sqrt(-std::log(0.05) / M_PI) / a_;
    impulse_density_ = number_of_impulses_per_kernel / (M_PI * kernel_radius_ * kernel_radius_);
  }
  float operator()(float x, float y, float z) const
  {
    static std::map<float, std::map<float, std::map<float, float> > > memo;
    if (memo.count(x) && memo[x].count(y) && memo[x][y].count(z)) return memo[x][y][z];

    x /= kernel_radius_, y /= kernel_radius_, z /= kernel_radius_;
    float int_x = std::floor(x), int_y = std::floor(y), int_z = std::floor(z);
    float frac_x = x - int_x, frac_y = y - int_y, frac_z = z - int_z;
    int i = int(int_x), j = int(int_y), k = int(int_z);
    float noise = 0.0;
    for (int di = -1; di <= +1; ++di) {
      for (int dj = -1; dj <= +1; ++dj) {
        for (int dk = -1; dk <= +1; ++dk) {
          noise += cell(i + di, j + dj, k + dk, frac_x - di, frac_y - dj, frac_z - dk);
        }
      }
    }
    return memo[x][y][z] = noise;
  }
  float cell(int i, int j, int k, float x, float y, float z) const
  {
    unsigned s = ((((unsigned(k) % period_) * period_) * period_) + ((unsigned(j) % period_) * period_) + (unsigned(i) % period_)) + random_offset_; // periodic noise
    
    if (s == 0) s = 1;
    pseudo_random_number_generator prng;
    prng.seed(s);
    double number_of_impulses_per_cell = impulse_density_ * kernel_radius_ * kernel_radius_;
    unsigned number_of_impulses = prng.poisson(number_of_impulses_per_cell);
    float noise = 0.0;
    for (unsigned i = 0; i < 1.5 * number_of_impulses; ++i) {
      float x_i = prng.uniform_0_1();
      float y_i = prng.uniform_0_1();
      float z_i = prng.uniform_0_1();
      float w_i = prng.uniform(-1.0, +1.0);
      float omega_0_i = prng.uniform(0.0, 2.0 * M_PI);
      float omega_1_i = prng.uniform(0.0, M_PI);
      float x_i_x = x - x_i;
      float y_i_y = y - y_i;
      float z_i_z = z - z_i;
      if (((x_i_x * x_i_x) + (y_i_y * y_i_y) + (z_i_z * z_i_z)) < 1.0) {
        //noise += w_i * gabor(K_, a_, F_0_, omega_0_, omega_1_, x_i_x * kernel_radius_, y_i_y * kernel_radius_, z_i_z * kernel_radius_); // anisotropic
        noise += w_i * gabor(K_, a_, F_0_, omega_0_i, omega_1_i, x_i_x * kernel_radius_, y_i_y * kernel_radius_, z_i_z * kernel_radius_); // isotropic
      }
    }
    return noise;
  }
  // TODO: Analytically determine 3d variance
  float variance() const
  {
    float integral_gabor_filter_squared = ((K_ * K_) / (4.0 * a_ * a_)) * (1.0 + std::exp(-(2.0 * M_PI * F_0_ * F_0_) / (a_ * a_)));
    return impulse_density_ * (1.0 / 3.0) * integral_gabor_filter_squared;
  }

  float derivative(float x, float y, float z, int dir) const
  {
    x /= kernel_radius_, y /= kernel_radius_, z /= kernel_radius_;
    float int_x = std::floor(x), int_y = std::floor(y), int_z = std::floor(z);
    float frac_x = x - int_x, frac_y = y - int_y, frac_z = z - int_z;
    int i = int(int_x), j = int(int_y), k = int(int_z);
    float noise_d = 0.0;
    for (int di = -1; di <= +1; ++di) {
      for (int dj = -1; dj <= +1; ++dj) {
        for (int dk = -1; dk <= +1; ++dk) {
          noise_d += cell_derivative(i + di, j + dj, k + dk, frac_x - di, frac_y - dj, frac_z - dk, dir);
        }
      }
    }
    return noise_d;
  }
  float cell_derivative(int i, int j, int k, float x, float y, float z, int dir) const
  {
    unsigned s = ((((unsigned(k) % period_) * period_) * period_) + ((unsigned(j) % period_) * period_) + (unsigned(i) % period_)) + random_offset_; // periodic noise
    
    if (s == 0) s = 1;
    pseudo_random_number_generator prng;
    prng.seed(s);
    double number_of_impulses_per_cell = impulse_density_ * kernel_radius_ * kernel_radius_;
    unsigned number_of_impulses = prng.poisson(number_of_impulses_per_cell);
    float noise = 0.0;
    for (unsigned i = 0; i < 1.5 * number_of_impulses; ++i) {
      float x_i = prng.uniform_0_1();
      float y_i = prng.uniform_0_1();
      float z_i = prng.uniform_0_1();
      float w_i = prng.uniform(-1.0, +1.0);
      float omega_0_i = prng.uniform(0.0, 2.0 * M_PI);
      float omega_1_i = prng.uniform(0.0, M_PI);
      float x_i_x = x - x_i;
      float y_i_y = y - y_i;
      float z_i_z = z - z_i;
      if (((x_i_x * x_i_x) + (y_i_y * y_i_y) + (z_i_z * z_i_z)) < 1.0) {
        //noise += w_i * gabor_derivative(K_, a_, F_0_, omega_0_, omega_1_, x_i_x * kernel_radius_, y_i_y * kernel_radius_, z_i_z * kernel_radius_, dir); // anisotropic
        noise += w_i * gabor_derivative(K_, a_, F_0_, omega_0_i, omega_1_i, x_i_x * kernel_radius_, y_i_y * kernel_radius_, z_i_z * kernel_radius_, dir); // isotropic
      }
    }
    return noise;
  }

private:
  float K_;
  float a_;
  float F_0_;
  float omega_0_;
  float omega_1_;
  float kernel_radius_;
  float impulse_density_;
  unsigned period_;
  unsigned random_offset_;
};

/*
//////////////////////////////////////////////////////////////////////////////////////////
// x derivative of noise
//////////////////////////////////////////////////////////////////////////////////////////
static float GNoiseDx(Vec3 p, noise3d* data) { 
  return data->derivative(p[0], p[1], p[2], 0);
}

//////////////////////////////////////////////////////////////////////////////////////////
// y derivative of noise
//////////////////////////////////////////////////////////////////////////////////////////
static inline float GNoiseDy(Vec3 p, noise3d* data) { 
  return data->derivative(p[0], p[1], p[2], 1);
}

//////////////////////////////////////////////////////////////////////////////////////////
// z derivative of noise
//////////////////////////////////////////////////////////////////////////////////////////
static inline float GNoiseDz(Vec3 p, noise3d* data) { 
  return data->derivative(p[0], p[1], p[2], 2);
}
*/

#define NOISE_TILE_SIZE 128

// warning - noiseTileSize has to be 128^3!
#define modFast128(x) ((x) & 127)
#define modFast64(x)  ((x) & 63)

//////////////////////////////////////////////////////////////////////////////////////////
// x derivative of noise
//////////////////////////////////////////////////////////////////////////////////////////
static inline float GNoiseDx(Vec3 p, noise3d* data) { 
  int i, f[3], c[3], mid[3], n = NOISE_TILE_SIZE;
  float w[3][3], t, result = 0;
  
  mid[0] = ceil(p[0] - 0.5); 
  t = mid[0] - (p[0] - 0.5);
  w[0][0] = -t;
  w[0][2] = (1.f - t);
  w[0][1] = 2.0f * t - 1.0f;
  
  mid[1] = ceil(p[1] - 0.5); 
  t = mid[1] - (p[1] - 0.5);
  w[1][0] = t * t / 2; 
  w[1][2] = (1 - t) * (1 - t) / 2;
  w[1][1] = 1 - w[1][0] - w[1][2];

  mid[2] = ceil(p[2] - 0.5); 
  t = mid[2] - (p[2] - 0.5);
  w[2][0] = t * t / 2; 
  w[2][2] = (1 - t) * (1 - t)/2; 
  w[2][1] = 1 - w[2][0] - w[2][2];
 
  // to optimize, explicitly unroll this loop
  for (int z = -1; z <=1; z++)
    for (int y = -1; y <=1; y++)
      for (int x = -1; x <=1; x++)
      {
        float weight = 1.0f;
        c[0] = modFast128(mid[0] + x);
        weight *= w[0][x+1];
        c[1] = modFast128(mid[1] + y);
        weight *= w[1][y+1];
        c[2] = modFast128(mid[2] + z);
        weight *= w[2][z+1];
        result += weight * (*data)(c[0], c[1], c[2]);
      }
 return result;
}

//////////////////////////////////////////////////////////////////////////////////////////
// y derivative of noise
//////////////////////////////////////////////////////////////////////////////////////////
static inline float GNoiseDy(Vec3 p, noise3d* data) { 
  int i, f[3], c[3], mid[3], n=NOISE_TILE_SIZE; 
  float w[3][3], t, result =0;
  
  mid[0] = ceil(p[0] - 0.5); 
  t = mid[0]-(p[0] - 0.5);
  w[0][0] = t * t / 2; 
  w[0][2] = (1 - t) * (1 - t) / 2;
  w[0][1] = 1 - w[0][0] - w[0][2];
  
  mid[1] = ceil(p[1] - 0.5); 
  t = mid[1]-(p[1] - 0.5);
  w[1][0] = -t;
  w[1][2] = (1.f - t);
  w[1][1] = 2.0f * t - 1.0f;

  mid[2] = ceil(p[2] - 0.5); 
  t = mid[2] - (p[2] - 0.5);
  w[2][0] = t * t / 2; 
  w[2][2] = (1 - t) * (1 - t)/2; 
  w[2][1] = 1 - w[2][0] - w[2][2];
  
  // to optimize, explicitly unroll this loop
  for (int z = -1; z <=1; z++)
    for (int y = -1; y <=1; y++)
      for (int x = -1; x <=1; x++)
      {
        float weight = 1.0f;
        c[0] = modFast128(mid[0] + x);
        weight *= w[0][x+1];
        c[1] = modFast128(mid[1] + y);
        weight *= w[1][y+1];
        c[2] = modFast128(mid[2] + z);
        weight *= w[2][z+1];
        result += weight * (*data)(c[0], c[1], c[2]);
      }

  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////
// z derivative of noise
//////////////////////////////////////////////////////////////////////////////////////////
static inline float GNoiseDz(Vec3 p, noise3d* data) { 
  int i, f[3], c[3], mid[3], n=NOISE_TILE_SIZE; 
  float w[3][3], t, result =0;

  mid[0] = ceil(p[0] - 0.5); 
  t = mid[0]-(p[0] - 0.5);
  w[0][0] = t * t / 2; 
  w[0][2] = (1 - t) * (1 - t) / 2;
  w[0][1] = 1 - w[0][0] - w[0][2];
  
  mid[1] = ceil(p[1] - 0.5); 
  t = mid[1]-(p[1] - 0.5);
  w[1][0] = t * t / 2; 
  w[1][2] = (1 - t) * (1 - t) / 2;
  w[1][1] = 1 - w[1][0] - w[1][2];

  mid[2] = ceil(p[2] - 0.5); 
  t = mid[2] - (p[2] - 0.5);
  w[2][0] = -t;
  w[2][2] = (1.f - t);
  w[2][1] = 2.0f * t - 1.0f;

  // to optimize, explicitly unroll this loop
  for (int z = -1; z <=1; z++)
    for (int y = -1; y <=1; y++)
      for (int x = -1; x <=1; x++)
      {
        float weight = 1.0f;
        c[0] = modFast128(mid[0] + x);
        weight *= w[0][x+1];
        c[1] = modFast128(mid[1] + y);
        weight *= w[1][y+1];
        c[2] = modFast128(mid[2] + z);
        weight *= w[2][z+1];
        result += weight * (*data)(c[0], c[1], c[2]);
      }
  return result;
}
#endif