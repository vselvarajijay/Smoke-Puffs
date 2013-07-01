//
//  Fluid.h
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Smoke_Puffs_Fluid_h
#define Smoke_Puffs_Fluid_h

#include <sys/time.h>
#include <vector>

#import "TargetConditionals.h"

#if !TARGET_IPHONE_SIMULATOR
#import <arm_neon.h>
#endif  // !TARGET_IPHONE_SIMULATOR

#include <Eigen/Dense>

class Timer {
public:
  explicit Timer(const std::string& label);
  
  void BeginEvent();
  void EndEvent();
  void Print();
  
private:
  std::string label_;
  float total_time_;
  int total_events_;
  struct timeval this_event_start_;
};

class Fluid {
 public:
  Fluid(int width, int height);
  ~Fluid() {}
  
  void Advect(float dt);
  void ConfineVorticity(float dt);
  void Project();
  void FixBoundaries();

  void ApplyImpulses();
  void AddImpulse(float x0, float y0,
                  float vx, float vy);
  
  void GetLines(std::vector<float>* line_coords, float scale);
  void GetDensities(std::vector<float>* densities);
  
  void Step(double dt) {
    advection_timer_.BeginEvent();
    Advect(dt);
    advection_timer_.EndEvent();
    
    vorticity_timer_.BeginEvent();
    ConfineVorticity(dt);
    vorticity_timer_.EndEvent();
    
    impulse_timer_.BeginEvent();
    ApplyImpulses();
    impulse_timer_.EndEvent();
    
    projection_timer_.BeginEvent();
    Project();
    projection_timer_.EndEvent();
    
    density_timer_.BeginEvent();
    AdvectDensity(dt);
    density_timer_.EndEvent();
    
    impulse_timer_.Print();
    advection_timer_.Print();
    vorticity_timer_.Print();
    projection_timer_.Print();
    density_timer_.Print();
  }
  
  void AdvectDensity(float dt);
  void AdvectPoint(float dt, float x0, float y0, float* xf, float* yf);

  void set_smoke_radius(int radius) { smoke_radius_ = radius; }
  
 private:
  inline int Clip(int x, int lower, int upper) { return std::max(lower, std::min(upper, x)); }
  inline float Clip(float x, float lower, float upper) { return std::max(lower, std::min(upper, x)); }
  
  inline Eigen::Vector2f ClipPoint(const Eigen::Vector2f& src) {
    return Eigen::Vector2f(Clip(src[0], 2.0f, static_cast<float>(w_-2)-1e-8), Clip(src[1], 2.0f, static_cast<float>(h_-2)-1e-8));
  }
  
  inline int fidx(int axis, int x, int y) {
    int res = axis * w_ * h_ + x * h_ + y;
#ifdef DEBUG
    assert(res >= 0);
    assert(res < w_ * h_ * 2);
#endif
    return res;
  }
  inline int vidx(int x, int y) {
    int res = x * h_ + y;
#ifdef DEBUG
    assert(res >= 0);
    assert(res < w_ * h_);
#endif
    return res;
  }
  
  void InterpolateXVelocities(const std::vector<float>& xs,
                              const std::vector<float>& ys,
                              float* vx);
  
  void InterpolateYVelocities(const std::vector<float>& xs,
                              const std::vector<float>& ys,
                              float* result) {
    for (int i = 0; i < xs.size(); ++i) {
      result[i] = InterpolateYVelocity(Eigen::Vector2f(xs[i], ys[i]));
    }
  }
  
  void InterpolateVelocities(const std::vector<float>& xs,
                             const std::vector<float>& ys,
                             float* x_result,
                             float* y_result) {
    InterpolateXVelocities(xs, ys, x_result);
    InterpolateYVelocities(xs, ys, y_result);
  }
  
  inline float BilinearInterp(float x,
                              float y,
                              const float* source) {
#if !TARGET_IPHONE_SIMULATOR
    int sx = static_cast<int>(x);
    int sy = static_cast<int>(y);
    float fx = x - sx;
    float fy = y - sy;
    int here = vidx(sx,sy);
    float lowx = source[here] + fy*(source[here+1] - source[here]);
    float highx = source[here+h_] + fy*(source[here+h_+1] - source[here+h_]);
    return lowx + fx*(highx - lowx);
#else
    int sx = static_cast<int>(x);
    int sy = static_cast<int>(y);
    float fx = x - sx;
    float fy = y - sy;
    int here = vidx(sx,sy);
    float lowx = source[here] + fy*(source[here+1] - source[here]);
    float highx = source[here+h_] + fy*(source[here+h_+1] - source[here+h_]);
    return lowx + fx*(highx - lowx);
#endif  // !TARGET_IPHONE_SIMULATOR
  }
  
  inline float MonotoneCubicInterp(float fm1,
                                   float f0,
                                   float f1,
                                   float f2,
                                   float t) {    
    float delta_k = f1 - f0;
    float d_k = 0.5f * (f1 - fm1);
    float d_k1 = 0.5f * (f2 - f0);
    
    if (delta_k > 0.0f) {
      d_k = std::max(d_k, 0.0f);
      d_k1 = std::max(d_k1, 0.0f);
    } else if (delta_k < 0.0f) {
      d_k = std::min(d_k, 0.0f);
      d_k1 = std::min(d_k1, 0.0f);
    } else {
      d_k = 0.0f;
      d_k1 = 0.0f;
    }
    
    float a3 = d_k + d_k1 - 2.0f*delta_k;
    float a2 = 3.0f*delta_k - 2.0f*d_k - d_k1;
    float a1 = d_k;
    float a0 = f0;
    
    float result = a3*t*t*t + a2*t*t + a1*t + a0;
//    if (f0 < f1) {
//      assert(f0 <= result);
//      assert(result <= f1);
//    } else {
//      assert(f1 <= result);
//      assert(result <= f0);
//    }
    return result;
  }
  
  inline float BicubicInterp(float x,
                             float y,
                             const float* source) {
    int sx = static_cast<int>(x);
    int sy = static_cast<int>(y);
    float fx = x - sx;
    float fy = y - sy;
    int here = vidx(sx,sy);
    return MonotoneCubicInterp(MonotoneCubicInterp(source[here-h_-1], source[here-h_], source[here-h_+1], source[here-h_+2], fy),
                               MonotoneCubicInterp(source[here-1], source[here], source[here+1], source[here+2], fy),
                               MonotoneCubicInterp(source[here+h_-1], source[here+h_], source[here+h_+1], source[here+h_+2], fy),
                               MonotoneCubicInterp(source[here+2*h_-1], source[here+2*h_], source[here+2*h_+1], source[here+2*h_+2], fy),
                               fx
                               );
  }
  
  inline float InterpolateXVelocity(const Eigen::Vector2f& source) {
    return BilinearInterp(source[0], source[1] - 0.5f, &fluxes_[0]);
  }
  inline float InterpolateYVelocity(const Eigen::Vector2f& source) {
    return BilinearInterp(source[0] - 0.5f, source[1], &fluxes_[0] + w_*h_);
  }
  inline float InterpolateDensity(const Eigen::Vector2f& source) {
    return BilinearInterp(source[0] - 0.5f, source[1] - 0.5f, &densities_[0]);
  }
  inline Eigen::Vector2f InterpolateVelocity(const Eigen::Vector2f& source) {
    return Eigen::Vector2f(InterpolateXVelocity(source), InterpolateYVelocity(source));
  }
  inline Eigen::Vector2f CellCenterVelocity(int x, int y) {
    return InterpolateVelocity(Eigen::Vector2f(x + 0.5f, y + 0.5f));
  }
  
  Timer impulse_timer_;
  Timer advection_timer_;
  Timer vorticity_timer_;
  Timer projection_timer_;
  Timer density_timer_;

  std::vector<float> fluxes_;
  std::vector<float> densities_;
  std::vector<float> cell_centered_vels_;
  std::vector<float> vorticities_;
  std::vector<float> fluid_;
  std::vector<float> fluid_face_;

  std::vector<Eigen::Vector2f> pending_impulse_origins_;
  std::vector<Eigen::Vector2f> pending_impulse_velocities_;
  
  int w_;
  int h_;
  
  int smoke_radius_;

};

#endif
