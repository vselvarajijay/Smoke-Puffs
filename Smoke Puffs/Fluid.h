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
  void Project();
  void FixBoundaries();

  void ApplyImpulses();
  void AddImpulse(float x0, float y0,
                  float vx, float vy);
  
  void GetLines(std::vector<float>* line_coords, float scale);
  void GetDensities(std::vector<float>* densities);
  
  void Step(double dt) {
    impulse_timer_.BeginEvent();
    ApplyImpulses();
    impulse_timer_.EndEvent();

    advection_timer_.BeginEvent();
    Advect(dt);
    advection_timer_.EndEvent();
    
    projection_timer_.BeginEvent();
    Project();
    projection_timer_.EndEvent();
    
    density_timer_.BeginEvent();
    AdvectDensity(dt);
    density_timer_.EndEvent();
    
    impulse_timer_.Print();
    advection_timer_.Print();
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
    return Eigen::Vector2f(Clip(src[0], 1.0f, static_cast<float>(w_-1)-1e-6), Clip(src[1], 1.0f, static_cast<float>(h_-1)-1e-6));
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
  inline float InterpolateXVelocity(const Eigen::Vector2f& source) {
    int sx = static_cast<int>(source[0]);
    int sy = static_cast<int>(source[1] - 0.5f);
    float fx = source[0] - sx;
    float fy = source[1] - 0.5f - sy;
    float flx = (1.0f - fx);
    float fhx = fx;
    float fly = (1.0f - fy);
    float fhy = fy;
    float result = 0.0f;
    float count = 1e-20f;
    float weight = 0.0f;
    int here = fidx(0,sx,sy);
    weight = flx*fly*fluid_face_[here];
    result += weight*fluxes_[here];
    count += weight;
    weight = flx*fhy*fluid_face_[here+1];
    result += weight*fluxes_[here+1];
    count += weight;
    weight = fhx*fly*fluid_face_[here+h_];
    result += weight*fluxes_[here+h_];
    count += weight;
    weight = flx*fly*fluid_face_[here+h_+1];
    result += weight*fluxes_[here+h_+1];
    count += weight;
    return result / count;
  }
  inline float InterpolateYVelocity(const Eigen::Vector2f& source) {
    int sx = static_cast<int>(source[0] - 0.5f);
    int sy = static_cast<int>(source[1]);
    float fx = source[0] - 0.5f - sx;
    float fy = source[1] - sy;
    float flx = (1.0f - fx);
    float fhx = fx;
    float fly = (1.0f - fy);
    float fhy = fy;
    float result = 0.0f;
    float count = 1e-20f;
    float weight = 0.0f;
    int here = fidx(1,sx,sy);
    weight = flx*fly*fluid_face_[here];
    result += weight*fluxes_[here];
    count += weight;
    weight = flx*fhy*fluid_face_[here+1];
    result += weight*fluxes_[here+1];
    count += weight;
    weight = fhx*fly*fluid_face_[here+h_];
    result += weight*fluxes_[here+h_];
    count += weight;
    weight = flx*fly*fluid_face_[here+h_+1];
    result += weight*fluxes_[here+h_+1];
    count += weight;
    return result / count;
  }
  inline float InterpolateDensity(const Eigen::Vector2f& source) {
    int sx = static_cast<int>(source[0] - 0.5f);
    int sy = static_cast<int>(source[1] - 0.5f);
    float fx = source[0] - 0.5f - sx;
    float fy = source[1] - 0.5f - sy;
    float flx = (1.0f - fx);
    float fhx = fx;
    float fly = (1.0f - fy);
    float fhy = fy;
    float result = 0.0f;
    float count = 1e-20f;
    float weight = 0.0f;
    weight = flx*fly*fluid_[vidx(sx,sy)];
    result += weight*densities_[vidx(sx,sy)];
    count += weight;
    weight = flx*fhy*fluid_[vidx(sx,sy+1)];
    result += weight*densities_[vidx(sx,sy+1)];
    count += weight;
    weight = fhx*fly*fluid_[vidx(sx+1,sy)];
    result += weight*densities_[vidx(sx+1,sy)];
    count += weight;
    weight = fhx*fhy*fluid_[vidx(sx+1,sy+1)];
    result += weight*densities_[vidx(sx+1,sy+1)];
    count += weight;
    return result / count;
  }
  inline Eigen::Vector2f InterpolateVelocity(const Eigen::Vector2f& source) {
    return Eigen::Vector2f(InterpolateXVelocity(source), InterpolateYVelocity(source));
  }
  inline Eigen::Vector2f CellCenterVelocity(int x, int y) {
    return InterpolateVelocity(Eigen::Vector2f(x + 0.5f, y + 0.5f));
  }
  inline void SplatCenterVelocity(int x, int y,
                                  const Eigen::Vector2f& vel,
                                  std::vector<float>* new_fluxes) {
    if (x > 1) (*new_fluxes)[fidx(0, x, y)] += 0.5 * vel[0];
    if (x < w_-2) (*new_fluxes)[fidx(0, x+1, y)] += 0.5 * vel[0];
    if (y > 1) (*new_fluxes)[fidx(1, x, y)] += 0.5 * vel[1];
    if (y < h_-2) (*new_fluxes)[fidx(1, x, y+1)] += 0.5 * vel[1];
  }
  
  Timer impulse_timer_;
  Timer advection_timer_;
  Timer projection_timer_;
  Timer density_timer_;

  std::vector<float> fluxes_;
  std::vector<float> densities_;
  std::vector<float> fluid_;
  std::vector<float> fluid_face_;

  std::vector<Eigen::Vector2f> pending_impulse_origins_;
  std::vector<Eigen::Vector2f> pending_impulse_velocities_;
  
  int w_;
  int h_;
  
  int smoke_radius_;

};

#endif
