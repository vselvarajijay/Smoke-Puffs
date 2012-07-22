//
//  Fluid.h
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Smoke_Puffs_Fluid_h
#define Smoke_Puffs_Fluid_h

#include "Eigen/Dense"
#include "Eigen/SparseCore"

class Fluid {
 public:
  Fluid(int width, int height);
  ~Fluid() {}

  //void BuildMatrices();

  void Advect(float dt);
  void Project();
  void FixBoundaries();

  void ApplyImpulses();
  void AddImpulse(float x0, float y0,
                  float dx, float dy);
  
  void GetLines(std::vector<float>* line_coords, float scale);
  void GetDensities(std::vector<float>* densities);
  
  void Step(double dt) {
    ApplyImpulses();
    Advect(dt);
    Project();
    AdvectDensity(dt);
  }
  
  void AdvectDensity(float dt);
  void AdvectPoint(float dt, float x0, float y0, float* xf, float* yf);

 private:
  inline int Clip(int x, int lower, int upper) { return std::max(lower, std::min(upper, x)); }
  inline float Clip(float x, float lower, float upper) { return std::max(lower, std::min(upper, x)); }
  
  inline Eigen::Vector2f ClipPoint(const Eigen::Vector2f& src) {
    return Eigen::Vector2f(Clip(src[0], 0.0f, w_+1.0f), Clip(src[1], 0.0f, h_+1.0f));
  }
  
  inline int fidx(int axis, int x, int y) { return axis * w_ * h_ + x * h_ + y; }
  inline int vidx(int x, int y) { return x * h_ + y; }
  inline float InterpolateXVelocity(const Eigen::Vector2f& source) {
    int sx = static_cast<int>(source[0]);
    int sy = static_cast<int>(source[1] - 0.5f);
    float fx = source[0] - sx;
    float fy = source[1] - 0.5f - sy;
    int lx = std::max(1, sx);
    int ly = std::max(1, sy);
    int hx = std::min(w_-1, sx+1);
    int hy = std::min(h_-1, sy+1);
    return (1.0f-fx)*(1.0-fy)*fluxes_[fidx(0,lx,ly)] + (1.0-fx)*fy*fluxes_[fidx(0,lx,hy)] +
    fx*(1.0-fy)*fluxes_[fidx(0,hx,ly)] + fx*fy*fluxes_[fidx(0,hx,hy)];
  }
  inline float InterpolateYVelocity(const Eigen::Vector2f& source) {
    int sx = static_cast<int>(source[0] - 0.5f);
    int sy = static_cast<int>(source[1]);
    float fx = source[0]  - 0.5f - sx;
    float fy = source[1] - sy;
    int lx = std::max(1, sx);
    int ly = std::max(1, sy);
    int hx = std::min(w_-1, sx+1);
    int hy = std::min(h_-1, sy+1);
    return (1.0f-fx)*(1.0-fy)*fluxes_[fidx(1,lx,ly)] + (1.0-fx)*fy*fluxes_[fidx(1,lx,hy)] +
    fx*(1.0-fy)*fluxes_[fidx(1,hx,ly)] + fx*fy*fluxes_[fidx(1,hx,hy)];
  }
  inline float InterpolateDensity(const Eigen::Vector2f& source) {
    int sx = static_cast<int>(source[0] - 0.5f);
    int sy = static_cast<int>(source[1] - 0.5f);
    float fx = source[0] - 0.5f - sx;
    float fy = source[1] - 0.5f - sy;
    int lx = std::max(1, sx);
    int ly = std::max(1, sy);
    int hx = std::min(w_-1, sx+1);
    int hy = std::min(h_-1, sy+1);
    return (1.0f-fx)*(1.0-fy)*densities_[vidx(lx,ly)] + (1.0-fx)*fy*densities_[vidx(lx,hy)] +
    fx*(1.0-fy)*densities_[vidx(hx,ly)] + fx*fy*densities_[vidx(hx,hy)];
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
    /*
    int x = static_cast<int>(floor(pos[0]));
    int y = static_cast<int>(floor(pos[1]));
    float fx = pos[0] - x;
    float fy = pos[1] - y;
     */
    if (x > 0) (*new_fluxes)[fidx(0, x, y)] += 0.5 * vel[0];
    if (x < w_-1) (*new_fluxes)[fidx(0, x+1, y)] += 0.5 * vel[0];
    if (y > 0) (*new_fluxes)[fidx(1, x, y)] += 0.5 * vel[1];
    if (y < h_-1) (*new_fluxes)[fidx(1, x, y+1)] += 0.5 * vel[1];
  }

  std::vector<float> fluxes_;
  std::vector<float> densities_;

  std::vector<Eigen::Vector2f> pending_impulse_origins_;
  std::vector<Eigen::Vector2f> pending_impulse_deltas_;
  
  int w_;
  int h_;
};

#endif
