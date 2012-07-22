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
  
  void Step(double dt) {
    //ApplyImpulses();
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
  
  inline int fidx(int axis, int x, int y) { int res = axis * w_ * h_ + x * h_ + y;
    assert(res < fluxes_.size());
    assert(res >= 0);
    return res; }
  inline int vidx(int x, int y) { int res = x * h_ + y;
    assert(res < w_ * h_);
    assert(res >= 0);
    return res;}
  inline Eigen::Vector2f CellCenterVelocity(int x, int y) {
    x = std::max(x, 1);
    y = std::max(y, 1);
    int hx = std::min(x+1, w_-1);
    int hy = std::min(y+1, h_-1);
    return Eigen::Vector2f(0.5f * fluxes_[fidx(0,x,y)] + 0.5f * fluxes_[fidx(0,hx,y)],
                           0.5f * fluxes_[fidx(1,x,y)] + 0.5f * fluxes_[fidx(1,x,hy)]);
  }
  inline Eigen::Vector2f InterpolateVelocity(const Eigen::Vector2f& pos) {
    int x = static_cast<int>(floor(pos[0]));
    int y = static_cast<int>(floor(pos[1]));
    float fx = pos[0] - x;
    float fy = pos[1] - y;
    int lx = std::max(x, 1);
    int ly = std::max(y, 1);
    int hx = std::min(x+1, w_-1);
    int hy = std::min(y+1, h_-1);
    return Eigen::Vector2f((1.0f-fx) * fluxes_[fidx(0,lx,ly)] + fx * fluxes_[fidx(0,hx,ly)],
                           (1.0f-fy) * fluxes_[fidx(1,lx,ly)] + fy * fluxes_[fidx(1,lx,hy)]);
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
  std::vector<float> density_;

  std::vector<Eigen::Vector2f> pending_impulse_origins_;
  std::vector<Eigen::Vector2f> pending_impulse_deltas_;
  
  int w_;
  int h_;
};

#endif
