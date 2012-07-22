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
  //void Diffuse(double dt, double nu);
  void Project();
  //void FixBoundaries();

  void ApplyImpulses();
  void AddImpulse(float x0, float y0,
                  float dx, float dy);
  
  void GetLines(std::vector<float>* line_coords, float scale);
  
  void Step(double dt) {
    //ApplyImpulses();
    Advect(dt);
    Project();
  }
  
  void AdvectPoint(float dt, float x0, float y0, float* xf, float* yf);

 private:
  inline int fidx(int axis, int x, int y) { return axis * w_ * h_ + ((x+w_)%w_) * h_ + (y+h_)%h_;}
  inline int vidx(int x, int y) { return ((x+w_)%w_) * h_ + (y+h_)%h_; }
  inline Eigen::Vector2f Interpolate(const Eigen::Vector2f& pos) {
    int x = static_cast<int>(floor(pos[0]));
    int y = static_cast<int>(floor(pos[1]));
    float fx = pos[0] - x;
    float fy = pos[1] - y;
    return Eigen::Vector2f((1.0f-fx) * fluxes_[fidx(0,x,y)] + fx * fluxes_[fidx(0,x+1,y)],
                           (1.0f-fy) * fluxes_[fidx(1,x,y)] + fy * fluxes_[fidx(1,x,y+1)]);
  }

  void SetVelocitiesFromFluxes();
  void SetFluxesFromVelocities();
  
  //Eigen::SparseMatrix If_;
  //Eigen::SparseMatrix Iv_;
  //Eigen::SparseMatrix A_;
  //Eigen::SparseMatrix P_;

  Eigen::VectorXf fluxes_;
  std::vector<Eigen::Vector2f> vels_;
  std::vector<float> density_;

  std::vector<Eigen::Vector2f> pending_impulse_origins_;
  std::vector<Eigen::Vector2f> pending_impulse_deltas_;
  
  int w_;
  int h_;
};

#endif
