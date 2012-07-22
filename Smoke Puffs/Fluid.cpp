//
//  Fluid.cpp
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "Fluid.h"

#include <cmath>
#include <iostream>
#include <vector>

typedef Eigen::Triplet<double> Triplet;

Fluid::Fluid(int width, int height) {
  w_ = width;
  h_ = height;
  
  fluxes_.resize(w_*h_*2);
  vels_.clear();
  vels_.resize(w_*h_);
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      vels_[vidx(x,y)] = Eigen::Vector2f(0.0, 0.0);
    }
  }
  
  for (int x = 10; x < 20; ++x) {
    for (int y = 10; y < 20; ++y) {
      //vels_[vidx(x,y)] = Eigen::Vector2f(static_cast<float>(x) / w_, static_cast<float>(y) / h_);
      vels_[vidx(x,y)] = Eigen::Vector2f(12.0, 0.0);
    }
  }
  
  SetFluxesFromVelocities();
  
  for (int i = 0; i < 10; ++i) {
    Project();
  }
}

/*
void Fluid::BuildMatrices() {
  std::vector<Triplet> If_triplets;
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      If_triplets.push_back(Triplet(index(0, x, y),
                                    index(0, x-1, y),
                                    0.5));
      If_triplets.push_back(Triplet(index(0, x, y),
                                    index(0, x, y),
                                    0.5));
      If_triplets.push_back(Triplet(index(1, x, y),
                                    index(1, x, y-1),
                                    0.5));
      If_triplets.push_back(Triplet(index(1, x, y),
                                    index(1, x, y),
                                    0.5));

      Iv_triplets.push_back(Triplet(index(0, x, y),
                                    index(0, x+1, y),
                                    0.5));
      Iv_triplets.push_back(Triplet(index(0, x, y),
                                    index(0, x, y),
                                    0.5));
      Iv_triplets.push_back(Triplet(index(1, x, y),
                                    index(1, x, y+1),
                                    0.5));
      Iv_triplets.push_back(Triplet(index(1, x, y),
                                    index(1, x, y),
                                    0.5));


    }
  }
}
 */

void Fluid::Advect(float dt) {
  std::vector<Eigen::Vector2f> sources(w_ * h_);

  std::vector<Eigen::Vector2f> new_vels(w_ * h_);
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      const int here = vidx(x,y);
      //const Eigen::Vector2f mid_source = Eigen::Vector2f(x + 0.5, y + 0.5) - dt*0.5*vels_[here];
      const Eigen::Vector2f source = Eigen::Vector2f(x + 0.5, y + 0.5) - dt*vels_[here];
      const float rx = floorf(source[0] - 0.5);
      const float ry = floorf(source[1] - 0.5);
      const int lx = static_cast<int>(rx);
      const int ly = static_cast<int>(ry);
      const double fx = source[0] - rx - 0.5;
      const double fy = source[1] - ry - 0.5;
      new_vels[here] = ((1.0f-fy) * (vels_[vidx(lx,ly)]*(1.0f-fx) + vels_[vidx(lx+1, ly)]*fx) +
                        fy * (vels_[vidx(lx,ly+1)]*(1.0f-fx) + vels_[vidx(lx+1, ly+1)]*fx));
    }
  }

  vels_.swap(new_vels);
  
  SetFluxesFromVelocities();
}

void Fluid::Project() {
  std::vector<float> pressure(w_ * h_);
  std::vector<float> div(w_ * h_);
  const float omega = 1.2f;
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      pressure[vidx(x,y)] = fluxes_[fidx(0,x,y)] + fluxes_[fidx(1,x,y)] - fluxes_[fidx(0,x+1,y)] - fluxes_[fidx(1,x,y+1)];
      div[vidx(x,y)] = fluxes_[fidx(0,x,y)] + fluxes_[fidx(1,x,y)] - fluxes_[fidx(0,x+1,y)] - fluxes_[fidx(1,x,y+1)];
    }
  }
  
  for (int k = 0; k < 20; ++k) {
    float err = 0.0f;
    for (int x = 0; x < w_; ++x) {
      for (int y = 0; y < h_; ++y) {
        float diff = pressure[vidx(x,y)];
        const float sigma = -pressure[vidx(x-1,y)] - pressure[vidx(x+1,y)] - pressure[vidx(x,y-1)] - pressure[vidx(x,y+1)];
        pressure[vidx(x,y)] = (1.0f - omega) * pressure[vidx(x,y)] + omega * 0.25f * (div[vidx(x,y)] - sigma);
        diff -= pressure[vidx(x,y)];
        err += diff*diff;
      }
    }
  }
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      fluxes_[fidx(0,x,y)] -= pressure[vidx(x,y)] - pressure[vidx(x-1,y)];
      fluxes_[fidx(1,x,y)] -= pressure[vidx(x,y)] - pressure[vidx(x,y-1)];
    }
  }
  
  SetVelocitiesFromFluxes();
}

float Wrap(float x, float min, float max) {
  while (x < min) {
    x += (max - min);
  }
  while (x > max) {
    x += (min - max);
  }
  return x;
}

void Fluid::AdvectPoint(float dt, float x0, float y0, float* xf, float* yf) {
  Eigen::Vector2f source(x0, y0);
  Eigen::Vector2f result = source + dt * Interpolate((source + 0.5 * dt * Interpolate(source)));
  *xf = Wrap(result[0], 0.0, w_);
  *yf = Wrap(result[1], 0.0, h_);
}

void Fluid::GetLines(std::vector<float>* line_coords, float scale) {
  line_coords->clear();
  line_coords->reserve(w_*h_*4);
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      line_coords->push_back(x + 0.5); 
      line_coords->push_back(y + 0.5); 
      line_coords->push_back(x + 0.5 + scale*vels_[vidx(x,y)][0]);
      line_coords->push_back(y + 0.5 + scale*vels_[vidx(x,y)][1]);
    }
  }
}

void Fluid::AddImpulse(float x0, float y0, float dx, float dy) {
  pending_impulse_origins_.push_back(Eigen::Vector2f(x0, y0));
  pending_impulse_deltas_.push_back(Eigen::Vector2f(dx*50.0, dy*50.0));
}

void Fluid::ApplyImpulses() {
  return;
  assert(pending_impulse_origins_.size() == pending_impulse_deltas_.size());
  
  for (size_t i = 0; i < pending_impulse_origins_.size(); ++i) {
    const Eigen::Vector2f& origin = pending_impulse_origins_[i];
    const Eigen::Vector2f& delta = pending_impulse_deltas_[i];
    
    int x = static_cast<int>(origin[0]);
    int y = static_cast<int>(origin[1]);
    
    float fx = origin[0] - x;
    float fy = origin[1] - y;
    
    fluxes_[fidx(0,x,y)] += (1.0f-fx) * delta[0];
    fluxes_[fidx(0,x+1,y)] += fx * delta[0];
    fluxes_[fidx(1,x,y)] += (1.0f-fy) * delta[1];
    fluxes_[fidx(1,x,y+1)] += fy * delta[1];
  }
  
  pending_impulse_origins_.clear();
  pending_impulse_deltas_.clear();
  
  SetVelocitiesFromFluxes();
}

void Fluid::SetVelocitiesFromFluxes() {
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      vels_[vidx(x,y)][0] = 0.5 * (fluxes_[fidx(0, x, y)] + fluxes_[fidx(0, x+1, y)]);
      vels_[vidx(x,y)][1] = 0.5 * (fluxes_[fidx(1, x, y)] + fluxes_[fidx(1, x, y+1)]);
    }
  }
}

void Fluid::SetFluxesFromVelocities() {
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      fluxes_[fidx(0,x,y)] = 0.5 * (vels_[vidx(x-1,y)][0] + vels_[vidx(x,y)][0]);
      fluxes_[fidx(1,x,y)] = 0.5 * (vels_[vidx(x,y-1)][1] + vels_[vidx(x,y)][1]);
    }
  }
}

