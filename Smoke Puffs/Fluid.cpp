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

Fluid::Fluid(int width, int height) {
  w_ = width;
  h_ = height;
  
  fluxes_.resize(w_*h_*2);
  densities_.resize(w_*h_);
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      fluxes_[fidx(0,x,y)] = 0.0f;
      fluxes_[fidx(1,x,y)] = 0.0f;
      densities_[vidx(x,y)] = 0.0f;
    }
  }
  
  for (int i = 0; i < 10; ++i) {
    Project();
  }
  
  smoke_radius_ = 8;
}

void Fluid::Advect(float dt) {
  std::vector<float> new_fluxes(w_ * h_ * 2);
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      new_fluxes[fidx(0,x,y)] = 0.0f;
      new_fluxes[fidx(1,x,y)] = 0.0f;
    }
  }
  
  for (int x = 1; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      const Eigen::Vector2f mid_source = Eigen::Vector2f(x, y + 0.5f) - dt*0.5*InterpolateVelocity(Eigen::Vector2f(x, y + 0.5f));
      const Eigen::Vector2f source = Eigen::Vector2f(x, y + 0.5) - dt*InterpolateVelocity(ClipPoint(mid_source));
      new_fluxes[fidx(0,x,y)] = 0.996f*InterpolateXVelocity(ClipPoint(source));
    }
  }
  for (int x = 0; x < w_; ++x) {
    for (int y = 1; y < h_; ++y) {
      const Eigen::Vector2f mid_source = Eigen::Vector2f(x + 0.5f, y) - dt*0.5*InterpolateVelocity(Eigen::Vector2f(x + 0.5f, y));
      const Eigen::Vector2f source = Eigen::Vector2f(x + 0.5f, y) - dt*InterpolateVelocity(ClipPoint(mid_source));
      new_fluxes[fidx(1,x,y)] = 0.996f*InterpolateYVelocity(ClipPoint(source));
    }
  }
  fluxes_.swap(new_fluxes);
}

void Fluid::Project() {
  std::vector<float> pressure(w_ * h_);
  std::vector<float> div(w_ * h_);
  // For faster SOR convergence.
  const float omega = 1.2f;
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      int here = vidx(x, y);
      pressure[here] = 0.0f;
      if (x > 0) pressure[here] += fluxes_[fidx(0,x,y)];
      if (x < w_-1) pressure[here] -= fluxes_[fidx(0,x+1,y)];
      if (y > 0) pressure[here] += fluxes_[fidx(1,x,y)];
      if (y < h_-1) pressure[here] -= fluxes_[fidx(1,x,y+1)];
      div[here] = pressure[here];
    }
  }
  
  const int MAX_ITERS = 20;
  for (int k = 0; k < MAX_ITERS; ++k) {
    for (int x = 0; x < w_; ++x) {
      for (int y = 0; y < h_; ++y) {
        int here = vidx(x,y);
        float sigma = 0.0f;
        float count = 0.0f;
        if (x > 0) {
          sigma -= pressure[here-h_];
          count += 1.0f;
        }
        if (x < w_-1) {
          sigma -= pressure[here+h_];
          count += 1.0f;
        }
        if (y > 0) {
          sigma -= pressure[here-1];
          count += 1.0f;
        }
        if (y < h_-1) {
          sigma -= pressure[here+1];
          count += 1.0f;
        }
        pressure[here] = (1.0f - omega) * pressure[here] + omega / count * (div[here] - sigma);
      }
    }
  }
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      if (x > 0) fluxes_[fidx(0,x,y)] -= pressure[vidx(x,y)] - pressure[vidx(x-1,y)];
      if (y > 0) fluxes_[fidx(1,x,y)] -= pressure[vidx(x,y)] - pressure[vidx(x,y-1)];
    }
  }
}

void Fluid::AdvectPoint(float dt, float x0, float y0, float* xf, float* yf) {
  Eigen::Vector2f source(x0, y0);
  Eigen::Vector2f result = ClipPoint(source + dt * InterpolateVelocity(ClipPoint(source + 0.5 * dt * InterpolateVelocity(source))));
  *xf = result[0];
  *yf = result[1];
}

void Fluid::GetLines(std::vector<float>* line_coords, float scale) {
  line_coords->clear();
  line_coords->reserve(w_*h_*4);
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      line_coords->push_back(x + 0.5); 
      line_coords->push_back(y + 0.5); 
      line_coords->push_back(x + 0.5 + scale*CellCenterVelocity(x, y)[0]);
      line_coords->push_back(y + 0.5 + scale*CellCenterVelocity(x, y)[1]);
    }
  }
}

void Fluid::GetDensities(std::vector<float>* densities) {
  densities->clear();
  densities->reserve(densities_.size());
  for (int y = 0; y < h_; ++y) {
    for (int x= 0; x < w_; ++x) {
      densities->push_back(densities_[vidx(x,y)]);
    }
  }
}


void Fluid::AddImpulse(float x0, float y0, float vx, float vy) {
  pending_impulse_origins_.push_back(Eigen::Vector2f(x0, y0));
  pending_impulse_velocities_.push_back(Eigen::Vector2f(vx*50.0, vy*50.0));
}

void Fluid::ApplyImpulses() {
  static int call_count = 0;
  if (call_count < 3) {
    pending_impulse_origins_.clear();
    pending_impulse_velocities_.clear();
    ++call_count;
    return;
  }
  
  assert(pending_impulse_origins_.size() == pending_impulse_velocities_.size());
  
  for (size_t i = 0; i < pending_impulse_origins_.size(); ++i) {
    const Eigen::Vector2f& origin = pending_impulse_origins_[i];
    const Eigen::Vector2f& delta = pending_impulse_velocities_[i];
    
    int x = static_cast<int>(origin[0]);
    int y = static_cast<int>(origin[1]);
    
    float fx = origin[0] - x;
    float fy = origin[1] - y;
    
    if (x > 1 && x < w_-2 && y > 1 && y < h_-2) {
      fluxes_[fidx(0,x,y)] += (1.0f-fx) * delta[0];
      fluxes_[fidx(0,x+1,y)] += fx * delta[0];
      fluxes_[fidx(1,x,y)] += (1.0f-fy) * delta[1];
      fluxes_[fidx(1,x,y+1)] += fy * delta[1];
    }
    
    for (int k = -smoke_radius_; k <= smoke_radius_; ++k) {
      int xk = x + k;
      if (xk < 1 || xk >= w_ - 1) continue;
      for (int j = -smoke_radius_; j <= smoke_radius_; ++j) {
        int yj = y + j;
        if (k*k + j*j > smoke_radius_*smoke_radius_) continue;
        if (yj < 1 || yj >= h_ - 1) continue;
        densities_[vidx(xk,yj)] = std::max(densities_[vidx(xk, yj)], std::max(1.0f, 3.0f / (1 + sqrtf(k*k + j*j))));
      }
    }
  }
  
  pending_impulse_origins_.clear();
  pending_impulse_velocities_.clear();
}

void Fluid::AdvectDensity(float dt) {
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      const Eigen::Vector2f mid_source = Eigen::Vector2f(x + 0.5f, y + 0.5f) - dt*0.5*InterpolateVelocity(Eigen::Vector2f(x + 0.5f, y + 0.5f));
      const Eigen::Vector2f source = Eigen::Vector2f(x + 0.5f, y + 0.5f) - dt*InterpolateVelocity(ClipPoint(mid_source));
      densities_[vidx(x,y)] = 0.985f * InterpolateDensity(ClipPoint(source));
    }
  }
}