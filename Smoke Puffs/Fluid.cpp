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
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      fluxes_[fidx(0,x,y)] = 0.0f;
      fluxes_[fidx(1,x,y)] = 0.0f;
    }
  }
  
  for (int x = 20; x < 30; ++x) {
    for (int y = 20; y < 30; ++y) {
      SplatCenterVelocity(x, y, Eigen::Vector2f(12.0, 2.0), &fluxes_);
    }
  }
  
  for (int i = 0; i < 10; ++i) {
    Project();
  }
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
      int sx = static_cast<int>(source[0]);
      int sy = static_cast<int>(source[1] - 0.5f);
      float fx = source[0] - sx;
      float fy = source[1] - 0.5f - sy;
      int lx = std::max(1, sx);
      int ly = std::max(1, sy);
      int hx = std::min(w_-1, sx+1);
      int hy = std::min(h_-1, sy+1);
      new_fluxes[fidx(0,x,y)] = (1.0f-fx)*(1.0-fy)*fluxes_[fidx(0,lx,ly)] + (1.0-fx)*fy*fluxes_[fidx(0,lx,hy)] +
                                 fx*(1.0-fy)*fluxes_[fidx(0,hx,ly)] + fx*fy*fluxes_[fidx(0,hx,hy)];
    }
  }
  for (int x = 0; x < w_; ++x) {
    for (int y = 1; y < h_; ++y) {
      const Eigen::Vector2f mid_source = Eigen::Vector2f(x + 0.5f, y) - dt*0.5*InterpolateVelocity(Eigen::Vector2f(x + 0.5f, y));
      const Eigen::Vector2f source = Eigen::Vector2f(x + 0.5f, y) - dt*InterpolateVelocity(ClipPoint(mid_source));
      int sx = static_cast<int>(source[0] - 0.5f);
      int sy = static_cast<int>(source[1]);
      float fx = source[0]  - 0.5f - sx;
      float fy = source[1] - sy;
      int lx = std::max(1, sx);
      int ly = std::max(1, sy);
      int hx = std::min(w_-1, sx+1);
      int hy = std::min(h_-1, sy+1);
      new_fluxes[fidx(1,x,y)] = (1.0f-fx)*(1.0-fy)*fluxes_[fidx(1,lx,ly)] + (1.0-fx)*fy*fluxes_[fidx(1,lx,hy)] +
                                 fx*(1.0-fy)*fluxes_[fidx(1,hx,ly)] + fx*fy*fluxes_[fidx(1,hx,hy)];
    }
  }
  fluxes_.swap(new_fluxes);
}

void Fluid::Project() {
  std::vector<float> pressure(w_ * h_);
  std::vector<float> div(w_ * h_);
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
  
  for (int k = 0; k < 20; ++k) {
    float err = 0.0f;
    for (int x = 0; x < w_; ++x) {
      for (int y = 0; y < h_; ++y) {
        int here = vidx(x,y);
        float diff = pressure[here];
        float sigma = 0.0f;
        float count = 0.0f;
        if (x > 0) {
          sigma -= pressure[vidx(x-1,y)];
          count += 1.0f;
        }
        if (x < w_-1) {
          sigma -= pressure[vidx(x+1,y)];
          count += 1.0f;
        }
        if (y > 0) {
          sigma -= pressure[vidx(x,y-1)];
          count += 1.0f;
        }
        if (y < h_-1) {
          sigma -= pressure[vidx(x,y+1)];
          count += 1.0f;
        }
        pressure[here] = (1.0f - omega) * pressure[here] + omega / count * (div[here] - sigma);
        diff -= pressure[here];
        err += diff*diff;
      }
    }
    //std::cerr << "err: " << std::sqrt(err) << std::endl;
  }
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      if (x > 0) fluxes_[fidx(0,x,y)] -= pressure[vidx(x,y)] - pressure[vidx(x-1,y)];
      if (y > 0) fluxes_[fidx(1,x,y)] -= pressure[vidx(x,y)] - pressure[vidx(x,y-1)];
    }
  }
  
  //FixBoundaries();
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

void Fluid::AddImpulse(float x0, float y0, float dx, float dy) {
  pending_impulse_origins_.push_back(Eigen::Vector2f(x0, y0));
  pending_impulse_deltas_.push_back(Eigen::Vector2f(dx*50.0, dy*50.0));
}

void Fluid::ApplyImpulses() {
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
}

void Fluid::AdvectDensity(float dt) {
  
}