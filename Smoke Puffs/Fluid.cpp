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

#include <Accelerate/Accelerate.h>

#import "TargetConditionals.h"

void Fluid::InterpolateXVelocities(const std::vector<float>& xs,
                            const std::vector<float>& ys,
                            float* vx) {
  for (int i = 0; i < xs.size(); ++i) {
    float x = xs[i];
    float y = ys[i] - 0.5f;
    int sx = static_cast<int>(x);
    int sy = static_cast<int>(y);
    float fx = x - sx;
    float fy = y - sy;
    float flx = (1.0f - fx);
    float fhx = fx;
    float fly = (1.0f - fy);
    float fhy = fy;
    float result = 0.0f;
    float count = 1e-20f;
    float weight = 0.0f;
    int here = fidx(0,sx,sy);
    weight = flx*fly;//*fluid_face_[here];
    result += weight*fluxes_[here];
    count += weight;
    weight = flx*fhy;//*fluid_face_[here+1];
    result += weight*fluxes_[here+1];
    count += weight;
    weight = fhx*fly;//*fluid_face_[here+h_];
    result += weight*fluxes_[here+h_];
    count += weight;
    weight = flx*fly;//*fluid_face_[here+h_+1];
    result += weight*fluxes_[here+h_+1];
    count += weight;
    vx[i] = result / count;
  }
}

//class FluxInterpolator {
//public:
//  FluxInterpolator(int block_size) : block_size_(block_size) {
//    x_.resize(block_size);
//    y_.resize(block_size);
//    xf_.resize(block_size);
//    yf_.resize(block_size);
//  }
//  
//  void Interpolate(float* in_flux,
//                   int in_stride,
//                   float* x_coords,
//                   int x_stride,
//                   float* y_coords,
//                   int y_stride,
//                   float* out_flux,
//                   int out_stride,
//                   int length) {
//    vDSP_vfixu16(x_coords, 1, &x_[0], 1, length);
//    vDSP_vfixu16(y_coords, 1, &y_[0], 1, length);
//    vDSP_vfrac(x_coords, 1, &xf_[0], 1, length);
//    vDSP_vfrac(y_coords, 1, &yf_[0], 1, length);
//  }
//  
//private:
//  int block_size_;
//  std::vector<float> base_;
//  std::vector<float> xf_;
//  std::vector<float> yf_;
//};


#if !TARGET_IPHONE_SIMULATOR
extern "C" void arm7_jacobi_iteration(float* pressure,
                                      float* div,
                                      float* inv_count,
                                      int w,
                                      int h,
                                      float* new_pressure);

void JacobiARM(float* pressure,
               float* div,
               float* inv_count,
               int w,
               int h,
               float* new_pressure) {
  float* prev_pressure_start = pressure;
  float* pressure_start = pressure + h;
  float* next_pressure_start = pressure + 2*h;
  float* new_pressure_start = new_pressure + h;
  float* div_start = div + h;
  float* inv_count_start = inv_count + h;
  int length = (h * (w-2))/2;
  
  __asm__
  (
   ".align 4                         \n"
   "L3_%=:                           \n\t"
   "vldr.f64   d2, [%[p], #-4]       \n\t"
   "vldr.f64   d3, [%[p], #+4]       \n\t"
   "vldr.f64   d0, [%[pp]]           \n\t"
   "vldr.f64   d1, [%[np]]           \n\t"
   "vldr.f64   d4, [%[d]]            \n\t"
   "vldr.f64   d5, [%[ic]]           \n\t"
   "subs.w     %[l], %[l], #1        \n\t"
   "vadd.f32   d7, d2, d3            \n\t"
   "vadd.f32   d6, d0, d1            \n\t"
   "vadd.f32   d6, d6, d4            \n\t"
   "vadd.f32   d7, d7, d6            \n\t"
   "vmul.f32   d7, d7, d5            \n\t"
   "vstr.f64   d7, [%[n]]            \n\t"
   "add        %[p], %[p], #8        \n\t"
   "add        %[pp], %[pp], #8      \n\t"
   "add        %[np], %[np], #8      \n\t"
   "add        %[ic], %[ic], #8      \n\t"
   "add        %[n], %[n], #8        \n\t"
   "add        %[d], %[d], #8        \n\t"
   "bne        L3_%=                 \n\t"
   
   : [p] "+r" (pressure_start), [pp] "+r" (prev_pressure_start), [np] "+r" (next_pressure_start), [d] "+r" (div_start), [n] "+r" (new_pressure_start), [l] "+r" (length), [ic] "+r" (inv_count_start)
   :
   : "memory", "cc", "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7"
   );
  
  pressure[0] = 0.0f;
}

#endif  // !TARGET_IPHONE_SIMULATOR

Timer::Timer(const std::string& label) : label_(label) {
  total_time_ = 0.0f;
  total_events_ = 0;
}
  
void Timer::BeginEvent() {
  gettimeofday(&this_event_start_, NULL);
}

void Timer::EndEvent() {
  struct timeval end_time;
  gettimeofday(&end_time, NULL);
  
  float diff = 1e3 * (end_time.tv_sec - this_event_start_.tv_sec);
  diff += 1e-3 * (end_time.tv_usec - this_event_start_.tv_usec);
  total_time_ += diff;
  total_events_ += 1;
}

void Timer::Print() {
  std::cerr << label_ << ": " << (total_time_ / static_cast<float>(total_events_)) << " ms avg. over " << total_events_ << " events" << std::endl;
}

Fluid::Fluid(int width, int height) : impulse_timer_("Impulse"), advection_timer_("Advection"), projection_timer_("Projection"), density_timer_("Density") {
  w_ = width + 2;
  h_ = height + 2;
  
  fluxes_.resize(w_*h_*2);
  densities_.resize(w_*h_);
  fluid_.resize(w_*h_);
  fluid_face_.resize(w_*h_*2);
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      fluxes_[fidx(0,x,y)] = 0.0f;
      fluxes_[fidx(1,x,y)] = 0.0f;
      densities_[vidx(x,y)] = 0.0f;
      
      if (x == 0 || x == w_ - 1 || y == 0 || y == h_ - 1) {
        fluid_[vidx(x,y)] = 0.0f;
      } else {
        fluid_[vidx(x,y)] = 1.0f;
      }
      
      if (x == 0) {
        fluid_face_[fidx(0,x,y)] = 0.0f;
      }
      if (y == 0) {
        fluid_face_[fidx(1,x,y)] = 0.0f;
      }
      if (x > 0 && y > 0) {
        fluid_face_[fidx(0,x,y)] = 1.0f;
        fluid_face_[fidx(1,x,y)] = 1.0f;
      }
    }
  }
  
  for (int x = 1; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      if (fluid_[vidx(x,y)] == 0.0f || fluid_[vidx(x-1,y)] == 0.0f) {
        fluid_face_[fidx(0,x,y)] = 0.0f;
      }
    }
  }
  for (int x = 0; x < w_; ++x) {
    for (int y = 1; y < h_; ++y) {
      if (fluid_[vidx(x,y)] == 0.0f || fluid_[vidx(x,y-1)] == 0.0f) {
        fluid_face_[fidx(1,x,y)] = 0.0f;
      }
    }
  }
  
  for (int i = 0; i < 10; ++i) {
    //Project();
  }
  
  smoke_radius_ = 8;
}

void Fluid::Advect(float dt) {
  std::vector<float> new_fluxes(w_ * h_ * 2);
  
  std::vector<float> source_points_x(h_);
  std::vector<float> source_points_y(h_);
  std::vector<float> source_points_frac_x(h_);
  std::vector<float> source_points_frac_y(h_);
  std::vector<float> interpolated_x(h_);
  std::vector<float> interpolated_y(h_);
  
  float half = 0.5f;
  float one = 1.0f;
  float mdt = -dt;
  float hmdt = -dt * 0.5f;
  float zero = 0.0f;
  float max_x = w_ - 1e-6;
  float max_y = h_ - 1e-6;
  float quarter = 0.25f;
  
//  for (int x = 1; x < w_-1; ++x) {
//    vDSP_vramp(&half, &one, &source_points_y[0], 1, h_);
//    float x_coord = x;
//    vDSP_vfill(&x_coord, &source_points_x[0], 1, h_);
//    vDSP_vadd(&fluxes_[0] + w_*h_ + x*h_, 1, &fluxes_[0] + w_*h_ + x*h_ + 1, 1, &interpolated_y[0], 1, h_);
//    vDSP_vadd(&fluxes_[0] + w_*h_ + x*h_ - h_, 1, &interpolated_y[0], 1, &interpolated_y[0], 1, h_);
//    vDSP_vasm(&fluxes_[0] + w_*h_ + x*h_ - h_ + 1, 1, &interpolated_y[0], 1, &quarter, &interpolated_y[0], 1, h_);
//    vDSP_vsma(&fluxes_[0] + x*h_, 1, &hmdt, &source_points_x[0], 1, &source_points_x[0], 1, h_);
//    vDSP_vsma(&interpolated_y[0], 1, &hmdt, &source_points_y[0], 1, &source_points_y[0], 1, h_);
//
//    vDSP_vclip(&source_points_x[0], 1, &zero, &max_x, &source_points_x[0], 1, h_);
//    vDSP_vclip(&source_points_y[0], 1, &zero, &max_y, &source_points_y[0], 1, h_);
//    
//    vDSP_vfrac(&source_points_x[0], 1, &source_points_frac_x[0], 1, h_);
//    vDSP_vfrac(&source_points_y[0], 1, &source_points_frac_y[0], 1, h_);
//    
//    vDSP_vsub(&source_points_x[0], 1, &source_points_frac_x[0], 1, &source_points_x[0], h_);
//    vDSP_vsub(&source_points_y[0], 1, &source_points_frac_y[0], 1, &source_points_y[0], h_);
//  }
  

//  std::vector<float> x_points(w_*h_);
//  std::vector<float> y_points(w_*h_);
//  std::vector<float> x_mid_results(w_*h_);
//  std::vector<float> y_mid_results(w_*h_);
//  std::vector<float> x_results(w_*h_);
//  std::vector<float> y_results(w_*h_);
//  
//  
//  for (int x = 0;x < w_; ++x) {
//    for (int y = 0; y < h_; ++y) {
//      x_points[vidx(x,y)] = x;
//      y_points[vidx(x,y)] = y + 0.5f;
//    }
//  }
//  InterpolateVelocities(x_points, y_points, &x_mid_results[0], &y_mid_results[0]);
//  for (int x = 0;x < w_; ++x) {
//    for (int y = 0; y < h_; ++y) {
//      int here = vidx(x,y);
//      x_mid_results[here] = Clip(x_points[here] - dt*0.5f*x_mid_results[here], 1.0f, (w_ - 1.0f) - 1e-6f);
//      y_mid_results[here] = Clip(y_points[here] - dt*0.5f*y_mid_results[here], 1.0f, (h_ - 1.0f) - 1e-6f);
//    }
//  }
//  InterpolateVelocities(x_mid_results, y_mid_results, &x_results[0], &y_results[0]);
//  for (int x = 0;x < w_; ++x) {
//    for (int y = 0; y < h_; ++y) {
//      int here = vidx(x,y);
//      x_results[here] = Clip(x_points[here] - dt*x_results[here], 1.0f, (w_ - 1.0f) - 1e-6f);
//      y_results[here] = Clip(y_points[here] - dt*y_results[here], 1.0f, (h_ - 1.0f) - 1e-6f);
//    }
//  }
//  InterpolateXVelocities(x_results, y_results, &new_fluxes[0]);
  
  
  
  for (int x = 1;x < w_-1; ++x) {
    for (int y = 1; y < h_-1; ++y) {
      const Eigen::Vector2f mid_source = Eigen::Vector2f(x, y + 0.5f) - dt*0.5*InterpolateVelocity(Eigen::Vector2f(x, y + 0.5f));
      const Eigen::Vector2f source = Eigen::Vector2f(x, y + 0.5) - dt*InterpolateVelocity(ClipPoint(mid_source));
      new_fluxes[fidx(0,x,y)] = 0.999f*InterpolateXVelocity(ClipPoint(source));
    }
  }
  for (int x = 1; x < w_-1; ++x) {
    for (int y = 1; y < h_-1; ++y) {
      const Eigen::Vector2f mid_source = Eigen::Vector2f(x + 0.5f, y) - dt*0.5*InterpolateVelocity(Eigen::Vector2f(x + 0.5f, y));
      const Eigen::Vector2f source = Eigen::Vector2f(x + 0.5f, y) - dt*InterpolateVelocity(ClipPoint(mid_source));
      new_fluxes[fidx(1,x,y)] = 0.999f*InterpolateYVelocity(ClipPoint(source));
    }
  }
  fluxes_.swap(new_fluxes);
}

void Fluid::Project() {
  std::vector<float> pressure(w_ * h_);
  std::vector<float> div(w_ * h_);
  
  for (int x = 0; x < w_; ++x) {
    for (int y = 0; y < h_; ++y) {
      int xidx = fidx(0, x, y);
      int yidx = fidx(1, x, y);
      
      fluxes_[xidx] *= fluid_face_[xidx];
      fluxes_[yidx] *= fluid_face_[yidx];
    }
  }
  
  for (int x = 1; x < w_-1; ++x) {
    for (int y = 1; y < h_-1; ++y) {
      int here = vidx(x, y);
      pressure[here] = 0.0f;
      pressure[here] += fluxes_[fidx(0,x,y)];
      pressure[here] -= fluxes_[fidx(0,x+1,y)];
      pressure[here] += fluxes_[fidx(1,x,y)];
      pressure[here] -= fluxes_[fidx(1,x,y+1)];
      div[here] = pressure[here];
      
      // initialized with div => converges Gauss-Seidel, but not Jacobi
      // initialized with zero => converges Jacobi (and maybe GS)
      pressure[here] = 0.0f;
    }
  }
  
  std::vector<float> inv_count(w_*h_);
  vDSP_vadd(&inv_count[0], 1, (&fluid_[0])+h_, 1, &inv_count[0], 1, w_*h_-h_);
  vDSP_vadd((&inv_count[0])+h_, 1, &fluid_[0], 1, (&inv_count[0])+h_, 1, w_*h_-h_);
  vDSP_vadd(&inv_count[0], 1, (&fluid_[0])+1, 1, &inv_count[0], 1, w_*h_-1);
  vDSP_vadd((&inv_count[0])+1, 1, &fluid_[0], 1, (&inv_count[0])+1, 1, w_*h_-1);
  float one = 1.0f;
  vDSP_vthr(&inv_count[0], 1, &one, &inv_count[0], 1, w_*h_);
  vDSP_svdiv(&one, &inv_count[0], 1, &inv_count[0], 1, w_*h_);
  vDSP_vmul(&inv_count[0], 1, &fluid_[0], 1, &inv_count[0], 1, w_*h_);
  
  // needed for SOR
//  float omega = 1.9f;
//  float one_minus_omega = 1.0f - omega;
//  vDSP_vsmul(&inv_count[0], 1, &omega, &inv_count[0], 1, w_*h_);
  
  //std::vector<float> sigma(w_*h_);
  std::vector<float> new_pressure(w_*h_);
  const int MAX_ITERS = 50;
  for (int k = 0; k < MAX_ITERS; ++k) {
    // Jacobi
#if !TARGET_IPHONE_SIMULATOR
    arm7_jacobi_iteration(&pressure[0],
                          &div[0],
                          &inv_count[0],
                          w_,
                          h_,
                          &new_pressure[0]);
#else
    for (int x = 1; x < w_-1; ++x) {
      float* pressure_start = &pressure[0] + x*h_;
      float* new_pressure_start = &new_pressure[0] + x*h_;
      
      vDSP_vadd(pressure_start + h_, 1, pressure_start + 1, 1, new_pressure_start, 1, h_);
      vDSP_vadd(pressure_start - h_, 1, new_pressure_start, 1, new_pressure_start, 1, h_);
      vDSP_vadd(pressure_start - 1, 1, new_pressure_start, 1, new_pressure_start, 1, h_);
      vDSP_vam(new_pressure_start, 1, &div[0] + x*h_, 1, &inv_count[0] + x*h_, 1, new_pressure_start, 1, h_);
    }
#endif  // !TARGET_IPHONE_SIMULATOR
    


    
    // Gauss-Seidel
//        if (k == MAX_ITERS - 1) {
//          std::copy(pressure.begin(), pressure.end(), new_pressure.begin());
//        }
//    for (int offset = 0; offset < 2; ++offset) {
//      for (int x = 1; x < w_-1; ++x) {
//        float* pressure_start = &pressure[0] + x*h_ + (offset+x)%2;
//        
//        vDSP_vadd(pressure_start + h_, 2, pressure_start + 1, 2, pressure_start, 2, h_/2);
//        vDSP_vadd(pressure_start - h_, 2, pressure_start, 2, pressure_start, 2, h_/2);
//        vDSP_vadd(pressure_start - 1, 2, pressure_start, 2, pressure_start, 2, h_/2);
//        vDSP_vam(pressure_start, 2, &div[0] + x*h_ + (offset+x)%2, 2, &inv_count[0] + x*h_ + (offset+x)%2, 2, pressure_start, 2, h_/2);
//      }
//    }
//    pressure.swap(new_pressure);

    
    // SOR
//    if (k == MAX_ITERS - 1) {
//      std::copy(pressure.begin(), pressure.end(), new_pressure.begin());
//    }
//    
//    int even_offset = 0;
//    int odd_offset = 1;
//    for (int x = 1; x < w_-1; ++x) {
//      vDSP_vadd(&pressure[0] + x*h_ + h_ + even_offset, 2, &pressure[0] + x*h_ + 1 + even_offset, 2, &sigma[0] + x*h_ + even_offset, 2, h_/2);
//      vDSP_vadd(&pressure[0] + x*h_ - h_ + even_offset, 2, &sigma[0] + x*h_ + even_offset, 2, &sigma[0] + x*h_ + even_offset, 2, h_/2);
//      vDSP_vadd(&pressure[0] + x*h_ - 1 + even_offset, 2, &sigma[0] + x*h_ + even_offset, 2, &sigma[0] + x*h_ + even_offset, 2, h_/2);
//      vDSP_vsmul(&pressure[0] + x*h_ + even_offset, 2, &one_minus_omega, &pressure[0] + x*h_ + even_offset, 2, h_/2);
//      vDSP_vam(&sigma[0] + x*h_ + even_offset, 2, &div[0] + x*h_ + even_offset, 2, &inv_count[0] + x*h_ + even_offset, 2, &sigma[0] + x*h_ + even_offset, 2, h_/2);
//      vDSP_vadd(&sigma[0] + x*h_ + even_offset, 2, &pressure[0] + x*h_ + even_offset, 2, &pressure[0] + x*h_ + even_offset, 2, h_/2);
//    }
//    for (int x = 1; x < w_-1; ++x) {
//      vDSP_vadd(&pressure[0] + x*h_ + h_ + odd_offset, 2, &pressure[0] + x*h_ + 1 + odd_offset, 2, &sigma[0] + x*h_ + odd_offset, 2, h_/2);
//      vDSP_vadd(&pressure[0] + x*h_ - h_ + odd_offset, 2, &sigma[0] + x*h_ + odd_offset, 2, &sigma[0] + x*h_ + odd_offset, 2, h_/2);
//      vDSP_vadd(&pressure[0] + x*h_ - 1 + odd_offset, 2, &sigma[0] + x*h_ + odd_offset, 2, &sigma[0] + x*h_ + odd_offset, 2, h_/2);
//      vDSP_vsmul(&pressure[0] + x*h_ + odd_offset, 2, &one_minus_omega, &pressure[0] + x*h_ + odd_offset, 2, h_/2);
//      vDSP_vam(&sigma[0] + x*h_ + odd_offset, 2, &div[0] + x*h_ + odd_offset, 2, &inv_count[0] + x*h_ + odd_offset, 2, &sigma[0] + x*h_ + odd_offset, 2, h_/2);
//      vDSP_vadd(&sigma[0] + x*h_ + odd_offset, 2, &pressure[0] + x*h_ + odd_offset, 2, &pressure[0] + x*h_ + odd_offset, 2, h_/2);
//    }
    


//    if (k == MAX_ITERS - 1) {
//      float norm = 0.0f;
//      float err = 0.0f;
//      for (int x = 1; x < w_-1; ++x) {
//        for (int y = 1; y < h_-1; ++y) {
//          int here = vidx(x, y);
//          norm += new_pressure[here]*new_pressure[here];
//          float diff = new_pressure[here] - pressure[here];
//          err += diff*diff;
//        }
//      }
//      std::cerr << "FINAL REL ERR: " << (err / norm) << std::endl;
//    }
    
    pressure.swap(new_pressure);
  }
  
  for (int x = 1; x < w_-1; ++x) {
    for (int y = 1; y < h_-1; ++y) {
      fluxes_[fidx(0,x,y)] -= fluid_[vidx(x,y)]*fluid_[vidx(x-1,y)]*(pressure[vidx(x,y)] - pressure[vidx(x-1,y)]);
      fluxes_[fidx(1,x,y)] -= fluid_[vidx(x,y)]*fluid_[vidx(x,y-1)]*(pressure[vidx(x,y)] - pressure[vidx(x,y-1)]);
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
  line_coords->reserve(w_*h_*8);
  for (int x = 1; x < w_-1; ++x) {
    for (int y = 1; y < h_-1; ++y) {
      if (x > 0 && y > 0) {
        line_coords->push_back(x - 1.0f + 0.5);
        line_coords->push_back(y - 1.0f + 0.5);
        line_coords->push_back(x - 1.0f + 0.5 + scale*CellCenterVelocity(x, y)[0]);
        line_coords->push_back(y - 1.0f + 0.5 + scale*CellCenterVelocity(x, y)[1]);
      } else {
        line_coords->push_back(0.0f);
        line_coords->push_back(0.0f);
        line_coords->push_back(0.0f);
        line_coords->push_back(0.0f);
      }
//      line_coords->push_back(x - 1.0f );
//      line_coords->push_back(y - 1.0f + 0.5f);
//      line_coords->push_back(x - 1.0f + fluxes_[fidx(0,x,y)]);
//      line_coords->push_back(y - 1.0f + 0.5f);
//      line_coords->push_back(x - 1.0f + 0.5f);
//      line_coords->push_back(y - 1.0f);
//      line_coords->push_back(x - 1.0f + 0.5f);
//      line_coords->push_back(y - 1.0f + fluxes_[fidx(1,x,y)]);
    }
  }
}

void Fluid::GetDensities(std::vector<float>* densities) {
  densities->clear();
  densities->reserve(densities_.size());
  for (int y = 1; y < h_-1; ++y) {
    for (int x= 1; x < w_-1; ++x) {
      densities->push_back(densities_[vidx(x,y)]);
    }
  }
}

void Fluid::AddImpulse(float x0, float y0, float vx, float vy) {
  pending_impulse_origins_.push_back(Eigen::Vector2f(x0-1.0f, y0-1.0f));
  pending_impulse_velocities_.push_back(Eigen::Vector2f(vx*100.0f, vy*100.0f));
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
    
    if (x > 2 && x < w_-3 && y > 2 && y < h_-3) {
      fluxes_[fidx(0,x,y)] += delta[0];
      fluxes_[fidx(0,x+1,y)] += delta[0];
      fluxes_[fidx(1,x,y)] += delta[1];
      fluxes_[fidx(1,x,y+1)] += delta[1];
    }
    
    for (int k = -smoke_radius_; k <= smoke_radius_; ++k) {
      int xk = x + k;
      if (xk < 2 || xk >= w_ - 2) continue;
      for (int j = -smoke_radius_; j <= smoke_radius_; ++j) {
        int yj = y + j;
        if (k*k + j*j > smoke_radius_*smoke_radius_) continue;
        if (yj < 2 || yj >= h_ - 2) continue;
        densities_[vidx(xk,yj)] = std::max(densities_[vidx(xk, yj)], std::max(1.0f, 1.0f / (1 + sqrtf(k*k + j*j))));
      }
    }
  }
  
  pending_impulse_origins_.clear();
  pending_impulse_velocities_.clear();
}

void Fluid::AdvectDensity(float dt) {
  std::vector<float> new_densities(densities_.size());
  for (int x = 1; x < w_-1; ++x) {
    for (int y = 1; y < h_-1; ++y) {
      const Eigen::Vector2f mid_source = Eigen::Vector2f(x + 0.5f, y + 0.5f) - dt*0.5*InterpolateVelocity(Eigen::Vector2f(x + 0.5f, y + 0.5f));
      const Eigen::Vector2f source = Eigen::Vector2f(x + 0.5f, y + 0.5f) - dt*InterpolateVelocity(ClipPoint(mid_source));
      new_densities[vidx(x,y)] = 0.993f * InterpolateDensity(ClipPoint(source));
    }
  }
  densities_.swap(new_densities);
}
