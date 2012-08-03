//
//  BlobTracker.cpp
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/24/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "BlobTracker.h"

BlobTracker::BlobTracker(int width, int height, int bpp, int bpr)
: detect_threshold_(10),
width_(width),
height_(height),
bpp_(bpp),
bpr_(bpr),
frames_seen_(0) {
  background_.clear();
  background_.resize(height * bpr);
}

void BlobTracker::FindBlobs(const uint8_t* buffer,
                            bool* green_found, int* green_x, int* green_y, bool* red_found, int* red_x, int* red_y) {
  int greens_found = 0;
  int greens_x = 0;
  int greens_y = 0;
  
  int reds_found = 0;
  int reds_x = 0;
  int reds_y = 0;
  
  *green_found = false;
  *red_found = false;
  
  if (frames_seen_ > 5) {
    for(int x = 0;x<width_;++x){
      for(int y = 0;y<height_;++y) {
        int rVal = buffer[y*bpr_ + x*bpp_ + 2] - background_[y*bpr_ + x*bpp_ + 2];
        int gVal = buffer[y*bpr_ + x*bpp_ + 1] - background_[y*bpr_ + x*bpp_ + 1];
        int bVal = buffer[y*bpr_ + x*bpp_ + 0] - background_[y*bpr_ + x*bpp_ + 0];
        
        for (int c = 0; c < 3; ++c) {
          background_[y*bpr_ + x*bpp_ + c] -= background_[y*bpr_ + x*bpp_ + c] >> 3;
          background_[y*bpr_ + x*bpp_ + c] += buffer[y*bpr_ + x*bpp_ + c] >> 3;
        }
        
        if(gVal>50 && rVal<45 && bVal<45){
          greens_found++;
          greens_x = greens_x + x;
          greens_y = greens_y + y;
        }
        
        if(gVal<50 && bVal<50 && rVal>100){
          reds_found++;
          reds_x = reds_x + x;
          reds_y = reds_y + y;
        }
      }                    
    }    
    
    const int detect_threshold = 0;
    
    if(greens_found>detect_threshold){
      *green_found = true;
      *green_x = greens_x/greens_found;
      *green_y = greens_y/greens_found;
    }
    
    if(reds_found>detect_threshold){
      *red_found = true;
      *red_x = reds_x/reds_found;
      *red_y = reds_y/reds_found;            
    }
  }
  
  ++frames_seen_;
}