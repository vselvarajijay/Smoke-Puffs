//
//  BlobTracker.h
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/24/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Smoke_Puffs_BlobTracker_h
#define Smoke_Puffs_BlobTracker_h

#include <stdint.h>
#include <cstddef>
#include <vector>

#include <Eigen/Dense>

class BlobTracker {
public:
  BlobTracker(int width, int height, int bpp, int bpr);
  
  // Assumes BGR byte ordering.
  void FindBlobs(const uint8_t* buffer,
                 bool* green_found, int* green_x, int* green_y, bool* red_found, int* red_x, int* red_y);
  
private:
  int detect_threshold_;
  int width_;
  int height_;
  int bpp_;
  int bpr_;
  std::vector<short> background_;
  std::vector<uint8_t> last_frame_;
  std::vector<uint8_t> detect_image_;
  int frames_seen_;
};

#endif
