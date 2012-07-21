//
//  ViewController.h
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import <UIKit/UIKit.h>
#import <GLKit/GLKit.h>
#import <AVFoundation/AVFoundation.h>
#import <CoreGraphics/CoreGraphics.h>
#import <CoreVideo/CoreVideo.h>
#import <CoreMedia/CoreMedia.h>



@interface ViewController : GLKViewController <AVCaptureVideoDataOutputSampleBufferDelegate> {
    cv::VideoCapture *_videoCapture;
    cv::Mat _lastFrame;
    
    AVCaptureSession *_captureSession;
    
}

@property (nonatomic, retain) AVCaptureSession *captureSession;

- (void)initCapture;


@end
