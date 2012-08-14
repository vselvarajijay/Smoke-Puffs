//
//  ViewController.m
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "ViewController.h"

#import <CoreGraphics/CoreGraphics.h>
#import <CoreVideo/CoreVideo.h>
#import <CoreMedia/CoreMedia.h>

#import "Fluid.h"

#define RENDER_LINES 0
#define SHOW_PUPPETS 0

// Uniform index.
enum
{
  UNIFORM_MODELVIEWPROJECTION_MATRIX,
  UNIFORM_SAMPLER,
  UNIFORM_XRES,
  UNIFORM_YRES,
  NUM_UNIFORMS
};
GLint uniforms[NUM_UNIFORMS];

@interface ViewController () {
  AVCaptureSession *_captureSession;
  
  GLuint _program;
  
  GLKMatrix4 _modelViewProjectionMatrix;
  float _rotation;
  
  GLuint _vertexArray;
  GLuint _vertexBuffer;
  
  int _width;
  int _height;
  
  Fluid* _fluid;
  
  float _dt;
  
  float _ball_x;
  float _ball_y;
  
  GLint sampler;
  GLuint density_texture;
  
  UIView* red_ball;
  UIView* green_ball;
  UIView* goal_left;
  UIView* goal_right;
  UIImageView* soccer_ball;
  UIImageView *androidView;
  UIImageView *cmView;
  int ball_size;
  int soccer_ball_size;
  int goal_height;
  int goal_width;
  int points_left;
  int points_right;
  UILabel* score_label;
  int ball_wait_frames;
  
  float _view_width;
  float _view_height;
  
  BOOL soccer_on;
  
  int green_x, green_y;
  int previous_green_x, previous_green_y;
  int red_x, red_y;
  int previous_red_x, previous_red_y;
}
@property (nonatomic, retain) AVCaptureSession *captureSession;

@property (strong, nonatomic) EAGLContext *context;
@property (strong, nonatomic) GLKBaseEffect *effect;

@property (nonatomic, assign) int width;
@property (nonatomic, assign) int height;

@property (nonatomic, assign) float view_width;
@property (nonatomic, assign) float view_height;

@property (nonatomic, assign) Fluid* fluid;
@property (nonatomic, assign) float dt;

@property (nonatomic, assign) float ball_x;
@property (nonatomic, assign) float ball_y;


- (void)setupGL;
- (void)tearDownGL;

- (void)initCapture;
- (BOOL)loadShaders;
- (BOOL)compileShader:(GLuint *)shader type:(GLenum)type file:(NSString *)file;
- (BOOL)linkProgram:(GLuint)prog;
- (BOOL)validateProgram:(GLuint)prog;
- (CGPoint)screenCoordsFromFluidCoords:(CGPoint)pt;
- (CGPoint)fluidCoordsFromScreenCoords:(CGPoint)pt;

-(void)setupSoccer;
-(void)tearDownSoccer;
-(void)pointScored;
@end

@implementation ViewController

@synthesize captureSession = _captureSession;

@synthesize context = _context;
@synthesize effect = _effect;
@synthesize view_width = _view_width;
@synthesize view_height = _view_height;
@synthesize width = _width;
@synthesize height = _height;
@synthesize fluid = _fluid;
@synthesize dt = _dt;
@synthesize ball_x = _ball_x;
@synthesize ball_y = _ball_y;

- (CGPoint)screenCoordsFromFluidCoords:(CGPoint)pt
{
  return CGPointMake(self.view_width / self.width * pt.x,
                     self.view_height - self.view_height / self.height * pt.y);
}

- (CGPoint)fluidCoordsFromScreenCoords:(CGPoint)pt
{
  return CGPointMake(self.width / self.view_width * pt.x,
                     (self.view_height - pt.y) * self.height / self.view_height);
}

- (void)pointScored
{
  score_label.text = [NSString stringWithFormat: @"%d : %d", points_left, points_right];
  self.ball_x = self.width / 2.0;
  self.ball_y = self.height / 2.0;
  ball_wait_frames = 30;
  [soccer_ball removeFromSuperview];
}

- (void)viewDidLoad
{
  [super viewDidLoad];
  
  self.context = [[EAGLContext alloc] initWithAPI:kEAGLRenderingAPIOpenGLES2];
  
  if (!self.context) {
    NSLog(@"Failed to create ES context");
  }
  
  GLKView *view = (GLKView *)self.view;
  view.context = self.context;
  view.drawableDepthFormat = GLKViewDrawableDepthFormat24;
  
  self.view_width = self.view.frame.size.height;
  self.view_height = self.view.frame.size.width;
  self.width = 160;//MIN(128, self.view_width/6);
  self.height = 120;//MIN(96, self.view_height/6);
  self.fluid = new Fluid(self.width, self.height);
  self.fluid->set_smoke_radius(self.width / 20.0f);
  self.dt = 0.25f;
  self.ball_x = self.width * 0.5f;
  self.ball_y = self.height * 0.5f;
  
  [self initCapture];
  
  [self setupGL];
  [self.view setMultipleTouchEnabled:YES];
  
  ball_size = self.view_width / 34.0f;
  soccer_ball_size = self.view_width / 12.8f;
  float character_size = self.view_width / 12.8f;
  
  UIImage *cm = [UIImage imageNamed: @"cm.jpeg"];  
  cmView = [[UIImageView alloc] initWithImage:cm];
  cmView.frame =   CGRectMake(0, 0, character_size, character_size);
  
  red_ball = [[UIView alloc] init];
  [red_ball setBackgroundColor:[UIColor redColor]];
  red_ball.frame = CGRectMake(0, 0, ball_size, ball_size);
  red_ball.layer.cornerRadius = ball_size/2;
  
#if SHOW_PUPPETS
  [self.view addSubview:cmView];
#endif
  
  UIImage *android = [UIImage imageNamed: @"android.png"];
  androidView = [[UIImageView alloc] initWithImage:android];
  androidView.frame = CGRectMake(0, 0, character_size, character_size);
  
  green_ball = [[UIView alloc] init];
  [green_ball setBackgroundColor:[UIColor greenColor]];
  green_ball.frame = CGRectMake(0, 0, ball_size, ball_size);
  green_ball.layer.cornerRadius = ball_size/2;
  
#if SHOW_PUPPETS
  [self.view addSubview:androidView];
#endif
  
  soccer_ball = [[UIImageView alloc] initWithImage: [UIImage imageWithContentsOfFile:[[NSBundle mainBundle] pathForResource:@"soccer_ball" ofType:@"png"]]];
  soccer_ball.frame = CGRectMake(0, 0, soccer_ball_size, soccer_ball_size);
  
  goal_width = self.view_width / 7.0f;
  goal_height = self.view_height / 2.0f;
  
  goal_left = [[UIView alloc] init];
  goal_left.frame = CGRectMake(0, (self.view_height-goal_height)/2, goal_width, goal_height);
  goal_left.layer.borderWidth = self.view_width / 102.0f;
  goal_left.layer.borderColor = [[UIColor blackColor] CGColor];
  
  
  goal_right = [[UIView alloc] init];
  goal_right.frame = CGRectMake(self.view_width - goal_width, (self.view_height-goal_height)/2, goal_width, goal_height);
  goal_right.layer.borderWidth = self.view_width / 102.0f;
  goal_right.layer.borderColor = [[UIColor blackColor] CGColor];
  
  score_label = [ [UILabel alloc ] initWithFrame:CGRectMake(self.view_width * 0.33f, 0.0f, self.view_width*0.4f, self.view_height/7.6f) ];
  score_label.textAlignment =  UITextAlignmentCenter;
  score_label.textColor = [UIColor orangeColor];
  score_label.backgroundColor = [UIColor clearColor];
  score_label.font = [UIFont fontWithName:@"Arial Rounded MT Bold" size:(self.view_width/13.2f)];
  
  ball_wait_frames = -1;
  
  green_x = self.view_width / 2.0f;
  green_y = self.view_height / 2.0f;
  red_x = self.view_width / 2.0f;
  red_y = self.view_height / 2.0f;
  previous_green_x = self.view_width / 2.0f;
  previous_green_y = self.view_height / 2.0f;
  previous_red_x = self.view_width / 2.0f;
  previous_red_y = self.view_height / 2.0f;
  
  soccer_on = NO;
  
}


- (void)setupSoccer {
  [self.view addSubview:soccer_ball];
  [self.view addSubview:goal_left];
  [self.view addSubview:goal_right];
  
  points_left = 0;
  points_right = 0;
  ball_wait_frames = -1;
  self.ball_x = self.width / 2.0;
  self.ball_y = self.height / 2.0;
  soccer_ball.center = [self screenCoordsFromFluidCoords:CGPointMake(self.ball_x, self.ball_y)];
  [self.view addSubview:score_label];
  score_label.text = [NSString stringWithFormat: @"%d : %d", points_left, points_right];
  soccer_on = YES;
}

- (void)tearDownSoccer {
  [soccer_ball removeFromSuperview];
  [goal_left removeFromSuperview];
  [goal_right removeFromSuperview];
  [score_label removeFromSuperview];
  
  soccer_on = NO;
}

- (void)initCapture {
	/*We setup the input*/
  
  NSArray *devices = [AVCaptureDevice devices];
  AVCaptureDeviceInput *captureInput = nil;
  for (AVCaptureDevice *device in devices) {
    // Get the front or the back camera
    if (device.position == AVCaptureDevicePositionFront) {
      captureInput = [AVCaptureDeviceInput deviceInputWithDevice:device error:nil];
      break;
    }
  }
  
	/*We setup the output*/
	AVCaptureVideoDataOutput *captureOutput = [[AVCaptureVideoDataOutput alloc] init];
	/*While a frame is processed in -captureOutput:didOutputSampleBuffer:fromConnection: delegate methods no other frames are added in the queue.
	 If you don't want this behaviour set the property to NO */
	captureOutput.alwaysDiscardsLateVideoFrames = YES; 
	
	/*We create a serial queue to handle the processing of our frames*/
	dispatch_queue_t queue;
	queue = dispatch_queue_create("cameraQueue", NULL);
	[captureOutput setSampleBufferDelegate:self queue:queue];
	dispatch_release(queue);
	// Set the video output to store frame in BGRA (It is supposed to be faster)
	NSString* key = (NSString*)kCVPixelBufferPixelFormatTypeKey; 
	NSNumber* value = [NSNumber numberWithUnsignedInt:kCVPixelFormatType_32BGRA]; 
	NSDictionary* videoSettings = [NSDictionary dictionaryWithObject:value forKey:key]; 
	[captureOutput setVideoSettings:videoSettings]; 
	/*And we create a capture session*/
	self.captureSession = [[AVCaptureSession alloc] init];
	/*We add input and output*/
  if (captureInput != nil) {
    [self.captureSession addInput:captureInput];
  }
	[self.captureSession addOutput:captureOutput];
  /* Low quality capture, 192x144. That's enough for the detection we're doing. */
  [self.captureSession setSessionPreset:AVCaptureSessionPresetLow];	
	[self.captureSession startRunning];
	
}

#pragma mark -
#pragma mark AVCaptureSession delegate
- (void)captureOutput:(AVCaptureOutput *)captureOutput 
didOutputSampleBuffer:(CMSampleBufferRef)sampleBuffer 
       fromConnection:(AVCaptureConnection *)connection 
{ 	
  CVImageBufferRef imageBuffer = CMSampleBufferGetImageBuffer(sampleBuffer); 
  /*Lock the image buffer*/
  CVPixelBufferLockBaseAddress(imageBuffer,0); 
  /*Get information about the image*/
  uint8_t *baseAddress = (uint8_t *)CVPixelBufferGetBaseAddress(imageBuffer);
  size_t bytesPerRow = CVPixelBufferGetBytesPerRow(imageBuffer); 
  size_t width = CVPixelBufferGetWidth(imageBuffer); 
  size_t height = CVPixelBufferGetHeight(imageBuffer);
  [self findBlobsInBuffer:baseAddress withWidth:width Height:height bytesPerPixel:4 bytesPerRow:bytesPerRow];
  
  /*We unlock the  image buffer*/
  CVPixelBufferUnlockBaseAddress(imageBuffer,0);
}

// Expects BGR[A] pixel ordering
-(void)findBlobsInBuffer:(uint8_t*)buffer withWidth:(size_t)width Height:(size_t)height bytesPerPixel:(size_t)bpp bytesPerRow:(size_t)bpr { 
  int greens_found = 0;
  int greens_x = 0;
  int greens_y = 0;
  
  int reds_found = 0;
  int reds_x = 0;
  int reds_y = 0;
  for(int x = 0;x<width;x++){
    for(int y = 0;y<height;y++) {
      size_t base = y*bpr + x*bpp;
      uint8_t rVal = buffer[base + 2];
      uint8_t gVal = buffer[base + 1];
      uint8_t bVal = buffer[base + 0];
      
      if(gVal>100 && rVal<90 && bVal<90){
        greens_found++;
        greens_x = greens_x + x;
        greens_y = greens_y + y;
      }
      
      if(gVal<100 && bVal<100 && rVal>200){
        reds_found++;
        reds_x = reds_x + x;
        reds_y = reds_y + y;
      }
    }                    
  }    
  
  const int detect_threshold = 0;
  if(greens_found>detect_threshold){
    int x = greens_x/greens_found;
    int y = greens_y/greens_found;       
    green_x = (int)(x * (self.view_width/width));
    green_y = self.view_height - (int)(y * (self.view_height/height));               
  }
  
  if(reds_found>detect_threshold){
    int x = reds_x/reds_found;
    int y = reds_y/reds_found;        
    red_x = (int)(x * (self.view_width/width));
    red_y = self.view_height - (int)(y * (self.view_height/height));               
  }
}

- (void)viewDidUnload
{    
  [super viewDidUnload];
  
  [self tearDownGL];
  
  if ([EAGLContext currentContext] == self.context) {
    [EAGLContext setCurrentContext:nil];
  }
  delete self.fluid;
	self.context = nil;
}

- (void)didReceiveMemoryWarning
{
  [super didReceiveMemoryWarning];
  // Release any cached data, images, etc. that aren't in use.
}

- (BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)interfaceOrientation
{
  return interfaceOrientation == UIInterfaceOrientationLandscapeRight;
}

- (void)setupGL
{
  [EAGLContext setCurrentContext:self.context];
  
  [self loadShaders];
  
  self.effect = [[GLKBaseEffect alloc] init];
}

- (void)tearDownGL
{
  [EAGLContext setCurrentContext:self.context];
  
  glDeleteBuffers(1, &_vertexBuffer);
  glDeleteVertexArraysOES(1, &_vertexArray);
  
  self.effect = nil;
  
  if (_program) {
    glDeleteProgram(_program);
    _program = 0;
  }
}

#pragma mark - GLKView and GLKViewController delegate methods

- (void)update
{
  self.fluid->Step(self.dt);
  
  if (soccer_on) {
    if (ball_wait_frames > 0) {
      --ball_wait_frames;
    } else {
      
      if (ball_wait_frames == 0) {
        --ball_wait_frames;
        [self.view addSubview:soccer_ball];
      }
      
      self.fluid->AdvectPoint(self.dt, self.ball_x, self.ball_y, &_ball_x, &_ball_y);
      float ball_r = soccer_ball_size / 2.0f * self.width / self.view_width;
      self.ball_x = std::max(ball_r, self.ball_x);
      self.ball_x = std::min(self.width - ball_r, self.ball_x);
      self.ball_y = std::max(ball_r, self.ball_y);
      self.ball_y = std::min(self.height - ball_r, self.ball_y);
      
      float ball_screen_r = soccer_ball_size / 2.0f - 4.0f;
      CGPoint ball_screen = [self screenCoordsFromFluidCoords:CGPointMake(self.ball_x, self.ball_y)];
      if (ball_screen.y < ((self.view_height + goal_height)/2 + ball_screen_r) && ball_screen.y > ((self.view_height - goal_height)/2 - ball_screen_r)) {
        if (ball_screen.x < (ball_screen_r + goal_width)) {
          ++points_right;
          [self pointScored];
        }
        if (ball_screen.x > (self.view_width - ball_screen_r - goal_width)) {
          ++points_left;
          [self pointScored];
        }
      }
      soccer_ball.center = ball_screen;
    }
  }
  
  static const float motion_smoothness = 0.75f;
  
  green_x = motion_smoothness * previous_green_x + (1.0f - motion_smoothness) * green_x;
  green_y = motion_smoothness * previous_green_y + (1.0f - motion_smoothness) * green_y;
  
  red_x = motion_smoothness * previous_red_x + (1.0f - motion_smoothness) * red_x;
  red_y = motion_smoothness * previous_red_y + (1.0f - motion_smoothness) * red_y;
  
  red_ball.center = CGPointMake(red_x , red_y);
  cmView.center = CGPointMake(red_x , red_y);
  
  green_ball.center = CGPointMake(green_x , green_y);
  androidView.center = CGPointMake(green_x , green_y);
  
  CGPoint green_pt = [self fluidCoordsFromScreenCoords:CGPointMake(green_x, green_y)];
  CGPoint previous_green_pt = [self fluidCoordsFromScreenCoords:CGPointMake(previous_green_x, previous_green_y)];
  
#if SHOW_PUPPETS
  self.fluid->AddImpulse( green_pt.x,
                         green_pt.y,
                         green_pt.x - previous_green_pt.x,
                         green_pt.y - previous_green_pt.y);
#endif
  
  
  CGPoint red_pt = [self fluidCoordsFromScreenCoords:CGPointMake(red_x, red_y)];
  CGPoint previous_red_pt = [self fluidCoordsFromScreenCoords:CGPointMake(previous_red_x, previous_red_y)];
  
#if SHOW_PUPPETS
  self.fluid->AddImpulse( red_pt.x,
                         red_pt.y,
                         red_pt.x - previous_red_pt.x,
                         red_pt.y - previous_red_pt.y);
#endif
  
  previous_green_x = green_x;
  previous_green_y = green_y;
  
  previous_red_x = red_x;
  previous_red_y = red_y;
}

- (void)glkView:(GLKView *)view drawInRect:(CGRect)rect
{
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glClearColor(0.3f, 0.3f, 0.6f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  
  self.effect.transform.projectionMatrix = GLKMatrix4MakeOrtho(0, self.width, 0, self.height, 1, -1);
  self.effect.transform.modelviewMatrix = GLKMatrix4Identity;
  
  // Use the program object
  glUseProgram ( _program );
  
  float drawrect[] = {
    0.0, 0.0, 0.0, 0.0,
    self.width, 0.0, 1.0, 0.0,
    0.0, self.height, 0.0, 1.0,
    self.width, self.height, 1.0, 1.0
  };
  GLushort indices[] = {0, 1, 2, 2, 3, 1};
  
  // Load the vertex position
  glVertexAttribPointer ( GLKVertexAttribPosition, 2, GL_FLOAT, 
                         GL_FALSE, 4 * sizeof(GLfloat), drawrect );
  // Load the texture coordinate
  glVertexAttribPointer ( GLKVertexAttribTexCoord0, 2, GL_FLOAT,
                         GL_FALSE, 4 * sizeof(GLfloat), &drawrect[2] );
  
  glEnableVertexAttribArray ( GLKVertexAttribPosition );
  glEnableVertexAttribArray (GLKVertexAttribTexCoord0 );
  
  // Bind the texture
  glActiveTexture ( GL_TEXTURE0 );
  glBindTexture ( GL_TEXTURE_2D, density_texture );
  std::vector<float> density;
  self.fluid->GetDensities(&density);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, self.width, self.height, GL_LUMINANCE, GL_FLOAT, &density[0]);
  
  // Set the sampler texture unit to 0
  glUniform1i ( uniforms[UNIFORM_SAMPLER], 0 );
  glUniform1f ( uniforms[UNIFORM_XRES], self.width );
  glUniform1f ( uniforms[UNIFORM_YRES], self.height );
  glUniformMatrix4fv(uniforms[UNIFORM_MODELVIEWPROJECTION_MATRIX], 1, GL_FALSE, self.effect.transform.projectionMatrix.m);
  
  glDrawElements ( GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, indices );
  glDisableVertexAttribArray(GLKVertexAttribPosition);
  glDisableVertexAttribArray(GLKVertexAttribTexCoord0);
  
  // Render velocity fields as lines with GLKit
#if RENDER_LINES
  [self.effect prepareToDraw];
  self.effect.constantColor = GLKVector4Make(1.0f, 0.0f, 0.0f, 1.0f);
  std::vector<float> lines;
  self.fluid->GetLines(&lines, 0.3);
  glEnableVertexAttribArray(GLKVertexAttribPosition);
  glVertexAttribPointer(GLKVertexAttribPosition, 2, GL_FLOAT, GL_FALSE, 0, &(lines[0]));
  glDrawArrays(GL_LINES, 0, lines.size()/2);
  glDisableVertexAttribArray(GLKVertexAttribPosition);
#endif
}

#pragma mark -  OpenGL ES 2 shader compilation

- (BOOL)loadShaders
{
  GLuint vertShader, fragShader;
  NSString *vertShaderPathname, *fragShaderPathname;
  
  // Create shader program.
  _program = glCreateProgram();
  
  // Create and compile vertex shader.
  vertShaderPathname = [[NSBundle mainBundle] pathForResource:@"Shader" ofType:@"vsh"];
  if (![self compileShader:&vertShader type:GL_VERTEX_SHADER file:vertShaderPathname]) {
    NSLog(@"Failed to compile vertex shader");
    return NO;
  }
  
  // Create and compile fragment shader.
  fragShaderPathname = [[NSBundle mainBundle] pathForResource:@"Shader" ofType:@"fsh"];
  if (![self compileShader:&fragShader type:GL_FRAGMENT_SHADER file:fragShaderPathname]) {
    NSLog(@"Failed to compile fragment shader");
    return NO;
  }
  
  // Attach vertex shader to program.
  glAttachShader(_program, vertShader);
  
  // Attach fragment shader to program.
  glAttachShader(_program, fragShader);
  
  // Bind attribute locations.
  // This needs to be done prior to linking.
  glBindAttribLocation(_program, GLKVertexAttribPosition, "position");
  glBindAttribLocation(_program, GLKVertexAttribTexCoord0, "a_texCoord");
  
  // Link program.
  if (![self linkProgram:_program]) {
    NSLog(@"Failed to link program: %d", _program);
    
    if (vertShader) {
      glDeleteShader(vertShader);
      vertShader = 0;
    }
    if (fragShader) {
      glDeleteShader(fragShader);
      fragShader = 0;
    }
    if (_program) {
      glDeleteProgram(_program);
      _program = 0;
    }
    
    return NO;
  }
  
  // Get uniform locations.
  uniforms[UNIFORM_MODELVIEWPROJECTION_MATRIX] = glGetUniformLocation(_program, "modelViewProjectionMatrix");
  uniforms[UNIFORM_SAMPLER] = glGetUniformLocation(_program, "s_texture");
  uniforms[UNIFORM_XRES] = glGetUniformLocation(_program, "xres");
  uniforms[UNIFORM_YRES] = glGetUniformLocation(_program, "yres");
  
  
  // Release vertex and fragment shaders.
  if (vertShader) {
    glDetachShader(_program, vertShader);
    glDeleteShader(vertShader);
  }
  if (fragShader) {
    glDetachShader(_program, fragShader);
    glDeleteShader(fragShader);
  }
  
  return YES;
}

- (BOOL)compileShader:(GLuint *)shader type:(GLenum)type file:(NSString *)file
{
  GLint status;
  const GLchar *source;
  
  source = (GLchar *)[[NSString stringWithContentsOfFile:file encoding:NSUTF8StringEncoding error:nil] UTF8String];
  if (!source) {
    NSLog(@"Failed to load vertex shader");
    return NO;
  }
  
  *shader = glCreateShader(type);
  glShaderSource(*shader, 1, &source, NULL);
  glCompileShader(*shader);
  
#if defined(DEBUG)
  GLint logLength;
  glGetShaderiv(*shader, GL_INFO_LOG_LENGTH, &logLength);
  if (logLength > 0) {
    GLchar *log = (GLchar *)malloc(logLength);
    glGetShaderInfoLog(*shader, logLength, &logLength, log);
    NSLog(@"Shader compile log:\n%s", log);
    free(log);
  }
#endif
  
  glGetShaderiv(*shader, GL_COMPILE_STATUS, &status);
  if (status == 0) {
    glDeleteShader(*shader);
    return NO;
  }
  
  // Generate a texture object
  glGenTextures ( 1, &density_texture );
  
  // Bind the texture object
  glBindTexture ( GL_TEXTURE_2D, density_texture );
  
  // Set the filtering mode
  glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
  glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
  std::vector<float> density;
  self.fluid->GetDensities(&density);
  // Load the texture
  glTexImage2D ( GL_TEXTURE_2D, 0, GL_LUMINANCE, self.width, self.height, 0, GL_LUMINANCE, GL_FLOAT, &density[0] );
  
  return YES;
}

- (BOOL)linkProgram:(GLuint)prog
{
  GLint status;
  glLinkProgram(prog);
  
#if defined(DEBUG)
  GLint logLength;
  glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLength);
  if (logLength > 0) {
    GLchar *log = (GLchar *)malloc(logLength);
    glGetProgramInfoLog(prog, logLength, &logLength, log);
    NSLog(@"Program link log:\n%s", log);
    free(log);
  }
#endif
  
  glGetProgramiv(prog, GL_LINK_STATUS, &status);
  if (status == 0) {
    return NO;
  }
  
  return YES;
}

- (BOOL)validateProgram:(GLuint)prog
{
  GLint logLength, status;
  
  glValidateProgram(prog);
  glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLength);
  if (logLength > 0) {
    GLchar *log = (GLchar *)malloc(logLength);
    glGetProgramInfoLog(prog, logLength, &logLength, log);
    NSLog(@"Program validate log:\n%s", log);
    free(log);
  }
  
  glGetProgramiv(prog, GL_VALIDATE_STATUS, &status);
  if (status == 0) {
    return NO;
  }
  
  return YES;
}

-(void)touchesBegan:(NSSet *)touches withEvent:(UIEvent *)event {
  return;
  if ([touches anyObject] != nil) {
    if (soccer_on) {
      [self tearDownSoccer];
    } else {
      [self setupSoccer];
    }
  }
}

-(void)touchesMoved:(NSSet *)touches withEvent:(UIEvent *)event {
  [touches enumerateObjectsUsingBlock:^(id obj, BOOL *stop) {
    UITouch* touch = obj;
    CGPoint touch_point = [touch locationInView:self.view];
    CGPoint previous_touch_point = [touch previousLocationInView:self.view];
    
    CGPoint fluid_pos = [self fluidCoordsFromScreenCoords:touch_point];
    CGPoint previous_fluid_pos = [self fluidCoordsFromScreenCoords:previous_touch_point];
    self.fluid->AddImpulse(fluid_pos.x, fluid_pos.y, fluid_pos.x - previous_fluid_pos.x, fluid_pos.y - previous_fluid_pos.y);
  }];
}


@end
