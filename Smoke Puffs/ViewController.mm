//
//  ViewController.m
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "ViewController.h"

#import <QuartzCore/QuartzCore.h>

#import "Fluid.h"

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

// Uniform index.
enum
{
    UNIFORM_MODELVIEWPROJECTION_MATRIX,
    UNIFORM_NORMAL_MATRIX,
    NUM_UNIFORMS
};
GLint uniforms[NUM_UNIFORMS];

// Attribute index.
enum
{
    ATTRIB_VERTEX,
    ATTRIB_NORMAL,
    NUM_ATTRIBUTES
};

@interface ViewController () {
  GLuint _program;
  
  GLKMatrix4 _modelViewProjectionMatrix;
  GLKMatrix3 _normalMatrix;
  float _rotation;
  
  GLuint _vertexArray;
  GLuint _vertexBuffer;
  
  int _width;
  int _height;
  
  Fluid* _fluid;
  
  float _dt;
  
  float _ball_x;
  float _ball_y;
}
@property (strong, nonatomic) EAGLContext *context;
@property (strong, nonatomic) GLKBaseEffect *effect;

@property (nonatomic, assign) int width;
@property (nonatomic, assign) int height;

@property (nonatomic, assign) Fluid* fluid;
@property (nonatomic, assign) float dt;

@property (nonatomic, assign) float ball_x;
@property (nonatomic, assign) float ball_y;


- (void)setupGL;
- (void)tearDownGL;

- (BOOL)loadShaders;
- (BOOL)compileShader:(GLuint *)shader type:(GLenum)type file:(NSString *)file;
- (BOOL)linkProgram:(GLuint)prog;
- (BOOL)validateProgram:(GLuint)prog;
@end

@implementation ViewController

@synthesize context = _context;
@synthesize effect = _effect;
@synthesize width = _width;
@synthesize height = _height;
@synthesize fluid = _fluid;
@synthesize dt = _dt;
@synthesize ball_x = _ball_x;
@synthesize ball_y = _ball_y;

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
  
  self.width = 60;
  self.height = 45;
  self.fluid = new Fluid(self.width, self.height);
  self.dt = 0.25f;
  self.ball_x = self.width * 0.5f;
  self.ball_y = self.height * 0.5f;
  
  [self setupGL];
  [self.view setMultipleTouchEnabled:YES];
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
    if ([[UIDevice currentDevice] userInterfaceIdiom] == UIUserInterfaceIdiomPhone) {
        return (interfaceOrientation != UIInterfaceOrientationPortraitUpsideDown);
    } else {
        return YES;
    }
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
  self.fluid->AdvectPoint(self.dt, self.ball_x, self.ball_y, &_ball_x, &_ball_y);
  
  NSArray *subviews = [self.view subviews];
  for (UIView *view in subviews) {
    [view removeFromSuperview];
  }
  
  UIView *touchView = [[UIView alloc] init];
  [touchView setBackgroundColor:[UIColor redColor]];
  touchView.frame = CGRectMake(self.ball_x * 1024 / self.width, 768 - self.ball_y * 768 / self.height, 30, 30);
  touchView.layer.cornerRadius = 15;
  [self.view addSubview:touchView];
}

- (void)glkView:(GLKView *)view drawInRect:(CGRect)rect
{
  std::vector<float> lines;
  self.fluid->GetLines(&lines, 3.0);
  
  glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  self.effect.transform.projectionMatrix = GLKMatrix4MakeOrtho(0, self.width, 0, self.height, 1, -1);
  self.effect.transform.modelviewMatrix = GLKMatrix4Identity;
    
  // Render the object with GLKit
  [self.effect prepareToDraw];

  glEnableVertexAttribArray(GLKVertexAttribPosition);
  glVertexAttribPointer(GLKVertexAttribPosition, 2, GL_FLOAT, GL_FALSE, 0, &(lines[0]));
  glDrawArrays(GL_LINES, 0, lines.size()/2);
  glDisableVertexAttribArray(GLKVertexAttribPosition);
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
    glBindAttribLocation(_program, ATTRIB_VERTEX, "position");
    glBindAttribLocation(_program, ATTRIB_NORMAL, "normal");
    
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
    uniforms[UNIFORM_NORMAL_MATRIX] = glGetUniformLocation(_program, "normalMatrix");
    
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

-(void)touchesMoved:(NSSet *)touches withEvent:(UIEvent *)event {  
  // Remove old red circles on screen
//  NSArray *subviews = [self.view subviews];
//  for (UIView *view in subviews) {
//    [view removeFromSuperview];
//  }
  
  // Enumerate over all the touches and draw a red dot on the screen where the touches were
  [touches enumerateObjectsUsingBlock:^(id obj, BOOL *stop) {
    // Get a single touch and it's location
    UITouch *touch = obj;
    CGPoint touchPoint = [touch locationInView:self.view];
    CGPoint previousTouchPoint = [touch previousLocationInView:self.view];
    
    // Draw a red circle where the touch occurred
//    UIView *touchView = [[UIView alloc] init];
//    [touchView setBackgroundColor:[UIColor redColor]];
//    touchView.frame = CGRectMake(touchPoint.x, touchPoint.y, 30, 30);
//    touchView.layer.cornerRadius = 15;
//    [self.view addSubview:touchView];
//    
//    UIView *previousTouchView = [[UIView alloc] init];
//    [previousTouchView setBackgroundColor:[UIColor blueColor]];
//    previousTouchView.frame = CGRectMake(previousTouchPoint.x, previousTouchPoint.y, 30, 30);
//    previousTouchView.layer.cornerRadius = 15;
//    [self.view addSubview:previousTouchView];
    
    self.fluid->AddImpulse(touchPoint.x * self.width / 1024,
                           self.height - touchPoint.y * self.height / 768,
                           (touchPoint.x - previousTouchPoint.x) * self.width / 1024,
                           -(touchPoint.y - previousTouchPoint.y) * self.height / 768);
    NSLog(@"%g %g %g %g\n", touchPoint.x, touchPoint.y, previousTouchPoint.x, previousTouchPoint.y);
  }];
}


@end
