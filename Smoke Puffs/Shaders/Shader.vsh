//
//  Shader.vsh
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

attribute vec4 position;
uniform mat4 modelViewProjectionMatrix;

attribute vec2 a_texCoord;
varying vec2 v_texCoord;

void main()
{
  gl_Position = modelViewProjectionMatrix * position;
  v_texCoord = a_texCoord;
}
