//
//  Shader.fsh
//  Smoke Puffs
//
//  Created by Matt Stanton on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

precision mediump float;
varying vec2 v_texCoord;
uniform sampler2D s_texture;
uniform float xres;
uniform float yres;

vec4 texture2DBilinear( sampler2D textureSampler, vec2 uv )
{
  // in vertex shaders you should use texture2DLod instead of texture2D
  vec4 tl = texture2D(textureSampler, uv);
  vec4 tr = texture2D(textureSampler, uv + vec2(1.0/xres, 0));
  vec4 bl = texture2D(textureSampler, uv + vec2(0, 1.0/yres));
  vec4 br = texture2D(textureSampler, uv + vec2(1.0/xres , 1.0/yres));
  vec2 f = fract( uv.xy * vec2(xres, yres) ); // get the decimal part
  vec4 tA = mix( tl, tr, f.x ); // will interpolate the red dot in the image
  vec4 tB = mix( bl, br, f.x ); // will interpolate the blue dot in the image
  return mix( tA, tB, f.y ); // will interpolate the green dot in the image
}

void main()
{
  //vec4 color = vec4(v_texCoord.x, v_texCoord.y, 0.0, 1.0);
  //vec4 color = texture2DBilinear(s_texture, v_texCoord);
  vec4 color = texture2D(s_texture, v_texCoord);
  gl_FragColor = vec4(1.0, 1.0, 1.0, color.r);
}