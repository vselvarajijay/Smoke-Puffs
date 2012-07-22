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

void main()
{
  //gl_FragColor = vec4(v_texCoord.x, v_texCoord.y, 0.0, 1.0);
  gl_FragColor = texture2D(s_texture, v_texCoord);
}
