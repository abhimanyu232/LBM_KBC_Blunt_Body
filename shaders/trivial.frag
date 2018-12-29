/** 
 *  @file
 *  @author Fabian Bösch
 *  trivial fragment shader
 */

// varyings: vertex shader output, input to rasterizer
// ---------------------------------------------------
// shaded color for this vertex
varying vec4 color_out;


// uniforms: global data
// ---------------------

void main(void) 
{
    gl_FragColor = color_out;
}
