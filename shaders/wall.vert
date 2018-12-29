/** 
 *  @file
 *  @author Fabian Bösch
 *  vertex shader for walls
 */
 
// attributes: per vertex data
// ---------------------------
// vertex coordinates
attribute vec4 position;
// wall flag
attribute float wall;


// varyings: vertex shader output, input to rasterizer
// ---------------------------------------------------
// shaded color for this vertex
varying vec4 color_out;


void main(void)
{
	color_out = vec4(0.0, 0.0, 0.0, wall);
	gl_Position = position;
}
