/** 
 *  @file
 *  @author Fabian Bösch
 *  vertex shader for velocity magnitude
 */
 
// attributes: per vertex data
// ---------------------------
// vertex coordinates
attribute vec4 position;
// velocity in x-direction
attribute float u;
// velocity in y-direction
attribute float v;


// varyings: vertex shader output, input to rasterizer
// ---------------------------------------------------
// shaded color for this vertex
varying vec4 color_out;


// uniforms: global data
// ---------------------
// maximum velocity magnitude squared (for scaling)
uniform float max_magnitude2;
// maximum velocity magnitude squared (for scaling)
uniform float min_magnitude2;
// alpha value
uniform float alpha;

void main(void)
{
	float mag2_scaled =  clamp(((u*u + v*v)-min_magnitude2)/(max_magnitude2-min_magnitude2), 0.0, 1.0);
	color_out = vec4(get_color(mag2_scaled),alpha);
	gl_Position = position;
}
