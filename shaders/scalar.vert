/** 
 *  @file
 *  @author Fabian Bösch
 *  vertex shader for velocity
 */
 
// attributes: per vertex data
// ---------------------------
// vertex coordinates
attribute vec4 position;
// scalar field
attribute float s;


// varyings: vertex shader output, input to rasterizer
// ---------------------------------------------------
// shaded color for this vertex
varying vec4 color_out;


// uniforms: global data
// ---------------------
// maximum scalar field value
uniform float max_s;
// minimum scalar field value
uniform float min_s;
// alpha value
uniform float alpha;

void main(void)
{
	float s_scaled = clamp((s-min_s)/(max_s-min_s), 0.0, 1.0);
	color_out = vec4(get_color(s_scaled),alpha);
	gl_Position = position;
}
