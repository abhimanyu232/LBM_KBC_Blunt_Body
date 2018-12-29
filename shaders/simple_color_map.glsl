/** 
 *  @file
 *  @author Fabian Bösch
 *  simple color map
 */

// returns the appropriate value from the jet color function.
vec3 get_color(float value) 
{
	if (value > 0.5) return vec3(0.0,0.0,1.0);
	else return vec3(0.0);
}
