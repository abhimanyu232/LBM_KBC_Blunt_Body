/** 
 *  @file
 *  @author Fabian Bösch
 *  jet color map
 */

// returns the appropriate value from the jet color function.
vec3 get_color(float value) 
{
     float four_value = 4.0 * value;
     float red   = min(four_value - 1.5, -four_value + 4.5);
     float green = min(four_value - 0.5, -four_value + 3.5);
     float blue  = min(four_value + 0.5, -four_value + 2.5);
     return clamp( vec3(red, green, blue), 0.0, 1.0 );
}
