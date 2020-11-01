package utils;

import core.Simulation;
import ec.util.MersenneTwisterFast;

/**
 * This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
 * 
 * Represents a Levy distributed random variable. This is calculated using the method outlined in
 * Jacobs's "Stochastic Processes for Physicists". 2010, Cambridge University Press.  
 * 
 * @author Mark N. Read
 *
 */
public class LevyDistribution 
{

	private static double bounded_uniform(double low, double high)
	{
		MersenneTwisterFast rng = Simulation.instance.random;
		// returns a double in interval (0,1). IE, neither zero nor one will be returned. 		
		double x = rng.nextDouble(false, false);
		
		double range = high - low;
		x *= range;	// scale onto the required range of values
		x += low;	// translate x onto the values requested
		
		return x;
	}
	
	/**
	 * Samples a Levy distribution wherein the power law decay can be adjusted between 1/x and 1/x^3.
	 * 
	 * This method is based on that found in section 9.2.2 of 
	 * Jacobs's "Stochastic Processes for Physicists". 2010, Cambridge University Press.  
	 * 
	 * Note that this sampling method can return negative values. Values are symmetrical around zero.
	 * 
	 * @param mu must lie between 1 and 3. Corresponds to 1/x and 1/x^3
	 * @return
	 */
	public static double sample(double mu)
	{
		double X = bounded_uniform(-Math.PI/2.0, Math.PI/2.0);
		// uses Mersenne Twister random number generator to retrieve a value between (0,1) (does not include 0 or 1
		// themselves)
		double Y = -Math.log(Simulation.instance.random.nextDouble(false, false));
		double alpha = mu - 1.0;
		// there's a lot going on here, written over several lines to aid clarity.  		
		double Z = 	(	Math.sin(alpha * X) 
						/
						Math.pow( Math.cos(X) , 1.0 / alpha ) 
					) 
					* 
					Math.pow( 
								Math.cos((1.0-alpha) * X) / Y,
								(1.0 - alpha) / alpha)
					;
		return Z;
	}
	
	/** 
	 * Same as above, but ensures all values are positive. Negative values are simply negated, as the Levy distribution
	 * represented is symmetrical around zero.  
	 * @param mu
	 * @return
	 */
	public static double sample_positive(double mu, double scale)
	{
		double l = sample(mu) * scale;	
		if (l < 0.0)	
		{	return -1.0 * l;	}
		return l;
	}
	
	/** Default value case, scale=1 */
	public static double sample_positive(double mu)
	{	return sample_positive(mu, 1.0);		}
}
