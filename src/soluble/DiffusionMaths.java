package soluble;


/**
 * Handles the calculation of diffusion gradients in 3D using the heat kernel. 
 * 
 * Easy to adjust to 2D (or another other dimensionality; not currently implemented) - the coefficient in the
 * heat kernel calculation needs to be adusted to accommodate the new dimensionality. 
 * 
 * @author Mark N. Read, 2018
 *
 */
public class DiffusionMaths 
{
	/**
	 * Calculates the concentration of substance at a location, given that it was secreted at a certain
	 * point in space and time in the past.
	 * 
	 * This assumes no decay of diffusion factor. 
	 * 
	 * The quantity secreted can be multiplied by the returned quantity downstream, which accommodates variable
	 * secretion rates. 
	 * 
	 * 
	 * https://en.wikipedia.org/wiki/Heat_kernel
	 * 
	 * 
	 * @param x0 Location where secretion took place
	 * @param y0
	 * @param z0
	 * @param x Location of detection at current time
	 * @param y
	 * @param z
	 * @param dt Time difference from secretion event to now
	 * @param diffusionCoeff
	 * @return
	 */
	public static double heatKernel3D(double x0, double y0, double z0,  // Secretion source
										 double x, double y, double z,  // Deteection point
										 double dt, // Time since secretion event
										 double diffusionCoeff)
	{
		final double dx = x - x0;
		final double dy = y - y0;
		final double dz = z - z0;
		final double dotProd = dx*dx + dy*dy + dz*dz;
	
		return heatKernel3D(dotProd, dt, diffusionCoeff);
	}
	
	
	public static double heatKernel3D(double dist2,  // Squared Eucliean distance of source from perception point 
									double dt,  // Time since secretion event 
									double diffusionCoeff)
	{
		// Power exponent here is # of dimensions, then square rooted = #d/2.
		final double coeff = 1. / Math.pow((4 * Math.PI * diffusionCoeff * dt), (3./2.));
		double diffusion = Math.exp(-dist2 / (4 * diffusionCoeff * dt)) ;
		return coeff * diffusion;
	}
}
