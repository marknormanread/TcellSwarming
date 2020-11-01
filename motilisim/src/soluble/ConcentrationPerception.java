package soluble;

import sim.util.Double3D;



public interface ConcentrationPerception 
{
	/** Find maximum absolute concentration in this perception */
	public double maxConcentration();
	
	/** Returns a Double3D representing the direction of the gradient relative to the orientation of the concentrations
	 *  supplied. Ie, if the concentrations represent points around the cell, then the vector will be relative to the
	 *  cell's orientation. 
	 *  
	 *  The vector is not a unit vector.  
	 *  It may have magnitude zero. 
	 *  Its magnitude corresponds to the gradient differential as measured around the cell.
	 */
	public Double3D gradientDirection();
	
	/** 
	 * Maximum concentration differential measured across any of the cell's three principle axes (up/down, etc).
	 * This is akin to the maximum gradient the cell can perceive. 
	 */
	public double maxConcentrationDifferential();
	
	/** If the time at which this peception was taken is in the past, the it is out of date; expired */
	public boolean expired();
}
