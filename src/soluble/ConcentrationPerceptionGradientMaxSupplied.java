package soluble;

import core.Simulation;
import sim.util.Double3D;

/**
 * Allows max concentration and gradient to be directly supplied. Max is assumed to be front of cell  
 * 
 * @author Mark N. Read, 2018
 *
 */
public class ConcentrationPerceptionGradientMaxSupplied implements ConcentrationPerception 
{
	final double sampleTime;
	
	private double maxConcentration;
	private Double3D gradient;
	
	public ConcentrationPerceptionGradientMaxSupplied(double maxConc, Double3D gradient, double sampleTime)
	{
		this.sampleTime = sampleTime;
		this.maxConcentration = maxConc;
		this.gradient = gradient;
	}
	
	/** Find maximum absolute concentration in this perception */
	public double maxConcentration()
	{
		return this.maxConcentration;
	}
	
	/** Returns a Double3D representing the direction of the gradient relative to the orientation of the concentrations
	 *  supplied. Ie, if the concentrations represent points around the cell, then the vector will be relative to the
	 *  cell's orientation. 
	 *  
	 *  The vector is not a unit vector.  
	 *  It may have magnitude zero. 
	 *  Its magnitude corresponds to the gradient differential as measured around the cell.
	 */
	public Double3D gradientDirection()
	{
		return this.gradient;
	}
	
	public double maxConcentrationDifferential()
	{
		return Double.NaN;  // Cannot be implemented with how this class works. 
	}
	
	/** If the time at which this peception was taken is in the past, the it is out of date; expired */
	public boolean expired()
	{	return Simulation.instance.schedule.getTime() > this.sampleTime;	}
}
