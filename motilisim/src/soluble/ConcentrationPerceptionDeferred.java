package soluble;

import core.Simulation;
import environment.Compartment3D;
import sim.util.Double3D;

public class ConcentrationPerceptionDeferred implements ConcentrationPerception 
{
	final private Compartment3D compartment;
	final private double sampleTime;
	
	// Locations around the cell where chemoattractant should be measured, if asked. 
	final private Double3D frontLoc;
	final private Double3D backLoc;
	final private Double3D leftLoc;
	final private Double3D rightLoc;	
	final private Double3D upLoc;
	final private Double3D downLoc;
	
	// Calculated on the fly as needed. 
	private double frontConc = Double.NaN;
	private double backConc = Double.NaN;
	private double leftConc = Double.NaN;
	private double rightConc = Double.NaN;
	private double upConc = Double.NaN;
	private double downConc = Double.NaN;
	
	private Double3D gradient;
	private double maxGradientDifferential = Double.NaN;
	
	public ConcentrationPerceptionDeferred(Compartment3D compartment, 
			Double3D front, Double3D back, Double3D left, Double3D right, Double3D up, Double3D down,
			double sampleTime)
	{
		this.compartment = compartment;
		
		frontLoc = front;
		backLoc = back;
		leftLoc = left;
		rightLoc = right;
		upLoc = up;
		downLoc = down;
		
		this.sampleTime = sampleTime;
	}
	
	/** Find maximum absolute concentration in this perception */
	public double maxConcentration()
	{		
		if (Double.isNaN(frontConc))
			frontConc = this.compartment.perceiveAttractantAtPoint(frontLoc);
		
		return frontConc;
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
		if (this.gradient == null)
		{	
			sampleAllConcentrations();
			calculateGradient();
		} 
		return this.gradient;
	}
	
	/** 
	 * Maximum gradient differential measured across any of the cell's four principle axes (up/down, etc)
	 */
	public double maxConcentrationDifferential()
	{
		if (Double.isNaN(this.maxGradientDifferential))
		{	
			sampleAllConcentrations();
			calculateGradient();
		}
		return this.maxGradientDifferential;
	}
	
	/** Sample all concentrations around the cell */
	private void sampleAllConcentrations()
	{
		if (Double.isNaN(frontConc))
			frontConc = this.compartment.perceiveAttractantAtPoint(frontLoc);
		if (Double.isNaN(backConc))  
			backConc = this.compartment.perceiveAttractantAtPoint(backLoc);

		if (Double.isNaN(leftConc))  
			leftConc = this.compartment.perceiveAttractantAtPoint(leftLoc);
		if (Double.isNaN(rightConc))  
			rightConc = this.compartment.perceiveAttractantAtPoint(rightLoc);		
		
		if (Double.isNaN(upConc))  
			upConc = this.compartment.perceiveAttractantAtPoint(upLoc);
		if (Double.isNaN(downConc))  
			downConc = this.compartment.perceiveAttractantAtPoint(downLoc);
	}
	
	private void calculateGradient()
	{
		// Differential of gradients across orthogonal axes (relative to the cell, cell faces 'x'). 
		final double dx = this.frontConc - this.backConc;	
		final double dy = this.leftConc - this.rightConc; 	
		final double dz = this.upConc - this.downConc;
		this.gradient = new Double3D(dx, dy, dz);
		this.maxGradientDifferential = Math.max(Math.max(dx, dy), dz);
	}
	
	/** If the time at which this peception was taken is in the past, the it is out of date; expired */
	public boolean expired()
	{	return Simulation.instance.schedule.getTime() > this.sampleTime;		}
}
