package core;

import java.util.Vector;

import sim.engine.Schedule;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.Double3D;
import soluble.AttractantSecretor;

/**
 * 
 * @author Mark N. Read, 2018
 *
 */
public class StationaryChemokineSource extends Cell implements AttractantSecretor
{
	public static double radius = 20;
	
	public Double3D location;
	
	public double attractantSecretionRate;  // Molecules per minute.

	// Where and when the cell has been secreting
	private Vector<SecretionRecord> secretionHistory = new Vector<SecretionRecord>();
	private double secretionMemory = Double.POSITIVE_INFINITY;  // Period of recent secretion events maintained. 
	
	public static enum SecretionType 
	{
		EPOCH,  // Single secretion event at start of simulation.
		CONTINUOUS
	}
	private SecretionType secretionMode;
	
	public StationaryChemokineSource(Double3D location, SecretionType secretionMode, double chemokineSecretionRate)
	{
		this.location = location;
		this.secretionMode = secretionMode;
		this.attractantSecretionRate = chemokineSecretionRate;
	}
	
	@Override
	public void step(SimState state) 
	{
		secrete((Simulation)state);	
	}
	
	public void secrete(Simulation sim)
	{
		final double secretionRate = attractantSecretionRate * Simulation.timeSlice_min;
		
		secretionHistory.add(new SecretionRecord(
				location,  // Location where secretion takes place. Back of cell.
				sim.schedule.getTime(),  // Time of secretion
				secretionRate));  // Quantity of secretion
		trimSecretionHistory();
	}
	
	/** For computational efficiency. Very old secretion histories have negligible effect in present time. */
	private void trimSecretionHistory()
	{
		if (secretionMemory == Double.POSITIVE_INFINITY)
			return;  // Never trim the history. 
		
		double time = Simulation.instance.schedule.getTime();
		Vector<SecretionRecord> delete = new Vector<SecretionRecord>();
		for (SecretionRecord h : this.secretionHistory)
			if (h.getTime() + secretionMemory < time)  // Time is expressed in minutes
				delete.add(h);
		this.secretionHistory.removeAll(delete);
	}

	@Override
	public Vector<SecretionRecord> getSecretionHistory() 
	{	return secretionHistory;	}

	@Override
	public boolean isSecretingAttractant() 
	{	
		if (secretionMode == SecretionType.EPOCH)
			return Simulation.instance.schedule.getTime() == 0.;
	
		// Continuous secretion case. 
		return true;	
	}

	@Override
	public double getAttractantSecretionRate() 
	{
		return attractantSecretionRate;
	}

	@Override
	public double getDiameter() 
	{	return 2. * getRadius();		}

	@Override
	public double getRadius() 
	{	return radius;	}

	@Override
	public String getType() {
		return "Stationary chemoattractant source";
	}
}
