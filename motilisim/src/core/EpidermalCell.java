package core;

import java.util.Vector;

import sim.engine.SimState;
import sim.util.Double3D;
import soluble.AttractantSecretor;

@SuppressWarnings("serial")
public class EpidermalCell extends Cell implements AttractantSecretor
{
	public final static double diameter = 10.0;

	private final boolean lysed;
	private static double attractantSecretionRate = 1e6;

	// Where and when the cell has been secreting
	private final Vector<SecretionRecord> secretionHistory = new Vector<SecretionRecord>();

	public EpidermalCell(Double3D location)
	{
		this.lysed = false;
		this.location = location;
	}

	public EpidermalCell(boolean lysed, Double3D location)
	{
		super(Simulation.instance.schedule);
		this.lysed = lysed;
		this.location = location;
	}

	public double getDiameter()
	{	return diameter;	}

	public double getRadius()
	{ 	return diameter / 2.0;	}


	@Override
	public void step(SimState state)
	{
		if(lysed)	{ 	secrete((Simulation)state); 	}
		unSchedule.stop();  // Remove from the schedule; cell lysed.
	}

	public void secrete(Simulation sim)
	{
		secretionHistory.add(new SecretionRecord(
				this.location,  // Location where secretion takes place. Back of cell.
				sim.schedule.getTime(),  // Time of secretion
				attractantSecretionRate));  // Quantity of secretion; the entire cell.
	}

	public boolean isSecretingAttractant()
	{	return lysed;	}

	public double getAttractantSecretionRate()
	{	return attractantSecretionRate;	}

	public static double getMaxAttractantSecretionRate()
	{	return attractantSecretionRate;	}

	public String getType()
	{	return "epidermal cell";	}

	public Vector<SecretionRecord> getSecretionHistory()
	{	return secretionHistory;	}
}
