package reenact.neutrophil;

import sim.util.Double3D;

public class NeutrophilStateCollection 
{
	public double time;
	public Double3D location;
	public boolean ltb4Stimulated;
	public boolean attractantStimulated;
	public boolean inSwarm;
	
	public NeutrophilStateCollection(double time, Double3D loc, boolean attractant, boolean ltb4, boolean inSwarm)
	{
		this.time = time;
		this.location = loc;
		this.attractantStimulated = attractant;
		this.ltb4Stimulated = ltb4;
		this.inSwarm = inSwarm;
	}
}