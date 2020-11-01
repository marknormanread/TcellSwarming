package reenact;

import sim.util.Double3D;

public class CellStateCollection 
{
	public double time;
	public Double3D location;
	public boolean chemotactic;
	public int timePointID;
	public int trackID;
	
	public CellStateCollection(Double3D loc, boolean chemotactic, double time, int timePointID, int trackID)
	{
		this.time = time;
		this.location = loc;
		this.chemotactic = chemotactic;
		this.timePointID = timePointID;
		this.trackID = trackID;
	}
	
	@Override
	public String toString()
	{
		return "trackID: " + trackID + "; timePointID: " + timePointID + "; time:" + time + 
				"; location: " + location + "; " + chemotactic;
	}
}
