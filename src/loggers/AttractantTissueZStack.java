package loggers;

import java.io.File;

import sim.util.Double3D;
import core.Simulation;
import environment.BoundedCylinder;

/**
 * Class takes a z-stack snapshot of attractant. 
 * @author Mark N. Read, 2017
 *
 */
public class AttractantTissueZStack extends TissueZStack
{
	public AttractantTissueZStack(String dir)
	{	this.dir = dir;		}
	
	@Override
	double getConcentration(Double3D location)
	{
		return ((BoundedCylinder)Simulation.instance.space).perceiveAttractantAtPoint(location);
	}
	
	@Override
	File getImageFile()
	{
		double time = Simulation.instance.schedule.getTime();
		File directory = new File(this.dir + "/AttractantStacks/");
		if (!directory.exists())
			directory.mkdirs();
		
		return new File(directory + "/Attractant-ZStack-" + time + ".png");
	}
}
