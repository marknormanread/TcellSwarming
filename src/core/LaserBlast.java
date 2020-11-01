package core;

import environment.BoundedCubeNeutrophil;
import environment.Compartment3D;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.Double3D;
import sim.util.IntBag;

/**
 * Class implements a laser blast that damages the epidermal skin cells, causing them
 * to secrete damage signals. 
 * 
 * 
 * @author Mark N. Read
 *
 */
@SuppressWarnings("serial")
public class LaserBlast implements Steppable
{
	private double blastStartTime;
	private double blastRadius = 20.0;
	public Double3D location;
	
	public LaserBlast(double startTime)
	{	
		this.blastStartTime = startTime;		
	}
	
	@Override
	public void step(SimState state) 
	{
		/* Will try to find a location to place the lysed epidermal cell, near the centre of the
		 * tissue space (though on the surface). */ 		
		IntBag xPos = new IntBag();  // Candidate locations will be placed here. 
		IntBag yPos = new IntBag();
		
		// Find all locations in 2D within the laser blast that might be suitable. 
		Compartment3D.radial2DLocations (
			(int)Math.round(Simulation.space.extremeWidth / 2.0), 
			(int)Math.round(Simulation.space.extremeHeight / 2.0), 
			(int)Simulation.space.extremeWidth, (int)Simulation.space.extremeHeight, 
			blastRadius, xPos, yPos);		 				 
		
		/* Try to find a location that is not within a follicle, and is within the laser blast radius.
		 * If we run out of possible locations, because the blast radius is smaller than a follicle
		 * placed dead centre in the tissue volume, we place the epidermal cell there anyway.
		 * This may be interpreted as the laser burning a hair. */ 
		boolean foundLocation = false;
		while(foundLocation == false)
		{
			int i = Simulation.instance.random.nextInt(xPos.size());
			int j = Simulation.instance.random.nextInt(yPos.size());
			double x = xPos.get(i);
			double y = yPos.get(j);	
			this.location = new Double3D(x,y,0);
			if (! ((BoundedCubeNeutrophil)Simulation.space).insideFollicle(this.location, blastRadius)) {
			 	foundLocation = true;	
			} else {
				xPos.remove(i);	 // Remove this location, since it is not suitable. 
				yPos.remove(j);				
			}
			// Safety, in case there is a hair follicle greater in size than the blast radius. 
			if (xPos.size() == 0 || yPos.size() == 0)	
			{	foundLocation = true;	}
		}
		((BoundedCubeNeutrophil)Simulation.space).placeEpidermalCell
			(new EpidermalCell(true, this.location), this.location);
	}	
}