package environment;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Vector;

import core.Cell;
import core.EpidermalCell;
import core.HairFollicle;
import core.SecretionRecord;
import core.Simulation;
import core.SimulationNeutrophil;
import sim.field.continuous.Continuous3D;
import sim.field.grid.DoubleGrid3D;
import sim.util.Double3D;
import sim.util.IntBag;
import soluble.DiffusionMaths;
import soluble.LTB4;
import soluble.LTB4Secretor;

/**
 * Specialization of the BoundedCube environment to include soluble factors pertinent to netrophil swarming. 
 * 
 * @author Mark N. Read, 2018
 *
 */
public class BoundedCubeNeutrophil extends BoundedCube
{
	public DoubleGrid3D attractantField;  	// discrete field for damage signal. 
	public DoubleGrid3D ltb4Field;			// discrete field for ltb4
	public Continuous3D sampleField;		// stores points at which space is to be sampled for concentrations of
											// soluble factor. 
	
	public ArrayList<EpidermalCell> lysedEpidermals = new ArrayList<EpidermalCell>();
	
	/* Store a list of cells secreting the respective soluble factors, and the time at which the list was compiled.
	 * This is for efficiency purposes, to prevent the list being needlessly recompiled in a single timestep. */
	private SecretorsTimepoint<LTB4Secretor> ltb4Secretors = new SecretorsTimepoint<LTB4Secretor>(null, -1.0);
	
	
	public BoundedCubeNeutrophil(boolean topBounded, boolean bottomBounded)
	{
		super(topBounded, bottomBounded);

		createHairFollicles();	
		findAttractantSecretors();  // Initialise lists of secretors. 
		findLTB4Secretors();		
	}
		
	public void placeEpidermalCell(EpidermalCell cell, Double3D location)
	{
		System.out.println("Placing epidermal cell at location " + location);
		this.cellField.setObjectLocation(cell, location);
		this.lysedEpidermals.add(cell);
	}
	
	protected boolean restrictedSpace(Double3D location, double radius)
	{
		return insideFollicle(location, radius);
	}
	
	/**
	 * Method places hair follicles in the environment. 
	 */
	private void createHairFollicles()
	{
		double totalArea = width * height;  // How much skin area is being represented?
		// How much area of hair follicle should there be?
		double requiredHairCoveredArea = totalArea * SimulationNeutrophil.follicleRatio;	
		
		while(requiredHairCoveredArea > 0.0)		
		{			
			HairFollicle h = new HairFollicle();
			Double3D loc;								// the proposed location of the new follicle. 
			do{				
				loc = new Double3D(
						Math.round(width * Simulation.instance.random.nextDouble()),
						Math.round(height * Simulation.instance.random.nextDouble()),
						restrictedField.getLength()
						);
			} while(! sufficientlySeparatedFollicle(h, loc));
			h.setLocation(loc);
			cellField.setObjectLocation(h, loc);
			// Set value of restrictedField locations comprising the hair follicle to 1. 
			// This is done by finding the circle within follicle radius at each z-layer in the tissue in turn. 
			IntBag xPos = new IntBag();  // coordinates of locations placed in these two. 
			IntBag yPos = new IntBag();			
			radial2DLocations((int)Math.round(loc.x), (int)Math.round(loc.y),  // The center of the hair follicle. 
						restrictedField.getWidth(), restrictedField.getHeight(),  // Size of the grid.
						(h.getDiameter()/2.0),  // The radius of the hair follicle.
						xPos, yPos);  // Coordinates within radius placed in here.
			// Set the locations' values to 1. 
			for(int z = 0; z < depth; z++)
				for (int i = 0; i < xPos.numObjs; i++)								
					restrictedField.set(xPos.get(i), yPos.get(i), z, 1);				
			
			// Area of circle is pi r^2. Implicit conversion from diameter to radius. 
			double hArea = Math.PI * ((0.5 * h.getDiameter()) * (0.5 * h.getDiameter()));  
			requiredHairCoveredArea -= hArea;
			
			follicles.add(h);  // Store follicle for use later.
		}		
	}
	
	/**
	 * Returns true if the supplied sphere, consisting of the supplied radius and location, resides within a hair 
	 * follicle. By setting the radius to zero, this method can indicate if a particular point in space resides within
	 * a follicle.  
	 */
	public boolean insideFollicle(Double3D location, double radius)
	{
		// Scan through every follicle checking if the location lies within its boundary.
		for (HairFollicle hf : follicles)
		{
			Double3D hfLoc = hf.getCentreLocation();
			double dx = location.x - hfLoc.x;  // Differences in x and y coordinates. 
			double dy = location.y - hfLoc.y;  // Note no need to consider z, since follicles are cylinders along z. 
			double dist2 = dx*dx + dy*dy;
			double minDist2 = (hf.getRadius() + radius) * (hf.getRadius() + radius);
			if (dist2 < minDist2)  // Comparing squares more efficient than taking square roots. 
				return true;				
		}
		return false;
	}
	
	/** Check that the proposed hair follicle can reside at the proposed location. This entails checking that the area 
	 *  it will cover does not infringe on the area of any other existing follicle. 
	 * 
	 * @param proposedFollicle The new follicle that we are attempting to place
	 * @param loc The proposed location of the new follicle.  
	 */
	private boolean sufficientlySeparatedFollicle(HairFollicle proposedFollicle, Double3D loc)
	{
		// Check for sufficient distance between each follicle, and the proposed location of a new one, `loc`.
		for(HairFollicle h : follicles)
		{			
			// Minimum allowable separation is the sum of both hairs' radius, and 10% extra.  
			double minSep = 1.1* (h.getRadius() + proposedFollicle.getRadius());
			double minSep2 = minSep * minSep;
			Double3D hLoc = cellField.getObjectLocation(h);
			double actualSeparation2 = ((loc.x - hLoc.x) * (loc.x - hLoc.x)) + ((loc.y - hLoc.y) * (loc.y - hLoc.y));
			if(actualSeparation2 < minSep2)
				return false;
		}		
		return true;
	}
	
	/**
	 * Will return a collection of all the hair follicles which the supplied cell will collide with
	 * were it to be placed at the specified location. 
	 */
	private Collection<HairFollicle> collidingFollicles(Double3D cellLoc, Cell cell)
	{
		Vector<HairFollicle> colliders = new Vector<HairFollicle>();
		for (HairFollicle hf : follicles)
		{
			Double3D hfLoc = cellField.getObjectLocation(hf);
			// Collision if the actual distance between follicle and cell is less than this. 
			double minDistance = hf.getRadius() + cell.getRadius();
			double minDistance2 = minDistance*minDistance;
			double xDiff = cellLoc.x - hfLoc.x;
			double yDiff = cellLoc.y - hfLoc.y;
			double distance2 = xDiff*xDiff + yDiff*yDiff;
			if (distance2 < minDistance2)	// Comparing distances squared is more efficient than taking sqrt. 
				colliders.add(hf);
		}		
		return colliders;
	}
	
	
	/** Returns a vector containing all the entities currently secreting LTB4.  
	 */
	private LTB4Secretor[] findLTB4Secretors()
	{	
		// Compile a new list of secretors only if time has elapsed since last time. 
		if (ltb4Secretors.expired())
		{	
			ArrayList<LTB4Secretor> secretors = new ArrayList<LTB4Secretor>();
		
			for (Object o : cellField.allObjects)
				if (o instanceof LTB4Secretor && ((LTB4Secretor)o).getSecretionHistory().size() > 0)
					secretors.add((LTB4Secretor)o);
			
			double timeNow = Simulation.instance.schedule.getTime();
			LTB4Secretor[] array = new LTB4Secretor[secretors.size()];
			ltb4Secretors = new SecretorsTimepoint<LTB4Secretor>(secretors.toArray(array), timeNow);
		}
		
		return ltb4Secretors.getSecretors();
	}
			
	/**
	 * Method returns the concentration of LTB4 at the supplied `location` based on the PDE heat equation.
	 * 
	 * The equation is too complex to write here in the comments, see the accompanying simulation documentation for 
	 * specifics.  
	 */
	public double perceiveLTB4AtPoint(Double3D location)
	{
		// get all the LTB4-secreting neutrophils.
		LTB4Secretor[] secretors = findLTB4Secretors();
		double concentration = 0.0;  // Cumulative calculation.
		for (LTB4Secretor s : secretors)
		{
			Vector<SecretionRecord> history = s.getSecretionHistory();
			for (SecretionRecord r : history)
			{
				double dt = Simulation.instance.schedule.getTime() - r.getTime();	// Time passed since secretion.
				// Safety. Updating neutrophils is not an atomic operation. 
				// Some may have recorded secreting at this timestep before another evaluates the gradient.  
				// This ensures that only secretion records from the  past are considered.
				if (dt != 0.0)
				{
					double quantity = r.getQuantity() * 
							DiffusionMaths.heatKernel3D(
							r.getLocation().x, r.getLocation().y, r.getLocation().z, 
							location.x, location.y, location.z, 
							dt, LTB4.diffusion);
					
					// Apply decay, if relevant.
					assert Double.isFinite(LTB4.diffusion);
					if (LTB4.decayConstant != 0.)
						quantity *= Math.exp(-LTB4.decayConstant * dt);
					
					concentration += quantity; 
				}
			}
		}		
		return concentration;
	}
}
