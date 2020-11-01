package environment;

import java.util.ArrayList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.Cell;
import core.HairFollicle;
import core.Simulation;
import sim.field.continuous.Continuous3D;
import sim.field.grid.IntGrid3D;
import sim.util.Double3D;
import sim.util.IntBag;
import soluble.Attractant;
import soluble.HeatKernelLookupTableArrayList;

/**
 * Compartment represents a cubic (or rectangular) unbounded imaging volume (cells can leave and re-enter).
 * 
 * 
 * @author Mark N. Read, 2017
 *
 */
public class BoundedCube extends Compartment3D
{
	public static double width;
	public static double height;
	public static double depth;
	public static double side_bufferSize;
	public static double bottom_bufferSize;  // Allows top and bottom bounded cube
	public static double top_bufferSize;  // Allows top and bottom bounded cube
	
	public static double follicleRatio = 0.0;  // Looking top down, how much of the skin area comprises hair?
	
	public IntGrid3D restrictedField;  // Stores restrictions on cell movement. 0 = no restriction, 1 = blocked.
	
	private boolean bottomBounded;
	private boolean topBounded;
	private double volumeSimulated;
	private double volumeImaged;
	
	
	public BoundedCube(boolean topBounded, boolean bottomBounded)
	{
		this.topBounded = topBounded;
		this.bottomBounded = bottomBounded;
		if (bottomBounded)
			bottom_bufferSize = 0;
		if (topBounded)
			top_bufferSize = 0;
		
		cellField = new Continuous3D(
				/* Discretization, dividing space into regions for maintaining a map of objects' locations. 
				 * This can have a considerable detrimental effect on computational efficiency if poorly chosen. */
				15.0,    
				width, height, depth);
		cellField.clear();		
		restrictedField = new IntGrid3D((int)width, (int)height, (int)depth, 0);
		createHairFollicles();

		// Volume simulated can be larger (for unbounded environments) than that simulated. 
		// Allows realistic cell entry/departure
		double Vw = width + (2.0 * side_bufferSize * width);
		double Vh = height + (2.0 * side_bufferSize * height);
		double Vd = depth + (bottom_bufferSize * depth) + (top_bufferSize * depth);
		volumeSimulated = Vw * Vh * Vd;  // total volume, including imaging volume and buffer zone.
		volumeImaged = width * height * depth;  // volume of the imaging volume.
		
		attractantLookup = new HeatKernelLookupTableArrayList
				(6.0, Simulation.instance.timeSlice_min, Attractant.diffusion);
	}
	
	/**
	 * Populate the space with cells, placing each cell randomly in the volume such that cells do not overlap. 
	 */
	@Override
	public void placeCellRandomly(Cell cell)
	{
		Double3D loc;				// will try to find a location for the neutrophil.
		// local copies of parameter values, to make equations below more readable. 
		final double w = width;
		final double h = height;
		final double d = depth;
		final double bf = side_bufferSize;
		final double dbf = bottom_bufferSize;
		final double tbf = top_bufferSize;
		do
		{
			// Select a random location.  
			// Cells can be placed in the tissue volume, and similar sized volumes all around it (except 
			// above it, because that breaches the skin).
			final double lw = (Simulation.instance.random.nextDouble() * (w + (2.0 * w * bf))) - (w * bf);
			final double lh = (Simulation.instance.random.nextDouble() * (h + (2.0 * h * bf))) - (h * bf);
			final double ld = (Simulation.instance.random.nextDouble() * (d + (d * dbf) + (d * tbf))) - (d * tbf);
			
			loc = new Double3D(lw, lh, ld);
		} while (! isOccupiableSpaceStartup(loc, cell));
		this.cellField.setObjectLocation(cell, loc);
		cell.setCurrentLocation(loc);
	}
	
	
	/** 
	 * Concrete classes queried on whether the sphere placed at location with given radius breaches any other 
	 * environmental objects.
	 * Currently used only at startup. Collision detection is used during runtime. 
	 * Overridden by subclasses. 
	 */
	@Override
	protected boolean restrictedSpace(Double3D location, double radius)
	{ 
		if (bottomBounded && (location.z > depth))
			return true;
		if (topBounded && (location.z < 0.0))			
			return true;
		try{
			// Currently just checks for hair follicles
			if (restrictedField.get((int)location.x, (int)location.y, (int)location.z) == 1)
				return true;
		} catch (ArrayIndexOutOfBoundsException e) {
			// Do nothing. In unbounded space cells can be placed outside of areas that the restricted field 
			// represents. It only represents bounded space in the imaging volume.
		}
		return false;
	}	
	
	@Override
	public boolean insideImagingVolume(double x, double y, double z)
	{
		if (x < 0)			return false;
		if (x > width)		return false;
		if (y < 0)			return false;
		if (y > height) 	return false;
		if (z < 0) 			return false;
		if (z > depth) 		return false;
		
		return true;
	}
	
	public boolean insideRecordingVolume(double x, double y, double z, Cell cell)
	{
		return true;  // Not currently implemented, everything is in the recording volume.
	}

	/**
	 * Detects collisions with boundaries of the cube (as defined in its creation). 
	 * @param cell cell to make a move
	 * @param move the move to be made
	 */
	@Override
	public MoveResults boundaryCollisionDetection(Cell cell, Double3D move)
	{
		if (topBounded && move.z < 0)  // Moving towards "top" (boundary, z = 0)
			return boundaryCollisionDetectionTop(cell, move);
		else if (bottomBounded && move.z > 0)  // Cell moving towards "bottom" (boundary, z > 0)
		{
			return boundaryCollisionDetectionBottom(cell, move, depth);
		}
		// No collision detected. 
		final Double3D cl = cell.getCurrentLocation();  
		Double3D proposed_location = new Double3D(cl.x + move.x, cl.y + move.y, cl.z + move.z);
		return new MoveResults(proposed_location, new Double3D(0,0,0), move, new ArrayList<Cell>());
	}

	/**
	 * Method places hair follicles in the environment. 
	 */
	private void createHairFollicles()
	{
		// How much skin area is being represented?
		double totalArea = width * height;  
		// How much area of hair follicle should there be?
		double requiredHairCoveredArea = totalArea * follicleRatio;	
		
		while(requiredHairCoveredArea > 0.0)		
		{			
			HairFollicle h = new HairFollicle();
			Double3D loc;  // The proposed location of the new follicle. 
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
			IntBag xPos = new IntBag();  // Coordinates of locations placed in these two. 
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
		
	@Override
	public double volumeImaged() 
	{	return volumeImaged;	}

	@Override
	public double volumeSimulated() 
	{	return volumeSimulated;		}
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;	
		n = (Node) xPath.compile("/params/Environment/BoundedCube").evaluate(params, XPathConstants.NODE);;
		if (n != null)
		{
			n = (Node) xPath.compile("/params/Environment/BoundedCube/width").evaluate(params, XPathConstants.NODE);
			width = Double.parseDouble(n.getTextContent());
			n = (Node) xPath.compile("/params/Environment/BoundedCube/height").evaluate(params, XPathConstants.NODE);
			if(n.getTextContent() != null && n.getTextContent().length() != 0)
				height = Double.parseDouble(n.getTextContent());
			else 
				height = width;
			n = (Node) xPath.compile("/params/Environment/BoundedCube/depth").evaluate(params, XPathConstants.NODE);
			depth = Double.parseDouble(n.getTextContent());
			n = (Node) xPath.compile("/params/Environment/BoundedCube/bufferSize")
					.evaluate(params, XPathConstants.NODE);
			side_bufferSize = Double.parseDouble(n.getTextContent());
			n = (Node) xPath.compile("/params/Environment/BoundedCube/bufferSize")
					.evaluate(params, XPathConstants.NODE);
			side_bufferSize = Double.parseDouble(n.getTextContent());
			
			n = (Node) xPath.compile("/params/Environment/BoundedCube/follicleRatio")
					.evaluate(params, XPathConstants.NODE);
			if (n != null && n.getTextContent().length() != 0)
				follicleRatio = Double.parseDouble(n.getTextContent());
		} else
			// deprecated use of parameters
			loadParameters_backwardsCompatibility(params);
		
		bottom_bufferSize = side_bufferSize;  // both of these may be set to zero depending on cube boundings
		top_bufferSize = side_bufferSize;  		
		
		// needed for MASON's GUI. Stored in Compartment3D because the GUI doesn't know at launch time what compartment
		// is being used (and thus which one to check).
		extremeWidth = width;
		extremeHeight = height;
		extremeDepth = depth;
	}
	
	/** Deprecated parameter structure, but provided for backwards compatibility with older parameters files */
	private static void loadParameters_backwardsCompatibility(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;	
		n = (Node) xPath.compile("/params/Simulation/tissueWidth").evaluate(params, XPathConstants.NODE);
		width = Integer.parseInt(n.getTextContent());
		height = width;
		n = (Node) xPath.compile("/params/Simulation/tissueDepth").evaluate(params, XPathConstants.NODE);
		depth = Integer.parseInt(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Simulation/bufferSize").evaluate(params, XPathConstants.NODE);
		side_bufferSize = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Simulation/Neutrophils/follicleRatio").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() != 0)
			follicleRatio = Double.parseDouble(n.getTextContent());
		else
			follicleRatio = 0.;
	}
}
