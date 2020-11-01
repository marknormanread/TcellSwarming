package environment;

import java.util.ArrayList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.Cell;
import core.Simulation;
import core.StationaryChemokineSource;
import core.TCell;
import core.SimulationChemotacticAgent.BolusChemokineSourceType;
import jogamp.opengl.SystemUtil;
import loggers.CellLogger;
import sim.field.continuous.Continuous3D;
import sim.util.Double2D;
import sim.util.Double3D;
import soluble.Attractant;
import soluble.HeatKernelLookupTable2DArray;
import soluble.StationaryDiffusionTimeIntegration;
import utils.Utils;


/**
 * Created to represent Biro's entire-well imaging experiments. 
 * 
 * @author Mark N. Read
 *
 */
public class BoundedCylinder extends Compartment3D
{
	public static double radius;  // Dimensions of well
	public static double depth;
	
	public static boolean bolusPresent = false; // Default 
	public static double bolusRadius;
	private static double bolusSpeed;
	public static BolusCylinder bolus;  // Does nothing, used for visualisation
	private static boolean bolusTCells = false; // Are there T cells in the bolus at instigation?
	public static Double3D centre;  // Location of cylinder cross-section centre, used to bound the space.
	
	private double volumeSimulated;
	private double volumeImaged;	
	
	public static enum BolusChemokineSourceType {
		NONE,
		EPOCH,  // Single secretion event at start of simulation.
		CONTINUOUS
	}
	public static BolusChemokineSourceType bolusChemokineSourceMode = BolusChemokineSourceType.NONE;
	public static double bolusChemokineSecretionRate = 0;

	public BoundedCylinder()
	{
		cellField = new Continuous3D(
				/* Discretisation, dividing space into regions for maintaining a map of objects' locations. 
				 * This can have a considerable detrimental effect on computational efficiency if poorly chosen. */
				15.0,    
				2*radius, 2*radius, depth);
		cellField.clear();	
		
		volumeSimulated = Math.PI * radius * radius * depth;
		volumeImaged = volumeSimulated;
		if (bolusPresent)
		{
			bolus = new BolusCylinder();
			cellField.setObjectLocation(bolus, centre);
		}
		
		double max_distance = BoundedCylinder.radius * 2;  // Diameter. 
		double max_time = 960;
		
		switch (secretionMotilityMode)
		{
		case MOTILE: 
			attractantLookup = new HeatKernelLookupTable2DArray
						(6.0, Simulation.instance.timeSlice_min, Attractant.diffusion, max_distance, max_time);
			break;
		default:
			// Note that this assumes stationary secreters that secrete CONTINUOUSLY.
			attractantLookup = new StationaryDiffusionTimeIntegration
						(6.0, Simulation.instance.timeSlice_min, Attractant.diffusion, max_distance, max_time);
			break;
		}
		
		// Place secretor in bolus centre, if specified in parameters file. 
		if (bolusChemokineSourceMode != BolusChemokineSourceType.NONE)
		{
			// The following functionality relies on having a bolus, currently only implemented in BoundedCylinder.
			assert(Simulation.space instanceof BoundedCylinder); 
			assert(bolusPresent);		
			StationaryChemokineSource chemokineSource = null;
			
			if (bolusChemokineSourceMode==BolusChemokineSourceType.EPOCH) 
			{
				chemokineSource = new StationaryChemokineSource
						(centre, StationaryChemokineSource.SecretionType.EPOCH, bolusChemokineSecretionRate);
				Simulation.instance.schedule.scheduleOnce(0., Simulation.compartmentOrdering, chemokineSource);								

			}
			if (bolusChemokineSourceMode==BolusChemokineSourceType.CONTINUOUS)	
			{
				chemokineSource = new StationaryChemokineSource
						(centre, StationaryChemokineSource.SecretionType.CONTINUOUS, bolusChemokineSecretionRate);
				Simulation.instance.schedule.scheduleRepeating(0., Simulation.compartmentOrdering, 
						                   					   chemokineSource, 
						                   					   Simulation.timeSlice_min);
			}
			cellField.setObjectLocation(chemokineSource, centre);
		}
			
	}
	
	/** The supplied proposed move by the supplied cell may be curtailed. */
	@Override
	public MoveResults boundaryCollisionDetection(Cell cell, Double3D move) {
		MoveResults mr = null;
		if (bolusPresent)
		{
			// New location can be ignored (more curtailing to come), no bounce and no collisions. 
			// Only move needs update.
			mr = bolusInfiltrationDetection(cell, move);
			move = mr.executedMove;  
		}
		// Top and bottom boundary detection
		if (move.z < 0) // Moving towards "top" (boundary, z = 0)
		{
			mr = boundaryCollisionDetectionTop(cell, move);
		} else {  // Cell moving towards "bottom" (boundary, z > 0)
			mr = boundaryCollisionDetectionBottom(cell, move, depth);
		}
		// Top and bottom collision detection may have curtailed desired cell move, feed that for further restriction
		// into cylinder collision detection. 
		MoveResults cmr = cylinderCollisionDetection(cell, mr.executedMove);  
		cmr.bounce = mr.bounce.add(cmr.bounce);		
		// At this point the MoveResults cmr contains a move that could place a cell on the cylinder boundaries, but
		// not beyond them. 
		return cmr;
	}	
	
	/**
	 * Operates in a 2D plane, conceptualising the cylinder as a circle. Z-boundary breaches handled elsewhere. 
	 * Assuming the cell is inside the cylinder (circle), a check is performed on whether the proposed move would take
	 * the cell outside of the circle. If not, nothing else needs to be done. 
	 * 
	 * If it does, then the point of contact with the circle needs to be calculated. The equations used for this return
	 * two points of contact, having traced the line to infinity. Dot product is used to calculate which point of 
	 * contact the cell was actually moving towards. The proposed move is curtailed accordingly to place the cell on 
	 * the circle boundary. This is all scaled back into 3D. 
	 * @param cell
	 * @param move
	 * @return
	 */
	private MoveResults cylinderCollisionDetection(Cell cell, Double3D move)
	{
		Double3D location = cell.getCurrentLocation();
		Double3D proposedLocation = location.add(move);
		// Work only with centroids (calculating point of contact of sphere with cylinder when the sphere is moving
		// across the cylinder's normal vector is complicated and a level of detail we don't need). 
		// If a cell's xy distance from centre is greater than the radius, then a collision has occurred.
		// The bounce vector (normal of contact) points straight into the centre. 
		// Interested only in location w.r.t. x and y axis; z ignored because it's a cylinder. 
		double distance_from_centre = utils.Utils.length(new Double3D(proposedLocation.x-centre.x, 
																	  proposedLocation.y-centre.y, 
																	  0.));
		final boolean breach = distance_from_centre >= BoundedCylinder.radius;
		Double3D bounce = new Double3D(0., 0., 0.);
		Double3D executedMove = move;
		Double3D newLocation = new Double3D(location.x + move.x, location.y + move.y, location.z + move.z);
		if (breach)
		{
			final Double2D centre_xy = new Double2D(centre.x, centre.y);
			final Double2D start = new Double2D(location.x, location.y);
			final Double2D end = new Double2D(start.x + move.x, start.y + move.y);
			
			// Tracing line to infinity, where are the two points of contact?
			// These intersection points can be NaN if start == end (or sufficiently close given precision). 
			final Double2D[] intersections = circleLineIntersectionPoints(radius, centre_xy, start, end);
			final Double2D intersection1 = intersections[0];
			final Double2D intersection2 = intersections[1];
			
			double move_proportion_executed = 0.0;  // Default case.
			if (Double.isFinite(intersection1.getX()))  // Avoid NaN if start ~ end. Precision limitations. 
			{
				// The vectors that would bring the cell to each point of contact. 
				final Double2D to_int1 = intersection1.subtract(start);
				final Double2D to_int2 = intersection2.subtract(start);
				final Double2D move_xy = end.subtract(start); 
				
				// The positive dot product indicates which intersection point the move is pointing towards
				final double dp1 = Utils.dotProduct(to_int1, move_xy);
				final double dp2 = Utils.dotProduct(to_int2, move_xy);
				
				Double2D contact_move = (dp1 >= 0.) ? to_int1 : to_int2;
				if (move_xy.length() > 0.0)  // Avoid div by zero.
					// Can use this to calculate cell resting point in 3D
					move_proportion_executed = collisionBleed * (contact_move.length() / move_xy.length());	
			}
			
			executedMove = new Double3D(move.x * move_proportion_executed, 
										move.y * move_proportion_executed, 
										move.z * move_proportion_executed);
			newLocation = new Double3D(location.x + executedMove.x, 
									   location.y + executedMove.y, 
									   location.z + executedMove.z);
			// Normal of contact with cylinder (circle) is vector of point of contact to centre
			bounce = new Double3D(centre.x - newLocation.x, centre.y-newLocation.y, 0.);
		}
		return new MoveResults(newLocation, bounce, executedMove, new ArrayList<Cell>());
	}
	
	/** 
	 * Cells within the bolus head for the centre and are slowed down. 
	 * This logic is folded within the collision detection framework.  
	 */
	private MoveResults bolusInfiltrationDetection(Cell cell, Double3D move)
	{
		final Double3D presentLocation = cell.getCurrentLocation();
		if (! insideBolus(presentLocation))
		{
			// No alterations made to the proposed move.
			return new MoveResults(null, null, move, null);
		}
		// At this point cell is within the bolus, and will not leave again. 
		Double3D vectorToCentre = vectorToCentre(presentLocation).normalize();
		vectorToCentre = vectorToCentre.resize(bolusSpeed * Simulation.instance.timeSlice_min);
		final MoveResults mr = new MoveResults(null, null, vectorToCentre, null);
		return mr;
	}
	
	private double sgn(double x)
	{
		if (x < 0)
			return -1.;
		else 
			return 1;
	}
	
	/**
	 * Given vector defined by given start and end points, this will return the coordinates of where the vector 
	 * extended to infinity intersects with the line. 
	 * 
	 * There is no relation between direction of vector and the two intersection points, you will need to check both
	 * points to see which one the vector points at and away from.
	 * 
	 * This implementation does not check that the line in fact does intersect with the circle, in that case the 
	 * behaviour is undefined - you must check beforehand that the line does in fact intersect the circle. In the 
	 * present anticipated usage in this simulation the start points (cell locations) should all be inside the circle,
	 * so the assumption holds. See the URL below for use of the determinate in differentiating the possibilities of
	 * {no intersections, tangent, two intersections}.
	 * 
	 * Implementation of this (accessed 2017-10-13):
	 * http://mathworld.wolfram.com/Circle-LineIntersection.html
	 * @param radius
	 * @param centre
	 * @param start
	 * @param end
	 * @return Intersection points can be NaN if start and end are the same (accounting for precision limitations). 
	 */
	private Double2D[] circleLineIntersectionPoints(double radius, Double2D centre, Double2D start, Double2D end)
	{
		// Transform everything such that the centre lies at 0,0. 
		Double2D s = start.subtract(centre);
		Double2D e = end.subtract(centre);
		double dx = e.x - s.x;
		double dy = e.y - s.y;
		double dr = Math.sqrt(dx*dx + dy*dy);  // Length of vector
		double D = s.x * e.y - s.y * e.x;
		
		// Squared values used in calculations
		double r2 = radius * radius;
		double dr2 = dr * dr;
		double D2 = D * D;
		
		// Coordinates of the two intersection points, i1 and i2
		double i1x = (D * dy + sgn(dy) * dx * Math.sqrt(r2 * dr2 - D2)) / dr2;
		double i2x = (D * dy - sgn(dy) * dx * Math.sqrt(r2 * dr2 - D2)) / dr2;
		double i1y = (-D * dx + Math.abs(dy) * Math.sqrt(r2 * dr2 - D2)) / dr2;
		double i2y = (-D * dx - Math.abs(dy) * Math.sqrt(r2 * dr2 - D2)) / dr2;
		
		// Translate out of circle-centre-origin space back into original space
		Double2D intersection1 = new Double2D(centre.x + i1x, centre.y + i1y);
		Double2D intersection2 = new Double2D(centre.x + i2x, centre.y + i2y);
		// Return two intersection points as an array
		return new Double2D[]{intersection1, intersection2};
	}
	
	/**
	 * Populate the space with cells upon simulation instigation, placing each cell randomly in the volume such that 
	 * cells do not overlap. 
	 */
	@Override
	public void placeCellRandomly(Cell cell) 
	{
		Double3D loc;  // Will try to find a location for the neutrophil.
		do
		{
			// select a random location.  
			// cells can be placed in the tissue volume, and similar sized volumes all around it (except 
			// above it, because that breaches the skin).
			final double x = Simulation.instance.random.nextDouble() * radius*2;
			final double y = Simulation.instance.random.nextDouble() * radius*2;
			final double z = Simulation.instance.random.nextDouble() * depth;
			
			loc = new Double3D(x, y, z);
		} while (! isOccupiableSpaceStartup(loc, cell) || 
				(!bolusTCells && insideBolus(loc)) /* can disable exclusion from bolus in instigation */);
		
		if (insideBolus(loc))  // make cells inside the bolus invisible
		{
			((TCell)cell).recordable = false;
			((TCell)cell).motile = false;
		}
		
		this.cellField.setObjectLocation(cell, loc);
		cell.setCurrentLocation(loc);
	}	
	
	@Override
	public boolean insideBolus(Double3D location)
	{		
		if (bolusPresent == false)	
			return false;
		// Interested only in location w.r.t. x and y axis; z ignored because bolus a cylinder. 
		double distance_from_centre = utils.Utils.length(new Double3D(location.x-centre.x, location.y-centre.y, 0.));
		return distance_from_centre < bolusRadius;  
	}
	
	@Override
	public boolean insideImagingVolume(double x, double y, double z) 
	{
		// If cell is not within restricted space, then it is inside the imaging volume. For this bounded cylinder this
		// should always be true; cell's can't escape, and the entire cylinder is imaged. 
		return !restrictedSpace(new Double3D(x, y, z), 0.);
	}
	
	public boolean insideRecordingVolume(double x, double y, double z, Cell cell)
	{
		// Any part of the cell can lie within the recording slice, hence use of radius. 
		return Math.abs(z - CellLogger.sliceZ) <= CellLogger.sliceDepth + cell.getRadius();
	} 
	
	@Override
	public double volumeImaged() 
	{	return volumeImaged;	}

	@Override
	public double volumeSimulated() 
	{	return volumeSimulated;		}
	
	/** 
	 * Concrete classes queried on whether the sphere placed at location with given radius breaches any other objects.
	 * Overridden by subclasses. 
	 */
	@Override
	protected boolean restrictedSpace(Double3D location, double radius)
	{ 
		// Interested only in location w.r.t. x and y axis; z ignored because it's a cylinder. 
		double distance_from_centre = utils.Utils.length(new Double3D(location.x - centre.x, location.y - centre.y, 0.));
		// Is the distance from the centre beyond the radius?
		final boolean xy = distance_from_centre > BoundedCylinder.radius;
		// Is the z coordinate above the top (z=0, so z<0 not permitted) or below the bottom?
		final boolean d = location.z < 0 || depth < location.z;
		// Any condition being true means the location is in restricted space.
		return xy || d;   
	}
		
	public Double3D vectorToCentre(Double3D location)
	{
		return new Double3D(centre.x-location.x, centre.y-location.y, 0);
	}
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;	
		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/radius").evaluate(params, XPathConstants.NODE);
		radius = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/depth").evaluate(params, XPathConstants.NODE);
		depth = Double.parseDouble(n.getTextContent());
		
		extremeWidth = 2* radius;
		extremeHeight = extremeWidth;
		extremeDepth = depth;

		centre = new Double3D(radius, radius, depth/2.);
		
		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/bolus/present")
				.evaluate(params, XPathConstants.NODE);
		if (n != null)
			bolusPresent = Boolean.parseBoolean(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/bolus/bolusTCells")
				.evaluate(params, XPathConstants.NODE);
		if (n != null)
			bolusTCells = Boolean.parseBoolean(n.getTextContent());		
		
		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/bolus/radius")
				.evaluate(params, XPathConstants.NODE);
		if (n != null)
			bolusRadius = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/bolus/bolusSpeed")
				.evaluate(params, XPathConstants.NODE);
		if (n != null)
			bolusSpeed = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/bolus/bolusChemokineSource").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().equals("NONE"))
			bolusChemokineSourceMode = BolusChemokineSourceType.NONE;
		if (n != null && n.getTextContent().equals("EPOCH"))
			bolusChemokineSourceMode = BolusChemokineSourceType.EPOCH;
		if (n != null && n.getTextContent().equals("CONTINUOUS")) 
			bolusChemokineSourceMode = BolusChemokineSourceType.CONTINUOUS;

		n = (Node) xPath.compile("/params/Environment/BoundedCylinder/bolus/bolusChemokineSecretionRate").evaluate(params, XPathConstants.NODE);
		if (n != null && ! "".equals(n.getTextContent()))
			bolusChemokineSecretionRate = Double.parseDouble(n.getTextContent());
	}
	
}