package environment;

import java.util.ArrayList;
import java.util.Vector;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.Cell;
import core.HairFollicle;
import core.SecretionRecord;
import core.Simulation;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.field.continuous.Continuous3D;
import sim.util.Bag;
import sim.util.Double3D;
import sim.util.IntBag;
import soluble.Attractant;
import soluble.AttractantSecretor;
import soluble.HeatKernelLookup;
import soluble.Secretor;
import utils.Utils;

/**
 *
 *  This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Represents a 3D tissue volume.
 * @author Mark N. Read
 *
 */
@SuppressWarnings("serial")
public abstract class Compartment3D implements Steppable
{
	/*
	 * Note that these fields are deliberately public - they have to be for the GUI to
	 * latch onto them for portrayal.
	 */
	public Continuous3D cellField;  // stores cells in a continuous space

	/* Store a list of cells secreting the respective soluble factors, and the time at which the list was compiled.
	 * This is for efficiency purposes, to prevent the list being needlessly recompiled in a single timestep. */
	private SecretorsTimepoint<AttractantSecretor> attractantSecretors =
			new SecretorsTimepoint<AttractantSecretor>(null, Double.NEGATIVE_INFINITY);

	// Storing all hair follicles here in case we need to access them later (e.g., the GUI does, and collision
	// detection may also). Here only because the GUI needs it.
	public ArrayList<HairFollicle> follicles = new ArrayList<HairFollicle>();

	/* Needed for visualization, not elegant I know, but necessary given how MASON setup works.
	 * Values are set when reading in parameters from concrete classes.
	 * These represent extreme dimensions (not necessarily cubes). */
	public static double extremeWidth;
	public static double extremeHeight;
	public static double extremeDepth;

	// Can disable inter-cellular collision detection.
	protected static boolean intercellularCollisionDetection = false;
	/* This shortens the length of a move that would place a cell next to something it has collided with.
	 * Was intended as a safeguard against mathematical precision errors, but not clear that it is really needed.
	 * Set to near but less than 1. */
	public static double collisionBleed = 1.0;

	protected HeatKernelLookup attractantLookup; 			
	
	protected static enum SecretionMotilityMode {
		// Secretors can move, means the Heat Kernel must be re-calculated for each location at which secretion
		// took place. 
		MOTILE,  
		// Secretors are assumed to have performed all secretion events at the first spatial location where 
		// secretion occurred. 
		STATIONARY
	}
	protected static SecretionMotilityMode secretionMotilityMode = SecretionMotilityMode.MOTILE;
	
	@Override
	public void step(SimState arg0)
	{}

	/**
	 * Populate the space with cells, placing each cell randomly in the volume such that cells do not overlap.
	 */
	public abstract void placeCellRandomly(Cell cell);

	public abstract boolean insideImagingVolume(double x, double y, double z);
	public abstract boolean insideRecordingVolume(double x, double y, double z, Cell cell);

	public abstract double volumeImaged();
	public abstract double volumeSimulated();

	/** False by default, overridden by classes simulating bolus. */
	public boolean insideBolus(Double3D location)	{ return false;	}

	/** Returns true if the supplied location resides within a cell. */
	public boolean insideCell(Double3D location)
	{
		for (Object obj : cellField.getAllObjects())
		{
			final Cell c = (Cell)obj;
			final Double3D nLoc = c.getCurrentLocation();
			if (Utils.displacement(nLoc, location) < c.getRadius())
				return true;
		}
		return false;
	}

	/**
	 * Checks whether the particular cell would collide with any others, were it to be placed at the supplied location.
	 *
	 * @param cellLoc
	 * @param cell
	 * @return	True if there is a collision.
	 */
	protected boolean cellularCollisions(Double3D cellLoc, Cell cell)
	{
		// This is a conservative estimate, true as long as the biggest cell's diameter isn't more than
		// twice the diameter as the smallest.
		double distance =  2 * cell.getDiameter();
		boolean collision = false;  // The default case, we will look for the exception below.
		Bag neighbours = cellField.getNeighborsExactlyWithinDistance(cellLoc, distance, false, true, true, null);
		// remove `cell` if it is already placed in the tissue space (clearly it will collide with itself).
		neighbours.remove(cell);
		for (Object other : neighbours)
		{
			double minDistance = ((Cell)other).getRadius() + cell.getRadius();
			// minimum distance at which there is no collision.
			Double3D otherLoc = cellField.getObjectLocation(other);
			// actual distance between the two cells.
			distance = Math.sqrt( ((cellLoc.x - otherLoc.x) * (cellLoc.x - otherLoc.x))
								+ ((cellLoc.y - otherLoc.y) * (cellLoc.y - otherLoc.y))
								+ ((cellLoc.z - otherLoc.z) * (cellLoc.z - otherLoc.z)) );
			if (distance <= minDistance)
			{
				collision = true;
				break;
			}
		}
		return collision;
	}

	public MoveResults moveCellCollisionDetection(Cell cell, Double3D move)
	{
		// Check for collisions with other cells first. A cell will generally encounter other cell (if it is to do
		// so at all before it encounters the barrier (because all other cells are within the boundaries).
		MoveResults cell_mr = cellularCollisionDetection(cell, move);
		// Check for collisions with spatial boundary. The move may already have been curtailed by cellular collision
		// detection. This can curtail it further if needed.		
		MoveResults boundary_mr = boundaryCollisionDetection(cell, cell_mr.executedMove);
		cell_mr.colliders.addAll(boundary_mr.colliders);
		
		// One last sanity check.		
		if (! restrictedSpace(boundary_mr.newLocation, cell.getRadius()))
		{
			cellField.setObjectLocation(cell, boundary_mr.newLocation);		
			return new MoveResults(boundary_mr.newLocation,
								   cell_mr.bounce.add(boundary_mr.bounce),  // add the bounces
								   boundary_mr.executedMove,
								   cell_mr.colliders);
		} else {
			/* If, for whatever reason, the new move would place the cell in non-occupiable space, don't move.
			 * This can happen using the Ballistic motility mode, for reasons I've not been able to pin down. I think
			 * it has to do with precision error near hard boundaries. Ballistic cells will not change their
			 * orientation on their own fruition, and a stuck cell will remain such. This is not a problem if cells
			 * are performing a random walk, as far as I can tell. Returning a random bounce vector here breaks the
			 * problem in Ballistic mode. */
			return new MoveResults(cell.getCurrentLocation(),  // New location
					new Double3D(Simulation.instance.random.nextDouble(),  // Bounce
							     Simulation.instance.random.nextDouble(), 
							     Simulation.instance.random.nextDouble()),
					new Double3D(0,0,0), new ArrayList<Cell>());  // Executed move, list of colliders
		}
	}

	/**
	 * Moves a cell in space, providing collision avoidance with other cells. A cell must be
	 * specified, and its preferred movement (movement, not preferred location!). The preferred movement is used to
	 * calculate the eventual movement, but is adjusted to ensure that no spatial co-occupancies occur.
	 *
	 * Returns an array containing two Double3Ds.
	 * The first is the new location is returned, as it may be useful to store this locally.
	 * The second is a 'bounce' vector. This is a sum of the normals of contacts that 'cell' makes with neighbours.
	 * The bounce vector can be used to allow cells to slide over one another, if it is summed with the next movement
	 * vector. Bounce vector is of zero length if no contact occurred.
	 *
	 * Note that his algorithm does not work so well for what seems like simultaneous collisions with two things. The
	 * reason is that with cells that don't squash, at this level of precision, one object is always hit before another,
	 * and hence the cell will make very small bounces between the two, rather than hitting both simultaneously.
	 */
	public MoveResults cellularCollisionDetection(Cell cell, Double3D move)
	{
		/* This algorithm is based on that found in
		 * http://www.gamasutra.com/view/feature/131424/pool_hall_lessons_fast_accurate_.php?page=2
		 * Accessed on the 13th June 2014.
		 *
		 * Some notation for what is happening here. The 'cell's current location is labelled C.
		 * Cell is to make a move, initially proposed to be 'move', but shortened to avoid collisions as needed.
		 * The location that the move will take C to is S.
		 *
		 * 'Other's current location is labelled O.
		 * Hence, vector CO is the distance from the centre of 'cell', and the centre of 'other'.
		 * Vector CS is the distance from the centre of 'cell' to its resting place at the end of the move, S.
		 *
		 * In calculating if there is a collision between 'cell' and 'other', the point along CS closest to O needs to
		 * be calculated. This point is R.
		 *
		 * Vector RO is the distance between the centre of 'other', and where 'cell'  will be closest to O.
		 * Point T lies along CS. It is the nearest point along CR where 'cell' touches 'other'. This is where 'cell'
		 * should stop to avoid a collision. T is only calculated if there will be a collision.
		 * Vectors ending in 2 indicate that they are the squared distance. This is done to avoid taking square roots,
		 * which are computationally expensive.
		 */
		double bouncex = 0.0;
		double bouncey = 0.0;
		double bouncez = 0.0;
		double CSx = move.x;
		double CSy = move.y;
		double CSz = move.z;
		double move_length = move.length();
		final double cell_radius = cell.getRadius();


		Double3D C = cellField.getObjectLocation(cell);
		// Record the cells with which there was an actual collision.
		ArrayList<Cell> colliders = new ArrayList<Cell>();		
		
		if (intercellularCollisionDetection)
		{
			// Find cells other cells that the length of 'move' could bring 'cell' to collide with (omnidirectional).
//			Bag potentialColliders = cellField.allObjects;
			// Note that if cells vary considerably in their radius, then this should be altered. 
			final double potential_collision_length = move_length + (3 * cell_radius);
			Bag potentialColliders = cellField.getNeighborsExactlyWithinDistance(C, potential_collision_length);
			potentialColliders.remove(cell);
			for (Object o : potentialColliders)
			{
				Cell other = (Cell)o;
//				if (!(o == cell))  // Cell can't collide with itself.
//				{
				Double3D O = cellField.getObjectLocation(other);
				// Dealing with hair follicles. This effectively places a sphere representing the follicle at the exact
				// same z-coordinate as the 'cell'.
				if (other instanceof HairFollicle)
					O = new Double3D(O.x, O.y, C.z);
				/* Check if current movement will bring 'cell' within touching range of 'other', ignoring direction. */
				double COx = O.x - C.x;
				double COy = O.y - C.y;
				double COz = O.z - C.z;
				double touchRadii = cell.getRadius() + other.getRadius();
				double touchRadii2 = touchRadii*touchRadii;
				double CO2 = COx*COx + COy*COy + COz*COz;  // CO squared.
				double CO = Math.sqrt( COx*COx + COy*COy + COz*COz );
				double CS = Math.sqrt( CSx*CSx + CSy*CSy + CSz*CSz );

				// If true, 'other' is within touching distance of 'cell'.
				// This is simply whether they are within range, not accounting for direction.
				if(CS > (CO - touchRadii))
				{
					/* Check if the movement vector is in the direction of the 'other' cell. Ignore if it is not. (it
					 * might be moving away from 'other', rather than towards it). */
					double dotP = COx*CSx + COy*CSy + COz*CSz;	// dot product.
					if (dotP > 0.0)		// if true, 'cell' is moving towards 'other'.
					{
						/* Check whether the CS vector will actually bring 'cell' within touching distance of 'other'.
						 * Is it here that point R must be calculated. */
						double CSunitx = CSx / CS;
						double CSunity = CSy / CS;
						double CSunitz = CSz / CS;
						// The geometry here is a right angle triangle, between C, O and R. The rightangle is CRO.
						// we have three equations - note that dp means dot product:
						// 1. cos(a) = adjacent / hypotenuse. Angle a is RCO.
						// 2. dp(v,w) = |v|*|w|*cos(a).
						// 3. dp(v,w) = vx*wx + vy*wy + vz*wz
						//
						// combining 1 and 2 gives:
						// 4. dp(v,w) = |v|*|w|* (|adj| / |hyp|)
						//
						// plug some variables into 4. Construct a right angle triangle. Hyp = CO, Adj = CR (what we are
						// trying to find). Opposite = RO, though we don't use it in this calculation.
						// dp(CS,CO) = |CS| * |CO| * |CR| / |CO|      the COs cancel each other out
						// dp(CS,CO) = |CS| * |CR|
						//
						// We can use the unit vector of CS, which has the same direction as CS, but length 1. Hence,
						// dp(CSunit,CO) = |CR|
						//
						// Hence, dot product of CS unit and CO gives the distance along 'CS' closest to O (i.e., CR).
						double CR = COx*CSunitx + COy*CSunity + COz*CSunitz; // this is the dot product, gives opposite.
						double CR2 = CR*CR;									 // gives opposite squared.
						double RO2 = CO2 - CR2;
						// if true, then cells get within contact range of one another.
						if (RO2 <= touchRadii2) // square minSeparation to avoid square roots.
						{
							/* The cells will colide, here we find point T, and vector CT. */
							colliders.add(other);
							// This is also done using a right angle triangle, between O, R and T. Right angle is ORT.
							// hyp = TO. The other two sides are TR and OR. We wish to find TR.
							double TO2 = touchRadii2;
							double TR2 = TO2 - RO2;
							double TR = Math.sqrt(TR2);
							double CT = CR - TR;
							double shrinkFactor = CT / CS;
							double CTx = CSx * shrinkFactor;
							double CTy = CSy * shrinkFactor;
							double CTz = CSz * shrinkFactor;
							// Set point S to point T, ready for consideration of other potentially colliding cells.
							CSx = CTx; 	CSy = CTy;  CSz = CTz;
							// Calculate vector OT, since this is the normal to the contact between 'cell' and 'other'.
							// This can be used to create a 'bounce' vector, which can be used to have cells slide over
							// one another.
							double Tx = C.x+CTx;	// these are coordinates, not magnitudes.
							double Ty = C.y+CTy;
							double Tz = C.z+CTz;
							// These are directional vector components.
							// They constitute the normal of the contact between 'cell' and 'other'.
							// Unit vectors so that mutliple bounces count equally.
							double bx = Tx-O.x;
							double by = Ty-O.y;
							double bz = Tz-O.z;
							double blen = Math.sqrt(bx*bx + by*by + bz*bz);
							bouncex += bx / blen;
							bouncey += by / blen;
							bouncez += bz / blen;
						}
					}
				}
//				}
			}
		}
		Double3D newLoc = new Double3D(C.x + CSx, C.y + CSy, C.z + CSz);
		Double3D executedMove = new Double3D(CSx, CSy, CSz);
		Double3D bounce = new Double3D(bouncex,bouncey,bouncez);
		if (bounce.lengthSq() > 0.0)
			bounce = Utils.unitVector(bounce);
		return new MoveResults(newLoc, bounce, executedMove, colliders);  // return information.
	}

	/**
	 * Collision with boundary detection for cells wanting to move past z=0 (defined as top in this geometry)
	 *
	 * This is a utility function that may be used in concrete classes.
	 */
	public MoveResults boundaryCollisionDetectionTop(Cell cell, Double3D move)
	{
		Double3D bounce = new Double3D(0, 0, 0);  // default case
		Double3D executedMove = move;  // default case
		/* Because the top (and bottom) boundaries lie perpendicular to the z-axis, we can simply check for whether the
		 * proposed move would breach these z-coordinates, and adjust the magnitude of the move vector accordingly,
		 * such that the cell touches but does not breach the boundary. */
		final Double3D cl = cell.getCurrentLocation();
		Double3D proposed_location = new Double3D(cl.x + move.x, cl.y + move.y, cl.z + move.z);
		if (proposed_location.z <= 0.)
		{
			// Cell centroid breaches the boundary.
			/* Because boundary = z = 0, these absolute numbers can be used. Move.z is negative, as it is taking a cell
			 * from permissible space (positive z values) into restricted space (negative z values) */
			double move_dz = move.z;
			double allowable_dz = -(cl.z);  // Hence, must make this negative
			double allowable_move_proportion = 0;  // Default case if move_dz = 0
			if (move_dz > 0)  // Avoid div by zero				
				allowable_move_proportion = collisionBleed * (allowable_dz / move_dz);  // This will be positive
			
			executedMove = new Double3D( move.x * allowable_move_proportion,
										 move.y * allowable_move_proportion,
										 move.z * allowable_move_proportion);
			proposed_location = new Double3D(cl.x + executedMove.x,
											 cl.y + executedMove.y,
											 cl.z + executedMove.z);
			// Calculate the bounce vector
			bounce = new Double3D(0., 0., 1); // Unit length. Points into permissible space
		}
		return new MoveResults(proposed_location, bounce, executedMove, new ArrayList<Cell>());
	}

	/**
	 * Collision with boundary detection for cells wanting to move past z=defined (defined as bottom in this geometry).
	 *
	 * This is a utility function that may be used in concrete classes.
	 */
	public MoveResults boundaryCollisionDetectionBottom(Cell cell, Double3D move, double zBoundary)
	{
		Double3D bounce = new Double3D(0, 0, 0);  // Default case
		Double3D executedMove = move;  // default case
		/* Because the top (and bottom) boundaries lie perpendicular to the z-axis, we can simply check for whether the
		 * proposed move would breach these z-coordinates, and adjust the magnitude of the move vector accordingly,
		 * such that the cell touches but does not breach the boundary. */
		final Double3D cl = cell.getCurrentLocation();
		Double3D proposed_location = new Double3D(cl.x + move.x, cl.y + move.y, cl.z + move.z);
		if (proposed_location.z >= zBoundary)
		{
			// Cell centroid breaches the boundary.
			double move_dz = move.z;  // Positive if it is breaching the boundary (z = height)
			double allowable_dz = zBoundary - cl.z;  // Hence, must make this negative
			double allowable_move_proportion = 0;  // Default case if move_dz = 0. 
			if (move_dz > 0.0)  // Avoid div by zero.
				allowable_move_proportion = collisionBleed * (allowable_dz / move_dz);  // This will be positive
			
			executedMove = new Double3D(move.x * allowable_move_proportion,
										move.y * allowable_move_proportion,
										move.z * allowable_move_proportion);
			proposed_location = new Double3D(cl.x + executedMove.x,
											 cl.y + executedMove.y,
											 cl.z + executedMove.z);
			// Calculate the bounce vector
			bounce = new Double3D(0., 0., -1); // Unit length. Points into permissible space
		}		
		return new MoveResults(proposed_location, bounce, executedMove, new ArrayList<Cell>());
	}

	public abstract MoveResults boundaryCollisionDetection(Cell cell, Double3D move);

	/**
	 * Used solely for grouping together information passed back to cell following an attempted move in space
	 * The new location is returned, as it may be useful to store this locally.
	 * The 'bounce' vector is a sum of the normals of contacts that 'cell' makes with neighbours.
	 * The bounce vector can be used to allow cells to slide over one another, if it is summed with the next movement
	 * vector. Bounce vector is of zero length if no contact occurred.
	 * 'Colliders' represents objects (usually other cells) that have been collided with.
	 */
	public static class MoveResults
	{
		public Double3D newLocation;
		public Double3D executedMove;
		public Double3D bounce;
		public ArrayList<Cell> colliders;

		public MoveResults(Double3D newLoc, Double3D bounce, Double3D executedMove, ArrayList<Cell> colliders)
		{
			this.newLocation = newLoc;
			this.executedMove = executedMove;
			this.bounce = bounce;
			this.colliders = colliders;
			assert Double.isFinite(this.newLocation.x);  // Check nothing went wrong.
		}
	}

	/**
	 * Concrete classes queried on whether the sphere placed at location with given radius breaches any other
	 * environmental objects.
	 * Currently used only at startup. Collision detection is used during runtime.
	 * Overridden by subclasses.
	 */
	protected abstract boolean restrictedSpace(Double3D location, double radius);

	/**
	 * Returns true if the specified location (in continuous space) can be occupied be supplied cell.
	 * Currently used only at startup.
	 * Collision detection is used during runtime.
	 */
	public boolean isOccupiableSpaceStartup(Double3D location, Cell cell)
	{
		boolean collision = cellularCollisions(location, cell);
		boolean restricted = restrictedSpace(location, cell.getRadius());	// Subclasses indicate restricted space.

		// Location cannot be occupied if it is restricted, or if it will result in a collision.
		return ! (restricted || collision);
	}

	/**
	 * Will return all the locations within a 2D discrete grid that lie within `distance` of the coordinate
	 * (originX, originY).
	 * The grid has dimensions (0,0) to (fieldWidth, fieldHeight).
	 * A location is considered within the circle if its center is within distance of the specified point.
	 * `xPos` and `yPos` must be supplied, and the locations are returned in these positions.
	 * The origin is always included.
	 *
	 * This method was necessary because the MASON methods in 3D don't work properly. The MASON code for these
	 * methods looks pretty hacked, so I have created my own.
	 */
	public static void radial2DLocations(int originX, int originY, int fieldWidth, int fieldHeight, double distance,
										IntBag xPos, IntBag yPos)
	{
		if(xPos == null)	{ xPos = new IntBag();	}
		if(yPos == null)	{ yPos = new IntBag();	}
		xPos.clear();
		yPos.clear();

		// All the locations to be returned must lie within a square, centered on (originX, originY).
		// This defines the boundaries of that square.
		int squareWMin = (int) Math.round(originX - distance);
		int squareWMax = (int) Math.round(originX + distance);
		int squareHMin = (int) Math.round(originY - distance);
		int squareHMax = (int) Math.round(originY + distance);
		// Ensure the boundaries do not lie outside of the grid.
		if (squareWMin < 0) 			{ squareWMin = 0; }
		if (squareWMax >= fieldWidth) 	{ squareWMax = fieldWidth - 1;  }
		if (squareHMin < 0) 			{ squareHMin = 0; }
		if (squareHMax >= fieldHeight) 	{ squareHMax = fieldHeight - 1; }

		// Go through each location in the grid, calculate its distance from the specified origin,
		// and stores those within range.
		for (int x = squareWMin; x <= squareWMax; x++)
			for (int y = squareHMin; y <= squareHMax; y++)
			{
				double hyp = Math.sqrt( ((x - originX) * (x - originX))
									  + ((y - originY) * (y - originY)));
				if(hyp <= distance)
				{
					xPos.add(x);
					yPos.add(y);
				}
			}
	}

	/** Performs the move operation for a cell. Assumes checks on validity of move have already taken place */
	public void moveCell(Cell cell, Double3D newLocation)
	{	cellField.setObjectLocation(cell, newLocation);	}

	public Double3D getCellLocation(Cell cell)
	{	return cellField.getObjectLocation(cell);	}

	/**
	 * Returns a Bag of those cells that lie within a radius of the specified cell.
	 * The specified cell is not included in the Bag.
	 */
	public Bag cellsInLocale(Cell cell, double radius)
	{
		Double3D loc = cellField.getObjectLocation(cell);
		Bag neighbours = cellField.getNeighborsExactlyWithinDistance(loc, radius);
		neighbours.remove(cell);
		return neighbours;
	}

	/** Returns a vector containing all the entities currently secreting attractant.
	 */
	protected AttractantSecretor[] findAttractantSecretors()
	{
		// This can be called many times in a single time step (once per cell performing query).
		// Hence, recompile list only when the time advances.
		if (attractantSecretors.expired())
		{
			ArrayList<AttractantSecretor> secretors = new ArrayList<AttractantSecretor>();

			for (Object o : cellField.allObjects)
			{
				// Check whether there is any record of having secreted
				if (o instanceof AttractantSecretor && ((AttractantSecretor)o).getSecretionHistory().size() > 0)
					secretors.add((AttractantSecretor)o);
			}
			double timeNow = Simulation.instance.schedule.getTime();
			AttractantSecretor[] array = new AttractantSecretor[secretors.size()];
			attractantSecretors = new SecretorsTimepoint<AttractantSecretor>(secretors.toArray(array),
																			 timeNow);
		}
		return attractantSecretors.getSecretors();
	}

	public double perceiveAttractantAtPoint(Double3D location)
	{
		switch (secretionMotilityMode)
		{
		case MOTILE: 
			return perceiveAttractantAtPointMotileSecretors(location);
		default:
			return perceiveAttractionAtPointStationarySecretors(location);
		}
	}
	
	/**
	 * Returns the concentration of attractant at the supplied `location` based on the PDE heat equation.
	 *
	 * The equation is too complex to write here in the comments, see the accompanying simulation documentation for
	 * specifics.
	 */
	public double perceiveAttractantAtPointMotileSecretors(Double3D location)
	{
		// Get all the attractant-secreting cells
		Secretor[] secretors = findAttractantSecretors();
		double concentration = 0.0;  // Cumulative calculation.
		for (Secretor s : secretors)
		{
			Vector<SecretionRecord> history = s.getSecretionHistory();
			for (SecretionRecord r : history)
			{
				double dt = Simulation.instance.schedule.getTime() - r.getTime();  // Time since secretion event
				// Safety. Updating cells is not an atomic operation.
				// Some may have recorded secreting at this timestep before another evaluates the gradient.
				// This ensures that only secretion records from the past are considered.
				if (dt != 0.0)
				{
					double distance2 = r.getLocation().subtract(location).lengthSq();
					double quantity = r.getQuantity() * attractantLookup.get(distance2, dt);

					// Apply decay, if relevant.
					assert Double.isFinite(Attractant.decayConstant);
					if (Attractant.decayConstant != 0)
						quantity *= Math.exp(-Attractant.decayConstant * dt);

					concentration += quantity;
				}
			}
		}
		return concentration;
	}
	
	/**
	 * Returns the concentration of attractant at the supplied `location`. 
	 * Makes several assumptions:
	 * 1. Secretors don't move. All secretions are calculated as having occurred at the first secretion event. 
	 * 2. Secretion rates are constant, that of the first secretion event.
	 * 3. Attractant does not decay. 
	 * @param location
	 * @return
	 */
	public double perceiveAttractionAtPointStationarySecretors(Double3D location)
	{
		// Get all the attractant-secreting cells
		Secretor[] secretors = findAttractantSecretors();
		double concentration = 0.0;  // Cumulative calculation.
		for (Secretor s : secretors)
		{	
			Vector<SecretionRecord> history = s.getSecretionHistory();
			final SecretionRecord firstSecretion = history.firstElement();
			final SecretionRecord lastSecretion = history.lastElement();
			final double secretionDuration = lastSecretion.getTime() - firstSecretion.getTime();
			
			// Safety. Updating cells is not an atomic operation.
			// Some may have recorded secreting at this timestep before another evaluates the gradient.
			// This ensures that only secretion records from the past are considered.
			if (secretionDuration > 0.)
			{
				final double distance2 = lastSecretion.getLocation().subtract(location).lengthSq();
				double quantity = firstSecretion.getQuantity() 
								* attractantLookup.get(distance2, secretionDuration);
				concentration += quantity;
			}
		}
		return concentration;
	}

	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;
		n = (Node) xPath.compile("/params/Environment/intercellularCollisionDetection")
				.evaluate(params, XPathConstants.NODE);
		if (n != null)
			intercellularCollisionDetection = Boolean.parseBoolean(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Environment/secretorMotilityMode")
				.evaluate(params, XPathConstants.NODE);
		if (n != null)
			if (n.getTextContent().equals("Stationary"))
				secretionMotilityMode = SecretionMotilityMode.STATIONARY;
			else 
				secretionMotilityMode = SecretionMotilityMode.MOTILE;
			
	}

	/**
	 * Encapsulates a list of cells and a time point at which the list was compiled.
	 *
	 * @author Mark N. Read, 2017
	 *
	 * @param <T> Type of the cells to be stored.
	 */
	protected class SecretorsTimepoint<T>
	{
		private T[] secretors;
		private double time;

		public SecretorsTimepoint(T[] secretors, double t)
		{
			this.secretors = secretors;
			this.time = t;
		}

		/** Returns true if the current list of secretors is no longer valid, because time has elapsed since it was
		 *  created.
		 */
		public boolean expired()
		{	return this.time != Simulation.instance.schedule.getTime();		}

		public T[] getSecretors()
		{	return secretors;		}
	}
}
