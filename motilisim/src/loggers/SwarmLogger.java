package loggers;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.Cell;
import core.MigratoryCell;
import core.Neutrophil;
import core.Simulation;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.Bag;
import sim.util.Double3D;
import utils.Utils;

/**
 * Encapuslates information pertaining to tracking of cells at a swarm level.
 * 
 * If the GUI is to colour by swarm membership, then this must be computed every step.  
 * What is written to file system is currently only the last step.  
 * 
 * Includes the criteria for what constitutes a swarm. 
 * At its minimum, a swarm must contain at least one cell with the given number of neighbours within the
 * given distance. 
 *  
 * There are helper wrappers for cells, which wrap both a cell and the connections it has to other cells; 
 * a connection is a neighbouring cell that lies within a threshold distance. 
 * There is also a swarm wrapper, which groups together cells into swarms and provides statistics thereon. 
 * 
 * @author Mark N. Read, 2018
 *
 */
@SuppressWarnings("serial")
public class SwarmLogger implements Steppable
{
	/**
	 * A simple wrapper for a cell and the connections it has to nearby neighbouring cells. 
	 */
	private static class CellConnectionWrapper 
	{
		public MigratoryCell cell;
		public LinkedList<CellConnectionWrapper> connections = new LinkedList<CellConnectionWrapper>();
		public LinkedList<Double> distances = new LinkedList<Double>();
		
		public CellConnectionWrapper(MigratoryCell cell)
		{  this.cell = cell;  }
		
		public void addConnection(CellConnectionWrapper cell, double distance)
		{
			connections.add(cell);
			distances.add(distance);
		}
		
		public int numberConnections()
		{	return connections.size();	}		
	}
	
	/**
	 * Convenience class that groups together cells as a swarm, and contains information about the swarm. 
	 */
	private static class Swarm
	{
		private LinkedList<CellConnectionWrapper> constituents;
		private double diameter = -1.0;  // Initial value is nonsensical. Calculated upon first request.		
		private Double3D centre = null;
		
		public Swarm(LinkedList<CellConnectionWrapper> members)
		{	this.constituents = members;	}
		
		public int getNumberConstituents()
		{	return constituents.size();		}
		
		public LinkedList<CellConnectionWrapper> getConstituents()
		{	return (LinkedList<CellConnectionWrapper>) constituents.clone();	}
		
		/**
		 * Returns the largest diameter of the swarm. 
		 * This is the biggest distance between two cells' centroids in the swarm. 
		 * Calculated once, upon request. 
		 */
		public double getDiameter()
		{
			if (diameter != -1.0)
				return diameter;
			// For each cell, find the distance to every other cell in the swarm. Store the biggest. 
			for (int i = 0; i < constituents.size(); i++)
			{
				Double3D cell1Loc = constituents.get(i).cell.getCurrentLocation();
				for (int j = i; j < constituents.size(); j++)
				{
					Double3D cell2Loc = constituents.get(j).cell.getCurrentLocation();
					double distance = Utils.displacement(cell1Loc, cell2Loc);
					
					if (distance > diameter)
						diameter = distance;
				}
			}
			return diameter;				
		}
		
		/**
		 * Calculates the centre of the swarm, based on the median position of its constituents. 
		 */
		public Double3D getSwarmCentre()
		{
			if (centre != null)
				return centre;
			double[] x = new double[constituents.size()];
			double[] y = new double[constituents.size()];
			double[] z = new double[constituents.size()];
			// Scan through constituents, getting median value
			for (int i = 0; i < constituents.size(); i++)	
			{
				Double3D loc = constituents.get(i).cell.getCurrentLocation();
				x[i] = loc.x;
				y[i] = loc.y;
				z[i] = loc.z;
			}
			double centreX = Utils.median(x);
			double centreY = Utils.median(y);
			double centreZ = Utils.median(z);
			
			return new Double3D(centreX, centreY, centreZ);
		}
	}

	
	// Store reference to all the simulation cells here. Cells are not created or removed during simulation.
	private LinkedList<MigratoryCell> cells = null;
	// Stores all the cells that are currently part of any swarm. 
	// Not currently used, but could facilitate future analysis. 
	private LinkedList<CellConnectionWrapper> constituents = null; 
	private LinkedList<Swarm> swarms = null;
	
	// Requirements for what constitutes a swarm. 
	// At its minimum, a swarm must contain at least one cell with the given number of neighbours within the
	// given distance. 
	private static int swarmConstituentNumberCriteria = 0;  // Non-optional parameter.
	private static double swarmConstituentProximityCriteria = 0;  // Non-optional parameter.
	
	public void step(SimState state)
	{
		// Collect together all the cells in the simulation that are currently recordable. 
		Bag all = Simulation.space.cellField.getAllObjects();
		cells = new LinkedList<MigratoryCell>();
		for (Object o : all)
			if (o instanceof MigratoryCell)
				if (CellLogger.trackable((MigratoryCell)o))
					cells.add((MigratoryCell)o);

		// Identify which neutrophils comprise swarms.
		// Cell in a swarm if it has e.g. 4 neighbours within 1.2 cell diameters range. 
		identifySwarms(swarmConstituentNumberCriteria, 
				       swarmConstituentProximityCriteria*Neutrophil.diameter);	
		
		// Probably more efficient to set them all to false and then reset those that are in a swarm to 
		// true, than to do complete arithmetic involving set difference. 
		for (Cell n : cells)
			if (n instanceof Neutrophil)
				((Neutrophil) n).setInSwarm(false);
		for (CellConnectionWrapper n : constituents)
			if (n.cell instanceof Neutrophil) 
				((Neutrophil) n.cell).setInSwarm(true);	
	}
	
	/**
	 * Based on the supplied thresholds, this method identifies those cells that constitute 'core' members of 
	 * the swarm, and from there the entire swarm. 
	 * Multiple swarms can be identified. 
	 * 
	 * A core swarm member has more than `minNum` cells within `threshold` distance. 
	 * This is calculated based on a cell's location, which is its centre. 
	 * @param minNum
	 * @param threshold
	 */
	private void identifySwarms(int minNum, double threshold)
	{
		// Core cells are swarm members with above-threshold number of neighbours within threshold distance. 
		// This is a subset of overall swarm constituents, which includes all cells connected to core cells.
		// Hence, you can never have a swarm of size 1. 
		// At minimum it must be 1 core, and its connections. 
		LinkedList<CellConnectionWrapper> core = new LinkedList<CellConnectionWrapper>();
		swarms = new LinkedList<Swarm>();
		// Wrappers will be created for all cells, and stored here. 
		LinkedList<CellConnectionWrapper> allCCW = new LinkedList<CellConnectionWrapper>();
		// Create a connection wrapper for every cell.
		for (MigratoryCell c : cells) 
			allCCW.add(new CellConnectionWrapper(c));
		// Calculate and store the distances between all supplied cells.
		for (int i=0 ; i<cells.size(); i++)		
		{
			final Double3D ci = cells.get(i).getCurrentLocation();
			final CellConnectionWrapper ccw_i = allCCW.get(i);  // Convenience, repeatedly referred to.
			
			for (int j=i+1; j<cells.size(); j++)  
			{
				final Double3D cj = cells.get(j).getCurrentLocation();
				final double dist = Utils.displacement(ci, cj);
				
				// Log connection if it is within the threshold distance 
				if (dist <= threshold) 	
				{
					ccw_i.addConnection(allCCW.get(j), dist);
					allCCW.get(j).addConnection(allCCW.get(i), dist);  // Relationships are bidirectional.
				}
			}
			// If cell 'i' has sufficient connections, log it as a swarm constituent. 
			if (ccw_i.numberConnections() >= minNum)
				core.add(ccw_i);
		}

		// Reset constituents and will re-compile them now. 		
		constituents = new LinkedList<CellConnectionWrapper>();
		// Starting with the core (highly connected) cells, assign cells into swarms.  
		while(core.size() > 0)  // Every 'core' cell must be assigned into a swarm. Iteratively removed.
		{
			LinkedList<CellConnectionWrapper> swarm = new LinkedList<CellConnectionWrapper>();
			// Assigns 'core' and their connected cells into swarms. 
			recursionSwarmIdentifier(swarm, core.pop());
			core.removeAll(swarm);
			constituents.addAll(swarm);		
			swarms.add(new Swarm(swarm));
		}	
	}
	
	/**
	 * Recursive method, which given a seed cell `n`, follows all its connections and dynamically constructs 
	 * a list of cells comprising the swarm. 
	 * Method ensures that no cell is counted twice. 
	 * 
	 * This method counts as a swarm any cell that is within `threshold` distance of a `constituent` cell. 
	 * Hence, swarms of size 1 cannot exist. 
	 */
	private void recursionSwarmIdentifier(LinkedList<CellConnectionWrapper> swarm, 
										  CellConnectionWrapper seedCell)
	{
		swarm.add(seedCell);		
		// Recursion for the connected cells that have not yet been allocated.
		for (CellConnectionWrapper cell : seedCell.connections)		
			if (! swarm.contains(cell))  // Don't revisit cells already allocated to this swarm.
				recursionSwarmIdentifier(swarm, cell);
	}
	
	/**
	 * Writes statistics on swarms at the current time in the simulation to the file system. 
	 * Files will be created under the specified directory. 
	 */
	public void writeSwarmData(String dir)
	{	
		System.out.println("Writing swarm data to the filesystem.");
		final double currentTime = Simulation.instance.schedule.getTime();
		try 
		{
			BufferedWriter centreOut = setupSwarmCentrePositionFile(dir + "/_SwarmCentrePosition.csv");
			BufferedWriter diameterOut = setupSwarmDiameterFile(dir + "/_SwarmDiameter.csv");			
			for (int id = 0; id < swarms.size(); id++)
			{
				Swarm s = swarms.get(id);
				Double3D loc = s.getSwarmCentre();
				centreOut.write(loc.x + ", " + loc.y + ", " + loc.z + ", " + currentTime + ", " + id + "\n");
				diameterOut.write(s.getDiameter() + ", " + currentTime + ", " + id + "\n");
			}			
			centreOut.close();
			diameterOut.close();
		} catch (IOException ex) {
			System.out.println("ERROR: exception when writing to filesystem, " + ex.toString());
		}
		System.out.println("Completed writing swarm data to the filesystem.");
	}
	
	/*
	 * Sets up a file for writing data on the swarms' centre positions to the filesystem.  
	 * Note that the returned BufferedWriter must be closed elsewhere, it is not performed here. 
	 */
	private BufferedWriter setupSwarmCentrePositionFile(String fileName) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter out = new BufferedWriter(fw);
		out.write("Position_X, Position_Y, Position_Z, Time, ID\n");
		return out;
	}	

	/*
	 * Sets up a file for writing data on the swarms' diameter to the filesystem. 
	 * Note that the returned BufferedWriter must be closed elsewhere, it is not performed here. 
	 */
	private BufferedWriter setupSwarmDiameterFile(String fileName) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter out = new BufferedWriter(fw);
		out.write("Diameter, Time, ID\n");
		return out;
	}	
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{	
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		// These are required parameters. 
		n = (Node) xPath.compile("/params/Swarm/swarmConstituentNumberCriteria")
				.evaluate(params, XPathConstants.NODE);				
		swarmConstituentNumberCriteria = Integer.parseInt(n.getTextContent());  
		assert swarmConstituentNumberCriteria > 0;  // Must be set. 
		
		n = (Node) xPath.compile("/params/Swarm/swarmConstituentProximityCriteria")
				.evaluate(params, XPathConstants.NODE);
		swarmConstituentProximityCriteria = Double.parseDouble(n.getTextContent());
		assert swarmConstituentProximityCriteria > 0;  // Must be set.
	}
}