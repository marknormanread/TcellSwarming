package reenact.neutrophil;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

import loggers.SwarmLogger;
import reenact.ReenactorPlayspeedController;
import sim.engine.SimState;
import sim.field.continuous.Continuous3D;
import sim.util.Double3D;
import au.com.bytecode.opencsv.CSVReader;
import core.EpidermalCell;
import core.HairFollicle;
import core.Simulation;
import environment.Compartment3D;


/**
 * Class is designed to re-enact a pre-run simulation. It reads in the cell positions, and renders their movements
 * 
 * @author Mark N. Read, 2018
 *
 */
@SuppressWarnings("serial")
public class ReenactorNeutrophil extends SimState
{	
	public static boolean trackNeutrophils = true;
	public static SwarmLogger swarmLogger = null;
	//public static NeutrophilLogger neutroLogger = null;  // TODO can this be deleted?
	
	public Continuous3D cellField;   // Stores cells in a continuous space
	// Storing all hair follicles here in case we need to access them later (e.g., the GUI does, and collision 
	// detection may also).
	public ArrayList<HairFollicle> follicles = new ArrayList<HairFollicle>();

	public ReenactorNeutrophil(long seed)
	{	super(seed);	}
	
	private void setupEnvironment()
	{		
		cellField = new Continuous3D(				
				1.0,  // Discretization, dividing space into regions for maintaining a map of objects' locations.  
				Compartment3D.extremeWidth, 
				Compartment3D.extremeHeight, 
				Compartment3D.extremeDepth);
		cellField.clear();

		// used to slow play-back to a desirable speed. 
		ReenactorPlayspeedController delay = new ReenactorPlayspeedController(250);
		schedule.scheduleRepeating(0.0, 0, delay, Simulation.timeSlice_min);
		
		HashMap<Integer, Vector<NeutrophilStateCollection>> cells = readCellLocations("_CellState.csv");
		for (Integer id : cells.keySet())
		{
			Vector<NeutrophilStateCollection> states = cells.get(id);
			ReenactableNeutrophil neutro = new ReenactableNeutrophil(states, this);
			schedule.scheduleRepeating(0.0, 1, neutro, Simulation.timeSlice_min);
		}
		
		Vector<FollicleLocation> follicles = readFollicleLocations("_FolliclePosition.csv");
		for (FollicleLocation f : follicles)
		{
			HairFollicle follicle = new HairFollicle(f.diameter, f.location);
			cellField.setObjectLocation(follicle, f.location);
			this.follicles.add(follicle);
		}
		
		Vector<EpidermalCell> burns = readBurnLocations("_BurnPosition.csv");
		for (EpidermalCell ec: burns)
		{
			System.out.println("found burn");
			cellField.setObjectLocation(ec, ec.getCurrentLocation());
		}
	}
	
	/**
	 * Sets up the simulation. 
	 */
	public void start()
	{
		System.out.println("starting sim.");
		super.start();
		setupEnvironment();
	}
	
	/*
	 * Reads the supplied CSV file of cell locations in time, and return a hashmap containing the results. 
	 * HashMap is indexed by an Integer key which uniquely identifies each cell. Values constitute a Vector of 
	 * TimeLocation objects, which store locations in 3D space and a time for that cell. 
	 */
	private HashMap<Integer, Vector<NeutrophilStateCollection>> readCellLocations(String filename)
	{
		HashMap<Integer, Vector<NeutrophilStateCollection>> cells = 
					new HashMap<Integer, Vector<NeutrophilStateCollection>>();
		try{
			CSVReader reader = new CSVReader(new FileReader(filename));
		    String [] line;
		    
		    // Throw away the header line. 
		    reader.readNext();

		    while ((line = reader.readNext()) != null) {
		    	// Read in data from the file. 
		    	// Rows are : "Parent,Time,Position X,Position Y,Position Z,Swarm,Attractant,LTB4
		    	Integer cell = new Integer(line[0]);
		    	double time = new Double(line[1]);		    	
		    	double posX = new Double(line[2]);
		    	double posY = new Double(line[3]);
		    	double posZ = new Double(line[4]);
		    	boolean swarm = new Boolean(line[5]);
		    	boolean attractant = new Boolean(line[6]);
		    	boolean ltb4 = new Boolean(line[7]);
		    	
		    	Double3D pos = new Double3D(posX, posY, posZ);
		    	// If the cell has not previously been encountered, create the constructs for it. 
		    	if (! cells.containsKey(cell))
		    		cells.put(cell, new Vector<NeutrophilStateCollection>());
		    	// log cell position.
		    	NeutrophilStateCollection state = new NeutrophilStateCollection(time, pos, attractant, ltb4, swarm);
		    	cells.get(cell).add(state);
		    }
		    reader.close();
		} catch(Exception e) {
			System.out.println("Exception thrown reading smiulation data. " + e);
		}	
		return cells;
	}
	
	private Vector<FollicleLocation> readFollicleLocations(String filename)
	{
		Vector<FollicleLocation> follicles = new Vector<FollicleLocation>();
		try
		{
			CSVReader reader = new CSVReader(new FileReader(filename));
		    String [] line;
		    // Throw away the header line. 
		    reader.readNext();
		    
		    while ((line = reader.readNext()) != null) 
		    {
		    	double posX = new Double(line[0]);
		    	double posY = new Double(line[1]);
		    	double posZ = new Double(line[2]);
		    	double diameter = new Double(line[3]);		    	
		    	
		    	Double3D location = new Double3D(posX,posY,posZ);
		    	follicles.add(new FollicleLocation(diameter, location));
		    }
		    reader.close();
		} catch(Exception e) {
			System.out.println("Exception thrown reading smiulation data. " + e);
		}
		return follicles;
	}
	
	
	private Vector<EpidermalCell> readBurnLocations(String filename)
	{
		Vector<EpidermalCell> burns = new Vector<EpidermalCell>();
		try
		{
			CSVReader reader = new CSVReader(new FileReader(filename));
		    String [] line;
		    // Throw away the header line. 
		    reader.readNext();
		    
		    while ((line = reader.readNext()) != null) 
		    {
		    	double posX = new Double(line[0]);
		    	double posY = new Double(line[1]);
		    	double posZ = new Double(line[2]);
		    	Double3D location = new Double3D(posX, posY, posZ);
		    	EpidermalCell burn = new EpidermalCell(true, location);
		    	cellField.setObjectLocation(burn, burn.getCurrentLocation());  // These cells don't move.
		    }
		    reader.close();
		} catch(Exception e) {
			System.out.println("Exception thrown reading smiulation data. " + e);
		}
		return burns;
	}
	
	/**
	 * Tears down the simulation, can be used for writing IO. 
	 */
	public void finish()
	{}
}
