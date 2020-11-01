package reenact;

import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import core.Simulation;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.Double3D;

@SuppressWarnings("serial")
public class ReenactTCell extends Simulation
{
	
	private String targetDataset;
	private long interframeDelay = 0;
			
	public ReenactTCell(String directory, long interframeDelay)
	{	
		super(directory, directory + "/parameters.xml", true); 
		
		this.targetDataset = directory;
		this.interframeDelay = interframeDelay;
	
		/* In reenacting, we rely on data written to the file system as the smallest time-slice of simulation stepping.
		 * This may not be the same as the time-slice of the original simulation - sometimes positional recording is 
		 * done less frequently. For the purposes of reenaction here, we need the simulation time slice to be the same
		 * as the time step at which data was recorded.  
		 */
		Simulation.timeSlice_min = Simulation.sampleTimeSlice_min;
	}
	
	/**
	 * Sets up the simulation. 
	 */
	public void start()
	{
		System.out.println("Starting reenaction.");
		super.start();	
		
		// Used to slow play-back to a desirable speed.
		if (interframeDelay > 0)
		{
			ReenactorPlayspeedController delay = new ReenactorPlayspeedController(250);
			schedule.scheduleRepeating(0.0, 0, delay, Simulation.timeSlice_min);
		}
		
		HashMap<Integer, HashMap<Double, CellStateCollection>> cells =
						readCellLocations(this.targetDataset + "/_Position.csv");
		
		for (Integer trackID : cells.keySet())
		{
			HashMap<Double, CellStateCollection> states = cells.get(trackID);
			Steppable s = new ReenactableCell(states);  // Schedules itself for execution.
		}		
	}
	
	/*
	 * Reads the supplied CSV file of cell locations in time, and return a hashmap containing the results. 
	 * HashMap is indexed by cell track ID, and contains a HashMap =(timePointID, cellState). 
	 * CellState contains, amongst other things, positional data for the cell at a given point in time.  
	 */
	private HashMap<Integer, HashMap<Double, CellStateCollection>> readCellLocations(String filename)
	{
		HashMap<Integer, HashMap<Double, CellStateCollection>> cells = 
					new HashMap<Integer, HashMap<Double, CellStateCollection>>();
		
		double maxTime = 0;  // Used to stop reenaction when no more data left. 
		try{
			Reader reader = Files.newBufferedReader(Paths.get(filename));
            CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT
                    .withFirstRecordAsHeader()
                    .withIgnoreHeaderCase()
                    .withTrim());
		
            
            for (CSVRecord csvRecord : csvParser)
            {
            	double posX = new Double(csvRecord.get("Position_X"));
		    	double posY = new Double(csvRecord.get("Position_Y"));
		    	double posZ = new Double(csvRecord.get("Position_Z"));
		    	
		    	boolean chemotactic = false;
		    	try {
		    		chemotactic = new Boolean(csvRecord.get("Chemotactic"));
		    	} catch (IllegalArgumentException e) {  /* column not suppled; ignore */ }
		    	
		    	int trackID = new Integer(csvRecord.get("TrackID"));
		    	int timePointID = new Integer(csvRecord.get("TimePointID"));
		    	double time = new Double(csvRecord.get("Time"));
		    	
		    	Double3D pos = new Double3D(posX, posY, posZ);
		    	// If the cell has not previously been encountered, create the constructs for it. 
		    	if (! cells.containsKey(trackID))
		    		cells.put(trackID, new HashMap<Double, CellStateCollection>());
		    	// Log cell position.
		    	CellStateCollection state = new CellStateCollection(pos, chemotactic, time, timePointID, trackID);
		    	cells.get(trackID).put(time, state);
		    	
		    	if (time > maxTime) 
		    		maxTime = time;
            }
            csvParser.close();
            
		} catch(Exception e) {
			System.out.println("Exception thrown reading smiulation data. " + e);
			e.printStackTrace();
		}	
		
		// Create anonymous class to handle stopping the reenaction when there's nothing more to reenact. 
		Steppable killReenaction = new Steppable()
		{
			public void step(SimState state)
			{
				Simulation.instance.finish();
			}
		};
		Simulation.instance.schedule.scheduleOnce(maxTime + 1., killReenaction);		
		return cells;
	}
	
	/**
	 * Tears down the simulation, can be used for writing IO. 
	 */
	public void finish()
	{}
	
	public void populateCells() {}  // Unused. 
	public String getDefaulParametersPath() {	return null;	}  // Unused.
}
