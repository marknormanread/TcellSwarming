package core;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import environment.BoundedCube;
import environment.BoundedCylinder;
import environment.Compartment3D;
import filesystem.FileSystemIO;
import loggers.AttractantTissueZStack;
import loggers.CellLogger;
import loggers.SeedLogger;
import loggers.TimeLogger;
import sim.engine.SimState;
import sim.engine.Steppable;

/**
 * This program is free software: you can redistribute it and/or modify
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
 * This is the simulation driver, it is the top level organising entity of a simulation execution.
 *
 * @author Mark N. Read
 *
 */
@SuppressWarnings("serial")
public abstract class Simulation extends SimState
{
	// java representation of the XML document holding parameters.
	public static String parametersPath = null;
	
	public static String defaultOutputDir = "results/test_reduced_infiltration/";

	public static Document parameters = null;
	public static String outputPath = defaultOutputDir;
	public static String[] cmdLineArgs = new String[0];  // Command line arguments. 

	public static long seed = 999;

	public static Compartment3D space;

	public static double endTime_min = 50.0;		// when to terminate the simulation.
	public static double timeSlice_min = 0.5;  		// duration covered by a stimulation time-step. In minutes.
	// times at which simulation state is sampled and recorded. Should match in vivo work.
	public static double sampleTimeSlice_min = 0.5;
	// how many iterations the simulation has been through. Starts at 1 to be consistent with IMARIS.
	public long timeIter = 1;

	public static boolean trackCells = true;
	public static CellLogger cellLogger = null;
	// make images of cross section of tissue volume to visualise soluble factors
	public static boolean zstacks = false;
	public static double snapper_timeslice_min = 5.0;
	
	public static boolean reenactment = false;  // Don't write data to filesystem if reenacting. 

	// cells reside outside the imaging volume, and enter it. Hence, simulation can be set up with cells
	// occupying this space (buffer), at similar density to those in the imaging volume (at launch time). This parameter
	// specifies the size of the buffer as a proportion of the imaging volume's size in each dimension. Hence, a value
	// of 0.0 has no buffer; a value of 1.0 has an entire imaging volume's worth of buffer in each dimension.
	public static double bufferSize = 0.0;	// must be >= 0.0.

	// will only have one simulation running at a time, so this instance can be used by other classes.
	public static Simulation instance;

	// ordering for simulation components added to the schedule.
	public static final int blastOrdering = 0;			// laser blast happens before everything else.
	public static final int compartmentOrdering = 1;	// compartments are always stepped before cells.
	public static final int cellOrdering = 2;			// the order for cells to be stepped by the schedule.
	public static final int loggerOrdering = 3;			// logger stepping order, comes after cells.
	public static final int snapperOrdering = 0;
	public static final int timeIterOrdering = 10;		// increments the iteration count of simulation time steps.


	// Represents different environments
	public static Environment environment = Environment.TOP_BOUNDED_CUBE;
	public static enum Environment {
		UNBOUNDED_CUBE,
		TOP_BOUNDED_CUBE,
		TOP_BOTTOM_BOUNDED_CUBE,
		WELL
	}

	public Simulation()
	{
		this(Simulation.outputPath, Simulation.parametersPath, false);
	}

	public Simulation(String outDir, String parametersPath, boolean reenactment)
	{
		this(outDir, parametersPath, reenactment, Simulation.seed);
	}
	
	public Simulation(String outDir, String parametersPath, boolean reenactment, long seed)
	{
		super(seed);
		instance = this;
		Simulation.outputPath = outDir;
		Simulation.parametersPath = parametersPath;
		Simulation.reenactment = reenactment;
		parameters = FileSystemIO.openXMLFile(parametersPath);
		Simulation.setupSimulationParameters(reenactment);
		random.setSeed(seed);
	}

	public abstract String getDefaulParametersPath();
	public abstract void populateCells();

	/**
	 * Subclasses can override this if a different compartment is needed.
	 */
	public Compartment3D initializeCompartment()
	{
		switch (environment)
		{
			case UNBOUNDED_CUBE: return new BoundedCube(false, false);
			case TOP_BOUNDED_CUBE: return new BoundedCube(true, false);
			case TOP_BOTTOM_BOUNDED_CUBE: return new BoundedCube(true, true);
			case WELL: return new BoundedCylinder();
		}
		throw new RuntimeException("No environment specified.");
	}

	/** Sets up the simulation. */
	public void start()
	{
		super.start();
		
		// ensure the experimental directory has been set up.
		File dir = new File(outputPath);
		if (!dir.exists())
		{
			System.out.println("creating directory : " + outputPath);
			boolean result = dir.mkdirs();
			if (!result)	System.out.println("ERROR: could not create directory " + outputPath);
		}


		final String outParam = outputPath + "/parameters.xml";
		// Don't overwrite the file if it is also the input parameters file.
		try {
			final Path srcPath = Paths.get(parametersPath);
			final Path dstPath = Paths.get(outParam);

			if (!Files.exists(srcPath))
				throw new RuntimeException("Fatal! Supply input parameters.xml file.\nSupplied path was "
						+ parametersPath);

			// Copy parameters file if the destination doesn't exist, or destination is not also source (rd/wr error).
			if(! Files.exists(dstPath) || ! Files.isSameFile(Paths.get(parametersPath), Paths.get(outParam)))
				FileSystemIO.copyFile(parametersPath, outParam);
		} catch (IOException e) {
			System.out.println("Error working with parameters files.");
			e.printStackTrace();
			throw new RuntimeException("Fatal. Don't mess around with parameter files, fix the problem.");
		}
		
		if (! reenactment) 
		{
			TimeLogger.writeTimeData_Sec(outputPath);
			SeedLogger.writeSeedData(outputPath);
		}
		
		// Anonymous inner class to increment the time iteration count.
		Steppable tih = new Steppable() {
			@Override
			public void step(SimState state)
			{
				((Simulation)state).timeIter ++;
			}
		};
		schedule.scheduleRepeating(0., timeIterOrdering, tih, sampleTimeSlice_min);

		space = initializeCompartment();
		schedule.scheduleRepeating(0., Simulation.compartmentOrdering, space, timeSlice_min);

		populateCells();
		
		schedule.scheduleRepeating(CellLogger.cellLoggingInterval, Simulation.loggerOrdering, cellLogger, 
								   CellLogger.cellLoggingInterval);
		
		if (zstacks)
			schedule.scheduleRepeating(0.0, loggerOrdering, new AttractantTissueZStack(outputPath),
					snapper_timeslice_min);

		timeIter += 1;  // must increment this too, else two records for iteration 1 are generated.
	}
	
    /**
     * Loads the parameters.xml config file and loads the default parameters for all classes in the simulation. Should be the first thing that is done when
     * running the simulation, with GUI or without. Abstract classes must be called before concrete classes.
     *
     * Reenactions don't require full parameters to be loaded.
     */
    public static void setupSimulationParameters(boolean reenaction)
    {
    	try
    	{
    		/* read in the default parameters for the various classes in the simulation */
            loadParameters(parameters);
            if (! reenaction)
            {
            	MigratoryCell.loadParameters(parameters);
            	CellLogger.loadParameters(parameters);
            }
            Compartment3D.loadParameters(parameters);
            loadEnvironmentParameters(parameters);

    	}
    	catch(Exception e)
    	{
			System.out.println("ERROR reading in parameters: " + e.toString());
			e.printStackTrace();
		}
    }

    public static void loadEnvironmentParameters(Document params) throws XPathExpressionException
    {
        if (environment == Environment.TOP_BOUNDED_CUBE ||  environment == Environment.TOP_BOTTOM_BOUNDED_CUBE
        		|| environment == Environment.UNBOUNDED_CUBE)
        	BoundedCube.loadParameters(params);
        else if (environment == Environment.WELL)
        	BoundedCylinder.loadParameters(params);
    }

    /**
     * Given the parameters.xml file (represented as a 'Document') this method loads the relevant default values for
     * the top level simulation.
     * @param params
     */
	private static void loadParameters(Document params)
	{
		try {
			System.out.println("params = " + params);
			XPath xPath =  XPathFactory.newInstance().newXPath();
			Node n;

			n = (Node) xPath.compile("/params/Simulation/endTime").evaluate(params, XPathConstants.NODE);
			endTime_min = Double.parseDouble(n.getTextContent());
			n = (Node) xPath.compile("/params/Simulation/timeSlice_s").evaluate(params, XPathConstants.NODE);
			timeSlice_min = Double.parseDouble(n.getTextContent()) / 60.;
			n = (Node) xPath.compile("/params/Simulation/sampleTimeSlice_s").evaluate(params, XPathConstants.NODE);
			if (n != null && !n.getTextContent().equals(""))
				sampleTimeSlice_min = Double.parseDouble(n.getTextContent()) / 60.;
			else
				sampleTimeSlice_min = timeSlice_min;  // Default case

			n = (Node) xPath.compile("/params/Simulation/snapper_timeslice_s").evaluate(params, XPathConstants.NODE);
			snapper_timeslice_min = sampleTimeSlice_min;  // Default case. 
			if (n != null && !n.getTextContent().equals(""))
			{
				snapper_timeslice_min = Double.parseDouble(n.getTextContent()) / 60.;
				System.out.println("Setting snapper_timeslice_min to" + snapper_timeslice_min);
			}
			
			n = (Node) xPath.compile("/params/Simulation/zstacks").evaluate(params, XPathConstants.NODE);
			zstacks = n != null ? Boolean.parseBoolean(n.getTextContent()) : false;

			environment = selectEnvironment(params);

		} catch(XPathExpressionException e) {
				System.out.println("ERROR reading in parameters: " + e.toString());
		}
	}

	public static Environment selectEnvironment(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;
		n = (Node) xPath.compile("/params/Simulation/environment").evaluate(params, XPathConstants.NODE);
		if (n != null)
		{
			System.out.println("environment =  " + n.getTextContent());
			if(n.getTextContent().equals("well"))
				return Environment.WELL;
			if(n.getTextContent().equals("topBottomCube"))
				return Environment.TOP_BOTTOM_BOUNDED_CUBE;
			if(n.getTextContent().equals("topCube"))
				return Environment.TOP_BOUNDED_CUBE;
			if(n.getTextContent().equals("unboundedCube"))
				return Environment.UNBOUNDED_CUBE;
		}
		throw new RuntimeException("Please specify environment in xml parameters file");
	}

	/** Read command line arguments */
	public static void readArgs(String[] args)
	{
		cmdLineArgs = args;  // Save for use elsewhere in program. 
		seed = -1;
		int i = 0;
		while (i < args.length)
		{
			String command = args[i];
			if (command.equals("-e"))		// the end time, as a double. If not supplied here, should be supplied
			{								// with the parameter file.
				i++;
				endTime_min = Double.parseDouble(args[i]);
			}
			else if (command.equals("-p"))	// location of the parameters file.
			{
				i++;
				parametersPath = args[i];
			}
			else if (command.equals("-o")) 	// where the simulation output files are to be written.
			{
				i++;
				outputPath = args[i];
			}
			else if (command.equals("-s"))	// seed
			{
				i++;
				seed = Long.parseLong(args[i]);
			}
			i++;
		}
//		if (seed == -1)
//			seed = System.currentTimeMillis();
	}

	/** Sets up and then executes the simulation. This is the primary driver loop. */
	public static void execute(Simulation state)
	{
		state.start();
		// This is the main driver loop.
		do
		{
			System.out.println("simulated time = " + state.schedule.getTime());
			if (!state.schedule.step(state))
				break;

		} while(state.schedule.getTime() < endTime_min);
		state.finish();

		System.exit(0);
	}
}
