package reenact;

import java.awt.Color;
import java.util.Arrays;
import java.util.List;

import loggers.Snapper;
import portrayal.GUI_TCell;
import portrayal.SimulationGUI;
import portrayal.TCellPortrayal;
import sim.display.Console;
import sim.display.GUIState;
import sim.engine.SimState;
import core.Simulation;
import core.TCell;

/**
 * Entry point into a re-enactment of a previously executed simulation. 
 * 
 * @author Mark N. Read
 *
 */
public class GUI_Reenact extends GUI_TCell
{
	private static String dataDirectory = "/Users/markread/dropbox_usyd/projects/biro/neutroswarm/imaris/simulation/20181213-well-bootstrap_reconstruction-chemotaxis-1h_20s/rep1/side_bolus";
	private static long interframeDelay = 0;  // Delay between subsequent time steps in milliseconds. 
	
	private GUI_Reenact(SimState sim)
	{	super(sim); 	}
	
	public void start()
	{
		System.out.println("ReenactorGUI - starting");
		super.start();
	}		
	
	@Override
	protected void setupPortrayals()
	{
		super.setupPortrayals();  // Do everything the super-class would do here. 
		// ReenactableCells are not caught by TCells. Whatever does this can't handle subclasses!
		cellPortrayal.setPortrayalForClass(ReenactableCell.class, new TCellPortrayal() );
	}
	
	public static void main (String[] args)
	{
		Simulation.cmdLineArgs = args;
		
		/* Use this to take snap shots for making a video. */
		SimulationGUI.imageVolume = true;
		/* There's a problem with simulation timing when reenacting. 
		 * Reenacted snapped images are recorded as being one timestep later than they should be (and are when 
		 * visualising directly in a normal non-reenacted simulation). 
		 * I'm not sure what the problem is, but I think it relates to the ordering in which cells are stepped, the 
		 * simulation time is stepped, the GUI is drawn, and cell locations are drawn. 
		 * It's like the GUI is updated, and then the simulation time is stepped in antipication of the next action 
		 * (but not drawn on the GUI).
		 * I haven't been able to correct this behaviour. 
		 * This is a hack, but it works.
		 * We set the snapper to record times as one timestep duration less than the simulation's schedule reports for
		 * reenactments. 
		 */
		Snapper.reenacting = true; 
		
		List<String> list = Arrays.asList(args);
		int di = list.indexOf("-d");  // -1 if item not in List
		if (di >= 0)  
			dataDirectory = list.get(di + 1);
		
		// Reenacting existing data, so write any snapshots of space to that same directory (not some other one).
		Simulation.defaultOutputDir = dataDirectory;
		Simulation.outputPath = dataDirectory;
		
		SimState sim = new ReenactTCell(dataDirectory, interframeDelay);
		GUIState gui = new GUI_Reenact(sim);
		// Not quite sure (yet) which this is needed, but if it is left as a slice, no cells are drawn. 
		SimulationGUI.imageOutsideRecorded = SimulationGUI.RecordingVolumeMode.VOLUME;
		
        Console c = new Console(gui);
        c.setWhenShouldEndTime(Simulation.endTime_min);
        c.setVisible(true);	
	}
}
