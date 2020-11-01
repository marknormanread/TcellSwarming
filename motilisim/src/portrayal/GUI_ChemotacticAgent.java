package portrayal;

import java.awt.Color;

import core.ChemotacticAgent;
import core.Simulation;
import core.SimulationChemotacticAgent;
import core.StationaryChemokineSource;
import environment.BoundedCylinder;
import environment.Compartment3D;
import loggers.Snapper;
import reenact.ReenactableCell;
import sim.display.Console;
import sim.engine.SimState;
import sim.portrayal3d.simple.CylinderPortrayal3D;
import sim.portrayal3d.simple.SpherePortrayal3D;
import sim.portrayal3d.simple.TransformedPortrayal3D;

public class GUI_ChemotacticAgent extends SimulationGUI
{
	public static Color bolusColour = new Color(0.5f, 0.0f, 0.5f, 0.4f);
	
	public GUI_ChemotacticAgent()
	{	
		super(new SimulationChemotacticAgent(
				"results/swarm_metric_v3/repulsion_epoch/test",
				"results/swarm_metric_v3/repulsion_epoch/parameters.xml",
				4));
	}

	public GUI_ChemotacticAgent(SimState sim)
	{	super(sim);	}
	
	public static String getString() {	return "ChemotacticAgents";	}
	
	public void start()
	{
		printStarting();
		super.start();
	}
		
	@Override
	protected void setupPortrayals()
	{
		Simulation simulation = (Simulation) state;

		cellPortrayal.setField(simulation.space.cellField);
		cellPortrayal.setPortrayalForClass(ChemotacticAgent.class, new ChemotacticAgentPortrayal() );
		cellPortrayal.setPortrayalForClass(StationaryChemokineSource.class, 
				new SpherePortrayal3D(Color.WHITE, StationaryChemokineSource.radius * SimulationGUI.cellScalar) );
		// ReenactableCells are not caught by TCells. Whatever does this can't handle subclasses!
		cellPortrayal.setPortrayalForClass(ReenactableCell.class, new ChemotacticAgentPortrayal() );

		// Draw the bolus, if there is one.
		if (simulation.space instanceof BoundedCylinder && BoundedCylinder.bolusPresent)
		{
			// Need to perform transformations on the resultant cylinder to make it the right size.
			TransformedPortrayal3D trans = new TransformedPortrayal3D(new CylinderPortrayal3D(bolusColour));
			// Turn the cylinder upright
			trans.rotateX(90.0);
			// Set its diameter, and have it extend the length of the tissue volume
			trans.scale(BoundedCylinder.bolusRadius*2, BoundedCylinder.bolusRadius*2, Compartment3D.extremeDepth);
			cellPortrayal.setPortrayalForObject(BoundedCylinder.bolus, trans);
		}

		// Redraw the scene.
		display.setBackdrop(Color.BLACK);
		display.createSceneGraph();
		display.reset();
		if (imageVolume)
		{
			Snapper snapper = new Snapper(display);
			simulation.schedule.scheduleRepeating(0. + Simulation.snapper_timeslice_min * 0.1,
					Simulation.snapperOrdering, snapper,
					Simulation.snapper_timeslice_min );
		}
	}
	
	public void printStarting()
	{	System.out.println("Chemotactic Agent Simulation - starting");		}
	
	public static void main(String[] args)
	{
        SimulationGUI sim = new GUI_ChemotacticAgent();
        Console c = new Console(sim);
        c.setVisible(true);
	}
}
