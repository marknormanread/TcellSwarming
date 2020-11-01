package reenact.neutrophil;

import core.Simulation;
import portrayal.GUI_Neutrophil;
import sim.display.Console;
import sim.display.GUIState;
import sim.engine.SimState;

public class GUI_ReenactNeutrophil extends GUI_Neutrophil
{	
	public GUI_ReenactNeutrophil(SimState sim)
	{
		super(sim);
	}
	
	public void start()
	{
		System.out.println("ReenactorGUI - starting");
		super.start();
		setupPortrayals();
	}		
	
	
	public static void main (String[] args)
	{
		Simulation.cmdLineArgs = args;
		SimState sim = new ReenactorNeutrophil(0);
		GUIState gui = new GUI_ReenactNeutrophil(sim);
        Console c = new Console(gui);
        c.setVisible(true);	
	}
}