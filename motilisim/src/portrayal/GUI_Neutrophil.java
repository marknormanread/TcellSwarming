package portrayal;

import java.awt.Color;

import loggers.Snapper;
import sim.display.Console;
import sim.engine.Schedule;
import sim.engine.SimState;
import sim.portrayal3d.continuous.ContinuousPortrayal3D;
import sim.portrayal3d.grid.ValueGridPortrayal3D;
import sim.portrayal3d.simple.CylinderPortrayal3D;
import sim.portrayal3d.simple.TransformedPortrayal3D;
import core.EpidermalCell;
import core.HairFollicle;
import core.Neutrophil;
import core.Simulation;
import core.SimulationChemotacticAgent;
import core.SimulationNeutrophil;
import environment.Compartment3D;

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
    
 * @author Mark N. Read
 *
 */
public class GUI_Neutrophil extends SimulationGUI
{
	public static Color hairFollicleColor = new Color(0.2f, 0.2f, 0.2f, 0.5f);
	
	public ContinuousPortrayal3D samplePortrayal;
	// These are inefficient cube-based portrayals that can be turned off. 
	public ValueGridPortrayal3D restrictedPortrayal;
	
	
	private GUI_Neutrophil()
	{	
		super(new SimulationNeutrophil(
			"results/ak_bugfix",
			"results/ak_bugfix/parameters.xml"));		
	}
	
	protected GUI_Neutrophil(SimState sim)  // Used by subclasses
	{	super(sim);		}
	
	
	public static String getString() {	return "NeutroSwarm";	}
	

	public void start()
	{
		System.out.println("NeutroSwarm3D - starting");
		super.start();
	}
	
	/** Called once, sets up how things will be drawn. */
	public void setupPortrayals()
	{		
		Simulation simulation = (Simulation) state;
		Compartment3D compartment = (Compartment3D)simulation.space;
		
		cellPortrayal.setField(simulation.space.cellField);
		cellPortrayal.setPortrayalForClass(Neutrophil.class, new NeutrophilPortrayal() );		
			
		// Lysed epidermal cell portrayal. This are cylinders. 
		TransformedPortrayal3D epidermalTrans = new TransformedPortrayal3D(
				new CylinderPortrayal3D(new Color((float)1.0, (float)0.0, (float)0.0, (float)0.5)));
		epidermalTrans.rotateX(90.0);	// Turn the cylinder upright
		// Scale the size of the cylinder. This is wide, but flat. 
		epidermalTrans.scale(EpidermalCell.diameter, EpidermalCell.diameter, 1.0);	
		// Move the cylinder out of the tissue volume slightly. 
		epidermalTrans.translate(0.0, 0.0, -0.5);		
		cellPortrayal.setPortrayalForClass(EpidermalCell.class, epidermalTrans);		
				
		// Follicles are drawn as cylinders.
		// Create a custom portrayal for each hair follicle. 
		for(HairFollicle h : compartment.follicles)
		{
			// Need to perform transformations on the resultant cylinder to make it the right size. 
			TransformedPortrayal3D trans = new TransformedPortrayal3D(new CylinderPortrayal3D(hairFollicleColor));
			// Turn the cylinder upright
			trans.rotateX(90.0);		
			// Set its diameter, and have it extend the length of the tissue volume
			trans.scale(h.getDiameter(), h.getDiameter(), Compartment3D.extremeDepth);
			// Move it so that it sits within the tissue volume. Because the follicle is set at depth 0, it needs to 
			// Be translated after scaling. 
			trans.translate(0.0, 0.0, -Compartment3D.extremeDepth / 2.0);
			cellPortrayal.setPortrayalForObject(h, trans);
		}
				
		// Redraw the scene. 
		display.setBackdrop(Color.BLACK);		
		display.createSceneGraph();
		display.reset();
		if (imageVolume)
		{
			Snapper snapper = new Snapper(display);
			// If this starts misbehaving, look at how T cell does it, where it appears to be working.
			// Need to bring the snapper out of step with simulation iteration. 
			// I havent been able to figure out why. 
			simulation.schedule.scheduleRepeating(Schedule.EPOCH, 100, snapper, Simulation.timeSlice_min);
		}		
	}
		
	public static void main(String[] args)
	{
        SimulationGUI sim = new GUI_Neutrophil();
        Console c = new Console(sim);
        c.setVisible(true);	
	}
}
