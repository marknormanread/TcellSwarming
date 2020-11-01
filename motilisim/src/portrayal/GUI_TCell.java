package portrayal;

import java.awt.Color;

import loggers.Snapper;
import reenact.ReenactableCell;
import sim.display.Console;
import sim.engine.Schedule;
import sim.engine.SimState;
import sim.portrayal3d.simple.CylinderPortrayal3D;
import sim.portrayal3d.simple.TransformedPortrayal3D;
import core.Simulation;
import core.SimulationTCell;
import core.TCell;
import environment.BoundedCylinder;
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
public class GUI_TCell extends SimulationGUI
{
//	public static Color bolusColour = new Color(0.5f, 0.0f, 0.5f, 0.4f);  # Dark maroon
	public static Color bolusColour = new Color(0.2f, 0.2f, 0.0f, 0.2f);

	public GUI_TCell()
	{	
		super(new SimulationTCell(
			"results/20200616-well-saturation/rep1",
			"results/20200616-well-saturation/parameters.xml"
			));	
	}

	public GUI_TCell(SimState sim)
	{	super(sim);	}

	public static String getString() {	return "TCells";	}

	public void start()
	{
		printStarting();
		super.start();
	}
	
	public void printStarting()
	{	System.out.println("TCell Simulation - starting");		}

	@Override
	protected void setupPortrayals()
	{
		Simulation simulation = (Simulation) state;

		cellPortrayal.setField(simulation.space.cellField);
		cellPortrayal.setPortrayalForClass(TCell.class, new TCellPortrayal() );
		// ReenactableCells are not caught by TCells. Whatever does this can't handle subclasses!
		cellPortrayal.setPortrayalForClass(ReenactableCell.class, new TCellPortrayal() );

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
		display.setBackdrop(Color.WHITE);
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

	public static void main(String[] args)
	{
        SimulationGUI sim = new GUI_TCell();
        Console c = new Console(sim);
        c.setVisible(true);
	}
}
