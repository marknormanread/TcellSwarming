package reenact;

import java.util.Collections;
import java.util.HashMap;

import sim.engine.SimState;
import sim.engine.Steppable;
import sim.engine.Stoppable;
import core.Cell;
import core.Simulation;
import core.TCell;


/**
 * Reenacts a cell, free of any cell-specific state, simply capturing location.
 * 
 * @author Mark N. Read
 *
 */
public class ReenactableCell extends TCell 
{
	// Indexed by time. Captures the cell's state (mainly location) at that time. 
	private HashMap<Double, CellStateCollection> states; 
	
	public ReenactableCell(HashMap<Double, CellStateCollection> states)
	{
		this.states = states; 
		
		// Schedule this cell to be stepped at each time point for which it has data. 
		for (double time : states.keySet())
		{
			
			Simulation.instance.schedule.scheduleOnce(time, Simulation.cellOrdering, this);
			
			// Removes cells from compartment after every iteration. 
			// This doesn't seem to be necessary now, but I leave it commented out. 
			// The snapper has been misbehaving quite a lot, and I can't figure out why. 
			// I think it comes down synchronization and completion inconsistencies between threads of the simulation and the GUI. 
			// if you use this, then remove the stuff below. 
			final RemoveFromCompartment remover = new RemoveFromCompartment(this);
			final double removalTime = time + (Simulation.instance.timeSlice_min * 0.9);
			Simulation.instance.schedule.scheduleOnce(removalTime, 1, remover);
		}
		
		// Removes the cell from the simulation's space just after the GUI has drawn the frame, but before the next 
		// time iteration. 
//		final double maxTime = Collections.max(states.keySet());
//		final double removalTime = maxTime + (Simulation.instance.timeSlice_min * 0.9);		
//		final RemoveFromCompartment remover = new RemoveFromCompartment(this);
//		Simulation.instance.schedule.scheduleOnce(removalTime, 1, remover);
	}	
	
	@Override
	public void step(SimState state) 
	{	
		CellStateCollection stateAtTime = states.get(Simulation.instance.schedule.getTime());
		Simulation.space.cellField.setObjectLocation(this, stateAtTime.location);
	}
	
	/** 
	 * Called by GUI. Currently defaults to false, can be implemented if this information is carried over in
	 * the recorded cell state files. 
	 */
	public boolean isChemotactic()
	{	
		CellStateCollection stateAtTime = states.get(Simulation.instance.schedule.getTime());
		if (stateAtTime != null)
			return stateAtTime.chemotactic;
		
		return false;
	}
	
	/** Used to remove cells from the simulation's spatial representation once they have no more data to show. */
	private class RemoveFromCompartment implements Steppable
	{
		private Cell removee;
		public RemoveFromCompartment(Cell removee)
		{
			this.removee = removee;
		}
		
		@Override
		public void step(SimState state)
		{
			Simulation.space.cellField.remove(removee);
		}
	}
}
