package reenact.neutrophil;

import java.util.Vector;

import core.Neutrophil;
import sim.engine.SimState;

public class ReenactableNeutrophil extends Neutrophil
{
	private Vector<NeutrophilStateCollection> states;
	private int progress;  // Iterates through recorded states in sequence 
	
	
	public ReenactableNeutrophil(Vector<NeutrophilStateCollection> states, ReenactorNeutrophil sim)
	{
		this.states = states;
		this.progress = 0;	
		// Place cell at its starting location.
		this.location = this.states.get(this.progress).location;
		sim.cellField.setObjectLocation(this, this.location);  // Tell environment where the cell is. 
	}
	
	@Override
	public void step(SimState state) 
	{	
		// Don't point to an index bigger the data available.
		if (progress+1 < this.states.size())
			progress++;
		
		// Where should cell be situated? Place it there. 
		this.location = this.states.get(progress).location;
		((ReenactorNeutrophil)state).cellField.setObjectLocation(this, this.location);
	}
	
	@Override
	public boolean isStateRandomWalk()
	{ 	// Random walk when not one of the others.
		return ! (isStateAttractantStimulated() || isStateLTB4Stimulated()); 
	}
	
	@Override
	public  boolean isStateAttractantStimulated()
	{	return states.get(progress).attractantStimulated;	}
	
	@Override
	public boolean isStateLTB4Stimulated()
	{ 	return states.get(progress).ltb4Stimulated;	}
	
	// TODO should this be re-introduced into Neutrophil?
//	@Override
//	public boolean inSwarm()
//	{	return states.get(progress).inSwarm;	} 
}
