package reenact;

import sim.engine.SimState;
import sim.engine.Steppable;

@SuppressWarnings("serial")
public class ReenactorPlayspeedController implements Steppable
{
	private long delayMS;  // Delay in miliseconds.
	public ReenactorPlayspeedController(long delay)
	{
		this.delayMS = delay;  // Delay in miliseconds. 
	}


	public void step(SimState state) 
	{
		try {
		    Thread.sleep(this.delayMS);
		} catch(InterruptedException ex) {
		    Thread.currentThread().interrupt();
		}
	}
	
}