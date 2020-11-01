package environment;

import sim.engine.SimState;
import core.Cell;

/** This class does nothing, it is simply used for drawing the bolus boundary in the GUI */
public class BolusCylinder extends Cell
{

	@Override
	public void step(SimState arg0) {	/* does nothing */	}

	/* used for hard shell spherical collision detection, hence return zero */
	@Override
	public double getDiameter() {	return 0;	}

	/* used for hard shell spherical collision detection, hence return zero */
	@Override
	public double getRadius() {	return 0;	}

	@Override
	public String getType() {
		return null;
	}

}
