package movement;

import utils.Quaternion;
import core.MigratoryCell;
import core.Simulation;

/**
 * Predominantly used for testing purposes, allows for selection of orientations that result in corkscrewing. This is 
 * accomplished by not randomly inverting the roll angle. Hence, if there is a non-zero mean roll anlge, it will always
 * be in the same rotational direction relative to the cell.
 *   
 * @author Mark N. Read, 2018
 *
 */
public class HomogeneousCRWCorkscrew_Orientation extends HomogeneousCRW_Orientation
{

	/* Singleton pattern */
	public static HomogeneousCRWCorkscrew_Orientation instance = new HomogeneousCRWCorkscrew_Orientation();	
	private HomogeneousCRWCorkscrew_Orientation()
	{}
	
	/** 
	 * Provides a new orientation of the cell, based on its current orientation (supplied as arg), and 
	 * distributions describing their propensity to chang direction.
	 * 
	 * Cells are rotated along their axis (direction in which they face, defined as the x-axis relative to the 
	 * cell), essentially rolling them, after which their pitch is changed (rotate around y-axis). 
	 * 
	 * @return The new orientation. 
	 */
	@Override
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell)
	{	
		/* Note the order that quaternion multiplication is done here - the rotations are expressed relative to the
		 * cell, hence orientation is multiplied by the rotation quaternion.  */
		// roll the cell along it's x-axis (axis in which it faces). 
		// This rolls the cell, but doesn't change it's heading or pitch. 
		double roll;
		if (rollRateMean < 0.0) {
			// if mean roll rate is negative, assume this indicates a uniform distribution
			// should be used. 
			roll = Simulation.instance.random.nextDouble() * 2.0 * Math.PI;
		} else {
			roll = (Simulation.instance.random.nextGaussian() * rollRateStd) + rollRateMean;
			// DO NOT RANDOMLY INVERT ROLL ANLGE, HENCE PROVIDING FOR CORKSCREWING. 
			roll *= Simulation.timeSlice_min;
		}
		// roll as a quaternion.
		Quaternion rotateQ = Quaternion.representRotation
				(roll, MigratoryCell.x_axis.x, MigratoryCell.x_axis.y, MigratoryCell.x_axis.z); 
		// multiply orientation by rotateQ, because rotateQ is calculated relative to cell, not in absolute space. 
		orientation = orientation.multiply(rotateQ).normalise();	// alter the cell's orientation. 
		
		// change cell pitch (roll along the y axis). Pitch can be changed in both positive and negative directions.
		double pitch = (Simulation.instance.random.nextGaussian() * pitchRateStd) + pitchRateMean;
		pitch *= Simulation.timeSlice_min;		// account for timestep.
		Quaternion pitchQ = Quaternion.representRotation
				(pitch, MigratoryCell.y_axis.x, MigratoryCell.y_axis.y, MigratoryCell.y_axis.z);
		// multiply orientation by rotateQ, because pitchQ is calculated relative to cell, not in absolute space.
		orientation = orientation.multiply(pitchQ).normalise();
		return orientation;
	}
	
}
