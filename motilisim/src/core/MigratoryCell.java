package core;


import java.io.IOException;
import java.util.ArrayList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import movement.Ballistic;
import movement.Bootstrap;
import movement.BootstrapDynamic;
import movement.BootstrapGradient;
import movement.Brownian;
import movement.HeterogeneousCRW_Orientation;
import movement.HeterogeneousCRW_Translation;
import movement.HomogeneousCRWCorkscrew_Orientation;
import movement.HomogeneousCRW_Orientation;
import movement.HomogeneousCRW_Translation;
import movement.InverseHeterogeneousCRW;
import movement.InverseHomogeneousCRW;
import movement.LevyFlight;
import movement.Orientation;
import movement.Translation;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.w3c.dom.Document;
import org.w3c.dom.Node;

import sim.engine.Schedule;
import sim.util.Double3D;
import soluble.ConcentrationPerception;
import utils.Quaternion;
import utils.Utils;


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
 * 
 * 
 * Cells that can be moved
 * @author Mark N. Read
 *
 */
public abstract class MigratoryCell extends Cell
{
	public enum TranslationParadigm 
	{
		BOOTSTRAP,  // resampling the jumps exhibited in a supplied input file.
		BOOTSTRAP_GRADIENT,  // always follows gradient (ensure there is one! 
		BOOTSTRAP_DYNAMIC,  // switches between non-directional and gradient. This is handled by a T cell subclass. 
		BALLISTIC,
		BROWNIAN,
		LEVY_VPL,	// Levy distributions with variable power law decay 
		HOMO_CRW,
		INVERSE_HOMO_CRW,
		HETERO_CRW,
		INDIVIDUALCORRELATED,
		INVERSE_HETERO_CRW
	}
	
	public enum OrientationParadigm
	{
		BOOTSTRAP,  // resampling the jumps exhibited in a supplied input file.
		BOOTSTRAP_GRADIENT,  // always follows gradient (ensure there is one!)
		BOOTSTRAP_DYNAMIC,  // switches between non-directional and gradient. This is handled by a T cell subclass.
		BALLISTIC,
		BROWNIAN,
		LEVY_VPL,	// Levy distributions with variable power law decay		
		HOMO_CRW,
		HOMO_CORKSCREW_CRW,
		INVERSE_HOMO_CRW,
		HETERO_CRW,
		INDIVIDUALCORRELATED,
		INVERSE_HETERO_CRW		
	}
	public static TranslationParadigm translationParadigm;
	public static OrientationParadigm orientationParadigm;
	
	// These are used by concrete subclasses to provide random walk movements. 
	protected Translation translationActuator;
	protected Orientation orientationActuator;
	
		
	// vector of movement resulting from knocking into other cells. Used to allow cells to slide over one another. 
	protected Double3D bounce;
	// Record of other cells that this cell is currently in contact with. 
	protected ArrayList<Cell> collidedCells = new ArrayList<Cell>();	 
		
	public enum CollisionResponse
	{
		RETURN,		// reverse path, back where cell came from (180 degree turn). 
		NORMAL,	 	// reorient along the normal of the contact
		SNOOKER,    // as snooker balls, bounce off. 
		SLIDE		// directionally persistent, but slide around the contacted cell 
	}
	public CollisionResponse collisionResponse = CollisionResponse.SLIDE;
		
	public MigratoryCell()
	{}

	public MigratoryCell(Schedule sched)
	{
		super(sched);
		
		bounce = new Double3D(0.0,0.0,0.0);		
		// assign a random orientation.		
		orientation = Quaternion.randomUniform();
		
		if (translationParadigm == TranslationParadigm.HOMO_CRW)
			translationActuator = HomogeneousCRW_Translation.instance;
		else if (translationParadigm == TranslationParadigm.HETERO_CRW)
			translationActuator = new HeterogeneousCRW_Translation();
		
		if (orientationParadigm == OrientationParadigm.HOMO_CRW)
			orientationActuator = HomogeneousCRW_Orientation.instance;
		else if (orientationParadigm == OrientationParadigm.HETERO_CRW)
			orientationActuator = new HeterogeneousCRW_Orientation();
		else if (orientationParadigm == OrientationParadigm.LEVY_VPL)
		{
			LevyFlight lf = new LevyFlight();
			orientationActuator = lf;
			translationActuator = lf;
		}
		else if (orientationParadigm == OrientationParadigm.INVERSE_HETERO_CRW)
		{
			InverseHeterogeneousCRW icto = new InverseHeterogeneousCRW();
			orientationActuator = icto;
			translationActuator = icto;
		}
		else if (orientationParadigm == OrientationParadigm.INVERSE_HOMO_CRW)
		{			
			orientationActuator = InverseHomogeneousCRW.instance;
			translationActuator = InverseHomogeneousCRW.instance;
		}
		else if (orientationParadigm == OrientationParadigm.BROWNIAN)
		{			
			orientationActuator = Brownian.instance;
			translationActuator = Brownian.instance;
		}
		else if (orientationParadigm == OrientationParadigm.BALLISTIC)
		{
			orientationActuator = Ballistic.instance;
			translationActuator = Ballistic.instance;
		}
		else if (orientationParadigm == OrientationParadigm.BOOTSTRAP)
		{
			Bootstrap bs = new Bootstrap();
			orientationActuator = bs;
			translationActuator = bs;
		}
		else if (orientationParadigm == OrientationParadigm.BOOTSTRAP_GRADIENT)
		{
			BootstrapGradient bs = new BootstrapGradient();
			orientationActuator = bs;
			translationActuator = bs;
		}
		else if (orientationParadigm == OrientationParadigm.BOOTSTRAP_DYNAMIC)
		{
			BootstrapDynamic bs = new BootstrapDynamic();
			orientationActuator = bs;
			translationActuator = bs;
		}
		else if (orientationParadigm == OrientationParadigm.HOMO_CORKSCREW_CRW)
		{
			translationActuator = HomogeneousCRW_Translation.instance;
			orientationActuator = HomogeneousCRWCorkscrew_Orientation.instance;
		}
	}
	
	public Double3D getChemotacticGradient()
	{	return null;	}
	
	public abstract ConcentrationPerception perceiveAttractant();
	
	
	/** Apply some rotation to the cell's current orientation in response to it having collided.
	 */
	protected void bounce()
	{	
		/* apply some rotation to the cell's current orientation in response to it having collided. Cells only
		 * move in the direction they are facing, so for cells to slide around obstacles, their orientation
		 * must be altered. The bounce represents the normal of the collision(s) that the cell had. A new
		 * orientation is selected that is part way between the desired heading, and the bounce heading. 
		 * 
		 * Note on application of quaternions. The calculation of the rotation is performed in absolute space (not
		 * relative to the cell). This is because the bounce vector is expressed in absolute space. Hence, the cell's 
		 * heading must be converted to absolute space, the angle between the two is calculated, and the rotation 
		 * is applied around the normal of the two vectors. */
		if (bounce.lengthSq() != 0.0)		
		{	
			switch(collisionResponse)
			{
			case RETURN:	collisionResponseReturn(); break;
			case NORMAL:	collisionResponseNormal(); break;
			case SNOOKER:	collisionResponseSnooker(); break;
			case SLIDE:		collisionResponseSlide(); break;
			}
		} 		
	}	
	
	private void collisionResponseReturn()
	{
		Double3D facing = orientation.transform(x_axis);
		orientation = Quaternion.faceVector(new Double3D(-facing.x, -facing.y, -facing.z));		
	}
	private void collisionResponseNormal()
	{
		orientation =  Quaternion.faceVector(bounce);
	}
	
	/** T cell persists in it's direction, but will slide around the obstacle cell */
	private void collisionResponseSlide()
	{
		/* apply some rotation to the cell's current orientation in response to it having collided. Cells only
		 * move in the direction they are facing, so for cells to slide around obstacles, their orientation
		 * must be altered. The bounce represents the normal of the collision(s) that the cell had. A new
		 * orientation is selected that is part way between the desired headin, and the bounce heading. 
		 * 
		 * Note on application of quaternions. The calculation of the rotation is performed in absolute space (not
		 * relative to the cell). This is because the bounce vector is expressed in absolute space. Hence, the cell's 
		 * heading must be converted to absolute space, the angle between the two is calculated, and the rotation 
		 * is applied around the normal of the two vectors. */
		Double3D facing = orientation.transform(x_axis);
		// the normal is perpendicular to the plane on which bounce and facing vectors lie. 
		Double3D normal = Utils.crossProduct(facing, bounce);		
		// this calculation takes vector direction into account. Ie, because facing points into the object that			
		// bounce points out of, 90 <= angle <= 180 (eqivalen in radians)
		double angle = Utils.angleBetweenVectors(facing, bounce);	// in radians.
		// want to slide along the object, so subtract a little less than 90 degrees. 90 exactly still results in
		// some collision, and cells grind along their obstacles more slowly. This is a balance between the two. 
		angle -= (Math.PI / 2.0) * 0.9;			 			
		// rotation applied to orientation to represent the cell bouncing off a collision.
		Quaternion bounceRot = Quaternion.representRotation(angle, normal.x, normal.y, normal.z);			
		/* apply the rotation. Note the reverse order of mutliplication. Mutliply orientation by the bounce
		 * rotation. This is because the bounce rotation is calculated in absolute space, rather than relative to 
		 * the cell. The bounce resulting from collision is given in absolute space, and the cell's heading
		 * in absolute space must be calculated before the rotation can be derived.  */ 
		Quaternion slideOrientation = bounceRot.multiply(orientation);	
		orientation = slideOrientation;
	}
	
	/**
	 * Alters cell's orientation such that it bounces off the contact cell as a snooker ball would hitting another 
	 * ball. 
	 */
	private void collisionResponseSnooker()
	{
		/* apply some rotation to the cell's current orientation in response to it having collided. Cells only
		 * move in the direction they are facing, so for cells to slide around obstacles, their orientation
		 * must be altered. The bounce represents the normal of the collision(s) that the cell had. 
		 * 
		 * Note on application of quaternions. The calculation of the rotation is performed in absolute space (not
		 * relative to the cell). This is because the bounce vector is expressed in absolute space. Hence, the cell's 
		 * heading must be converted to absolute space, the angle between the two is calculated, and the rotation 
		 * is applied around the normal of the two vectors. */
		Double3D facing = orientation.transform(x_axis);
		// the normal is perpendicular to the plane on which bounce and facing vectors lie. 
		Double3D normal = Utils.crossProduct(facing, bounce);		
		// this calculation takes vector direction into account. Ie, because facing points into the object that			
		// bounce points out of, 90 <= angle <= 180 (eqivalen in radians)
		double angle = Utils.angleBetweenVectors(facing, bounce);	// in radians.
		// 
		// the angle minus 90 degrees will allow the cell to 'slide' (persist) along the contacted obstacle. To have a
		// snooker ball reaction, double this angle.
		angle -= (Math.PI / 2.0); // subtract 90 degrees
		angle *= 2.0;
		// rotation applied to orientation to represent the cell bouncing off a collision.
		Quaternion bounceRot = Quaternion.representRotation(angle, normal.x, normal.y, normal.z);			
		/* apply the rotation. Note the reverse order of mutliplication. Mutliply orientation by the bounce
		 * rotation. This is because the bounce rotation is calculated in absolute space, rather than relative to 
		 * the cell. The bounce resulting from collision is given in absolute space, and the cell's heading
		 * in absolute space must be calculated before the rotation can be derived.  */ 
		orientation = bounceRot.multiply(orientation);			
	}	
	
	public double[] eulersAngles()
	{	return orientation.toEulerAngles();		}
	
	/** 
	 * Rotates the cell to face the supplied vector, which is expressed relative to the cell's current orientation. 
	 * Rotation is limited by the maximum rate of rotation, hence the rotation may not be achieved in one call. 
	 *  
	 * Method offers safety against vectors of zero length, which if supplied result in no action. 
	 */
	protected void rotateToVectorLimited(Double3D direction, double perceptionCone)
	{
		if (direction.lengthSq() != 0.0)
		{
			// Vector is normal to the plane between orientation of cell, and direction of gradient relative to cell.
			Double3D normal = Utils.unitVector(Utils.crossProduct(x_axis, direction));
			if (normal.lengthSq() == 0)
				// Can happen for very small gradients; precision errors. Can't proceed if no plane exists. 
				return;
			// Angle between orientation of cell and gradient direction relative to cell. 
			double angle = Utils.angleBetweenVectors(x_axis, direction);  // always positive, in radians. 
			// If the angle is zero, then there is no rotation to be applied: the two vectors are parallel. 
			if (angle == 0.0)	return;
			
			if (angle > perceptionCone)  // Limit based on a cell's angle of perception. In radians. 
				angle = perceptionCone;
			Quaternion rot = Quaternion.representRotation(angle, normal.x, normal.y, normal.z);
			// Multiply orientation by rot, because rot is given relative to cell, not in absolute space. 			
			orientation = orientation.multiply(rot);
		}
	}
	

	public static void loadParameters(Document params) throws XPathExpressionException, IOException
	{	
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		/* Set up motility paradigm */
		n = (Node) xPath.compile("/params/Motility/paradigm").evaluate(params, XPathConstants.NODE);				
		String selection = n.getTextContent();
		if (selection.equals("HomogeneousCRW"))
		{
			translationParadigm = TranslationParadigm.HOMO_CRW;
			HomogeneousCRW_Translation.loadParameters(params);
			orientationParadigm = OrientationParadigm.HOMO_CRW;
			HomogeneousCRW_Orientation.loadParameters(params);
		}
		else if (selection.equals("HeterogeneousCRW"))
		{
			translationParadigm = TranslationParadigm.HETERO_CRW;
			HeterogeneousCRW_Translation.loadParameters(params);
			orientationParadigm = OrientationParadigm.HETERO_CRW;
			HeterogeneousCRW_Orientation.loadParameters(params);
		}
		else if (selection.equals("LevyFlight"))
		{
			translationParadigm = TranslationParadigm.LEVY_VPL;
			LevyFlight.loadParameters(params);		
			orientationParadigm = OrientationParadigm.LEVY_VPL;
			LevyFlight.loadParameters(params);
		}
		else if (selection.equals("InverseHeterogeneousCRW"))
		{
			translationParadigm = TranslationParadigm.INVERSE_HETERO_CRW;
			orientationParadigm = OrientationParadigm.INVERSE_HETERO_CRW;			
			InverseHeterogeneousCRW.loadParameters(params);
		}
		else if (selection.equals("InverseHomogeneousCRW"))
		{
			translationParadigm = TranslationParadigm.INVERSE_HOMO_CRW;
			orientationParadigm = OrientationParadigm.INVERSE_HOMO_CRW;
			InverseHomogeneousCRW.loadParameters(params);
		}
		else if (selection.equals("Brownian"))
		{
			translationParadigm = TranslationParadigm.BROWNIAN;
			orientationParadigm = OrientationParadigm.BROWNIAN;
			Brownian.loadParameters(params);
		}
		else if (selection.equals("Ballistic"))
		{
			translationParadigm = TranslationParadigm.BALLISTIC;
			orientationParadigm = OrientationParadigm.BALLISTIC;
			Ballistic.loadParameters(params);			
		}
		else if (selection.equals("Bootstrap"))
		{
			translationParadigm = TranslationParadigm.BOOTSTRAP;
			orientationParadigm = OrientationParadigm.BOOTSTRAP;
			Bootstrap.loadParameters(params);
		}
		else if (selection.equals("BootstrapGradient"))
		{
			translationParadigm = TranslationParadigm.BOOTSTRAP_GRADIENT;
			orientationParadigm = OrientationParadigm.BOOTSTRAP_GRADIENT;
			BootstrapGradient.loadParameters(params);
		}
		else if (selection.equals("BootstrapDynamic"))
		{
			translationParadigm = TranslationParadigm.BOOTSTRAP_DYNAMIC;
			orientationParadigm = OrientationParadigm.BOOTSTRAP_DYNAMIC;
			Bootstrap.loadParameters(params);
			BootstrapGradient.loadParameters(params);
		}
		else if  (selection.equals("Corkscrew"))
		{
			translationParadigm = TranslationParadigm.HOMO_CRW;
			orientationParadigm = OrientationParadigm.HOMO_CORKSCREW_CRW;
			HomogeneousCRW_Translation.loadParameters(params);
			HomogeneousCRW_Translation.speedMean = 5.0;
			HomogeneousCRW_Translation.speedStD = 0.0;			
			HomogeneousCRW_Orientation.loadParameters(params);
			HomogeneousCRW_Orientation.rollRateStd = 0.0;
		}
		else throw new RuntimeException("Must select an appropriate translation paradigm.");
	}
}
