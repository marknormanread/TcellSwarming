package portrayal;

import java.awt.Color;

import javax.media.j3d.TransformGroup;

import core.Neutrophil;
import core.Simulation;
import sim.portrayal3d.simple.SpherePortrayal3D;
import sim.util.Double3D;

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
public class NeutrophilPortrayal extends SpherePortrayal3D
{
	// naive = green. 
	private static Color neutrophilColor = new Color(0.0f, 0.8f, 0.0f, 1.0f);
	// attractant  = blue
	private static Color attractantStimulatedColor = new Color(0.0f, 0.0f, 0.8f, 1.0f);
	// ltb4 simulated = red
	private static Color ltb4StimulatedColor = new Color(0.8f, 0.0f, 0.0f, 1.0f);
	private static Color inSwarmColor = new Color(0.5f, 0.5f, 0.0f, 1.0f);
	
	public NeutrophilPortrayal()
	{
		super(neutrophilColor, Neutrophil.diameter);
	}
	

	/**
	 * Overridden to provide a custom portrayal for neutrophils. 
	 */
	public TransformGroup getModel(Object obj, TransformGroup j3dModel)
	{
		Neutrophil neutrophil = (Neutrophil)obj;

		return colorByState(neutrophil, j3dModel);
	}
		
	
	/**
	 * Colours neutorphils in accordance to their state of activation/recognition of chemokine factors. 
	 */
	private TransformGroup colorByState(Neutrophil neutrophil, TransformGroup j3dModel)
	{
		if (neutrophil.isStateRandomWalk())
		{
			setAppearance(j3dModel, appearanceForColors(
					neutrophilColor, 		// ambient color
					null, 					// emissive color (black)
					neutrophilColor, 		// diffuse color
					null, 					// specular color (white)
					1.0f, 					// shininess, none. 
					1.0f));					// opacity
		}
		else if (neutrophil.isStateAttractantStimulated())
		{
			setAppearance(j3dModel, appearanceForColors(
					attractantStimulatedColor, 		// ambient color
					null, 							// emissive color (black)
					attractantStimulatedColor, 		// diffuse color
					null, 							// specular color (white)
					1.0f, 							// shininess, none. 
					1.0f));							// opacity
		}
		else if (neutrophil.isStateLTB4Stimulated())
		{
			setAppearance(j3dModel, appearanceForColors(
					ltb4StimulatedColor, 		// ambient color
					null, 						// emissive color (black)
					ltb4StimulatedColor, 		// diffuse color
					null, 						// specular color (white)
					1.0f, 						// shininess, none. 
					1.0f));						// opacity
		}
		if (neutrophil.isInSwarm())
		{
			setAppearance(j3dModel, appearanceForColors(
					inSwarmColor, 				// ambient color
					null, 						// emissive color (black)
					inSwarmColor, 				// diffuse color
					null, 						// specular color (white)
					1.0f, 						// shininess, none. 
					1.0f));						// opacity			
		}
		if (SimulationGUI.imageOutsideVolume == false && outsideImagingVolume(neutrophil))
		{
			setAppearance(j3dModel, appearanceForColors(
					neutrophilColor, 			// ambient color
					null, 						// emissive color (black)
					neutrophilColor, 			// diffuse color
					null, 						// specular color (white)
					1.0f, 						// shininess, none. 
					0.0f));						// opacity					
		}
		
		return super.getModel(neutrophil, j3dModel);		
	}
	
	private boolean outsideImagingVolume(Neutrophil n)
	{
		Double3D loc = n.getCurrentLocation();
		return !Simulation.space.insideImagingVolume(loc.x, loc.y, loc.z);
	}
	
}
