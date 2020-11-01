package portrayal;

import java.awt.Color;

import javax.media.j3d.TransformGroup;

import core.Simulation;
import core.TCell;
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
public class TCellPortrayal extends SpherePortrayal3D
{
//	private static Color tcellColor            = new Color(0.0f, 1.0f, 0.0f, 1.0f);  // Green.
//	private static Color tcellInvisibleColor   = new Color(0.0f, 1.0f, 0.0f, 0.0f);   
//	private static Color tcellChemotacticColor = new Color(1.0f, 0.0f, 1.0f, 1.0f);  // Purple. 
//	private static Color tcellBolus            = new Color(1.0f, 1.0f, 0.0f, 1.0f);

	
	// This looks awesome on black background, but lacks contrast. On white Grey isn't not as good. 
//	private static Color tcellColor            = new Color(0.0f, 0.0f, 0.75f, 1.0f);  // 	   
//	private static Color tcellChemotacticColor = new Color(0.0f, 1.0f, 1.0f, 1.0f);  //
//	private static Color tcellBolus            = new Color(1.0f, 0.7f, 0.0f, 1.0f);  // Orange	
//	private static Color tcellInvisibleColor   = new Color(0.0f, 1.0f, 0.0f, 0.0f);
	
//	private static Color tcellColor            = new Color(0.0f, 0.0f, 0.75f, 1.0f);  // 	   
//	private static Color tcellChemotacticColor = new Color(0.0f, 0.8f, 0.8f, 1.0f);  //
//	private static Color tcellBolus            = new Color(1.0f, 0.7f, 0.0f, 1.0f);  // Orange	
//	private static Color tcellInvisibleColor   = new Color(0.0f, 1.0f, 0.0f, 0.0f);	
	
	private static Color tcellColor            = new Color(0.f, 0.8f, 0.0f, 0.f); //
	private static Color tcellChemotacticColor = hex2Rgb("#8a37ff");  // Pink 	
	private static Color tcellBolus            = hex2Rgb("#ff37ac");  // Purple	
	private static Color tcellInvisibleColor   = new Color(0.0f, 1.0f, 0.0f, 0.0f);
//	
	
	public static Color hex2Rgb(String colorStr) {
	    return new Color(
	            Integer.valueOf( colorStr.substring( 1, 3 ), 16 ),
	            Integer.valueOf( colorStr.substring( 3, 5 ), 16 ),
	            Integer.valueOf( colorStr.substring( 5, 7 ), 16 ) );
	}
	
	public TCellPortrayal()
	{		
		super(tcellColor, SimulationGUI.cellScalar * TCell.diameter);
	}
	
	/**
	 * Overridden to provide a custom portrayal for T cells. 
	 */
	public TransformGroup getModel(Object obj, TransformGroup j3dModel)
	{
		TCell tcell = (TCell)obj;

		return colorByState(tcell, j3dModel);
	}
	
	/**
	 * Colours cells in accordance to their state of activation/recognition of chemokine factors. 
	 */
	private TransformGroup colorByState(TCell tcell, TransformGroup j3dModel)
	{
		// Can add conditional queries on T cell state here, and change color with state.
		Color col = tcellColor;
		
		if ( /* If cell outside of imaging volume, and user has asked such cells not to be displayed */
			 (SimulationGUI.imageOutsideVolume == false && outsideImagingVolume(tcell))
			 /* If only a subset of the imaging volume (slice) should actually be displayed, and this cell qualifies */
			 || (SimulationGUI.imageOutsideRecorded == SimulationGUI.RecordingVolumeMode.SLICE 
				 && outsideRecordingVolume(tcell)) 
		     /* If for whatever reason this cell should not be displayed */
			 || tcell.recordable == false)
		{	// Make cell completely transparent. 
			setAppearance(j3dModel, appearanceForColors(
					col,   // Ambient color
					null,  // Emissive color (black)
					col,   // Diffuse color 
					null,  // Specular color (white)
					1.0f,  // Shininess, none. 
					0.0f));  // Opacity					
		} else {
			if (Simulation.space.insideBolus(tcell.getCurrentLocation()))
				col = tcellBolus;
			else if (tcell.isChemotactic())
				col = tcellChemotacticColor;
			
			
			setAppearance(j3dModel, appearanceForColors(
					col,   // Ambient color
					null,  // Emissive color (black)
					col,   // Diffuse color
					null,  // Specular color (white)
					1.0f,  // Shininess, none. 
					1.0f));  // Opacity
		}
		return super.getModel(tcell, j3dModel);		
	}
	
	private boolean outsideRecordingVolume(TCell tcell)
	{
		Double3D loc = tcell.getCurrentLocation();
		return !Simulation.space.insideRecordingVolume(loc.x, loc.y, loc.z, tcell);
	}
	
	private boolean outsideImagingVolume(TCell n)
	{
		Double3D loc = n.getCurrentLocation();
		return !Simulation.space.insideImagingVolume(loc.x, loc.y, loc.z);
	}
}
