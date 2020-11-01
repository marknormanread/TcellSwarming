package portrayal;

import java.awt.Color;

import javax.media.j3d.TransformGroup;

import core.ChemotacticAgent;
import core.Simulation;
import core.TCell;
import sim.portrayal3d.simple.SpherePortrayal3D;
import sim.util.Double3D;

public class ChemotacticAgentPortrayal extends SpherePortrayal3D
{
	private static Color cellColor            = new Color(0.0f, 1.0f, 0.0f, 1.0f);  // Green.
	private static Color cellInvisibleColor   = new Color(0.0f, 1.0f, 0.0f, 0.0f);   
	private static Color cellChemotacticColor = new Color(1.0f, 0.0f, 1.0f, 1.0f);  // Purple. 
	private static Color cellBolus            = new Color(1.0f, 1.0f, 0.0f, 1.0f);
	
	public ChemotacticAgentPortrayal()
	{		
		super(cellColor, SimulationGUI.cellScalar * ChemotacticAgent.diameter);
	}
	
	/**
	 * Overridden to provide a custom portrayal for T cells. 
	 */
	public TransformGroup getModel(Object obj, TransformGroup j3dModel)
	{
		ChemotacticAgent tcell = (ChemotacticAgent)obj;

		return colorByState(tcell, j3dModel);
	}
	
	/**
	 * Colours cells in accordance to their state of activation/recognition of chemokine factors. 
	 */
	private TransformGroup colorByState(ChemotacticAgent cell, TransformGroup j3dModel)
	{
		// Can add conditional queries on T cell state here, and change color with state.
		Color col = cellColor;
		
		if ( /* If cell outside of imaging volume, and user has asked such cells not to be displayed */
			 (SimulationGUI.imageOutsideVolume == false && outsideImagingVolume(cell))
			 /* If only a subset of the imaging volume (slice) should actually be displayed, and this cell qualifies */
			 || (SimulationGUI.imageOutsideRecorded == SimulationGUI.RecordingVolumeMode.SLICE 
				 && outsideRecordingVolume(cell)) 
		     /* If for whatever reason this cell should not be displayed */
			 || cell.recordable == false)
		{	// Make cell completely transparent. 
			setAppearance(j3dModel, appearanceForColors(
					col,   // Ambient color
					null,  // Emissive color (black)
					col,   // Diffuse color
					null,  // Specular color (white)
					1.0f,  // Shininess, none. 
					0.0f));  // Opacity					
		} else {
			if (Simulation.space.insideBolus(cell.getCurrentLocation()))
				col = cellBolus;
			
			else if (cell.isChemotactic())
				col = cellChemotacticColor;
			
			setAppearance(j3dModel, appearanceForColors(
					col,   // Ambient color
					null,  // Emissive color (black)
					col,   // Diffuse color
					null,  // Specular color (white)
					1.0f,  // Shininess, none. 
					1.0f));  // Opacity
		}
		return super.getModel(cell, j3dModel);		
	}
	
	private boolean outsideRecordingVolume(ChemotacticAgent cell)
	{
		Double3D loc = cell.getCurrentLocation();
		return !Simulation.space.insideRecordingVolume(loc.x, loc.y, loc.z, cell);
	}
	
	private boolean outsideImagingVolume(ChemotacticAgent cell)
	{
		Double3D loc = cell.getCurrentLocation();
		return !Simulation.space.insideImagingVolume(loc.x, loc.y, loc.z);
	}
}
