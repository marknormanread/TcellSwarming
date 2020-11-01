package soluble;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

public abstract class LTB4 
{
	/*
	 * For the diffusion coefficient, small values result in less cytokine diffusion.  
	 * 
	 * The decay rate does not change the shape of the cyokine diffusion distribution, it simply
	 * adjusts the absolute quantities around at each point in space. Higher values result in
	 * higher cytokine concentrations.   
	 */	
	public static double diffusion = Double.NaN;
	// Convenient for human specification, but converted to decay constant when loading parameters. 
	public static double halflife = Double.NaN;  
	public static double decayConstant = Double.NaN;
	
	private LTB4()  // Cannot be instantiated
	{}

	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;	
		n = (Node) xPath.compile("/params/LTB4/halflife").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() != 0)  // Check parameter is present and set.
		{
			halflife = Double.parseDouble(n.getTextContent());
			if (Double.isInfinite(halflife))
				decayConstant = 0;  // Disables decay. 
			else
				decayConstant = Math.log(2) / halflife;  // Computationally convenient expression of half life
		}
		
		n = (Node) xPath.compile("/params/LTB4/diffusion").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() != 0)  // Check parameter is present and set.
			diffusion = Double.parseDouble(n.getTextContent());
	}	
	
}
