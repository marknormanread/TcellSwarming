package filesystem;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpression;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Comment;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

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
public class FileSystemIO 
{
	
	/**
	 * Opens an XML file from the filesystem as a Document, so that information can be extracted from it.  
	 * 
	 * @param file - the location of the file, as a string. 
	 */
	public static Document openXMLFile(String location)
	{
		DocumentBuilderFactory docBuilderFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder docBuilder;
		
		Document doc = null;
		
		try {
			docBuilder = docBuilderFactory.newDocumentBuilder();
			doc = docBuilder.parse(new File(location));

			//	normalize text representation
			doc.getDocumentElement().normalize();
			
			// remove any nodes in the DOM document that contain only whitespace. 
			removeWhitespaceNodes(doc);
			removeCommentNodes(doc);
			
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
		} catch (SAXException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return doc;
	}
	
	/**
	 * Code obtained (and slightly modified, but not much) from http://www.roseindia.net/java/beginners/CopyFile.shtml
	 * @param source
	 * @param dest
	 */
	public static void copyFile(String source, String dest){
	    System.out.println("copying file " + source + " to " + dest);
		try{
	      File sf = new File(source);
	      File df = new File(dest);
	      InputStream in = new FileInputStream(sf);
	      
	      OutputStream out = new FileOutputStream(df);

	      byte[] buf = new byte[1024];
	      int len;
	      while ((len = in.read(buf)) > 0){
	        out.write(buf, 0, len);
	      }
	      in.close();
	      out.close();
	      System.out.println("File " + source + " copied to " + dest);
	    }
	    catch(FileNotFoundException ex){
	      System.out.println(ex.getMessage() + " in the specified directory.");
	      System.exit(0);
	    }
	    catch(IOException e){
	      System.out.println(e.getMessage());      
	    }
	  }
	
	/**
	 * Method takes an XML document represented as a DOM document and writes it to the file system. 
	 * 
	 * @param document the Document representation of the XML file
	 * @param location the directory path to the location where the XML file is to be writen to 
	 * @param filename the file name of the file that the XML is to be written in.
	 */
	public static void writeXMLFile(Document document, String location, String filename)
	{
		File directories = new File(location);
		if(directories.exists() == false)									// if the directory structure does not already exist, then create it. 
		{
			try{
				directories.mkdirs();
			} catch (Exception e) {
				System.out.println("FileSystemIO: Failed to create the directory structure in which to write the XML file. " + e.getStackTrace());
			}
		}
		File outputFile = new File(location + File.separator + filename);
		if(outputFile.exists())
			throw new RuntimeException("FileSystemIO: designated output file already exists!");
		

		try {
			TransformerFactory tFactory = TransformerFactory.newInstance();	// transformer converts source to an output. 
			Transformer transformer = tFactory.newTransformer();

			DOMSource source = new DOMSource(document); 					// the source for the transformation.
			StreamResult result = new StreamResult(outputFile); 			// where the transformation is to be sent to.
			transformer.transform(source, result); 							// perform the transformation from source to result. 

		} catch (TransformerConfigurationException tce) {
			// Error generated by the parser
			System.out.println("* Transformer Factory error");
			System.out.println("  " + tce.getMessage());

			// Use the contained exception, if any
			Throwable x = tce;
			if (tce.getException() != null)
				x = tce.getException();
			x.printStackTrace();
		} catch (TransformerException te) {
			// Error generated by the parser
			System.out.println("* Transformation error");
			System.out.println("  " + te.getMessage());

			// Use the contained exception, if any
			Throwable x = te;
			if (te.getException() != null)
				x = te.getException();
			x.printStackTrace();

		}
		
	}

	/**
	 * Code to remove nodes from the DOM tree that contain only whitespace charaters. Adapted from code obtained from 
	 * 
	 * http://stackoverflow.com/questions/978810/how-to-strip-whitespace-only-text-nodes-from-a-dom-before-serialization
	 * 
	 */
	private static Document removeWhitespaceNodes(Document doc)
	{
		XPathFactory xpathFactory = XPathFactory.newInstance();
		try {
			//	 XPath to find empty text nodes.
			XPathExpression xpathExp = xpathFactory.newXPath().compile(
			        "//text()[normalize-space(.) = '']");  
			NodeList emptyTextNodes = (NodeList) 
		        xpathExp.evaluate(doc, XPathConstants.NODESET);

			//	Remove each empty text node from document.
			for (int i = 0; i < emptyTextNodes.getLength(); i++) 
			{	
			    Node emptyTextNode = emptyTextNodes.item(i);
			    emptyTextNode.getParentNode().removeChild(emptyTextNode);
			}
		} catch (Exception e)
		{ 
			System.out.println("exception raised whilst removing whitespace nodes from document : "+ e.getStackTrace());
		}
		
		return doc;
	}
	
	private static Document removeCommentNodes(Document doc)
	{
		Element e = doc.getDocumentElement();
		removeCommentNodesRecursion(e);
		
		return doc;
	}
	
	private static void removeCommentNodesRecursion(Node n)
	{
		if(n instanceof Comment)
		{
			n.getParentNode().removeChild(n);
		} else 	{
			NodeList ns = n.getChildNodes();
			if(ns == null) 	return;
			for(int i = 0; i < ns.getLength(); i++)
			{
				Node child = ns.item(i);
				removeCommentNodesRecursion(child);
			}
		}
	}
	
	public static void printDOMRecursive(int indent, Node n)
	{
		System.out.println();
		for(int i = 0; i < indent; i++)
			System.out.print("  ");
		
		System.out.print(n.getNodeName() + "  __  " + n.getTextContent());
		
		NodeList ns = n.getChildNodes();
		if (ns == null) return;
		for(int i = 0; i < ns.getLength(); i++)
		{
			Node child = ns.item(i);
			printDOMRecursive(indent + 1, child);
		}
	}
}
