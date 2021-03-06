package bio.uima;

import org.apache.uima.UIMAFramework;
import org.apache.uima.analysis_engine.AnalysisEngine;
import org.apache.uima.analysis_engine.AnalysisEngineDescription;
import org.apache.uima.jcas.JCas;
import org.apache.uima.util.XMLInputSource;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class BioUima {

	/**
	 * The main entry point for the BioUIMA application.
	 * 
	 * TODO: 
	 * 
	 * Complete the implementation of the ProteinSequenceAnnotator.
	 * There are some helpful comments in that src file to get you 
	 * started.  Don't worry about the biology, you should be able to
	 * complete the implementation just from reading the comments.  If
	 * you need help there is plenty of information on the web about how 
	 * to translate DNA sequences to protein sequences.
	 * 
	 * The second annotator should read the protein sequences computed
	 * by the first annotator and compute a minimum cost alignment using
	 * the Levenshtein distance where the cost for insertion and deletion
	 * are 1 and the cost for substitution is 2.  Note that this is 
	 * equivalent to computing the minimum edit distance between two strings.
	 * 
	 * See: http://en.wikipedia.org/wiki/Levenshtein_distance
	 * 
	 * The DNA sequences you need to work with are in the 'data' directory.
	 * There are 5 DNA sequences in one file.  You are to perform the
	 * translation and alignment on all possible pairs of sequences.
	 * 
	 * We would like you to implement this using multi-Sofa, that is
	 * you should create multiple views on the CAS to distinguish DNA from
	 * protein sequences.  You can print results to STDOUT or 
	 * serialize the CAS.  Don't worry about upper/lower case.
	 * We'll be impressed if you use the UIMA type system, implement your
	 * own FlowController and take advantage of other UIMA concepts
	 * that are applicable or describe in detail how your application
	 * could take advantage of things like UIMA-AS.
	 * 
	 * If you have trouble with the first annotator, one alternative 
	 * would be to just compute a simple reversal of the original DNA 
	 * sequences and store those in a separate CAS view.  Perform an alignment 
	 * on these reversed strings rather than the protein strings.
     *
	 * @param args
	 */
	
	private static final String INPUT_XML_FILE = "desc/BioDescriptor.xml";
	private static final String INPUT_DNA_FILE = "/home/farhang/workspace/BioUIMA/data/dna.txt";
	
	public static void main(String[] args) throws Exception {
		// create the AE
		XMLInputSource input = new XMLInputSource(INPUT_XML_FILE);
		AnalysisEngineDescription desc = UIMAFramework.getXMLParser().parseAnalysisEngineDescription(input);
		AnalysisEngine ae = UIMAFramework.produceAnalysisEngine(desc);

		// init a CAS
		JCas jCas = ae.newJCas();
		String fileText = Utils.readInputFile(INPUT_DNA_FILE);
		jCas.setDocumentText(fileText);

		// process the CAS
		ae.process(jCas);

		// print results to stdout
		JCas clusterView = jCas.getView("clusters");
		JCas alignmentView = jCas.getView("alignments");
		String[] alignments = alignmentView.getDocumentText().split(" ");
		String[] clusters = clusterView.getDocumentText().split(" ");
		System.out.println("Cluster\tProtein alignment");
		Utils.printTwoStringArrays(clusters, alignments);
	}
}