package bio.uima;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;

import java.util.HashMap;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/*
 * This annotator should read DNA sequences from the CAS
 * and translate each sequence to protein sequences.  
 * 
 * See http://en.wikipedia.org/wiki/DNA_codon_table
 *     http://en.wikipedia.org/wiki/Open_reading_frame
 */
public class ProteinSequenceAnnotator extends JCasAnnotator_ImplBase  {

	@Override
	public void process(JCas cas) throws AnalysisEngineProcessException {
		// read DNA sequences from the CAS
		Pattern DNAPattern = Pattern.compile("\\b[ACGT]+\\b");
		String sequencesString = cas.getDocumentText();
		Matcher matcher = DNAPattern.matcher(sequencesString);
		String[] proteins = new String[5];
		HashMap<String, String> codonTable = ProteinSequenceAnnotator.initMap();
		
		// Find DNA sequences, create annotations, and translate to protein sequences
		int iii = 0;
		int pos = 0;
		while(matcher.find(pos)) {
			DNASequence annotation = new DNASequence(cas);
			annotation.setBegin(matcher.start());
			annotation.setEnd(matcher.end());
			annotation.setValue(sequencesString.substring(matcher.start(), matcher.end()));
			proteins[iii++] = dnaToProtein(sequencesString.substring(matcher.start(), matcher.end()), codonTable);
			annotation.addToIndexes();
			pos = matcher.end();
		}
		
		// TODO: implement this
		// translate each DNA sequence to a protein sequence considering
		// all possible open reading frames and store in a CAS view
		try {
			JCas orf1 = cas.createView("protein");
			String proteinString = ProteinSequenceAnnotator.combineStringArray(proteins);
			orf1.setDocumentText(proteinString);
		} catch (CASException e) {
			e.printStackTrace();
		}
	}
	
	public static String combineStringArray(String[] array) {
		String returnString = "";
		for(int iii = 0; iii < array.length; iii++) {
			returnString += (array[iii] + " ");
		}
		return returnString;
	}

	// This method should generate the protein sequence for the specified
	// DNA sequence.  You man use regex.
	//
	// A codon is just a sequence of three characters.  Each codon maps to 
	// a protein.  For example the DNA codon 'TGT' maps to 'C' (cysteine).
	// So the sequence 'TGTTGTTGTTGT' should translate to 'CCCC'.
	//
	// See: http://en.wikipedia.org/wiki/DNA_codon_table
	private String dnaToProtein(String seq, HashMap<String, String> codonTable) {
		// TODO: consider all possible open reading frames
		String proteinSeq = "";
		for(int iii = 0; iii <= (seq.length() - 3); iii+=3) {
			String currentCodon = seq.substring(iii, iii + 3);
			proteinSeq += codonTable.get(currentCodon);
		}
		return proteinSeq;
	}

	/*
	 * @return HashMap<String, String> 	: hash map to map codon to protein
	 * Reads codon table from a file and creates a hash map.
	 */
	private static HashMap<String, String> initMap() {
		HashMap<String, String> codonMap = new HashMap<String, String>();
		final String codonFile = "/home/farhang/workspace/BioUIMA/resources/codonTable.txt";
		BufferedReader br = null;
		try {
			String currentLine;
			br = new BufferedReader(new FileReader(codonFile));
			while ((currentLine = br.readLine()) != null) {
				String[] array = currentLine.split("\t");
				codonMap.put(array[0], array[1]);
			}
		} catch	(IOException e) {
			e.printStackTrace();
			System.out.println("File not found: " + codonFile);
		} finally {
				try {
					if(br != null)	br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
		}
		return codonMap;
	}
}
