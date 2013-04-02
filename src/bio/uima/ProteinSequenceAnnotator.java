package bio.uima;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;

import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.ArrayList;

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
//		String[] proteins = new String[5];
		ArrayList<String> proteins = new ArrayList<String>();
		HashMap<String, String> codonTable = ProteinSequenceAnnotator.initMap();

		// Find DNA sequences
		int pos = 0;
		while(matcher.find(pos)) {
			DNASequence annotation = new DNASequence(cas);
			annotation.setBegin(matcher.start());
			annotation.setEnd(matcher.end());
			annotation.setValue(sequencesString.substring(matcher.start(), matcher.end()));
//			proteins[iii++] = dnaToProtein(sequencesString.substring(matcher.start(), matcher.end()), codonTable);
			String[] newProteins = dnaToProtein(sequencesString.substring(matcher.start(), matcher.end()), codonTable);
			addNewProteins(proteins, newProteins);
			annotation.addToIndexes();
			pos = matcher.end();
		}
		// TODO: implement this
		// translate each DNA sequence to a protein sequence considering
		// all possible open reading frames and store in a CAS view
		try {
			JCas orf1 = cas.createView("protein");
			System.out.println(proteins.size());
			String proteinString = ProteinSequenceAnnotator.combineStringArray(proteins.toArray(new String[proteins.size()]));
			orf1.setDocumentText(proteinString);
		} catch (CASException e) {
			e.printStackTrace();
		}
	}
	
	private void addNewProteins(ArrayList<String> proteins, String[] newProteins) {
		for(String temp : newProteins) {
			proteins.add(temp);
		}
	}

	public static String combineStringArray(String[] proteins) {
		String returnString = "";
		for(int iii = 0; iii < proteins.length; iii++) {
			returnString += (proteins[iii] + " ");
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
	
	private String[] dnaToProtein(String seq, HashMap<String, String> codonTable) {
		// TODO: consider all possible open reading frames
    	final int NUM_FRAMES = 3;
    	String[] proteins = new String[NUM_FRAMES];
		int length = seq.length();
		for(int frame = 0; frame < NUM_FRAMES; frame++) {
			StringBuilder proteinBuilder = new StringBuilder();
			for(int iii = frame; iii <= (length + frame - NUM_FRAMES); iii+=3) {
				StringBuilder codon = new StringBuilder();
				codon.append(seq.charAt(iii));
				codon.append(seq.charAt((iii + 1) % length));
				codon.append(seq.charAt((iii + 2) % length));
				proteinBuilder.append(codonTable.get(codon.toString()));
			}
			if(frame > 0)	proteinBuilder.deleteCharAt(proteinBuilder.length() - 1);
			proteins[frame] = proteinBuilder.toString();
		}
		return proteins;
	}
	
//	private String dnaToProtein(String seq, HashMap<String, String> codonTable) {
//		// TODO: consider all possible open reading frames
//		String proteinSeq = "";
//		for(int iii = 0; iii <= (seq.length() - 3); iii+=3) {
//			String currentCodon = seq.substring(iii, iii + 3);
//			proteinSeq += codonTable.get(currentCodon);
//		}
//		return proteinSeq;
//	}

	/*
	 * @return HashMap<String, String> 	: hash map to map codon to protein
	 * Reads codon table from a file and creates a hash map.
	 */
	private static HashMap<String, String> initMap() {
		final String codonFile = "/home/farhang/workspace/BioUIMA/resources/codonTable.txt";
    	String inputFileString = BioUima.readInputFile(codonFile);
    	HashMap<String, String> codonMap = new HashMap<String, String>();
    	String[] lines = inputFileString.split("\n");
    	for(String line : lines) {
    		String[] temp = line.split("\t");
    		codonMap.put(temp[0], temp[1]);
    	}
    	return codonMap;
    }
}
