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
	
	private static final int NUM_FRAMES = 3;

	@Override
	public void process(JCas cas) throws AnalysisEngineProcessException {
		// read DNA sequences from the CAS
		Pattern DNAPattern = Pattern.compile("\\b[ACGT]+\\b");
		String sequencesString = cas.getDocumentText();
		Matcher matcher = DNAPattern.matcher(sequencesString);
		
		ArrayList<String> proteins = new ArrayList<String>();
		HashMap<String, String> codonTable = Utils.initCodonTable();		// Codon table for translation
		String[] inputDNA = new String[5];									// input DNA sequences
		
		/*
		 * Find DNA sequences in input file using regular expressions.
		 * For every match, create a new Annotation, store features, and
		 * translate to protein sequences.
		 */
		int pos = 0;
		int dnaNumber = 0;
		while(matcher.find(pos)) {
			DNASequence annotation = new DNASequence(cas);
			annotation.setBegin(matcher.start());
			annotation.setEnd(matcher.end());
			annotation.setValue(sequencesString.substring(matcher.start(), matcher.end()));
			inputDNA[dnaNumber++] = sequencesString.substring(matcher.start(), matcher.end());
			String[] newProteins = dnaToProtein(sequencesString.substring(matcher.start(), matcher.end()), codonTable);
			Utils.addArrayToList(proteins, newProteins);
			annotation.addToIndexes();
			pos = matcher.end();
		}
		
		try {
			JCas orf1 = cas.createView("protein");							// create a view for proteins
			String proteinString = Utils.combineStringArray(proteins.toArray(new String[proteins.size()]));
			orf1.setDocumentText(proteinString);
			JCas orf2 = cas.createView("dna");								// create separate view for DNA
			String dnaString = Utils.combineStringArray(inputDNA);
			orf2.setDocumentText(dnaString);
		} catch (CASException e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * Returns an array of Strings which represent all the open frame
	 * translations of a particular DNA sequence to protein sequences
	 * 
	 * @param	dnaSequence		The input DNA sequence to be translated to protein
	 * @param	codonTable		HashMap that stores the codon table for translation
	 * @return					All open frame translations of a DNA sequence
	 */
	private String[] dnaToProtein(String dnaSequence, HashMap<String, String> codonTable) {
    	String[] proteins = new String[NUM_FRAMES];
		int length = dnaSequence.length();
		for(int frame = 0; frame < NUM_FRAMES; frame++) {
			StringBuilder proteinBuilder = new StringBuilder();
			for(int iii = frame; iii <= (length + frame - NUM_FRAMES); iii+=3) {
				StringBuilder codon = new StringBuilder();
				codon.append(dnaSequence.charAt(iii));
				codon.append(dnaSequence.charAt((iii + 1) % length));
				codon.append(dnaSequence.charAt((iii + 2) % length));
				proteinBuilder.append(codonTable.get(codon.toString()));
			}
			if(frame > 0)	proteinBuilder.deleteCharAt(proteinBuilder.length() - 1);
			proteins[frame] = proteinBuilder.toString();
		}
		return proteins;
	}
}
