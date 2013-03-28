package bio.uima;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;

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
		String[] seqs = cas.getDocumentText().split(" ");
			System.out.println();
			System.out.println("Protein Sequence Annotator -> process()");
		System.out.println(seqs[0]);
		System.out.println(seqs[1]);
		
		// TODO: implement this
		// translate each DNA sequence to a protein sequence considering
		// all possible open reading frames and store in a CAS view
		String orf1_protein1 = dnaToProtein(seqs[0]);
		String orf1_protein2 = dnaToProtein(seqs[1]); 
		try {
			JCas orf1 = cas.createView("orf1");
			orf1.setDocumentText(orf1_protein1+" "+orf1_protein2);
		} catch (CASException e) {
			e.printStackTrace();
		}
	}

	// This method should generate the protein sequence for the specified
	// DNA sequence.  You man use regex.
	//
	// A codon is just a sequence of three characters.  Each codon maps to 
	// a protein.  For example the DNA codon 'TGT' maps to 'C' (cysteine).
	// So the sequence 'TGTTGTTGTTGT' should translate to 'CCCC'.
	//
	// See: http://en.wikipedia.org/wiki/DNA_codon_table
	private String dnaToProtein(String seq) {
		// TODO: implement this
		//return new StringBuilder(seq).reverse().toString();
		return seq;
	}
}