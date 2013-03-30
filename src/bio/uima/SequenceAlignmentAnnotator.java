package bio.uima;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;

/*
 * This annotator should compute a minimum cost alignment
 */
public class SequenceAlignmentAnnotator extends JCasAnnotator_ImplBase  {

	@Override
	public void process(JCas cas) throws AnalysisEngineProcessException {
		// read protein sequences from the CAS 
		try {
			String[] seqs = cas.getView("orf1").getDocumentText().split(" ");
			System.out.println(seqs[0]);
			System.out.println(seqs[1]);
			System.out.println();
			// TODO: implement this	
			String[] alignment = computeLevenshteinDistance(seqs[0], seqs[1]);
			JCas alignmentCas = cas.createView("alignment");
			alignmentCas.setDocumentText(alignment[0] + " " + alignment[1]);
		} catch (CASException e) {
			e.printStackTrace();
		}
	}

	// This method should compute a minimum cost alignment
	// for the specified sequences.  The cost of an insertion
	// or delete is 1 and the cost of substitution is 2.
	// Return marked up strings where insertions/deletes are
	// represented with a '-'.  
	//
	// For example an alignment of:
	//
	// 'AAGT' and 'AGT' 
	//
	// could be 
	//
	// {'AAGT', 'A-GT'}
	//
	
	private String[] getAlignment(int[][] distances, char[][] directions, String seq1, String seq2) {
		StringBuilder myStr = new StringBuilder(seq1);
		int iii = seq1.length();
		int jjj = seq2.length();
		while(iii > 0 && jjj > 0) {
			char currentDirection = directions[iii - 1][jjj - 1];
			if(currentDirection == 'd') {
				myStr.setCharAt(iii - 1, seq2.charAt(jjj - 1));
				iii--;
				jjj--;
			} else if (currentDirection == 'u') {
				myStr.setCharAt(iii - 1, '-');
				iii--;
			} else {
//				myStr.setCharAt(iii - 1, '-');
				jjj--;
			}
		}
		while(iii > 0) {
			myStr.setCharAt(--iii, '-');
		}
		String[] result = {seq1, myStr.toString()};
		return result;
	}
	
	private String[] computeLevenshteinDistance(String seq1, String seq2) {
		// enforce that seq1 is longer than seq2
		// if not, swap seq1 and seq2
		// for ease of computation
		if(seq1.length() < seq2.length()) {
			String temp = seq2;
			seq2 = seq1;
			seq1 = temp;
		}
		
		int[][] distances = new int[seq1.length() + 1][seq2.length() + 1];		// stores minimum cost
		char[][] directions = new char[seq1.length()][seq2.length()];			// stores direction of previous minimum cost
		
		// Initialize value
		distances[0][0] = 0;
		for(int iii = 1; iii <= seq1.length(); iii++){
			distances[iii][0] = iii * 2;
		}
		for(int jjj = 1; jjj <= seq2.length(); jjj++) {
			distances[0][jjj] = jjj * 2;
		}
		
		for (int iii = 1; iii <= seq1.length(); iii++) {
			for (int jjj = 1; jjj <= seq2.length(); jjj++) {
				int diagonal = distances[iii - 1][jjj - 1];
				// Add cost 1 for substitution
				if (seq1.charAt(iii - 1) != seq2.charAt(jjj - 1)) {
					diagonal++;
				}
				if (diagonal <= Math.min(distances[iii - 1][jjj] + 2, distances[iii][jjj - 1] + 2)) {
					distances[iii][jjj] = diagonal;
					directions[iii - 1][jjj - 1] = 'd';
				} else if (distances[iii - 1][jjj] <= distances[iii][jjj - 1]) {
					distances[iii][jjj] = distances[iii - 1][jjj] + 2;		// add 2 to cost for insertion
					directions[iii - 1][jjj - 1] = 'u';
				} else {
					distances[iii][jjj] = distances[iii][jjj - 1] + 2;		// add 2 to cost for deletion
					directions[iii - 1][jjj - 1] = 'l';
				}
			}
		}
		return getAlignment(distances, directions, seq1, seq2);
	}
}
