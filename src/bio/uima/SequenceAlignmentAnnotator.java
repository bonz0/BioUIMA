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
			String[] proteins = cas.getView("protein").getDocumentText().split(" ");
			// TODO: implement this
			String[] alignment = getAllPairsAlignment(proteins);
			JCas alignmentCas = cas.createView("alignment");
			alignmentCas.setDocumentText(ProteinSequenceAnnotator.combineStringArray(alignment));
		} catch (CASException e) {
			e.printStackTrace();
		}
	}
	
	private String[] getAllPairsAlignment(String[] proteins) {
		String alignmentString = "";
		for(int iii = 0; iii < (proteins.length - 1); iii++) {
			for(int jjj = (iii + 1); jjj < proteins.length; jjj++) {
				String[] singlePairAlignment = computeLevenshteinDistance(proteins[iii], proteins[jjj]);
				alignmentString += (singlePairAlignment[0] + " " + singlePairAlignment[1]);
			}
		}
		return alignmentString.split(" ");
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
		
		// Initialize values
		int[][] distances = SequenceAlignmentAnnotator.initDistanceMatrix(seq1.length(), seq2.length());
//		int[][] distances = new int[seq1.length() + 1][seq2.length() + 1];		// stores minimum cost
		char[][] directions = new char[seq1.length()][seq2.length()];			// stores direction of previous minimum cost
			
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
	
	private static int[][] initDistanceMatrix(int length1, int length2) {
		int[][] distances = new int[length1 + 1][length2 + 1];
		distances[0][0] = 0;
		for(int iii = 1; iii <= length1; iii++){
			distances[iii][0] = iii * 2;
		}
		for(int jjj = 1; jjj <= length2; jjj++) {
			distances[0][jjj] = jjj * 2;
		}
		return distances;
	}
}