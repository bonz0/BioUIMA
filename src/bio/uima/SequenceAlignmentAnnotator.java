package bio.uima;

import java.util.HashMap;

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
	
	/*
	 * @param	proteins		: Array of strings
	 * @return Array of strings that contains alignment of all pairs of input String array
	 */
	private String[] getAllPairsAlignment(String[] proteins) {
		HashMap<Character, HashMap<Character, Integer>> distances = initDistancesMap();
		String alignmentString = "";
		for(int iii = 0; iii < (proteins.length - 1); iii++) {
			for(int jjj = (iii + 1); jjj < proteins.length; jjj++) {
				String[] singlePairAlignment = hirschberg(proteins[iii], proteins[jjj], distances);
				alignmentString += (singlePairAlignment[0] + " " + singlePairAlignment[1] + " ");
			}
		}
		return alignmentString.split(" ");
	}
	
	/*
	 * @param 	A			:	String to align
	 * @param	B			:	String to align
	 * @param	distances	:	HashMap that stores BLOSUM62 matrix
	 * @return	Array of two strings that are the optimal alignment of A and B
	 */
	private String[] getAlignment(String A, String B, HashMap<Character, HashMap<Character, Integer>> distances) {
		final int gapPenalty = -2;
		int[][] scores = computeNWScore(A, B, distances, gapPenalty);
		
		// Compute alignment
		String alignmentA = "";
		String alignmentB = "";
		int iii = A.length();
		int jjj = B.length();
		while (iii > 0 || jjj > 0)
		{
			if ((iii > 0 && jjj > 0) && (scores[iii][jjj] == scores[iii - 1][jjj - 1] + distances.get(A.charAt(iii - 1)).get(B.charAt(jjj - 1)))) {
				alignmentA = Character.toString(A.charAt(iii - 1)) + alignmentA;
				alignmentB = Character.toString(B.charAt(jjj - 1)) + alignmentB;
				iii--;
				jjj--;
			} else if ((iii > 0) && (scores[iii][jjj] == scores[iii-1][jjj] + gapPenalty)) {
				alignmentA = Character.toString(A.charAt(iii - 1)) + alignmentA;
				alignmentB = Character.toString('-') + alignmentB;
				iii--;
			} else if((jjj > 0) && (scores[iii][jjj] == scores[iii][jjj - 1] + gapPenalty)) {
				alignmentA = Character.toString('-') + alignmentA;
				alignmentB = Character.toString(B.charAt(jjj - 1)) + alignmentB;
				jjj--;
			}
		}
		String[] alignment = {alignmentA, alignmentB};
		return alignment;
	}
	
	/*
	 * @param	lengthA			: 	length of String A
	 * @param	lengthB			: 	length of String B
	 * @param	gapPenalty		:	indel penalty
	 * @return	scores			:	Initialized 2D int matrix to store alignment scores
	 */
	private int[][] initDistanceMatrix(int lengthA, int lengthB, int gapPenalty) {
		int[][] scores = new int[lengthA + 1][lengthB + 1];
		for(int iii = 0; iii <= lengthA; iii++){
			scores[iii][0] = iii * gapPenalty;
		}
		for(int jjj = 0; jjj <= lengthB; jjj++) {
			scores[0][jjj] = jjj * gapPenalty;
		}
		return scores;
	}
	
	/*
	 * @param	A			:		String A to compute alignment scores
	 * @param	B			:		String B to compute alignment scores
	 * @param	distances	:		HashMap that stores BLOSUM62 matrix
	 * @param	gapPenalty	:		indel penalty
	 * @return	scores		:		calculated alignment scores		
	 */
	private int[][] computeNWScore(String A, String B, HashMap<Character, HashMap<Character, Integer>> distances, int gapPenalty) {
		int[][] scores = initDistanceMatrix(A.length(), B.length(), gapPenalty);
		for(int iii = 1; iii <= A.length(); iii++) {
			HashMap<Character, Integer> tempMap = distances.get(A.charAt(iii - 1));
			for(int jjj = 1; jjj <= B.length(); jjj++) {
				int value = tempMap.get(B.charAt(jjj - 1));
				int match = scores[iii - 1][jjj - 1] + value;
				int delete = scores[iii - 1][jjj] + gapPenalty;
				int insert = scores[iii][jjj - 1] + gapPenalty;
				scores[iii][jjj] = max(match, delete, insert);
			}
		}
		return scores;
	}
	
	/*
	 * HIRSCHBERG
	 */
	private int[] computeLastNWScore(String A, String B, HashMap<Character, HashMap<Character, Integer>> distances, int gapPenalty) {
		int[][] scores = computeNWScore(A, B, distances, gapPenalty);
		int[] lastLine = new int[B.length() + 1];
		for(int jjj = 0; jjj <= B.length(); jjj++) {
			lastLine[jjj] = scores[A.length()][jjj];
		}
		return lastLine;
	}
	
	/*
	 * HIRSCHBERG
	 */
	private String[] hirschberg(String A, String B, HashMap<Character, HashMap<Character, Integer>> distances) {
		final int gapPenalty = -2;
		String Z = "";
		String W = "";
		if(A.length() == 0 || B.length() == 0) {
			if(A.length() == 0) {
				for(int jjj = 0; jjj < B.length(); jjj++) {
					Z = Z + "-";
					W = W + Character.toString(B.charAt(jjj));
				}
			} else {
				for(int iii = 0; iii < A.length(); iii++) {
					Z = Z + Character.toString(A.charAt(iii));
					W = W + "-";
				}
			}
		} else if (A.length() == 1 || B.length() == 1) {
			String[] alignments = getAlignment(A, B, distances);
			Z = Z + alignments[0];
			W = W + alignments[1];
		} else {
			int midA = (int) (A.length() / 2);
			int[] scoreL = computeLastNWScore(A.substring(0, midA), B, distances, gapPenalty);
			int[] scoreR = computeLastNWScore(reverse(A.substring(midA)), reverse(B), distances, gapPenalty);
			int midB = partition(scoreL, scoreR);
			String[] leftAlignments = hirschberg(A.substring(0, midA), B.substring(0, midB), distances);
			Z = Z + leftAlignments[0];
			W = W + leftAlignments[1];	
			String[] rightAlignments = hirschberg(A.substring(midA), B.substring(midB), distances);
			Z = Z + rightAlignments[0];
			W = W + rightAlignments[1];
		}
		String[] result = {Z, W};
		return result;
	}
	
	/*
	 * HIRSCHBERG
	 */
	private static int partition(int[] scoreL, int[] scoreR) {
		int maxSum = Integer.MIN_VALUE;
		int index = 0;
		for(int iii = 0; iii < scoreL.length; iii++) {
			int sum = scoreL[iii] + scoreR[scoreL.length - iii - 1];
			if(sum >= maxSum) {
				maxSum = sum;
				index = iii;
			}
		}
		return index;
	}
	
	/*
	 * HIRSCHBERG
	 */
	private static String reverse(String inputString) {
		if(inputString == null) return "";
		if(inputString.length() == 1) return inputString;
		return new StringBuilder(inputString).reverse().toString();
	}

	private static HashMap<Character, HashMap<Character, Integer>> initDistancesMap() {
		HashMap<Character, HashMap<Character, Integer>> distancesMap = new HashMap<Character, HashMap<Character, Integer>>(); 
		final String distancesFileName = "/home/farhang/workspace/BioUIMA/resources/distances.txt";
		String distancesFile = BioUima.readInputFile(distancesFileName);
		String[] distancesFileArray = distancesFile.split("\n");
		String[] firstLine = distancesFileArray[0].split("\t");
		for(int iii = 1; iii < distancesFileArray.length; iii++) {
			HashMap<Character, Integer> tempMap = new HashMap<Character, Integer>();
			String[] singleLineArray = distancesFileArray[iii].split("\t");
			for(int jjj = 0; jjj < singleLineArray.length; jjj++) {
				tempMap.put(firstLine[jjj].charAt(0), Integer.parseInt(singleLineArray[jjj]));
			}
			distancesMap.put(firstLine[iii - 1].charAt(0), tempMap);
		}
		return distancesMap;
	}

	/*
	 * Returns greatest of three integers
	 */
	private static int max(int x, int y, int z) {
		return Math.max(x, Math.max(y, z));
	}
}