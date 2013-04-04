package bio.uima;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;

public class SequenceAlignmentAnnotator extends JCasAnnotator_ImplBase  {
	
	private static int GAP_PENALTY = -2;													// Penalty for insertion or deletion
	private static HashMap<Character, HashMap<Character, Integer>> distances = null;		// BLOSUM62 matrix
	
	@Override
	public void process(JCas cas) throws AnalysisEngineProcessException {
		try {
			String[] proteins = cas.getView("protein").getDocumentText().split(" ");
			String[] alignment = getAllPairsAlignment(proteins);
			JCas alignmentCas = cas.createView("alignment");								// create view for alignments
			alignmentCas.setDocumentText(Utils.combineStringArray(alignment));
		} catch (CASException e) {
			e.printStackTrace();
		}
	}

	/*
	 * Aligns and returns all pairs of protein sequences which are input.
	 * 
	 * @param	proteins		Array of strings
	 * @return 					Array of strings that contains alignment of all pairs of input String array
	 */
	private String[] getAllPairsAlignment(String[] proteins) {
		if(distances == null) {	distances = Utils.initSimilarityMap(); }
		String alignmentString = "";
		for(int iii = 0; iii < (proteins.length - 1); iii++) {
			for(int jjj = (iii + 1); jjj < proteins.length; jjj++) {
				String[] singlePairAlignment = hirschberg(proteins[iii], proteins[jjj]);
				alignmentString += (singlePairAlignment[0] + " " + singlePairAlignment[1] + " ");
			}
		}
		return alignmentString.split(" ");
	}

	/*
	 * Aligns two protein sequences and returns the optimal alignment as an array of Strings
	 * 
	 * @param 	A				String to align
	 * @param	B				String to align
	 * @param	distances		HashMap that stores BLOSUM62 matrix
	 * @return					Array of two strings that are the optimal alignment of A and B
	 */
	private String[] getAlignment(String A, String B) {
		if(distances == null) {	distances = Utils.initSimilarityMap(); }
		int[][] scores = computeNWScore(A, B);
		
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
			} else if ((iii > 0) && (scores[iii][jjj] == scores[iii-1][jjj] + GAP_PENALTY)) {
				alignmentA = Character.toString(A.charAt(iii - 1)) + alignmentA;
				alignmentB = Character.toString('-') + alignmentB;
				iii--;
			} else if((jjj > 0) && (scores[iii][jjj] == scores[iii][jjj - 1] + GAP_PENALTY)) {
				alignmentA = Character.toString('-') + alignmentA;
				alignmentB = Character.toString(B.charAt(jjj - 1)) + alignmentB;
				jjj--;
			}
		}
		String[] alignment = {alignmentA, alignmentB};
		return alignment;
	}

	/*
	 * Computes and returns the Needleman-Wunsch score matrix which contains
	 * the optimal alignment scores for the two input sequences
	 * 
	 * @param	sequence1			Sequence to compute alignment scores
	 * @param	sequence2			Sequence to compute alignment scores
	 * @return						calculated alignment score matrix		
	 */
	private int[][] computeNWScore(String sequence1, String sequence2) {
		if(distances == null) { distances = Utils.initSimilarityMap(); }
		int[][] scores = Utils.initDistanceMatrix(sequence1.length(), sequence2.length(), GAP_PENALTY);
		for(int iii = 1; iii <= sequence1.length(); iii++) {
			HashMap<Character, Integer> tempMap = distances.get(sequence1.charAt(iii - 1));
			for(int jjj = 1; jjj <= sequence2.length(); jjj++) {
				int value = tempMap.get(sequence2.charAt(jjj - 1));
				int match = scores[iii - 1][jjj - 1] + value;
				int delete = scores[iii - 1][jjj] + GAP_PENALTY;
				int insert = scores[iii][jjj - 1] + GAP_PENALTY;
				scores[iii][jjj] = Utils.max(match, delete, insert);
			}
		}
		return scores;
	}
	
	/*
	 * Computes and returns the last row of the Needleman-Wunsch score matrix
	 * 
	 * @param	sequence1			Sequence to compute alignment scores
	 * @param	sequence2			Sequence to compute alignment scores
	 * @return						array representing the last row of score matrix
	 */
	private int[] computeLastNWScore(String sequence1, String sequence2) {
		int[][] scores = computeNWScore(sequence1, sequence2);
		int[] lastLine = new int[sequence2.length() + 1];
		for(int jjj = 0; jjj <= sequence2.length(); jjj++) {
			lastLine[jjj] = scores[sequence1.length()][jjj];
		}
		return lastLine;
	}
	
	/*
	 * Computes and returns the optimal alignment of two sequences using
	 * the Hirschberg's algorithm.
	 * 
	 * @param	sequence1			Sequence to compute alignment scores
	 * @param	sequence2			Sequence to compute alignment scores
	 * @return						Array of strings that holds the optimal alignment of input strings
	 */
	private String[] hirschberg(String sequence1, String sequence2) {
		if(distances == null) { distances = Utils.initSimilarityMap(); }
		String Z = "";
		String W = "";
		if(sequence1.length() == 0 || sequence2.length() == 0) {
			if(sequence1.length() == 0) {
				for(int jjj = 0; jjj < sequence2.length(); jjj++) {
					Z = Z + "-";
					W = W + Character.toString(sequence2.charAt(jjj));
				}
			} else {
				for(int iii = 0; iii < sequence1.length(); iii++) {
					Z = Z + Character.toString(sequence1.charAt(iii));
					W = W + "-";
				}
			}
		} else if (sequence1.length() == 1 || sequence2.length() == 1) {
			String[] alignments = getAlignment(sequence1, sequence2);
			Z = Z + alignments[0];
			W = W + alignments[1];
		} else {
			int midA = (int) (sequence1.length() / 2);
			int[] scoreL = computeLastNWScore(sequence1.substring(0, midA), sequence2);
			int[] scoreR = computeLastNWScore(reverse(sequence1.substring(midA)), reverse(sequence2));
			int midB = partition(scoreL, scoreR);
			String[] leftAlignments = hirschberg(sequence1.substring(0, midA), sequence2.substring(0, midB));
			Z = Z + leftAlignments[0];
			W = W + leftAlignments[1];	
			String[] rightAlignments = hirschberg(sequence1.substring(midA), sequence2.substring(midB));
			Z = Z + rightAlignments[0];
			W = W + rightAlignments[1];
		}
		String[] result = {Z, W};
		return result;
	}
	
	/*
	 * Returns the index of the two input score arrays that yield the maximum
	 * sum of the elements at that index.
	 * In other words, returns argmax iii (scoreL[iii] + reverse(scoreR[iii])
	 * 
	 * @param	scoreL			Score matrix for left side of partition
	 * @param	scoreR			Score matrix for right side of partition
	 * @return					The index that maximizes the sum of the two array elements
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
	 * Returns the reversed string
	 * 
	 * @param	inputString			String to reverse
	 * @return						inputString reversed
	 */
	private static String reverse(String inputString) {
		if(inputString == null) return "";
		if(inputString.length() == 1) return inputString;
		return new StringBuilder(inputString).reverse().toString();
	}
	
	/*
	 * Performs DBSCAN clustering algorithm on an array of protein
	 * sequences.
	 * 
	 * @param	proteins		Input protein sequences
	 * @param	epsilon			Neighborhood to check points
	 * @param	minPts			Minimum number of points that must be in epsilon radius
	 * @return					Array that holds which cluster each protein belongs to
	 */
	private int[] dbScan(String[] proteins, int epsilon, int minPts) {
		int cluster = 0;
		boolean[] visited = new boolean[proteins.length];
		int[] proteinClusters = new int[proteins.length];
		for(int iii = 0; iii < proteins.length; iii++) {
			if(visited[iii]) continue;
			visited[iii] = true;
			HashSet<Integer> neighbors = epsilonNeighbors(proteins, iii, epsilon);	
			if(neighbors.size() < minPts) {
				proteinClusters[iii] = -1;		// mark protein iii as noise
			} else {
				++cluster;
				expandCluster(proteins, iii, neighbors, proteinClusters, visited, cluster, epsilon, minPts);
			}
		}
		return proteinClusters;
	}
    
	/*
	 * Expands cluster if a core point is found which has >= minPts in its
	 * epsilon radius
	 * 
	 * @param	proteins			Dataset containing protein sequences
	 * @param	index				Index of current protein
	 * @param	neighbors			Neighbors of current protein
	 * @param	proteinClusters		Array representing protein clusters
	 * @param	visited				Array representing if a protein has been visited
	 * @param	cluster				Current cluster
	 * @parame	epsilon				Epsilon radius
	 * @param	minPts				Minimum number of points to check
	 */
    private void expandCluster(String[] proteins, int index, HashSet<Integer> neighbors,
    		int[] proteinClusters, boolean[] visited, int cluster, int epsilon, int minPts) {
    	proteinClusters[index] = cluster;
    	Queue<Integer> seeds = new LinkedList<Integer>();
    	seeds.addAll(neighbors);
    	while(!seeds.isEmpty()) {
    		int currentProtein = (Integer) seeds.poll();
    		if(!visited[currentProtein]) {
    			visited[currentProtein] = true;
    			HashSet<Integer> currentNeighbors = epsilonNeighbors(proteins, currentProtein, epsilon);
    			if(currentNeighbors.size() >= minPts) {
    				seeds.addAll(currentNeighbors);
    			}
    		}
    		if(proteinClusters[currentProtein] == 0) proteinClusters[currentProtein] = cluster;
    	}
    }
	
    /*
     * Calculates and returns the similarity of two
     * protein sequences based on the BLOSUM62 matrix.
     * 
     * @param	protein1			First protein
     * @param	protein2			Second protein
     * @return						Similarity score of two proteins
     */
	private int similarity(String protein1, String protein2) {
		int similarity = 0;
		int length = (protein1.length() > protein2.length()) ? protein1.length(): protein2.length();
		for(int iii = 0; iii < length; iii++) {
			char p1, p2;
			try {
				p1 = protein1.charAt(iii);
			} catch(StringIndexOutOfBoundsException e) {
				p1 = '-';
			}
			try {
				p2 = protein2.charAt(iii);
			} catch(StringIndexOutOfBoundsException e) {
				p2 = '-';
			}
			if(p1 == '-' || p2 == '-') {
				similarity += GAP_PENALTY;
			} else {
				similarity += distances.get(p1).get(p2);
			}
		}
		
		return similarity;
	}
	
	/*
	 * Returns indexes of all the proteins which are in the epsilon
	 * neighborhood of the protein at index == currentProtein
	 * 
	 * @param	proteins			Dataset containing all proteins
	 * @param	currentProtein		Index of current protein in the dataset
	 * @param	epsilon				Minimum similarity to check neighborhood
	 * @return						HashSet containing the neighbors
	 */
	private HashSet<Integer> epsilonNeighbors(String[] proteins, int currentProtein, int epsilon) {
		HashSet<Integer> neighbors = new HashSet<Integer>();
		String protein = proteins[currentProtein];
		for(int iii = 0; iii < proteins.length; iii++) {
			int score = similarity(proteins[iii], protein);
			if(score >= epsilon) neighbors.add(iii);
		}
		return neighbors;
	}
}