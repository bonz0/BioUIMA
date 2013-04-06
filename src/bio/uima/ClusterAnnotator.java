package bio.uima;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;

import org.apache.uima.analysis_component.JCasAnnotator_ImplBase;
import org.apache.uima.analysis_engine.AnalysisEngineProcessException;
import org.apache.uima.cas.CASException;
import org.apache.uima.jcas.JCas;

public class ClusterAnnotator extends JCasAnnotator_ImplBase {
	private static int GAP_PENALTY = -2;													// Penalty for insertion or deletion
	private static HashMap<Character, HashMap<Character, Integer>> distances = Utils.initSimilarityMap();		// BLOSUM62 matrix
	private static int currentEpsilon = 30;
	private static int minPoints = 3;
	
	@Override
	public void process(JCas cas) throws AnalysisEngineProcessException {
		// TODO Auto-generated method stub
		try {
			String[] alignments = cas.getView("alignments").getDocumentText().split(" ");
			int[] proteinClusters = dbScan(alignments, currentEpsilon, minPoints);
			String clusters = Utils.combineIntArray(proteinClusters);
			JCas cas1 = cas.createView("clusters");
			cas1.setDocumentText(clusters);
		} catch (CASException e) {
			e.printStackTrace();
		}
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