package bio.uima;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Utils {
	
	private static final String DISTANCES_FILE_NAME = "/home/farhang/workspace/BioUIMA/resources/distances.txt";
	private static final String CODON_FILE = "/home/farhang/workspace/BioUIMA/resources/codonTable.txt";

	/*
	 * Initializes the score matrix, prepping it for score calculation
	 * 
	 * @param	lengthA			: 	length of String A
	 * @param	lengthB			: 	length of String B
	 * @return						Initialized 2D int matrix to store alignment scores
	 */
	public static int[][] initDistanceMatrix(int lengthA, int lengthB, int penalty) {
		int[][] scores = new int[lengthA + 1][lengthB + 1];
		for(int iii = 0; iii <= lengthA; iii++){
			scores[iii][0] = iii * penalty;
		}
		for(int jjj = 0; jjj <= lengthB; jjj++) {
			scores[0][jjj] = jjj * penalty;
		}
		return scores;
	}
	
	/*
	 * Reads codon table from a file and creates a HashMap for quick lookup
	 *
	 * @return HashMap<String, String>		HashMap to map codons to protein
	 */
	public static HashMap<String, String> initCodonTable() {
    	String inputFileString = Utils.readInputFile(CODON_FILE);
    	HashMap<String, String> codonMap = new HashMap<String, String>();
    	String[] lines = inputFileString.split("\n");
    	for(String line : lines) {
    		String[] temp = line.split("\t");
    		codonMap.put(temp[0], temp[1]);
    	}
    	return codonMap;
    }
	
	/*
	 * Returns greatest of three integers
	 */
	public static int max(int x, int y, int z) {
		return Math.max(x, Math.max(y, z));
	}

	/*
	 * Prints an array of Strings
	 */
	public static void printStringArray(String[] array) {
		for(int iii = 0; iii < array.length; iii++) {
			System.out.println(iii + "->\t" + array[iii]);
		}
		System.out.println();
	}
	
	/*
	 * Prints an array of Strings
	 */
	public static void printTwoStringArrays(String[] array1, String[] array2) {
		for(int iii = 0; iii < array1.length; iii++) {
			System.out.println(array1[iii] + "\t" + array2[iii]);
		}
		System.out.println();
	}
	
	/*
	 * Adds an array of strings to an ArrayList of strings
	 */
	public static void addArrayToList(ArrayList<String> proteins, String[] newProteins) {
		for(String temp : newProteins) {
			proteins.add(temp);
		}
	}
	
	/*
	 * Combines an input string array to output one single string which is separated
	 * by a space for each element in the array
	 * 
	 * @param	proteins			input string array to combine
	 * @return						String of concatenation of input string array
	 */
	public static String combineStringArray(String[] proteins) {
		String returnString = "";
		for(int iii = 0; iii < proteins.length; iii++) {
			returnString += (proteins[iii] + " ");
		}
		return returnString;
	}
	
	/*
	 * Combines an input integer array to output one single string which is separated
	 * by a space for each element in the array
	 * 
	 * @param	proteins			input integer array to combine
	 * @return						String of concatenation of input integer array
	 */
	public static String combineIntArray(int[] proteins) {
		String returnString = "";
		for(int iii = 0; iii < proteins.length; iii++) {
			returnString += (proteins[iii] + " ");
		}
		return returnString;
	}
	
	/*
	 * Reads from a file and returns the string that represents
	 * the data in the file.
	 *
	 * @param inputFile
	 * @return	String read from input file
	 */
	public static String readInputFile(String inputFile) {
		BufferedReader br = null;
		String fileText = "";
		try {
			br = new BufferedReader(new FileReader(inputFile));
			String currentLine;
			while((currentLine = br.readLine()) != null) {
				fileText += (currentLine + "\n");
			}
		} catch(IOException e) {
			System.out.println("File not found: " + inputFile);
			System.exit(-1);
		} finally {
			try {
				if(br != null)	br.close();
			} catch(IOException e) {
				e.printStackTrace();
			}
		}
		return fileText;
	}
	
	/*
	 * Initializes the similarity map which stores similarity score for alignment of
	 * all possible proteins. This is used to calculate the scores for an alignment
	 * of proteins.
	 * 
	 * @return					HashMap of similarity of all pairs of characters of proteins
	 */
	public static HashMap<Character, HashMap<Character, Integer>> initSimilarityMap() {
		HashMap<Character, HashMap<Character, Integer>> distancesMap = new HashMap<Character, HashMap<Character, Integer>>(); 
		String distancesFile = Utils.readInputFile(DISTANCES_FILE_NAME);
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
}
