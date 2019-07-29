package partitioner;

import org.broadinstitute.sv.commandline.CommandLineProgram;
import org.broadinstitute.sv.util.CommandRunner;
import org.broadinstitute.sv.util.GenomeInterval;
import org.broadinstitute.sv.util.fasta.FastaReader;
import org.broadinstitute.sv.util.fasta.IndexedFastaFile;

import jaligner.Alignment;
import jaligner.NeedlemanWunschGotoh;
import jaligner.Sequence;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;

import java.io.*;
import java.util.*;

public class VNTRPartitioner extends CommandLineProgram {

	private static final String DEFAULT_MUSCLE_EXECUTABLE = "/humgen/cnp04/sandbox/bobh/muscle/muscle";
	private static final String USER = "yiming@handsake-vm.broadinstitute.org";
	private static final String MATRIX = "Z:\\sandbox\\bobh\\workspace\\gap_git\\gap\\FeatureDB\\src\\server\\edu\\mit\\broad\\seqserver\\utils\\align\\matrices\\GAP_DNA_IUPAC.res";


	private File mAlignmentOutputFile = null;
	private File mMuscleOutputFile = null;
	private File mSortedAlignmentOutputFile = null;
	private File mAlignmentScoresOutputFile = null;
	private String mMuscleExecutable = null;




	@Input(fullName="scanFile", shortName="S", doc="VNTR scan file", required=true)
	private File mScanFile = null;

	@Input(fullName="referenceFile", shortName="R", doc="Genome reference file", required=true)
	private File mRefFile = null;

	@Argument(fullName="identifier", shortName="I", doc="VNTR identifier", required=true)
	private String mIdentifier = null;

	@Argument(fullName="contextBefore", shortName="B", doc="Number of periods of context before", required=false)
	private Integer mContextBefore = null;

	@Argument(fullName="contextAfter", shortName="A", doc="Number of periods of context after", required=false)
	private Integer mContextAfter = null;

	public static void main(String[] args) throws Exception {
		run(new VNTRPartitioner(), args);
	}

	protected int run() throws IOException {



		IndexedFastaFile referenceFile = IndexedFastaFile.open(mRefFile);

		GenomeInterval vntrInterval = parseNewInterval();
		int modePeriod = computeModeFrequency(vntrInterval);
		
		String vntr = referenceFile.getSequence(vntrInterval);
		long mP = modePeriod;
		long len = vntr.length();
		long prod = mP*mP*len;
		if(prod > 10000000000L) return 0;

		// Output files
		mAlignmentOutputFile = new File("Z:\\sandbox\\yiming\\AlignmentFiles\\chr" + vntrInterval.getSequenceName() + "\\" + mIdentifier + ".fasta");
		mMuscleOutputFile = new File("Z:\\sandbox\\yiming\\MuscleOutputs\\chr" + vntrInterval.getSequenceName() + "\\" + mIdentifier + ".clw");
		mSortedAlignmentOutputFile = new File("Z:\\sandbox\\yiming\\SortedAlignmentFiles\\chr" + vntrInterval.getSequenceName() + "\\" + mIdentifier + ".txt");
		mAlignmentScoresOutputFile = new File("Z:\\sandbox\\yiming\\AlignmentScores\\chr" + vntrInterval.getSequenceName() + "\\" + mIdentifier + ".txt");

		String initial = vntr.substring(0, modePeriod);
		List<Integer> localMaxima = new ArrayList<Integer>();
		float[] scores = new float[vntr.length()-modePeriod];
		for(int i=0; i<vntr.length()-modePeriod; i++) {
			Sequence seq1 = new Sequence(initial);
			Sequence seq2 = new Sequence(vntr.substring(i, i+modePeriod));
			Alignment align = null;
			try {
				align = NeedlemanWunschGotoh.align(seq1, seq2, MatrixLoader.load(MATRIX), 2f, 0.5f);
			} catch (MatrixLoaderException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			//			scores[i] = FuzzySearch.ratio(initial, vntr.substring(i, i+modePeriod));
			scores[i] = align.getScore();
			System.out.println("index: " + i + ", score: " + scores[i]);
		}

		printScores(scores, mAlignmentScoresOutputFile);

		localMaxima.add(0);

		// TODO: might potentially be a problem when there is no remainder, not sure yet
		for(int i=1; i<vntr.length()-modePeriod-1; i++) {
			boolean isLocalMax = true;
			for(int j=0; j<10; j++) {
				if(i-j >= 0 && scores[i] < scores[i-j]) {
					isLocalMax = false;
					break;
				}
				if(i+j < scores.length && scores[i] < scores[i+j]) {
					isLocalMax = false;
					break;
				}
			}
			if(i+10 >= scores.length) {
				isLocalMax = false;
			}
			if(isLocalMax) {
				localMaxima.add(i);
			}
		}

		/**
		 * If local maxima of equal score occur at consecutive indices, take the floor of the average of the start and end
		 * Note: doesn't always work, as in the case of chr 20
		 * Might want to confirm those with very high scores (e.g. only 1 or 2 bp substitutions), adding both their starts and ends 
		 * to localMaxima, this would avoid this problem but there would still be edge cases that need care
		 */
		int counter= 0;
		List<Integer> maximaToRemove = new ArrayList<Integer>();
		while(counter < localMaxima.size()) {
			int start = counter;
			float score = scores[localMaxima.get(counter)];
			while(counter < localMaxima.size() && scores[localMaxima.get(counter)] == score && Math.abs(localMaxima.get(counter) - localMaxima.get(start)) < 10) {
				counter++;
			}
			if(counter > start+1) {
				int end = counter-1;
				int ind = (start + end)/2;
				for(int i=start; i<=end; i++) {
					if(i != ind) {
						maximaToRemove.add(localMaxima.get(i));
					}
				}
			}
		}

		for(int i=0; i<maximaToRemove.size(); i++) {
			localMaxima.remove(maximaToRemove.get(i));
		}

		BufferedWriter bw = new BufferedWriter(new PrintWriter(mAlignmentOutputFile));
		for(int i=0; i<localMaxima.size()-1; i++) {
			int start = localMaxima.get(i);
			int end = localMaxima.get(i+1);
			bw.write(">[" + start + "," + (end-1) + "]\n" + vntr.substring(start, end) +"\n");
		}
		int start = localMaxima.get(localMaxima.size()-1);
		int end = start + modePeriod - 1;
		bw.write(">[" + start + "," + end + "]\n" + vntr.substring(start, end+1) + "\n");
		bw.write(">[" + (end + 1) + "," + vntr.length() + "]\n" + vntr.substring(end+1));

		bw.close();

		runMuscle(new File(convertToUnix(mAlignmentOutputFile.toString())), new File(convertToUnix(mMuscleOutputFile.toString())));
		List<String> sortedIntervals = sortIntervals(mMuscleOutputFile);

		BufferedWriter bw1 = new BufferedWriter(new PrintWriter(mSortedAlignmentOutputFile));
		String consensus = computeConsensusString(sortedIntervals, sortedIntervals.get(0).length());
		bw1.write(consensus + "\n");

		for(int i=0; i<sortedIntervals.size(); i++) {
			bw1.write(sortedIntervals.get(i) + "\n");
		}

		bw1.close();


		return 0;
	}

	/**
	 * Prints alignment scores to output file
	 * @param scores alignment scores
	 * @param outputFile output file
	 */
	private static void printScores(float[] scores, File outputFile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		for(int i=0; i<scores.length; i++) {
			bw.write(i + " " + scores[i] + "\n");
		}
		bw.close();
	}

	/** 
	 * Computes the consensus string given a list of strings of equal length
	 * @param strings list of strings
	 * @param period length of each string
	 * @return consensus string
	 */
	private static String computeConsensusString(List<String> strings, int period) {

		StringBuilder builder = new StringBuilder();

		for(int i=0; i<period; i++) {
			Map<Character, Integer> letterFrequencies = new HashMap<Character, Integer>();
			for(String s : strings) {
				if(s.length() > i) {
					if(letterFrequencies.containsKey(Character.valueOf(s.charAt(i)))) {
						letterFrequencies.put(Character.valueOf(s.charAt(i)), letterFrequencies.get(Character.valueOf(s.charAt(i))) + 1);
					}
					else {
						letterFrequencies.put(Character.valueOf(s.charAt(i)), 1);
					}
				}
			}
			builder.append(Collections.max(letterFrequencies.entrySet(), Map.Entry.comparingByValue()).getKey());
		}

		return builder.toString();
	}

	/**
	 * Finds starting index of nucleotide sequence from a line of a muscle output file
	 * @param s line in muscle output file
	 * @return starting index
	 */
	private static int findStartingIndex(String s) {
		for(int i=0; i<s.length(); i++) {
			char c = s.charAt(i);
			if(c == '|') {
				return i+1;
			}
		}
		return -1;
	}

	/**
	 * Reads in muscle output file line by line, sorts intervals, and truncates interval labels
	 * @param inputFile
	 * @return
	 * @throws IOException
	 */
	private static List<String> sortIntervals(File inputFile) throws IOException {
		FastaReader fr = new FastaReader(inputFile);
		//		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		//
		//		// discard header and empty lines
		//		br.readLine();
		//		br.readLine();
		//		br.readLine();
		//
		//		List<String> lines = new ArrayList<String>();
		//		String line;
		//		while((line = br.readLine()).charAt(0) == '[') {
		//			lines.add(line);
		//		}
		//		Collections.sort(lines);
		//
		//		br.close();

		List<Line> lines = new ArrayList<Line>();
		String line;
		while(fr.hasNext()) {
			String seq = fr.next().toString();
			String head = fr.getHeaderLine();
			line = head + "|" + seq;
			lines.add(new Line(line));
		}

		Collections.sort(lines);

		List<String> vntrLines = new ArrayList<String>();
		for(Line l : lines) {
			int startIndex = findStartingIndex(l.getLine());
			vntrLines.add(l.getLine().substring(startIndex));
		}
		return vntrLines;
	}

	/**
	 * Invokes muscle in the command line
	 * @param inputFile: the file inputed for alignment
	 * @param outputFile: the file outputed by muscle
	 */
	private void runMuscle(File inputFile, File outputFile) {
		List<String> command = new ArrayList<String>();
		command.add("ssh");
		command.add("-q");
		command.add(USER);
		command.add((mMuscleExecutable != null) ? mMuscleExecutable : DEFAULT_MUSCLE_EXECUTABLE);
		command.add("-in");
		command.add(inputFile.toString().replace('\\', '/'));
		command.add("-out");
		command.add(outputFile.toString().replace('\\', '/'));

		String[] commandArray = command.toArray(new String[command.size()]);
		CommandRunner runner = new CommandRunner();
		runner.setMergeOutput(true);
		int status = runner.runCommand(commandArray, null, null);

		if (status != 0) {
			throw new RuntimeException("Error running muscle command " + 
					formatCommand(commandArray) + " : " +
					runner.getStandardOutputString());
		}

	}

	/**
	 * Converts pathname to Unix (Z: --> \humgen\cnp04)
	 * @param pathname Windows pathname
	 * @return Unix pathname
	 */
	private static String convertToUnix(String pathname) {
		pathname.replace('\\', '/');
		return "/humgen/cnp04" + pathname.substring(2);		
	}


	/**
	 * Formats error message
	 * @param command
	 * @return
	 */
	private String formatCommand(String[] command) {
		StringBuilder builder = new StringBuilder();
		for (int i = 1; i < command.length; i++) {
			if (i > 1) {
				builder.append(' ');
			}
			builder.append(command[i]);
		}
		return builder.toString();
	}

	/**
	 * Gets the GenomeInterval from the identifier and updated start and end coordinates from the new dotplot file.
	 * @return The desired GenomeInterval
	 * @throws IOException
	 */
	private GenomeInterval parseNewInterval() throws IOException {
		int counter = 0;
		int start = 0;
		int end = mIdentifier.length()-1;
		for(int i=0; i<mIdentifier.length(); i++) {
			char currChar = mIdentifier.charAt(i);
			if(currChar == '_') {
				counter++;
			}
			if(counter == 2 && start == 0) {
				start = i+1;
			}
			else if(counter == 3) {
				end = i;
				break;
			}
		}

		String id = mIdentifier.substring(start, end);
		File sitesFile = new File("Z:\\sandbox\\bobh\\projects\\vntrs\\bobh_scan5\\sites\\chr" + id + "\\" + mIdentifier + "\\" + mIdentifier + ".region_updated.dat");
		BufferedReader br = new BufferedReader(new FileReader(sitesFile));

		br.readLine(); // discard header
		String line = br.readLine();
		String[] lineSplit = line.split("\\s+");
		br.close();
		return new GenomeInterval(id, Integer.parseInt(lineSplit[3]), Integer.parseInt(lineSplit[4]));
	}



	/**
	 * Parses start coordinates from scan file line by lines and computes periods by taking differences of adjacent
	 * start coordinates on each line. Then computes most frequent such period.
	 * @return The mode period
	 */
	private int computeModeFrequency(GenomeInterval interval) throws IOException {
		Map<Integer, Integer> periodFrequencies = new HashMap<Integer, Integer>();
		BufferedReader br = new BufferedReader(new FileReader(mScanFile));

		br.readLine(); // discard header	
		String line;
		while((line = br.readLine()) != null) {
			String[] lineSplit = line.split("\\s+");
			GenomeInterval currInterval = new GenomeInterval(lineSplit[0], Integer.parseInt(lineSplit[1]), Integer.parseInt(lineSplit[2]));

			if (currInterval.overlaps(interval)) {
				String[] startCoords = lineSplit[5].split(",");
				for(int i=1; i<startCoords.length; i++) {
					int freq = Integer.parseInt(startCoords[i]) - Integer.parseInt(startCoords[i-1]);
					if(periodFrequencies.containsKey(freq)) {
						periodFrequencies.put(freq, periodFrequencies.get(freq)+1);
					}
					else {
						periodFrequencies.put(freq, 1);
					}
				}

			}
		}

		br.close();
		return Collections.max(periodFrequencies.entrySet(), Map.Entry.comparingByValue()).getKey();
	}

	private static class Line implements Comparable<Line>{
		private String line;
		
		public Line(String line) {
			this.line = line;
		}
		
		public String getLine() {
			return line;
		}
		
		@Override
	    public int compareTo(Line other){
	       int start1 = this.line.indexOf('[') + 1;
	       int end1 = this.line.indexOf(',');
	       int start2 = other.getLine().indexOf('[') + 1;
	       int end2 = other.getLine().indexOf(',');
	       
	       return Integer.parseInt(line.substring(start1, end1)) - Integer.parseInt(other.getLine().substring(start2, end2));
	    }
		
	}


}
