package partitioner;

import org.broadinstitute.sv.commandline.CommandLineProgram;
import org.broadinstitute.sv.util.CommandRunner;
import org.broadinstitute.sv.util.GenomeInterval;
import org.broadinstitute.sv.util.fasta.FastaReader;
import org.broadinstitute.sv.util.fasta.IndexedFastaFile;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Refined VNTR Partitioner which uses Smith-Waterman to align probe sequences in a sliding window fashion to an initial 
 * (supposedly consensus) sequence. Alignment scores are then calculated from these Smith-Waterman alignments and estimated start indices
 * for new periods are calculated from the local maxima indices.
 * @author yiming
 */
public class SmithWatermanVNTRPartitioner extends CommandLineProgram {

    private static final String DEFAULT_MUSCLE_EXECUTABLE = "/humgen/cnp04/sandbox/bobh/muscle/muscle";
    private static final int PADDING_DISTANCE = 20;
    private static final float MATCH = 2;
    private static final float MISMATCH = -1;
    private static final float OPEN = 2;
    private static final float EXTEND = 0.5f;
    private static final int CONTEXT_PERIODS_BEFORE = 20;
    private static final int CONTEXT_PERIODS_AFTER = 20;
    private static final int MAX_ALLOWED_PERIOD = 1000;
    private static final long MAX_ALLOWED_COST = 100_000_000_000L;
    private static char fileSeparator;
    private static Matrix matrix = null;




    private File mAlignmentOutputFile = null;
    private File mMuscleOutputFile = null;
    private File mSortedAlignmentOutputFile = null;
    private File mAlignmentScoresOutputFile = null;
    private File mAlignmentScoresOutputFile2 = null;
    private File mPeriodsOutputFile = null;
    private File mContextFile = null;
    private String mMuscleExecutable = null;




    @Input(fullName="scanFile", shortName="S", doc="VNTR scan file", required=true)
    private File mScanFile = null;

    @Input(fullName="referenceFile", shortName="R", doc="Genome reference file", required=true)
    private File mRefFile = null;

    @Argument(fullName="identifier", shortName="I", doc="VNTR identifier", required=true)
    private String mIdentifier = null;
    
    @Argument(fullName="outputDirectory", shortName="O", doc="Root output directory", required=true)
    private String mOutputDirectory = null;

    @Argument(fullName="updatedRegion", shortName="U", doc="Updated scan directory", required=true)
    private String mUpdatedRegion = null;
    
    @Argument(fullName="username", shortName="N", doc="SSH username", required=false)
    private String mUsername = null;


    public static void main(String[] args) throws Exception {
        run(new SmithWatermanVNTRPartitioner(), args);
    }

    protected int run() throws IOException {
        disableJAlignerLogging();
        fileSeparator = File.separatorChar;
        IndexedFastaFile referenceFile = IndexedFastaFile.open(mRefFile);

        GenomeInterval vntrInterval = parseNewInterval();
        int estimatedModePeriod = computeModeFrequency(vntrInterval);
        long estimatedCost = ((long) vntrInterval.getLength()) * estimatedModePeriod * estimatedModePeriod;
        if (estimatedModePeriod > MAX_ALLOWED_PERIOD || estimatedCost > MAX_ALLOWED_COST) {
            System.out.println("Cannot run large VNTR " + mIdentifier +
                               " with length " + vntrInterval.getLength() +
                               " and estimated mode period " + estimatedModePeriod);
            return 0;
        }

        GenomeInterval beforeInterval = new GenomeInterval(vntrInterval.getSequenceName(), vntrInterval.getStart() - (CONTEXT_PERIODS_BEFORE * estimatedModePeriod), vntrInterval.getStart()-1);
        GenomeInterval afterInterval = new GenomeInterval(vntrInterval.getSequenceName(), vntrInterval.getEnd()+1, vntrInterval.getEnd() + (CONTEXT_PERIODS_AFTER * estimatedModePeriod));

        String vntr = referenceFile.getSequence(vntrInterval);
        String contextBefore = referenceFile.getSequence(beforeInterval);
        String contextAfter = referenceFile.getSequence(afterInterval);
        

        // Output files
        mAlignmentOutputFile = new File(mOutputDirectory + fileSeparator + "AlignmentFile.fasta");
        mMuscleOutputFile = new File(mOutputDirectory + fileSeparator + "MuscleOutput.fasta");
        mSortedAlignmentOutputFile = new File(mOutputDirectory + fileSeparator + "SortedAlignmentFile.txt");
        mAlignmentScoresOutputFile = new File(mOutputDirectory + fileSeparator + "AlignmentScore.txt");
        mAlignmentScoresOutputFile2 = new File(mOutputDirectory + fileSeparator + "AlignmentScore1.txt");
        mPeriodsOutputFile = new File(mOutputDirectory + fileSeparator + "Periods.txt");
        mContextFile = new File(mOutputDirectory + fileSeparator + "Context.txt");
        matrix = MatrixGenerator.generate(MATCH, MISMATCH);
        
        // print context
        BufferedWriter bwContext = new BufferedWriter(new FileWriter(mContextFile));
        bwContext.write(contextBefore + "\n");
        bwContext.write(contextAfter + "\n");
        bwContext.close();
       
        float[][] scores = new float[vntr.length()-estimatedModePeriod][2];
        scores = slidingWindowAlignment(vntr, estimatedModePeriod);
        
        // print scores from first sliding window
        printScores(scores, mAlignmentScoresOutputFile);
        
        List<Integer> estimatedLocalMaxima = refineLocalMaxima(scores, estimatedModePeriod, vntr.length());
        Collections.sort(estimatedLocalMaxima);
        
        List<Integer> estimatedLocalMaximaNoRepeats = new ArrayList<Integer>();
        for(Integer max : estimatedLocalMaxima) {
            if(!estimatedLocalMaximaNoRepeats.contains(max)) {
                estimatedLocalMaximaNoRepeats.add(max);
            }
        }
        
        List<Integer> periods = new ArrayList<Integer>();
        for(int i=0; i<estimatedLocalMaximaNoRepeats.size()-1; i++) {
            periods.add(estimatedLocalMaximaNoRepeats.get(i+1) - estimatedLocalMaximaNoRepeats.get(i));
        }
        
        int modePeriod = getMostFrequentElement(periods);
        float[][] updatedScores = slidingWindowAlignment(vntr, modePeriod);
        List<Integer> refinedLocalMaxima = refineLocalMaxima(updatedScores, modePeriod, vntr.length());

        // print scores from second sliding window
        printScores(updatedScores, mAlignmentScoresOutputFile2);
        
        // print original and new periods
        BufferedWriter periodWriter = new BufferedWriter(new PrintWriter(mPeriodsOutputFile));
        periodWriter.write("VNTR_ID\tOLD_PERIOD\tNEW_PERIOD\n");
        periodWriter.write(mIdentifier + "\t" + estimatedModePeriod + "\t" + modePeriod);
        periodWriter.close();

        BufferedWriter bw = new BufferedWriter(new PrintWriter(mAlignmentOutputFile));
        for(int i=0; i<refinedLocalMaxima.size()-1; i++) {
            int start = refinedLocalMaxima.get(i);
            int end = refinedLocalMaxima.get(i+1);
            bw.write(">[" + start + "," + (end-1) + "]\n" + vntr.substring(start, end) +"\n");
        }
        int start = refinedLocalMaxima.get(refinedLocalMaxima.size()-1);
        int end = start + modePeriod - 1;
        if(end < vntr.length()) {
            bw.write(">[" + start + "," + end + "]\n" + vntr.substring(start, end+1) + "\n");
            bw.write(">[" + (end + 1) + "," + vntr.length() + "]\n" + vntr.substring(end+1));
        }
        else {
            bw.write(">[" + start + "," + vntr.length() + "]\n" + vntr.substring(start) + "\n");
        }

        bw.close();
        
        if(runningOnWindows()) {
            runMuscle(new File(convertToUnix(mAlignmentOutputFile.toString())), new File(convertToUnix(mMuscleOutputFile.toString())));
        }
        else {
            runMuscle(mAlignmentOutputFile, mMuscleOutputFile);
        }

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
    
    private static List<Integer> refineLocalMaxima(float[][] scores, int period, int length) {
     // TODO: might potentially be a problem when there is no remainder, not sure yet
        List<Integer> potentialLocalMaxima = new ArrayList<Integer>();
        for(int i=0; i<length-period-1; i++) {
            boolean isLocalMax = true;
            for(int j=0; j<10; j++) {
                if(i-j >= 0 && scores[i][0] < scores[i-j][0]) {
                    isLocalMax = false;
                    break;
                }
                if(i+j < scores[0].length && scores[i][0] < scores[i+j][0]) {
                    isLocalMax = false;
                    break;
                }
            }
            // get rid of local maxima at the end, this also might lead to a problem
            if(i+10 >= scores.length) {
                isLocalMax = false;
            }
            if(isLocalMax) {
                potentialLocalMaxima.add(i);
            }
        }
        

        /**
         * If local maxima of equal score occur at consecutive indices, take the floor of the average of the start and end
         * Note: doesn't always work, as in the case of chr 20
         * Might want to confirm those with very high scores (e.g. only 1 or 2 bp substitutions), adding both their starts and ends 
         * to localMaxima, this would avoid this problem but there would still be edge cases that need care
         */
        int counter= 0;
        List<Integer> estimatedLocalMaxima = new ArrayList<Integer>();
        while(counter < potentialLocalMaxima.size()) {
            int start = counter;
            float score = scores[potentialLocalMaxima.get(counter)][0];
            float estimatedStart = scores[potentialLocalMaxima.get(counter)][1];
            while(counter < potentialLocalMaxima.size() && scores[potentialLocalMaxima.get(counter)][0] == score && scores[potentialLocalMaxima.get(counter)][1] == estimatedStart) {
                counter++;
            }
            if(counter > start+1) {
                List<Integer> tiedLocalMax = new ArrayList<Integer>();
                for(int i=start; i<counter; i++) {
                    tiedLocalMax.add((int) scores[potentialLocalMaxima.get(i)][1]);
                }
                estimatedLocalMaxima.add(getMostFrequentElement(tiedLocalMax));
            }
        }
        
        Collections.sort(estimatedLocalMaxima);
        return estimatedLocalMaxima;
    }
    
    private static float[][] slidingWindowAlignment(String referenceString, int period) {
        String initial = referenceString.substring(0, period);
        float[][] scores = new float[referenceString.length()-period][2];
        Sequence seq1 = new Sequence(initial);
        for(int i=0; i<referenceString.length()-period; i++) {
            int start = (i - PADDING_DISTANCE >= 0) ? (i - PADDING_DISTANCE) : 0;
            int end = (i + period + PADDING_DISTANCE <= referenceString.length()) ? (i + period + PADDING_DISTANCE) : referenceString.length();
            Sequence seq2 = new Sequence(referenceString.substring(start, end));
            Alignment align = null;
            try {
                align = SmithWatermanGotoh.align(seq1, seq2, matrix, OPEN, EXTEND);
            } catch (Exception e) {
                e.printStackTrace();
            }
            scores[i][0] = align.getScore();
            int initialGapLength = align.getStart1();
            if(initialGapLength > 0) {
                scores[i][0] += OPEN + (initialGapLength * EXTEND);
            }
            int terminalGapLength = initial.length() - computeNonGapLength(align.getSequence1()) - initialGapLength;
            if(terminalGapLength > 0) {
                scores[i][0] += OPEN + (terminalGapLength * EXTEND);
            }

            scores[i][1] = start + align.getStart2();
        }
        
        return scores;
    }
    
    /**
     * Check if operating on windows
     * @return
     */
    private static boolean runningOnWindows() { 
        return (fileSeparator == '\\');
    }

    /**
     * Get most frequently occurring element in a list
     * @param list most frequent element
     */
    private static int getMostFrequentElement(List<Integer> list) {
        Map<Integer, Integer> frequenciesMap = new HashMap<Integer, Integer>();
        for(Integer i : list) {
            if(!frequenciesMap.containsKey(i)) {
                frequenciesMap.put(i, 1);
            }
            else {
                frequenciesMap.put(i, frequenciesMap.get(i) + 1);
            }
        }

        return Collections.max(frequenciesMap.entrySet(), Map.Entry.comparingByValue()).getKey();
    }

    /**
     * Prints alignment scores to output file
     * @param scores alignment scores
     * @param outputFile output file
     */
    private static void printScores(float[][] scores, File outputFile) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
        bw.write("OFFSET\tSCORE\tESTIMATED_START_INDEX\n");
        for(int i=0; i<scores.length; i++) {
            bw.write(i + "\t" + scores[i][0] + "\t" + scores[i][1] + "\n");
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

    private static int computeNonGapLength(char[] seq) {
        int length = 0;
        for(int i=0; i<seq.length; i++) {
            if(seq[i] != '-') {
                length++;
            }
        }
        return length;
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
        if(runningOnWindows()) {
            command.add("ssh");
            command.add("-q");
            command.add(mUsername);
        }
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
        File sitesFile = new File(mUpdatedRegion + fileSeparator + "sites" + fileSeparator + "chr" + id + fileSeparator + mIdentifier + fileSeparator + mIdentifier + ".region_updated.dat");
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
    
    /**
     * JAligner logs by default.
     * This method disables the logging.
     */
    private static void disableJAlignerLogging() {
        // Disable JAligner logging, which is on by default.
        Logger.getLogger(SmithWatermanGotoh.class.getName()).setLevel(Level.OFF);
    }


}
