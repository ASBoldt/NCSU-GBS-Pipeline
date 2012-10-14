// AlignmentUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.alignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import net.maizegenetics.pal.datatype.*;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.util.Utils;

import java.io.PrintWriter;

import java.util.StringTokenizer;


/**
 * Helper utilities for alignments.
 */
public class AlignmentUtils {

    static FormattedOutput format = FormattedOutput.getInstance();

    private AlignmentUtils() {
        // To prevent instantiation of this class.
    }

    public static byte getAlleleValue(byte[] data, int site, Pack2Alignment.NUM_BITS_PER_VALUE numBits) {
        return getAlleleValue(data, site, numBits, 0);
    }

    public static byte getAlleleValue(byte[] data, int site, Pack2Alignment.NUM_BITS_PER_VALUE numBits, int padding) {

        if (numBits == Pack2Alignment.NUM_BITS_PER_VALUE.one) {
            int index = site / 8;
            int shift = 7 - (site % 8);
            return (byte) ((data[index + padding] >>> shift) & 0x1);
        } else if (numBits == Pack2Alignment.NUM_BITS_PER_VALUE.two) {
            int index = site / 4;
            int shift = 6 - (site % 4) * 2;
            return (byte) ((data[index + padding] >>> shift) & 0x3);
        } else if (numBits == Pack2Alignment.NUM_BITS_PER_VALUE.four) {
            int index = site / 2;
            if (site % 2 == 0) {
                return (byte) ((data[index + padding] >>> 4) & 0xf);
            } else {
                return (byte) (data[index + padding] & 0xf);
            }
        } else if (numBits == Pack2Alignment.NUM_BITS_PER_VALUE.eight) {
            return data[site + padding];
        } else {
            throw new IllegalArgumentException("AlignmentUtils: getAlleleValue: Don't know how to handle bits per value: " + numBits);
        }

    }

    /**
     *  Report number of sequences, sites, and data type
     *  @note does not alter alignment state. If data type not defined
     *  in alignment will find a suitable instance for report but will
     *  not change alignment!
     */
    public static void report(Alignment a, PrintWriter out) {

        DataType dt = a.getDataType();
        out.println("Number of sequences: " + a.getSequenceCount());
        out.println("Number of sites: " + a.getSiteCount());
        if (dt != null) {
            out.println("Data type: " + dt.getDescription());
        }

    }

    /** print alignment (default format: INTERLEAVED) */
    public static void print(Alignment a, PrintWriter out) {
        printInterleaved(a, out);
    }

    /** print alignment (in plain format) */
    public static void printPlain(Alignment a, PrintWriter out) {
        printPlain(a, out, false);
    }

    /** print alignment (in plain format) */
    public static void printPlain(Alignment a, PrintWriter out, boolean relaxed) {
        // PHYLIP header line
        out.println("  " + a.getSequenceCount() + " " + a.getSiteCount());

        for (int s = 0; s < a.getSequenceCount(); s++) {
            format.displayLabel(out, a.getIdGroup().getIdentifier(s).getName(), (relaxed ? 20 : 10));
            out.print("     ");
            printNextSites(a, out, false, s, 0, a.getSiteCount());
            out.println();
        }
    }

    /** print alignment (in PHYLIP SEQUENTIAL format) */
    public static void printSequential(Alignment a, PrintWriter out) {
        // PHYLIP header line
        out.println("  " + a.getSequenceCount() + " " + a.getSiteCount() + "  S");

        // Print sequences
        for (int s = 0; s < a.getSequenceCount(); s++) {
            int n = 0;
            while (n < a.getSiteCount()) {
                if (n == 0) {
                    format.displayLabel(out,
                            a.getIdGroup().getIdentifier(s).getName(), 10);
                    out.print("     ");
                } else {
                    out.print("               ");
                }
                printNextSites(a, out, false, s, n, 50);
                out.println();
                n += 50;
            }
        }
    }

    /** print alignment (in PHYLIP 3.4 INTERLEAVED format) */
    public static void printInterleaved(Alignment a, PrintWriter out) {
        int n = 0;

        // PHYLIP header line
        out.println("  " + a.getSequenceCount() + " " + a.getSiteCount());

        // Print sequences
        while (n < a.getSiteCount()) {
            for (int s = 0; s < a.getSequenceCount(); s++) {
                if (n == 0) {
                    format.displayLabel(out,
                            a.getIdGroup().getIdentifier(s).getName(), 10);
                    out.print("     ");
                } else {
                    out.print("               ");
                }
                printNextSites(a, out, true, s, n, 50);
                out.println();
            }
            out.println();
            n += 50;
        }
    }

    /** Print alignment (in CLUSTAL W format) */
    public static void printCLUSTALW(Alignment a, PrintWriter out) {
        int n = 0;

        // CLUSTAL W header line
        out.println("CLUSTAL W multiple sequence alignment");
        out.println();

        // Print sequences
        while (n < a.getSiteCount()) {
            out.println();
            for (int s = 0; s < a.getSequenceCount(); s++) {
                format.displayLabel(out, a.getIdGroup().getIdentifier(s).getName(), 10);
                out.print("     ");

                printNextSites(a, out, false, s, n, 50);
                out.println();
            }
            // Blanks in status line are necessary for some parsers)
            out.println("               ");
            n += 50;
        }
    }

    public static void saveDelimitedAlignment(Alignment theAlignment, String delimit, String saveFile) {

        if ((saveFile == null) || (saveFile.length() == 0)) {
            return;
        }
        saveFile = Utils.addSuffixIfNeeded(saveFile, ".txt");
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {

            fw = new FileWriter(new File(saveFile));
            bw = new BufferedWriter(fw);

            bw.write("Taxa");
            int numSites = theAlignment.getSiteCount();
            for (int j = 0; j < numSites; j++) {
                bw.write(delimit);
                bw.write(String.valueOf(theAlignment.getPositionInLocus(j)));
            }
            bw.write("\n");

            for (int r = 0, n = theAlignment.getSequenceCount(); r < n; r++) {
                bw.write(theAlignment.getIdGroup().getIdentifier(r).getFullName());
                for (int i = 0; i < numSites; i++) {
                    bw.write(delimit);
                    bw.write(theAlignment.getBaseString(r, i));
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            System.out.println("AlignmentUtils: saveDelimitedAlignment: problem writing file: " + e.getMessage());
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    /**
     * Returns state indices for a sequence.
     */
    public static final void getAlignedSequenceIndices(Alignment a, int i, int[] indices, DataType dataType) {

        String sequence = a.getAlignedSequenceString(i);

        for (int j = 0; j < a.getSiteCount(); j++) {
            indices[j] = dataType.getState(sequence.charAt(j));
        }
    }

    /**
     * Unknown characters are given the state of -1
     */
    public static final int[][] getAlignedStates(Alignment base) {
        return getAlignedStates(base, -1);
    }

    public static final int[][] getAlignedStates(Alignment base, int unknownState) {
        int numberOfSites = base.getSiteCount();
        DataType dt = base.getDataType();
        int[][] sequences = new int[base.getSequenceCount()][base.getSiteCount()];
        for (int i = 0; i < sequences.length; i++) {
            for (int j = 0; j < sequences[i].length; j++) {
                char c = base.getBaseChar(i, j);
                if (dt.isUnknownChar(c)) {
                    sequences[i][j] = unknownState;
                } else {
                    sequences[i][j] = dt.getState(c);
                }
            }
        }
        return sequences;
    }

    private static interface CharacterMatrix {

        public int getSequenceCount();

        public int getSiteCount();

        public char getBase(int sequence, int site);
    }

    private static class AlignmentCM implements CharacterMatrix {

        Alignment a_;

        public AlignmentCM(Alignment a) {
            this.a_ = a;
        }

        public int getSequenceCount() {
            return a_.getSequenceCount();
        }

        public int getSiteCount() {
            return a_.getSiteCount();
        }

        public char getBase(int sequence, int site) {
            return a_.getBaseChar(sequence, site);
        }
    }

    private static class StringCM implements CharacterMatrix {

        String[] sequences_;

        public StringCM(String[] sequences) {
            this.sequences_ = sequences;
        }

        public int getSequenceCount() {
            return sequences_.length;
        }

        public int getSiteCount() {
            return sequences_[0].length();
        }

        public char getBase(int sequence, int site) {
            return sequences_[sequence].charAt(site);
        }
    }

    private static class CharacterCM implements CharacterMatrix {

        char[][] sequences_;

        public CharacterCM(char[][] sequences) {
            this.sequences_ = sequences;
        }

        public int getSequenceCount() {
            return sequences_.length;
        }

        public int getSiteCount() {
            return sequences_[0].length;
        }

        public char getBase(int sequence, int site) {
            return sequences_[sequence][site];
        }
    }

    /** count states
     *  @note Alignment state does not change! That is,
     *  does not create a data type for alignment, and does not
     *  set the frequencies internally for the alignment
     */
    public static double[] estimateFrequencies(Alignment a) {
        DataType dt = a.getDataType();
        if (dt == null) {
            dt = new TextDataType();
        }

        int numStates = dt.getNumStates();

        double[] frequency = new double[numStates];

        long[] stateCount = new long[numStates + 1];

        for (int i = 0; i < numStates + 1; i++) {
            stateCount[i] = 0;
        }

        for (int i = 0; i < a.getSequenceCount(); i++) {
            for (int j = 0; j < a.getSiteCount(); j++) {
                int state = dt.getState(a.getBaseChar(i, j));
                if (dt.isUnknownState(state)) {
                    state = stateCount.length - 1;
                }
                stateCount[state] += 1;
            }
        }

        // Compute frequencies suitable for RateMatrix (sum = 1.0)
        long sumStates = a.getSiteCount() * a.getSequenceCount() - stateCount[numStates];
        for (int i = 0; i < numStates; i++) {
            frequency[i] = (double) stateCount[i] / sumStates;
        }

        //a.setFrequency(frequency);

        return frequency;
    }

    public static final boolean isSiteRedundant(Alignment a, int site) {
        int numSeq = a.getSequenceCount();
        for (int i = 0; i < numSeq; i++) {
            if (!isGap(a, i, site)) {
                return false;
            }
        }
        return true;
    }

    public static final Alignment removeRedundantSites(Alignment a) {
        boolean[] keep = new boolean[a.getSiteCount()];
        int toKeep = 0;
        for (int i = 0; i < keep.length; i++) {
            keep[i] = !isSiteRedundant(a, i);
            if (keep[i]) {
                toKeep++;
            }
        }
        String[] newSeqs = new String[a.getSequenceCount()];
        int numberOfSites = a.getSiteCount();
        for (int i = 0; i < newSeqs.length; i++) {
            StringBuffer sb = new StringBuffer(toKeep);
            for (int j = 0; j < numberOfSites; j++) {
                if (keep[j]) {
                    sb.append(a.getBase(i, j));
                }
            }
            newSeqs[i] = sb.toString();
        }
        return SimpleAlignment.getInstance(a.getIdGroup(), newSeqs, a.getDataType());
    }

    /**
     * Returns true if the alignment has a gap at the site in the
     * sequence specified.
     */
    public static final boolean isGap(Alignment a, int seq, int site) {
        return a.getDataType().isGapChar(a.getBaseChar(seq, site));
    }

    public static final String cleanSequenceString(String s) {
        s = s.replace('?', DataType.UNKNOWN_CHARACTER);
        s = s.replace('N', DataType.UNKNOWN_CHARACTER);
        s = s.replace('n', DataType.UNKNOWN_CHARACTER);
        s = s.replace('X', DataType.UNKNOWN_CHARACTER);
        s = s.replace('x', DataType.UNKNOWN_CHARACTER);
        s = s.replace('-', DataType.GAP_CHARACTER);
        s = s.replace('~', DataType.GAP_CHARACTER);
        s = s.replace(' ', DataType.GAP_CHARACTER);
        s = s.replace('_', DataType.GAP_CHARACTER);
        return s;
    }

    public static final void cleanSequenceString(byte[] seq) {

        for (int i = 0; i < seq.length; i++) {
            if ((seq[i] == '?') || (seq[i] == 'N') || (seq[i] == 'n')
                    || (seq[i] == 'X') || (seq[i] == 'x')) {
                seq[i] = DataType.UNKNOWN_BYTE;
            }
            if ((seq[i] == '-') || (seq[i] == '~')
                    || (seq[i] == ' ') || (seq[i] == '_')) {
                seq[i] = DataType.GAP_BYTE;
            }
        }

    }

    /**
    @param startingCodonPosition - from {0,1,2}, representing codon position of first value in sequences...
    @note uses middle nucelotide of code to display info...
     */
    /*
    public static void getPositionMisalignmentInfo(Alignment a, PrintWriter out, int startingCodonPosition) {
    DataType dt = a.getDataType();

    for (int i = 0; i < a.getSequenceCount(); i++) {
    int codonPosition = startingCodonPosition;
    out.print(a.getIdGroup().getIdentifier(i) + ":");
    for (int j = 0; j < a.getSiteCount(); j++) {
    char c = a.getBaseChar(i, j);
    if (dt.isGapChar(c)) {
    out.print(c);
    } else {
    switch (codonPosition) {
    case 0: {
    out.print('[');
    break;
    }
    case 1: {
    out.print(c);
    break;
    }
    case 2: {
    out.print(']');
    break;
    }
    }
    codonPosition = (codonPosition + 1) % 3;
    }

    }
    out.print("\n");
    }
    }
     */
    /** Concatenates an array of alignments such that the resulting alignment is
     *		all of the sub alignments place along side each other
     */
    public static final Alignment concatAlignments(Alignment[] alignments, DataType dt) {
        int maxSequence = -1;
        Alignment maxAlignment = null;
        int length = 0;
        for (int i = 0; i < alignments.length; i++) {
            if (alignments[i].getSequenceCount() > maxSequence) {
                maxAlignment = alignments[i];
                maxSequence = alignments[i].getSequenceCount();
            }
            length += alignments[i].getSiteCount();
        }
        char[][] sequences = new char[maxSequence][length];
        for (int j = 0; j < sequences.length; j++) {
            int base = 0;
            for (int i = 0; i < alignments.length; i++) {
                if (alignments[i].getSequenceCount() <= j) {
                    for (int k = 0; k < alignments[i].getSiteCount(); k++) {
                        sequences[j][base + k] = DataType.GAP_CHARACTER;
                    }
                } else {
                    for (int k = 0; k < alignments[i].getSiteCount(); k++) {
                        sequences[j][base + k] = alignments[i].getBaseChar(j, k);
                    }
                }
                base += alignments[i].getSiteCount();
            }
        }
        SimpleAlignment sa = SimpleAlignment.getInstance(maxAlignment.getIdGroup(), sequences, dt);
        return sa;

    }

    /** Returns a particular sequence of an alignment as a char array */
    public static final char[] getSequenceCharArray(Alignment a, int sequence) {
        char[] cs = new char[a.getSiteCount()];
        for (int i = 0; i < cs.length; i++) {
            cs[i] = a.getBaseChar(sequence, i);
        }
        return cs;
    }

    /** Returns a particular sequence of an alignment as a String */
    public static final String getSequenceString(Alignment a, int sequence) {
        return new String(getSequenceCharArray(a, sequence));
    }

    /** Returns an alignment which follows the pattern of the input alignment
    except that all sites which do not contain states in dt (excluding the
    gap character) are removed. The Datatype of the returned alignment is dt
     */
    public static final Alignment getChangedDataType(Alignment a, DataType dt) {
        int numberOfSites = a.getSiteCount();
        boolean[] include = new boolean[numberOfSites];
        int goodSiteCount = 0;
        for (int i = 0; i < numberOfSites; i++) {
            include[i] = isGoodSite(a, dt, i);
            if (include[i]) {
                goodSiteCount++;
            }
        }
        //Yes, I'm aware it may be slightly faster to nest sequence
        // in site but it's easier to program this way
        String[] sequences = new String[a.getSequenceCount()];
        for (int i = 0; i < sequences.length; i++) {
            char[] seq = new char[goodSiteCount];
            int count = 0;
            for (int j = 0; j < numberOfSites; j++) {
                if (include[j]) {
                    seq[count] = a.getBaseChar(i, j);
                    count++;
                }
            }
            sequences[i] = new String(seq);
        }
        return SimpleAlignment.getInstance(SimpleIdGroup.getInstance(a.getIdGroup()), sequences, dt);
    }

    /** Tests the characters of an alignment to see if there are any characters that
    are not within a data type.
    @teturn the number of invalid characters
     */
    public static final int countUnknowns(Alignment a, DataType dt) {
        int count = 0;
        for (int i = 0; i < a.getSequenceCount(); i++) {
            for (int j = 0; j < a.getSiteCount(); j++) {
                if (dt.isUnknownState(dt.getState(a.getBaseChar(i, j)))) {
                    count++;
                }
            }
        }
        return count;
    }

    /*
    private static final void stripLeadingIncompleteCodon(int[] states, int unknownState) {
    int numberOfCodons = states.length / 3;
    final Nucleotides n = Nucleotides.DEFAULT_INSTANCE;
    for (int codon = 0; codon < numberOfCodons; codon++) {
    int unknownCount = 0;
    final int index = codon * 3;
    for (int i = 0; i < 3; i++) {
    if (n.isUnknownState(states[index + i])) {
    unknownCount++;
    }
    }
    if (unknownCount == 0) {
    return; //First codon is not incomplete
    } else if (unknownCount != 3) {
    //We have an incomplete codon on our hands!
    for (int i = 0; i < 3; i++) {
    states[index + i] = unknownState;
    }
    }
    }
    }

    // PRIVATE METHODS
    private static final void outputChar(PrintWriter out, char c, int number) {
    for (int i = 0; i < number; i++) {
    out.print(c);
    }
    }
     */
    private static void printNextSites(Alignment a, PrintWriter out, boolean chunked, int seq, int start, int num) {
        // Print next num characters
        for (int i = 0; (i < num) && (start + i < a.getSiteCount()); i++) {
            // Chunks of 10 characters
            if (i % 10 == 0 && i != 0 && chunked) {
                out.print(' ');
            }
            out.print(a.getBaseChar(seq, start + i));
        }
    }

    /**
     * Returns the gap creation costs between sequences x and y from site start to site finish.
     * @param a alignment
     * @param x first sequence
     * @param y second sequence
     * @param start first site to consider (inclusive)
     * @param finish last site to consider (inclusive)
     */
    /*
    private static int getNaturalGapCost(Alignment a, int x, int y,
    int start, int finish) {
    DataType dt = a.getDataType();
    int totalCost = 0;
    boolean inGap = false;

    // get gap creation costs
    for (int i = start; i <= finish; i++) {
    // if not a gap in one of them then consider column for x
    if (!(dt.isGapChar(a.getBaseChar(y, i)) && dt.isGapChar(a.getBaseChar(x, i)))) {
    // if gap in x then its the start of gap or already in gap
    if (isGap(a, x, i)) {
    // if not in gap then new gap
    if (!inGap) {
    totalCost += 1;
    inGap = true;
    } // else in gap and no extra cost
    } else {
    // else not null in x therefore not in gap
    inGap = false;
    }
    }
    }

    return totalCost;
    }
     */
//=================================================================
    private static final boolean isGoodSite(Alignment a, DataType dt, int site) {
        int numberOfSequences = a.getSequenceCount();
        for (int i = 0; i < numberOfSequences; i++) {
            char c = a.getBaseChar(i, site);
            if (c != DataType.GAP_CHARACTER && dt.isUnknownState(dt.getState(c))) {
                return false;
            }
        }
        return true;
    }

    /**
     * Fills a [length][numsequences] matrix with indices.
     * Each indices points to a position in the unaligned sequence, -1 means a gap.
     */
    private static int[][] getAlignmentIndices(Alignment a) {

        int[][] indices = new int[a.getSequenceCount()][a.getSiteCount()];
        DataType dataType = a.getDataType();

        for (int i = 0; i < a.getSequenceCount(); i++) {
            int seqcounter = 0;
            for (int j = 0; j < a.getSiteCount(); j++) {
                int index = dataType.getState(a.getBaseChar(i, j));
                if (index != dataType.getNumStates()) {
                    indices[i][j] = seqcounter;
                    seqcounter += 1;
                } else {
                    indices[i][j] = -1;
                }
            }
        }

        return indices;
    }

    /**
     * @return the total number of homology assignments between seqa and seqb in this alignment.
     */
    private static int getAlignmentMatches(int[][] indices, int seqa, int seqb) {

        int matches = 0;

        for (int i = 0; i < indices[seqa].length; i++) {
            if ((indices[seqa][i] != -1) && (indices[seqb][i] != -1)) {
                matches += 1;
            }
        }

        return matches;
    }

    /**
     * @return the number of homology matches between two sequences in
     * two alignments.
     */
    private static int getAlignmentMatches(int[][] indices1, int[][] indices2,
            int seqa1, int seqb1, int seqa2, int seqb2) {

        int matches = 0;
        int counter2 = 0;

        for (int i = 0; i < indices1[0].length; i++) {
            if ((indices1[seqa1][i] != -1) && (indices1[seqb1][i] != -1)) {
                try {
                    while (indices2[seqa2][counter2] != indices1[seqa1][i]) {

                        counter2 += 1;
                    }
                } catch (ArrayIndexOutOfBoundsException ae) {

                    for (int j = 0; j < indices1[0].length; j++) {
                        System.out.println("indice1[" + j + "]" + indices1[seqa1][j]);
                    }
                    for (int j = 0; j < indices2[0].length; j++) {
                        System.out.println("indice2[" + j + "]" + indices2[seqa2][j]);
                    }
                    System.out.println(
                            "indices1[" + seqa1 + "][" + i + "] = " + indices1[seqa1][i]
                            + "\tindices2[" + seqa2 + "][" + (counter2 - 1) + "] = "
                            + indices2[seqa2][counter2 - 1]);
                    System.out.println("counter2 = " + counter2);
                }

                if (indices1[seqb1][i] == indices2[seqb2][counter2]) {
                    matches += 1;
                }
            }
        }

        return matches;
    }

    public static int[][] parseQualityScores(String[] qualScores, int[] seqLengths) {

        int seqNum = qualScores.length;

        int[][] qScore = new int[seqNum][];
        for (int i = 0; i < seqLengths.length; i++) {
        }
        for (int seq = 0; seq < seqNum; seq++) {

            boolean hasHappenedBefore = false;

            qScore[seq] = new int[seqLengths[seq]];

            int site = 0;

            // if a particular quality score String does not exist, fill in the values with -1
            if (qualScores[seq] == null) {
                for (int i = 0; i < seqLengths.length; i++) {
                    qScore[seq][i] = -1;
                }
            } else {
                for (StringTokenizer st = new StringTokenizer(qualScores[seq]); st.hasMoreElements();) {

                    String s = ((String) st.nextElement()).trim();
                    int score = -9;

                    if (s.equals("-")) {
                        s = "-1";
                    }
                    try {
                        score = Integer.parseInt(s);
                    } catch (NumberFormatException nfe) {
                        // unfortunate... but don't tell me about it more than once per sequence
                        if (!hasHappenedBefore) {
                            System.err.println("The quality score - " + s + " could not be parsed as an int");
                            hasHappenedBefore = true;
                            nfe.printStackTrace();
                        }
                    }
                    qScore[seq][site++] = score;
                }
                if (qScore[seq].length != seqLengths[seq]) {
                    throw new IllegalStateException("The length of the sequence and the number of quality scores do not match.");
                }
            }
        }
        return qScore;
    }
}
