/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.plasmid;

/**
 *
 * @author fl262
 */
public class baseEncoder {
	public static final char[] bases={'A','C','G','T'};
	public static final int intSize = 16;
	public static final int byteSize = 4;
	public static final int byteLen = 19;
	public static final int seqLen = byteSize * byteLen;

	public static byte getByteFromSeq (String seq) {
		byte v = 0;
		byte b = 0;
		for (int i = 0; i < seq.length(); i++) {
			b = getTwoBitBase(seq.charAt(i));
			if (b == -1) {
				return -1;
			}
			v = (byte) ((v << 2) + b);
		}
		if (seq.length() == byteSize) {
			return v;
		}
		if (seq.length() > intSize) {
			return -1;
		}
		v = (byte) (v << (2 * (intSize - seq.length()))); //if shorter fill with AAAA
		return v;
	}
	public static byte[] getByteArrayFromSeq (String seq) {
		byte[] byteArray = new byte[byteLen];
		for (int i = 0; i < byteArray.length; i++) {
			byteArray[i] = getByteFromSeq(seq.substring(i * byteSize, (i + 1) * byteSize));
		}
		return byteArray;
	}
	public static int getFirstLowQualityPos(String quality, int minQual) {
		int qualInt = 0;
		for (int i = 0; i < quality.length(); i++) {
			qualInt = (int) quality.charAt(i) - 64;
			if (qualInt < minQual) {
				return i;
			}
		}
		return quality.length();
	}

	public static byte getTwoBitBase(char base) {
		byte b = -1;
		switch (base) {
			case 'A':
				return 0;  //00
			case 'C':
				return 1;  //01
			case 'G':
				return 2;  //10
			case 'T':
				return 3;  //11
			case 'a':
				return 0;
			case 'c':
				return 1;
			case 'g':
				return 2;
			case 't':
				return 3;
		}
		return b;
	}

	public static int getIntFromSeq(String seq) {
		int v = 0;
		byte b = 0;
		for (int i = 0; i < seq.length(); i++) {
			b = getTwoBitBase(seq.charAt(i));
			if (b == -1) {
				return -1;
			}
			v = (v << 2) + b;
		}
		if (seq.length() == intSize) {
			return v;
		}
		if (seq.length() > intSize) {
			return -1;
		}
		v = (v << (2 * (intSize - seq.length()))); //if shorter fill with AAAA
		return v;
	}
	public static String getSeqFromByte(byte byteSeq) {
		StringBuilder seq = new StringBuilder();
		byte mask = 3;
		for (int i = 0; i < byteSize; i++) {
			byte base = (byte) (byteSeq & mask);
			seq.append(bases[base]);
			byteSeq = (byte) (byteSeq >> 2);
		}
		seq.reverse();
		return seq.toString();
	}
	public static String getSeqFromByteArray (byte[] byteArray) {
		StringBuilder seq = new StringBuilder();
		for (int i = 0; i < byteArray.length; i++) {
			seq.append(getSeqFromByte(byteArray[i]));
		}
		return seq.toString();
	}
	public static byte seqDifferencesForSubset(int seq1, int seq2, int lengthOfComp, int maxDivergence) {
		int mask = 3;
		byte cnt = 0;
		int diff = seq1 ^ seq2;
		diff = diff >> (2 * (intSize - lengthOfComp));  //shift to 5' end of sequence
		for (int x = 0; x < lengthOfComp && cnt < maxDivergence; x++) {
			if ((diff & mask) > 0) {
				cnt++;
			}
			diff = diff >> 2;
		}
		return cnt;
	}
	public static boolean ifSeqSame(int seq1, int seq2, int lengthOfComp) {
		int mask = 3;
		int diff = seq1 ^ seq2;
		diff = diff >> (2 * (intSize - lengthOfComp));
		for (int i = 0; i < lengthOfComp; i++) {
			if ((diff & mask) > 0) {
				return false;
			}
			diff = diff >> 2;
		}
		return true;
	}
}
