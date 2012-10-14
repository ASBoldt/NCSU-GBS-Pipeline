/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
 *
 * @author Fei Lu
 */
public class TagPair {
	long[][] haplotype;
	int[] order;

	public TagPair (String infileS, String outfileS) {
		this(infileS);
		this.sortBySeq();
		this.writeTagPair(outfileS);
	}

	public TagPair (String infileS) {
		File f = new File(infileS);
        int nHap =(int)(f.length()/20);  //20 is the number of bytes 2 longs (2*8b) + 1 int (4b)
        haplotype = new long[2][nHap];
        order = new int[nHap];
		this.readTagPair(infileS, nHap);
	}

	public void readTagPair(String infileS, int nHap) {
		try {
			DataInputStream dis = new DataInputStream(new BufferedInputStream (new FileInputStream(infileS),65536));
			for(int i = 0; i < nHap; i++) {
               haplotype[0][i]=dis.readLong();
               haplotype[1][i]=dis.readLong();
               order[i] = dis.readInt();
            }
			System.out.println(nHap/2 + " TagPair are read");
		}
		catch (Exception e) {
			System.out.println("Error occurred while reading " + infileS + " " + e.toString());
		}
	}

	public void writeTagPair(String outfileS) {
		try {
			DataOutputStream dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(outfileS),65536));
			for(int i = 0; i < this.haplotype[0].length; i++) {
				dos.writeLong(haplotype[0][i]);
				dos.writeLong(haplotype[1][i]);
				dos.writeInt(order[i]);
            }
			dos.flush();
			dos.close();
			System.out.println(this.haplotype[0].length/2 + " TagPair are written");
		}
		catch (Exception e) {
			System.out.println("Error occurred while writing " + outfileS + " " + e.toString());
		}
	}

	public void sortBySeq () {
		GenericSorting.quickSort(0, haplotype[0].length, compSeq, swapper);
		System.out.println("TagPair is sorted by sequence");
	}

	public void sortByOrder () {
		GenericSorting.quickSort(0, haplotype[0].length, compOrder, swapper);
		System.out.println("TagPair is sorted by pair order");
	}

	Swapper swapper = new Swapper() {
		public void swap(int a, int b) {
			long t1, t2;
			t1 = haplotype[0][a]; haplotype[0][a] = haplotype[0][b];	haplotype[0][b] = t1;
			t2 = haplotype[1][a]; haplotype[1][a] = haplotype[1][b]; haplotype[1][b] = t2;
			int t3;
			t3=order[a]; order[a]=order[b]; order[b]=t3;
		}
	};

	IntComparator compSeq = new IntComparator() {
		public int compare(int a, int b) {
			if (haplotype[0][a] == haplotype[0][b]) return haplotype[1][a] == haplotype[1][b] ? 0 : (haplotype[1][a] < haplotype[1][b] ? -1 : 1);
			return haplotype[0][a] < haplotype[0][b] ? -1 : 1;
		}
	};

	IntComparator compOrder = new IntComparator() {
		public int compare(int a, int b) {
			return order[a] < order[b] ? -1 : 1;
		}
	};
}
