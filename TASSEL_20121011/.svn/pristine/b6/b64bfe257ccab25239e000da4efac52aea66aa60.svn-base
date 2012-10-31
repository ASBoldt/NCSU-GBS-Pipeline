/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.DataOutputStream;
import java.util.Comparator;

/**
 *
 * @author fl262
 */
public class CorrRecord implements Comparable<CorrRecord> {

	int queryIndex, hitIndex;
	float rsquare;

	public CorrRecord(int queryIndex, int hitIndex, float rsquare) {
		this.queryIndex = queryIndex;
		this.hitIndex = hitIndex;
		this.rsquare = rsquare;
	}

	public void writeBinary(DataOutputStream dos) {
		try {
			dos.writeInt(queryIndex);
			dos.writeInt(hitIndex);
			dos.writeFloat(rsquare);
		} catch (Exception e) {
			System.out.println("Error occurred while writing binary correlation record " + e.toString());
		}
	}

	public int compareTo(CorrRecord o) {
		if (rsquare - o.rsquare < 0) {
			return -1;
		} else if (rsquare - o.rsquare > 0) {
			return 1;
		} else {
			return 0;
		}
	}
}

class sortByIndex implements Comparator<CorrRecord> {
	public int compare(CorrRecord o1, CorrRecord o2) {
		if (o1.queryIndex == o2.queryIndex) {
			return (o1.hitIndex - o2.hitIndex);
		}
		return o1.queryIndex - o2.queryIndex;
	}
}

