/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author Fei
 */
public class Count_dis {
	cd[] hr;
	public Count_dis (String rcFileS, String outfileS) {
		ReadCounts rc = new ReadCounts(rcFileS, true);
		this.getCd(rc);
		this.write(outfileS);
	}

	public void write (String outfileS) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
			for (int i = 0; i < hr.length; i++) {
				bw.write(String.valueOf(hr[i].cnt)+"\t"+String.valueOf(hr[i].dis));
				bw.newLine();
			}
			bw.flush();
			bw.close();
		}
		catch (Exception e) {

		}
	}

	public void getCd (ReadCounts rc) {
		int[] counts = new int[rc.getSize()];
		for (int i = 0; i < counts.length; i++) {
			counts[i] = rc.hapcount[i];
		}
		Arrays.sort(counts);
		ArrayList<cd> cdList = new ArrayList();
		int temp = counts[0];
		int cnt = 1;
		for (int i = 1; i < counts.length; i++) {
			if (counts[i] == temp) {
				cnt++;
			}
			else {
				cdList.add(new cd(temp, cnt));
				temp = counts[i];
				cnt = 1;
			}
		}
		cdList.add(new cd(temp, cnt));
		hr = cdList.toArray(new cd[cdList.size()]);
	}

	public static void main (String[] args) {
		String CombinedTags = "M:/Generator/Collapsed_backup/parents/P1-U518_high.cnt";
		String outfileS = "M:/out.txt";
		new Count_dis (CombinedTags, outfileS);
	}

	class cd {
		int cnt, dis;
		cd (int cnt, int dis) {
			this.cnt = cnt;
			this.dis = dis;
		}

	}
}
