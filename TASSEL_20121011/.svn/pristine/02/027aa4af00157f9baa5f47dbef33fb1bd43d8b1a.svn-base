/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.genome.GBS.switchgrass;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.TreeSet;
import net.maizegenetics.genome.BaseEncoder;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;

/**
 *
 * @author fl262
 */
public class test {

	public test () {
		this.countCoverSite();
	}
	public static void main (String[] args) {
		new test();

	}
	public double getLd(double[] a, double[] b) {
		double ld2;
		double hapfre = this.getHapfre(a, b);
		double[] frea = this.getFrequency(a);
		double[] freb = this.getFrequency(b);
		ld2 = Math.pow(hapfre-frea[0]*freb[0], 2)/frea[0]/frea[1]/freb[0]/freb[1];
		return ld2;
	}
	public double[] getFrequency (double[] a) {
		double[] fre = new double[2];
		for (int i = 0; i < a.length; i++) {
			if (a[i] == 1) {
				fre[0]++;
			}
		}
		fre[0] = fre[0]/a.length;
		fre[1] = 1-fre[0];
		return fre;
	}
	public double getHapfre (double[] a, double[] b) {
		double hapfre = 0;
		for (int i = 0; i < a.length; i++) {
			if (a[i] == 1 && b[i] == 1) {
				hapfre++;
			}
		}
		return hapfre/a.length;
	}
	public double[] getGenoValue2 (int size, double maf) {
		double[] a = new double[size];
		double bound1 = maf*maf;
		double bound2 = Math.pow(1-maf, 2) + bound1;
		for (int i = 0; i < size; i++) {
			double r = Math.random();
			if (r >= 0 && r < bound1) {
				a[i] = 2;
			}
			else if (r < bound2) {
				a[i] = 0;
			}
			else {
				a[i] = 1;
			}
		}
		return a;
	}
	public double[] getGenoValue (int size, double maf) {
		double[] a = new double[size];
		for (int i = 0; i < size; i++) {
			if (Math.random() < maf) {
				a[i] = 1;
			}
			else {
				a[i] = 0;
			}
		}
		return a;
	}
    
    public void writeFasta () {
        String tagFileS = "M:/Generator/CombinedTags/CombinedTags.tc";
        String fastaFileS = "M:/CombinedTags.fasta";
        ReadCounts rc = new ReadCounts (tagFileS, true);
        int cnt = 0;
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (fastaFileS), 65536);
           
            int c = 1;
            for (int i = 0; i < rc.getSize(); i++) {
                long[] tag = rc.getRead(i);
                String seq = BaseEncoder.getSequenceFromLong(tag);
                if (seq.endsWith("A")) {
                    cnt++;
                }
                else {
                    bw.write(">"+String.valueOf(c));
                    bw.newLine();
                    bw.write("G"+seq);
                    bw.newLine();
                }
                c++;
            }
            bw.flush();
            bw.close();
            
        }
        catch (Exception e) {
            
        }
        System.out.println(rc.getSize());
        System.out.println(cnt);
    }
    
    public void countCoverSite () {
        String infileS = "M:/out.txt";
        int cnt = 0;
         TreeSet<String> tSet = new TreeSet();
        try {
            BufferedReader br = new BufferedReader (new FileReader(infileS), 65536);
            for (int i = 0; i < 15; i++) br.readLine();
            String temp;
           
           
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                //if (Integer.valueOf(tem[1]) != 0) continue;
                if (tem[2].equals("*")) continue;
                int chr = Integer.valueOf(tem[2]);
                if (chr > 10 || chr < 1) continue;
                if (!tem[5].equals("65M")) continue;
                if (!temp.contains("NM:i:0"))  continue;
                tSet.add(tem[0]);
                cnt++;
            }
        }
        catch (Exception e) {
            System.out.println (e.toString());
        }
        System.out.println(cnt);
        System.out.println(tSet.size());
    }

}
