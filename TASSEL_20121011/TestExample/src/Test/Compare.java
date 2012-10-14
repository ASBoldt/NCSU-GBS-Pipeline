/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Test;

import java.util.Arrays;
import java.util.Comparator;




/**
 *
 * @author fl262
 */
public class Compare {
	public static void main (String args[]) {
		example[] e = new example[3];
		e[0] = new example (2, "wet");
		e[1] = new example (3, "mid");
		e[2] = new example (1, "hot");
		Arrays.sort (e, new sortbyname());
		for (int i = 0; i < 3; i++) {
			e[i].pr();
		}
		Arrays.sort (e, new sortbyint());
		for (int i = 0; i < 3; i++) {
			e[i].pr();
		}
	}
}
class example implements Comparable <example> {
	int i1;
	String s1;
	example (int i1, String s1) {
		this.i1 = i1;
		this.s1 = s1;
	}
	public void pr () {
		System.out.println (i1+" "+s1);
	}
	public int compareTo(example o) {
		if (i1 < o.i1) {
			return -1;
		}
		else {
			return 1;
		}
	}
}
class sortbyint implements Comparator <example> {
	public int compare(example o1, example o2) {
		return o1.i1 - o2.i1;
	}
}
class sortbyname implements Comparator <example> {
	public int compare (example o1, example o2) {
		return o1.s1.compareTo(o2.s1);
	}
}
