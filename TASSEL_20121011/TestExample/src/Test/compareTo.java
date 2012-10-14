/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Test;

import java.util.Arrays;
/**
 *
 * @author fl262
 */
public class compareTo {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String str = "sdf";
        tarray[] ob = new tarray[3];
        ob[0] = new tarray (2, 3, "ere");
        ob[1] = new tarray (1, 5, "rer");
        ob[2] = new tarray (3, 4, "fe");
        Arrays.sort(ob);
        
        for (int i = 0; i < 3; i++) {
            System.out.println (ob[i].a+"\t"+ob[i].b+"\t"+ob[i].c);
        }
    }
}
class tarray implements Comparable <tarray> {
    int a;
    int b;
    String c;
    tarray (int a, int b, String c) {
        this.a = a;
        this.b = b;
        this.c = c;
    }

    public int compareTo (tarray o) {
        if (a < o.a) {
            return -1;
        }
        else {
            return 1;
        }
    }
}
