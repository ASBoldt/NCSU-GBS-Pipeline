/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Test;

/**
 *
 * @author fl262
 */
public class inheri {
	public static void main(String[] args) {
		c1 example = new c1 ();
		example.MeA();
		c2 e2 = new c2 (example);
		e2.MeA();
	}
}
interface pub {
	void MeA ();
}
class c1 implements pub {
	int a = 0;
	c1 () {
		a = 1;
	}
	public void MeA () {
		System.out.println ("c1");
	}
}
class c2 implements pub {
	c1 su;
	c2 (c1 g) {
		su = g;
	}
	public void MeA () {
		System.out.println (su.a);
	}
}

