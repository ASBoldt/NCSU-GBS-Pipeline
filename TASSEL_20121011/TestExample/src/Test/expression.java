/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Test;

import java.util.regex.*;

/**
 *
 * @author fl262
 */
public class expression {

	public static void main(String args[]) {
		String a = "abcdef";
		Pattern p = Pattern.compile("((.).).");
		Matcher m = p.matcher(a);

		a = a.replaceAll("((.).).", "$2");


		System.out.println(a);

/*		if (m.find(0)) {
			System.out.println("Yes");
			System.out.println(m.start());
			System.out.println(m.group());
			System.out.println("number is " + m.groupCount());

			System.out.println(m.group(0)); // the whole regex
			System.out.println(m.group(1));
			System.out.println(m.group(2));
			System.out.println(m.toString());
		}
		if (m.find()) {
			System.out.println("Yes");
			System.out.println(m.start());
			System.out.println(m.group());
			for (int i = 0; i <= m.groupCount(); i++) {
				System.out.println(m.group(i));
			}
			System.out.println(m.toString());
		}
 */
		
	}
}

