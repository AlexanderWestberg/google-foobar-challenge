package task2;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;

public class gen_perm {

	static ArrayList<String> perms = new ArrayList<String>();
	static Queue<String> que = new LinkedList<String>();

	
	public static void main(String[] args) {
        int[] n = {3, 1, 4, 1, 5, 9};
        
        System.out.println(Solution(n));
    }
	
	public static int Solution(int[] l) {
		String s = list_to_str(l);
		permute(s, 0, s.length()-1);
		
		if (perms.size() != 0) {
			return find_max();
		}
		else {
			for (String i: find_combinations(s)) {
				que.add(i);
			}
			less_perms();
			return find_max();
		}
	}
	
	public static void less_perms() {
			for (int j = 0; j < que.size(); j++) {
				String temp = que.remove();
				permute(temp, 0, temp.length()-1);
				String[] combs = find_combinations(temp);
				for (int k = 0; k < combs.length; k++) {
					que.add(combs[k]);
				} 
			}
			
			if (perms.size() == 0) {
				for (int i = 0; i < que.size(); i++) {
					
				}
			}
			
	}
	
	public static String[] find_combinations(String a) {
		String[] s = new String[a.length()];
		for (int i = 0; i < a.length(); i++) {
			s[i] = charRemoveAt(a, i);;
		}
		return s;
	}
	
	
	public static String charRemoveAt(String str, int p) {  
        return str.substring(0, p) + str.substring(p + 1);  
     } 
	
	
	public static int find_max () {
		int max = 0;
		for (String i : perms) {
			int a = Integer.parseInt(i); 
			if (a > max) {
				max = a;
			}
		}
		return max;
	}
	
	public static void is_div_by_3 (String i) {
		if (Integer.parseInt(i) % 3 == 0) {
			perms.add(i);
		}	
	
	}
	
	public static String list_to_str (int[] a) {
		String s = "";
		for (int i : a) {
			s += Integer.toString(i);
		}
		return s;
	} 
 
    /**
    * permutation function
    * @param str string to calculate permutation for
    * @param l starting index
    * @param r end index
    */
    public static void permute(String str, int l, int r)
    {
        if (l == r)
        	is_div_by_3(str);
        else
        {
            for (int i = l; i <= r; i++)
            {
                str = swap(str,l,i);
                permute(str, l+1, r);
                str = swap(str,l,i);
            }
        }
    }
 
    /**
    * Swap Characters at position
    * @param a string value
    * @param i position 1
    * @param j position 2
    * @return swapped string
    */
    public static String swap(String a, int i, int j)
    {
        char temp;
        char[] charArray = a.toCharArray();
        temp = charArray[i] ;
        charArray[i] = charArray[j];
        charArray[j] = temp;
        return String.valueOf(charArray);
    }

}
