
public class doomsday_fuel {
	
	//
	static int terminal[] = new int[10]; // See which states that are stable
	static int terminal_count = 0; // Count how many stable states there are
	
	// 
	static int[] denom = new int[10]; // A list for the sum of each state
	static int final_denom = 1; // The output denominator
	
	// 
	static int[] term_hash = new int[10]; // A "hashlist" that maps the terminal 
	//states to the indicies in the outputlist
	static int[] output; // An output list, where the last item is the denominator
 	
	
	
	public static void state_flow_non_inf(int[][] m) {
		int n_count = 0;
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m.length; j++) {
				 denom[i] += m[i][j];
			}
			if (denom[i] != 0) {
				final_denom *= denom[i];
				//terminal[i] = 0;
			} else if (denom[i] == 0) {
				//terminal[i] = 1;
				terminal_count ++;
				term_hash[i] = n_count; 
				n_count ++;
			}
		}
		output = new int[terminal_count + 1];
		output[output.length - 1] = final_denom;
	}
	
	public static void state_flow_inf(int[][] m) {
		
	}
	
	
	public static boolean inf_path(int[][] m) {
		for (int i = 1; i < m.length; i++) {
			if (m[i][0] != 0) {
				return true;
			} 	
		} 
		return false;
	}
	
	public static int[] ones() {
		int[] out = new int[10];
		
		for (int i = 0; i < out.length; i++) {
			out[i] = 1;
		}
		return out;
	}
	
	public static int[] solution(int[][] m) {
		
		int[] ev_state = ones();
		boolean paths = inf_path(m);
		
		if (paths == false) {
			state_flow_non_inf(m);
			
			for (int i = 0; i < m.length; i++) {
				for (int j = 0; j < m.length; j++) {
					if (denom[j] == 0 && m[i][j] > 0) {
						if (i == 0) {
						output[term_hash[j]] += m[i][j]*final_denom/denom[i];
						} else {
							output[term_hash[j]] += m[i][j]*ev_state[i];
						}
					} else if (denom[j] != 0 && m[i][j] > 0) {
						ev_state[j] = m[i][j];
					}
				}
			}
		} else if (paths == true) {
			
		}
		
		return  output;
	}
	
	public static void printList(int[] list) {
		for (int i : list) {
			System.out.print(i);
			System.out.print(" ");
			
			
		}
	}
	
	public static void main(String[] args) {
		int[][] one = {{0, 2, 1, 0, 0}, {0, 0, 0, 3, 4}, {0, 0, 0, 0, 0},
				{0, 0, 0, 0,0}, {0, 0, 0, 0, 0}};
		int[][] two  = {{0, 1, 0, 0, 0, 1}, {4, 0, 0, 3, 2, 0},
				{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
		
		//System.out.println(solution(one));
		//System.out.println(one.length);
		//System.out.println(solution(two));
		//printList(solution(one));
		printList(solution(two));
		
		
		
	}
}
