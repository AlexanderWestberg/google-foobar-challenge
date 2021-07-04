
public class lovely_lucky_lambs {
	
public static int fibonacci(int k) {
		
		int sum = 0;
		int n = 0;
		int[] list = new int[100];
		
		if (k > 2) {
			list[0] = 1;
			list[1] = 1;
			n += 1;
			sum += list[n] + list[n-1]; 
			
			while (true) {
				if (k-sum >= list[n] + list[n-1] ) {
					list[n+1] = list[n] + list[n-1];
					sum += list[n] + list[n-1]; 
					n++;	
				} else {
					break;
				}
			}
		} 
		else if (k == 2)  {
			list[0] = 1;
			list[1] = 1;
			n++;
			
		}
		else if (k == 1) {
			list[0] = 1;
		}
		
		return n+1;
	} 
	
	public static int exponetial(int k) {

		int sum = 0;
		int n = 0;
		int[] list = new int[100];
		
		if (k > 2) {
			list[0] = 1;
			list[1] = 2;
			n += 1;
			sum += sum += list[n] + list[n-1]; 
			
			while (true) {
				if (k-sum > list[n] + list[n-1] ) {
					if (k-sum > list[n]*2) {
						list[n+1] = list[n]*2;
						sum += list[n]*2; 
						n++;
					} else {
						list[n+1] = list[n] + list[n-1];
						sum += list[n] + list[n-1]; 
						n++;	
					}
					
				} else {
					break;
				}
			}
		} else if (k == 2 || k == 1 )  {
			list[0] = 1;
		}
		return n+1;
	}

	public static int solution(int total_lambs) {
		return fibonacci(total_lambs) - exponetial(total_lambs);
	}
}
