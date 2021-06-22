
public class primes {
	
	static String n_primes = find_primes(10005);
	
	public static String find_primes(int n) {
		String prime_string = "";
		int i = 2;
		while (prime_string.length() < n) {
				boolean is_div = false;
				for (int j = 2; j < i; j++) {
					if (i%j == 0) {
						is_div = true;
						break;
					}	
				}
				if (is_div == false) {
					prime_string = prime_string.concat(Integer.toString(i));
				}
				i++;
			}
		return prime_string;
	}
	
	public static String solution(int i) {

		return n_primes.substring(i, i+5);
    }
	
	
	public static void main(String[] args) {

		//System.out.println(find_primes(10005).length());
		//System.out.println(n_primes.length());
		System.out.println(solution(0));
		
	}
}
