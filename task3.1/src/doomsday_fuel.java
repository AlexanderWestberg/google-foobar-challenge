
public class doomsday_fuel {
	// Inverse Functions
	// ------------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------
	// ------------------------------------------------------------------------------------------
	// Dimension of input square matrix
    // Function to get cofactor of
    // mat[p][q] in temp[][]. n is
    // current dimension of mat[][]
    static void getCofactor(Fraction mat[][], Fraction temp[][],
                            int p, int q, int n)
    {
        int i = 0, j = 0;
        // Looping for each element of
        // the matrix
        for (int row = 0; row < n; row++)
        {
            for (int col = 0; col < n; col++)
            {
                // Copying into temporary matrix
                // only those element which are
                // not in given row and column
                if (row != p && col != q)
                {
                    temp[i][j++] = mat[row][col];
                    // Row is filled, so increase
                    // row index and reset col
                    // index
                    if (j == n - 1)
                    {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }
 
    /* Recursive function for finding determinant
    of matrix. n is current dimension of mat[][]. */
    static Fraction determinantOfMatrix(Fraction mat[][], int n)
    {
        Fraction D = new Fraction(0,0); // Initialize result
        
        // Base case : if matrix contains single
        // element
        if (n == 1)
            return mat[0][0];
 
        // To store cofactors
        Fraction temp[][] = new Fraction[mat.length][mat.length];
 
        // To store sign multiplier
        Fraction sign = new Fraction(1,1);
 
        // Iterate for each element of first row
        for (int f = 0; f < n; f++)
        {
            // Getting Cofactor of mat[0][f]
            getCofactor(mat, temp, 0, f, n);
            D  = add_frac(D, mul_frac(mul_frac(sign,  mat[0][f])
                 ,determinantOfMatrix(temp, n - 1)) );
            
            // terms are to be added with
            // alternate sign

             sign = mul_frac(sign ,  new Fraction(-1, 1));
        }
 
        return D;
    }
 
    /* function for displaying the matrix */

    
	static void adjoint(Fraction A[][],Fraction[][] adj)
	{
	    if (A.length == 1)
	    {
	        adj[0][0] = new Fraction(1,1);
	        return ;
	    }
	 
	    // temp is used to store cofactors of A[][]
	    int sign = 1;
	    Fraction [][]temp = new Fraction[A.length][A.length];
	 
	    for (int i = 0; i < A.length; i++)
	    {
	        for (int j = 0; j < A.length; j++)
	        {
	            // Get cofactor of A[i][j]
	            getCofactor(A, temp, i, j, A.length);
	 
	            // sign of adj[j][i] positive if sum of row
	            // and column indexes is even.
	            sign = ((i + j) % 2 == 0)? 1: -1;
	 
	            // Interchanging rows and columns to get the
	            // transpose of the cofactor matrix
	            
	            Fraction sign_frac; 
	            
	            if  (sign  == 1) {
	            	 sign_frac = new Fraction(1,1);
	            } else {
	            	 sign_frac = new Fraction(-1,1);
	            }
	            
	            adj[j][i] = mul_frac(sign_frac, determinantOfMatrix(temp, A.length-1));
	        }
	    }
	}
	
	static boolean inverse(Fraction A[][], Fraction[][] inverse)
	{
	    // Find determinant of A[][]
	    Fraction det = determinantOfMatrix(A, A.length);
	    if (det.num == 0)
	    {
	        System.out.print("Singular matrix, can't find its inverse");
	        return false;
	    }
	 
	    // Find adjoint
	    Fraction[][] adj = new Fraction[A.length][A.length];
	    adjoint(A, adj);
	 
	    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	    for (int i = 0; i < A.length; i++)
	        for (int j = 0; j < A.length; j++)
	             inverse[i][j] = div_frac(adj[i][j],det); 
	 
	    return true;
	}
    
    
    
 
    // Driver code
    public static void main(String[] args)
    {
    	          
        /*                 
    	int[][] one = {{0, 2, 1, 0, 0}, {0, 0, 0, 3, 4}, {0, 0, 0, 0, 0},
				{0, 0, 0, 0,0}, {0, 0, 0, 0, 0}};
    	
    	int[][] two  = {{0, 1, 0, 0, 0, 1}, {4, 0, 0, 3, 2, 0},
				{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
				
    	
    	int[][] three  = {{1, 2, 3, 0, 0, 0}, {4, 5, 6, 0, 0, 0},
				{7, 8, 9, 1, 0, 0}, {0, 0, 0, 0, 1, 2},
				{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
    	
    	int[][] four = {{0}};
    	
    	
    	int[][] five = {{0, 0, 12, 0, 15, 0, 0, 0, 1, 8},
    	                {0, 0, 60, 0, 0, 7, 13, 0, 0, 0},
    	                {0, 15, 0, 8, 7, 0, 0, 1, 9, 0},
    	                {23, 0, 0, 0, 0, 1, 0, 0, 0, 0},
    	                {37, 35, 0, 0, 0, 0, 3, 21, 0, 0},
    	                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    	                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    	                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    	                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    	                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    	
    
    	
    	 
    	 
    	
    	Fraction[][] two_frac = int_frac_matrix(five);
    	
    	getQnR(two_frac);
    	getF();
    	display(F);

    	
    	double[] a = {1.823440896E12/1.215627264E12   ,6.07813632E11/9.11720448E11    ,3.646881792E11/3.646881792E11    ,3.646881792E11/1.823440896E12    ,2.9175054336E12/3.646881792E12    };
    	
    	for (double i : simplify(a)) {
    		System.out.print(i+ "  ");
    	}
    	*/
    		
    }
    
    public static int[] solution (int[][] m) {
    	
    	if (m.length == 1) {
    		int[] output = {1,1};
    		return output;
    	
    	} else {
	    	Fraction[][] two_frac = int_frac_matrix(m);
	    	getQnR(two_frac);
	    	getF();
	    	Fraction[] f_simps = simplify_all_fracs(matrix_mul_frac(F,R)[0]);
	    	f_simps = find_common_divider(f_simps);
	    	f_simps = simplify_common_den(f_simps);
	    	int[] output = new int[f_simps.length +1];
	    	
	    	for (int i = 0; i < f_simps.length; i++) {
	    		output[i] = (int)f_simps[i].num;
	    	}
	    	output[f_simps.length] = (int)find_den(f_simps);
	    	
	    	
    	
    	return output;
    	}
    }
    
    public static double find_den(Fraction[] f) {
    	
    	
    	for (int i = 0; i < f.length; i++) {
    		if (f[i].den != 0) {
    			return f[i].den;
    		}
    	}
    	return 0;
    }
    
    
    	// Transition Matrix Functions
	   //--------------------------------------------------------------------
	   //--------------------------------------------------------------------
	   //--------------------------------------------------------------------
    

	public static Fraction[][] int_frac_matrix (int[][] int_m ) {
		
		
		Fraction[][] m = new Fraction[int_m.length][int_m.length];
		sum_for_each_state(int_m);
		
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m.length; j++) {
				m[i][j] = new Fraction(int_m[i][j], sums_of_states[i]);
			}
		}
		return m;
	}
	
	static int non_termina_states = 0;
	static int[] sums_of_states;
	static Fraction[][] R;
	static Fraction[][] Q;
	static Fraction[][] F;
	
	public static void getQnR(Fraction[][] m) { //---------------------- Assumes that the states are ordered
		
		
		Q = new Fraction[non_termina_states][non_termina_states];
		R = new Fraction[non_termina_states][m.length - non_termina_states];
		
		for  (int i = 0; i<non_termina_states; i++) {
			for  (int j = 0; j<non_termina_states; j++) {
				 Q[i][j] = m[i][j];
			}
		}
		
		int count_i = 0;
		int count_j = 0;
	
		for  (int i = 0; i<non_termina_states; i++) {
			count_j = 0;
			for  (int j = non_termina_states; j< m.length; j++) {
				 R[count_i][count_j] = m[i][j];
				 count_j ++;
			}
			count_i ++;
		}
		
		Q = simplify_matrix_fracs(Q);
		R = simplify_matrix_fracs(R);
	}
	
	public static void getF() {  
		Fraction[][] I = new Fraction[Q.length][Q.length];
		for (int i = 0; i<Q.length; i++) {
			for (int j = 0; j<Q.length; j++) {
				if (i == j) {
					I[i][j] = new Fraction(1,1);
				} else {
					I[i][j] = new Fraction(0,0);
				}
			}
		}
		
		Fraction[][] I_minus_Q = matrix_sub_frac(I,Q);
		F = new Fraction[I_minus_Q.length][I_minus_Q.length];
		
	    inverse(I_minus_Q, F);
	    //F = simplify_matrix_fracs(F);
	}
	
	
	public static void  sum_for_each_state(int[][] m) {
		
		
		sums_of_states = new int[m.length];
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m.length; j++) {
				sums_of_states[i] += m[i][j];
			}
			if (sums_of_states[i] != 0) {
				non_termina_states++;
			}
		}
	}
    
    
    	// Fraction Functions
	   //--------------------------------------------------------------------
	   //--------------------------------------------------------------------
	   //--------------------------------------------------------------------
    
    
	public static class Fraction{
		
		
	    public double num;
	    public double den;
	    public Fraction(double numerator, double denominator) {
	        num = numerator;
	        den = denominator;
	    }
	}
	
	public static Fraction div_frac (Fraction a, Fraction b ) {
		
	Fraction c; //initialize
	
	if ( a.num*b.den < 0 && a.den*b.num < 0) {
		c = new Fraction(-a.num*b.den, -a.den*b.num);
	} else {
		c = new Fraction(a.num*b.den, a.den*b.num);
	}
	
	return c;
}

	public static Fraction mul_frac (Fraction a, Fraction b ) {
			
			if (a.num*b.num == 0 || a.den*b.den == 0) {
				
				return new Fraction(0, 0);
				
			} else if ( a.num*b.num < 0 && a.den*b.den < 0) {
				Fraction c  = new Fraction(-a.num*b.num, -a.den*b.den);
				return c;
			} else {
			
				 Fraction c  = new Fraction(a.num*b.num, a.den * b.den);
				 return c;
			 }
	}
	
	
	public static Fraction add_frac(Fraction t1, Fraction t2) {
		
		if (t1.den < 0) {
			t1 = new Fraction(-t1.num, -t1.den);
		} else if (t2.den < 0) {
			t1 = new Fraction(-t2.num, -t2.den);
		}
		
		
		if (t1.den == 0 || t1.num == 0) {
			return t2;
		} else if (t2.den == 0 || t2.num == 0) {
			return t1;
		}  else {
		

			double numerator;
			double denominator;
			
			if (t1.den % t2.den == 0) {
				numerator = t1.num + t2.num*t1.den/t2.den; 
				denominator = t1.den;
			} else if (t2.den % t1.den == 0) {
				numerator = t2.num + t1.num*t2.den/t1.den; 
				denominator = t2.den;
			} else {
				numerator = t1.num*t2.den + t2.num*t1.den; 
				denominator = t1.den*t2.den;
	    	}
			Fraction out = new Fraction(numerator , denominator);
			return out;
		}
	}
	
	public static Fraction sub_frac(Fraction t1, Fraction t2) {
		
		// Change the sign of the fraction
		if (t1.den < 0) {
			t1 = new Fraction(-t1.num, -t1.den);
		} else if (t2.den < 0) {
			t1 = new Fraction(-t2.num, -t2.den);
		}
		
		
		if (t1.den == 0 || t1.den == 0) {
			return new Fraction(-t2.num, t2.den);
		} else if (t2.den == 0 || t2.num == 0) {
			return t1;
		}  else {
			
			
			double numerator;
			double denominator;
			
			if (t1.den % t2.den == 0) {
				numerator = t1.num - t2.num*t1.den/t2.den; 
				denominator = t1.den;
				
			} else if (t2.den % t1.den == 0) {
				numerator =   t1.num*t2.den/t1.den - t2.num;  
				denominator = t2.den;
				
			} else {
				numerator = t1.num*t2.den - t2.num*t1.den; 
				denominator = t1.den*t2.den;
	    	}
			Fraction out = new Fraction(numerator , denominator);
			return out;
			}
	}

	public static Fraction[][] matrix_sub_frac(Fraction[][] m1, Fraction[][] m2 ) {
		
		Fraction[][] sub_m = new Fraction[m1.length][m1[0].length];
		
		for (int i = 0; i < m1.length; i++) {
			for (int j = 0; j < m1.length; j++) {
				sub_m[i][j] = sub_frac(m1[i][j] , m2[i][j]);
			}
		}
		return sub_m;
	}

	public static Fraction[][] matrix_mul_frac(Fraction[][] m1, Fraction[][] m2 ) {
		
		int row = m1.length;
		int col = m2[0].length; 
		int len = m1[0].length;
		
		Fraction[][] temp = new Fraction[row][col];
		
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				for (int k = 0; k < len; k++) {
					if (temp[i][j] == null) {
						temp[i][j] = mul_frac(m1[i][k], m2[k][j]); 
					} else {
						temp[i][j] = add_frac(temp[i][j], mul_frac(m1[i][k], m2[k][j]) );
					}
				} 
			}
		}
		return temp;
	}

	public static void display(Fraction[][] m) {
		int row = m.length;
		int col = m[0].length; 
	
		
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				System.out.print(m[i][j].num + "/" + m[i][j].den);
				System.out.print("    ");
			}
			System.out.print("\n");
		}
			
	}
	
	//--------------------------------------------------------------------------------------
	
	public static Fraction[] mat_frac_to_array (Fraction[][] m ) {
		
		int row = m.length;
		int col = m[0].length;
		
		Fraction[] a = new Fraction[row*col];
		
		int count = 0;
		
		for (int i = 0; i< row; i++) {
			for (int j = 0; j< col; j++) {
				a[count] = m[i][j];
				count ++;
			}
		}
		
		return a;
	}
	

	public static Fraction simplify_frac(Fraction f) {
		double[] f_arr = {f.num, f.den};
		f_arr = simplify(f_arr);
		f.num = f_arr[0];
		f.den = f_arr[1];
		return f;
	}
	
	public static Fraction[] simplify_all_fracs(Fraction[] f) {
		for (int i  = 0; i<f.length; i++) {
			f[i] = simplify_frac(f[i]);
		}
		return f;
	}
	
	
	public static Fraction[] simplify_common_den(Fraction[] t) {
		
		double[] nums = new double[t.length*2];
		
		int count1 = 0;
		
		for (int i = 0; i < t.length; i++) {
			nums[count1] = t[i].num;
			count1 ++;
			nums[count1] = t[i].den;
			count1 ++;
		}
		
		nums = simplify(nums);
		int count2 = 0;
		
		for (int j = 0; j<t.length; j++) {
			 t[j].num = nums[count2];
			 count2 ++;
			 t[j].den = nums[count2];
			 count2 ++;
		}		
		return t;
	}
	
	public static Fraction[][] simplify_matrix_fracs(Fraction[][] m) {
		for (int i  = 0; i<m.length; i++) {
			for (int j  = 0; j<m[0].length; j++) {
				m[i][j] = simplify_frac(m[i][j]);
			}
		}
		return m;
	}
	
	
public static Fraction[] find_common_divider(Fraction[] t) {
		
		
	double[] numerators = new double[t.length];
	double[] denominators = new double[t.length];
		
		for (int i = 0; i<t.length; i++) {
			numerators[i] = t[i].num;
			denominators[i] = t[i].den;
		}
		
		for (int j = 0; j<t.length; j++) {
			for (int k = 0; k<t.length; k++) {
				if (j != k && numerators[k] != 0 && denominators[k] != 0) {
					t[j].num *= denominators[k];
					t[j].den *= denominators[k];
				}
			}
		}
		return t;
	}

	public static int find_min(int[] a) {
		
		
		int m = a[a.length-1];
		for (int i = 0; i < a.length; i++) {
			if (a[i] < m && a[i] != 0) {
				m = a[i]; 
			}
		}
		return m;
	}
	
public static double find_min_abs(double[] a) {
		
		
	double m = Math.abs(a[a.length-1]);
		for (int i = 0; i < a.length; i++) {
			a[i] = Math.abs(a[i]);
			if (a[i] < m && a[i] != 0) {
				m = a[i]; 
			}
		}
		return m;
	}
	
	
	public static double[] div_array(double[] a, double d) {
		
		
		for (int i  = 0; i<a.length;i++) {
			a[i] = a[i] / d;
		}
		return a;
	}
	
	
	public static boolean is_div(double[] a, double d) {
		
		
		boolean divisble = true;
		for (int i  = 0; i<a.length;i++) {
			if (a[i] % d != 0) {
				divisble = false;
				break;
			}
		}
		return divisble;
	}
	
	public static double[] simplify(double[] a) {
		
		
		double min = find_min_abs(a);
		for (double i = 2; i < min+1; i++) {
			if (min % i == 0) {
					while (is_div(a, i) == true) {
						a = div_array(a, i);
					}	
				}
			}
		return a;
	}
}
