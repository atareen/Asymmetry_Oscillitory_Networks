package Neutral_Evolution_Correct_labels;



// The purpose of this class is to solve N coupled differential equations using the RK4 method.
// Start Date: 09/22/2014
// Updated: //11/25/2015

// This program can currently solve N=2 to N = 4. It can also do arbitrary N but Keff2D constants will have to be defined.

// Currently does two kinase and two phosphatases.

public class EvolutionNProteins {

	
	public  int N_kinase = 1;
	public  int N_phosphatase = 1;
	public  int N = N_kinase + N_phosphatase; // number of differential equations
	public  int N_total = N_kinase + N_phosphatase;
  //public  double[] keff = {720,5,720,25,0.1,10,75,1,30};
	public  double[] keff = {720,5,720,0.1,25,1,10,0.1,75,1,30,0.1,5,2,3,2};
	public  double[][] keff2D = new double[N][N];
	public  double[] y = new double[N]; // y: (y1,y2,...,yn)
	public  double[] diffEq = new double[N];
	public  double deltaX = 0.001;
	public  int steps = 100;	// don't change
	public  double x = 0;
	public  int time = 0;
	public  double[] bgCoeffs = new double [2*N]; // background coefficients e.g. alpha, beta...
	public  int numberIterations = 100;
	public  double[] kStar_Solution = new double[numberIterations*steps];
	public  double[] pStar_Solution = new double[numberIterations*steps];
	public  double[][] y_solution = new double[N][numberIterations*steps];
	public  int M = 25;	// total string length
	public  double k0 = 100;
	public  double h = 2;
	public  double epsilon = -0.2;
	/*
	public  int[] In1 = new int[M];
	public  int[] Out1 = new int[M];
	public  int[] In2 = new int[M];
	public  int[] Out2 = new int[M];
	*/
	
	public int[]  In_A = {0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0};
	public int[] Out_A = {1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1};
	public int[] In_D = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1};
	public int[] Out_D = {0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0};
	
	/*
	public  void formStrings()
	{
		for(int i=0;i<M;i++)
		{
			if(i<25)
			{
				In1[i]  =  1; //random.nextInt(2);
			}
			if(i>11 & i<25)
			{
				In2[i] = 1;   //random.nextInt(2);
			}
			if(i>2 & i<23)
			{
				Out2[i] = 1;   //random.nextInt(2);
			}
			if(i>4 & i<25)
			{
				Out1[i] =  1; //random.nextInt(2);
			}
		}
	}
	*/
	// the idea here is to keep the strings for all four proteins the same such the system is essentially a two protein system before a given time.
	// 0-1 are K's and 2-3 are P's
	public void duplicateStrings()
	{
			for(int i=0;i<M;i++)
			{
				In_D[i]  = In_A[i];
				Out_D[i] = Out_A[i];
			}
	}
	
	public  double computeK(int[] In, int[] Out)
	{
		double I = 0;
		double k = 0;
		for(int i=0;i<M;i++)
		{
			I+=In[i]*Out[i];
		}
		double E0 = M*epsilon;
		double E = I*epsilon;
		k = Math.pow(k0/(1+Math.exp(E-E0)),h);
		return k;
		
	}
	
	
	public  void initialize()
	{
		// Initialize Variables
		
		for(int i=0;i<N;i++)
		{
			//Random random1 = new Random();
			y[i] = 0.5;//random1.nextDouble();
		}
/*
		//Initialize rate constants
		if(N==2)	// 1 K 1 P
		{
			//public  double[] keff = {720,5,720,25,0.1,10,75,1,30};
			keff2D[0][0] = 723.29; keff2D[0][1] = 723.29;
			keff2D[1][0] = 47.8;  keff2D[1][1] = 32.86;
		}

		if(N==3)	// 2 K 1 P
		{
			//public  double[] keff = {720,5,720,25,0.1,10,75,1,30};
			keff2D[0][0] = 720; keff2D[0][1] = 5;   keff2D[0][2] = 720;
			keff2D[1][0] = 25;  keff2D[1][1] = 0.1; keff2D[1][2] = 10; 
			keff2D[2][0] = 75;  keff2D[2][1] = 1;   keff2D[2][2] = 30;
		}
		
		if(N==4)			// 2 K 2 P
		{
			//public  double[] keff = {720,5,720,0.1,25,1,10,0.1,75,1,30,0.1,5,2,3,2};
			keff2D[0][0] = 720; keff2D[0][1] = 5; keff2D[0][2] = 720; keff2D[0][3] = 0.1;
			keff2D[1][0] = 25;  keff2D[1][1] = 1; keff2D[1][2] = 10;  keff2D[1][3] = 0.1;
			keff2D[2][0] = 75;  keff2D[2][1] = 1; keff2D[2][2] = 30;  keff2D[2][3] = 0.1;
			keff2D[3][0] = 5;   keff2D[3][1] = 2; keff2D[3][2] = 3;   keff2D[3][3] = 2;
		}

		if(N>4)
		{
			Random random = new Random();
			for(int j=0;j<N;j++)
			{
				for(int i=0;i<N;i++)
				{
					keff2D[j][i] = random.nextDouble();
				}
			}
		}
*/		
		// Initialize background coefficients
		for(int i=0;i<bgCoeffs.length;i++)
		{
			bgCoeffs[i] = 1;
		}
	}
	
	public  double[] update = new double[N];
	
	public  double dydx(double x, int i) 
	{
		//dydxUpdate(update);
		dydxUpdateAuto(update);
		return diffEq[i];
	}
	

	// define and update differential equations here (Automated)
	public  void dydxUpdateAuto(double[] y) 
	{	
		for(int i=0;i<diffEq.length;i++)
		{
			diffEq[i] = 0;
		}
		/*
		for(int j=0;j<N_total;j++)
		{
			// sum over kinase
			for(int i=0;i<N_kinase;i++)
			{		// might have to use mods here.
					diffEq[j] += keff2D[j][i]*(y[i]*y[i])*(1-y[j]);  	
			}
			// sum over phosphatase
			for(int i=N_kinase;i<N_total;i++)
			{		// might have to use mods here.
					diffEq[j] += -(keff2D[j][i]*(y[i]*y[i])*(y[j]));  
			}
			// background terms
			diffEq[j] += (bgCoeffs[0]*(1-y[j])-bgCoeffs[1]*(y[j])); 
		}
		*/
		diffEq[0] = keff2D[0][0]*(y[0]*y[0])*(1-y[0])-keff2D[0][1]*(y[1]*y[1])*(y[0])+bgCoeffs[0]*(1-y[0])-1*y[0];
		diffEq[1] = keff2D[1][0]*(y[0]*y[0])*(1-y[1])-keff2D[1][1]*(y[1]*y[1])*(y[1])+bgCoeffs[0]*(1-y[1])-bgCoeffs[1]*y[1];
	}

	public  void RKSolveN() {

		double x;

		double[] k1 = new double[N];
		double[] k2 = new double[N];
		double[] k3 = new double[N];
		double[] k4 = new double[N];
		
		for(int i=0;i<y.length;i++)
		{
			update[i] = y[i];
		}
		
		for (int i = 0; i < steps; i++) {
			x = i * deltaX;

			for (int j = 0; j < N; j++) {
				k1[j] = deltaX*dydx(x,j);
			}
			
			for (int j = 0; j < N; j++) {
				update[j] = y[j]+k1[j]/2.0;
			}
			
			for (int j = 0; j < N; j++) {
				k2[j] = deltaX*dydx(x+deltaX/2.0,j);
			}
			
			for (int j = 0; j < N; j++) {
				update[j] = y[j]+k2[j]/2.0;
			}
			
			for (int j = 0; j < N; j++) {
				k3[j] = deltaX*dydx(x+deltaX/2.0,j);
			}
			
			for (int j = 0; j < N; j++) {
				update[j] = y[j]+k3[j];
			}
			
			for (int j = 0; j < N; j++) {
				k4[j] = deltaX*dydx(x+deltaX,j);
			}
					
			for(int j=0;j<N;j++)
			{
				y[j] += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6.0;
			}
			
			for(int j=0;j<N;j++)
			{
				y_solution[j][time] = y[j];
			}
			/*
				System.out.print(time+" ");
				for(int j=0;j<N;j++)
				{	
					System.out.print(" "+y[j]);
				}
				System.out.println();
			*/
			time++;
		}
	}
}
