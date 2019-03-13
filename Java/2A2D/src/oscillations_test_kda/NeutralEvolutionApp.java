package oscillations_test_kda;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;

import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.display.DrawingFrame;
import org.opensourcephysics.display.PlottingPanel;

public class NeutralEvolutionApp{

	public static int evolutionTimeTotal = 100000;
	public static double[] PVR;
	
	public static double[][] essentiality_N_T;
	public static double[] pvrTest;
	
	// for finding PVR
	private static double[] maxValue(double[][] solution, int steps,int numberIterations,int N) 
	{
		double[] max = new double[N];
		// start the max at 5000th point in the array so the initial condition doesn't get counted as the peak or valley.
		for(int i=0;i<N;i++)
		{
			max[i]=solution[i][5000];
		}
		for(int j=0;j<N;j++)
		{
			for (int i = 5000; i < numberIterations*steps; i++) 
			{
				if (solution[j][i] > max[j]) {
					max[j] = solution[j][i];
				}
			}
		}
		return max;
	}
	
	private static double[] minValue(double[][] solution, int steps,int numberIterations,int N) 
	{
		double[] min = new double[N];
		// start the max at 5000th point in the array so the initial condition doesn't get counted as the peak or valley.
		for(int i=0;i<N;i++)
		{
			min[i]=solution[i][5000];
		}
		for(int j=0;j<N;j++)
		{
			for (int i = 5000; i < numberIterations*steps; i++) 
			{
				if (solution[j][i] < min[j]) {
					min[j] = solution[j][i];
				}
			}
		}
		return min;
	}
	
	// for finding PVR
		private static double[] maxValue(double[][] solution, int steps,int numberIterations,int N, int start, int end) 
		{
			double[] max = new double[N];
			// start the max at 5000th point in the array so the initial condition doesn't get counted as the peak or valley.
			for(int i=0;i<N;i++)
			{
				max[i]=solution[i][start];
			}
			for(int j=0;j<N;j++)
			{
				for (int i = start; i < end; i++) 
				{
					if (solution[j][i] > max[j]) {
						max[j] = solution[j][i];
					}
				}
			}
			return max;
		}
	
		private static double[] minValue(double[][] solution, int steps,int numberIterations,int N, int start, int end) 
		{
			double[] min = new double[N];
			// start the max at 5000th point in the array so the initial condition doesn't get counted as the peak or valley.
			for(int i=0;i<N;i++)
			{
				min[i]=solution[i][start];
			}
			for(int j=0;j<N;j++)
			{
				for (int i = start; i < end; i++) 
				{
					if (solution[j][i] < min[j]) {
						min[j] = solution[j][i];
					}
				}
			}
			return min;
		}	
	
	// mutates a string and returns the index and the mutated value in the form of an array
	public static int[] mutateString(int[] in)
	{
			// select a string between 0 and 24
			Random stringIndex = new Random();
			int[] returnArray = new int[2];
			int index = stringIndex.nextInt(25);
			//System.out.println(" index: "+index+" value: "+in[index]);
			if(in[index]==0)
			{
				in[index] = 1;
				//System.out.println(" Changing value to 1");
			}
			else
			{
				in[index] = 0;
				//System.out.println(" Changing value to 0");
			}
			// to the check the value of the changed array
			//System.out.println(" after index: "+index+" value: "+in[index]);
			returnArray[0] = index;
			returnArray[1] = in[index];
			return returnArray;
			
		// end mutateString Snippet
	}
	
	public static void main(String[] args)
	{
		//NewtonRasphonNoPlot NR = new NewtonRasphonNoPlot(1,1,1);		// this object determines stables points and eigenvalues.
		EvolutionNProteins NE = new EvolutionNProteins();		// this object has the differential equation solver.
		PVR = new double[NE.N];	// will contain peak to valley ratios for all proteins.
		essentiality_N_T = new double[NE.N][evolutionTimeTotal];
		pvrTest = new double[NE.N*2];	// will contain peak to valley ratios for all proteins.
		
		double[][] average_Seq = new double[NE.N*2][NE.M];
		double[][] k_ij_seq_average = new double[NE.N*NE.N][NE.M];
		int acceptedMutations = 0;
		int rejectedMutations = 0;
		
		int xPos = 0; int yPos = 0;
		
		// initialize diffeqs
		NE.initialize();
		// form sequences of In and Out Strings
		//NE.formStrings();
		// duplicates sequences for the extra two proteins.
		
		for(int time=0;time<evolutionTimeTotal;time++)
		{
			
			boolean graphicsEnabled = false;
			boolean moveAccepted=false;
			NE.time=0;
			
			// compute all rate constants
			NE.keff2D[0][0] = NE.computeK(NE.In1, NE.Out1); NE.keff2D[0][1] = NE.computeK(NE.In1, NE.Out2); NE.keff2D[0][2] = NE.computeK(NE.In1, NE.Out3); NE.keff2D[0][3] = NE.computeK(NE.In1, NE.Out4);
			NE.keff2D[1][0] = NE.computeK(NE.In2, NE.Out1); NE.keff2D[1][1] = NE.computeK(NE.In2, NE.Out2); NE.keff2D[1][2] = NE.computeK(NE.In2, NE.Out3); NE.keff2D[1][3] = NE.computeK(NE.In2, NE.Out4);
			NE.keff2D[2][0] = NE.computeK(NE.In3, NE.Out1); NE.keff2D[2][1] = NE.computeK(NE.In3, NE.Out2); NE.keff2D[2][2] = NE.computeK(NE.In3, NE.Out3); NE.keff2D[2][3] = NE.computeK(NE.In3, NE.Out4);
			NE.keff2D[3][0] = NE.computeK(NE.In4, NE.Out1); NE.keff2D[3][1] = NE.computeK(NE.In4, NE.Out2); NE.keff2D[3][2] = NE.computeK(NE.In4, NE.Out3); NE.keff2D[3][3] = NE.computeK(NE.In4, NE.Out4);
	
			// the temporary array in contains the latest values of the In1 string.
			int[] in1 = new int[NE.In1.length];
			int[] in2 = new int[NE.In2.length];
			int[] in3 = new int[NE.In3.length];
			int[] in4 = new int[NE.In4.length];
			int[] out1 = new int[NE.Out1.length];
			int[] out2 = new int[NE.Out2.length];
			int[] out3 = new int[NE.Out3.length];
			int[] out4 = new int[NE.Out4.length];
			
			for(int i=0;i<NE.M;i++)
			{
				in1[i] = NE.In1[i];
				in2[i] = NE.In2[i];
				in3[i] = NE.In3[i];
				in4[i] = NE.In4[i];
				out1[i] = NE.Out1[i];
				out2[i] = NE.Out2[i];
				out3[i] = NE.Out3[i];
				out4[i] = NE.Out4[i];
			}
			
			int sequence=0;
			if(time>0)
			{
				Random randomString = new Random();
				sequence = randomString.nextInt(8);
				if(sequence==0){mutateString(in1);}
				if(sequence==1){mutateString(in2);}
				if(sequence==2){mutateString(in3);}
				if(sequence==3){mutateString(in4);}
				if(sequence==4){mutateString(out1);}
				if(sequence==5){mutateString(out2);}
				if(sequence==6){mutateString(out3);}
				if(sequence==7){mutateString(out4);}
				
			}

			NE.keff2D[0][0] = NE.computeK(in1, out1); NE.keff2D[0][1] = NE.computeK(in1, out2); NE.keff2D[0][2] = NE.computeK(in1, out3); NE.keff2D[0][3] = NE.computeK(in1, out4);
			NE.keff2D[1][0] = NE.computeK(in2, out1); NE.keff2D[1][1] = NE.computeK(in2, out2); NE.keff2D[1][2] = NE.computeK(in2, out3); NE.keff2D[1][3] = NE.computeK(in2, out4);
			NE.keff2D[2][0] = NE.computeK(in3, out1); NE.keff2D[2][1] = NE.computeK(in3, out2); NE.keff2D[2][2] = NE.computeK(in3, out3); NE.keff2D[2][3] = NE.computeK(in3, out4);
			NE.keff2D[3][0] = NE.computeK(in4, out1); NE.keff2D[3][1] = NE.computeK(in4, out2); NE.keff2D[3][2] = NE.computeK(in4, out3); NE.keff2D[3][3] = NE.computeK(in4, out4);		
				 			
			//essentiality_N_T[0][time] = 0 or 1;
			for (int i = 0; i < NE.numberIterations; i++) 
			{
				NE.Essentiality(1, 1, 1, 1);
				NE.RKSolveN();			
			}			
			for(int i=0;i<NE.N;i++)
			{
				PVR[i] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N)[i]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N)[i];
				//System.out.println(maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N)[i]+" "+minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N)[i]);
			}
			
			int start1 = 0;
			int end1 = 5000;
			int start2 = 7500;
			int end2 = NE.numberIterations*NE.steps; // numberIterations*steps
			
			pvrTest[0] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[0]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[0];
			pvrTest[1] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[0]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[0];
			
			pvrTest[2] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[1]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[1];
			pvrTest[3] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[1]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[1];
			
			pvrTest[4] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[2]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[2];
			pvrTest[5] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[2]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[2];

			pvrTest[6] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[3]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start1,end1)[3];
			pvrTest[7] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[3]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N,start2,end2)[3];

			
			
			// if all PVR's > then a threshold, system oscillating, else roll back time.
			// all accepted (oscillator) results
			
			if( ( (PVR[0]-1)<1e-2) && ((PVR[1]-1)<1e-2) && ((PVR[2]-1)<1e-2) && ((PVR[3]-1)<1e-2) )
			{
				time--;				
				rejectedMutations++;
				
			}	// move rejected
			//else if(pvrTest[0]-pvrTest[1]<0.5 & pvrTest[2]-pvrTest[3]< 0.5)
			else if(pvrTest[0]-pvrTest[1]<0.5 & pvrTest[2]-pvrTest[3]< 0.5 & pvrTest[4]-pvrTest[5]<0.5 & pvrTest[6]-pvrTest[7]< 0.5)
			{			
				moveAccepted=true;
				acceptedMutations++;
				//System.out.print(time+" ");
				
				
				/*
				// check Essentiality K1
				EvolutionNProteins NE1 = new EvolutionNProteins();
				NE1.initialize();
				double[] PVR_Ess_K1 = new double[NE1.N];
				for (int i = 0; i < NE1.numberIterations; i++) 
				{					
					NE1.Essentiality(1e-5, 1, 1, 1);
					NE1.keff2D[0][0] = NE.computeK(in1, out1); NE1.keff2D[0][1] = NE.computeK(in1, out2); NE1.keff2D[0][2] = NE.computeK(in1, out3); NE1.keff2D[0][3] = NE.computeK(in1, out4);
					NE1.keff2D[1][0] = NE.computeK(in2, out1); NE1.keff2D[1][1] = NE.computeK(in2, out2); NE1.keff2D[1][2] = NE.computeK(in2, out3); NE1.keff2D[1][3] = NE.computeK(in2, out4);
					NE1.keff2D[2][0] = NE.computeK(in3, out1); NE1.keff2D[2][1] = NE.computeK(in3, out2); NE1.keff2D[2][2] = NE.computeK(in3, out3); NE1.keff2D[2][3] = NE.computeK(in3, out4);
					NE1.keff2D[3][0] = NE.computeK(in4, out1); NE1.keff2D[3][1] = NE.computeK(in4, out2); NE1.keff2D[3][2] = NE.computeK(in4, out3); NE1.keff2D[3][3] = NE.computeK(in4, out4);			
					NE1.RKSolveN();	
				}
				for(int i=0;i<NE1.N;i++)
				{
					PVR_Ess_K1[i] = maxValue(NE1.y_solution,NE1.steps,NE1.numberIterations,NE1.N)[i]/minValue(NE1.y_solution,NE1.steps,NE1.numberIterations,NE1.N)[i];
					// if stops oscillating, essentiality = 1
				}
				
				if( (PVR_Ess_K1[0]-1)<1e-2 & (PVR_Ess_K1[1]-1)<1e-2 & (PVR_Ess_K1[2]-1)<1e-2 & (PVR_Ess_K1[3]-1)<1e-2 )
				{
					essentiality_N_T[0][time] = 1;
					System.out.print(" 1 ");
				}
				else
				{
					essentiality_N_T[0][time] = 0;
					System.out.print(" 0 ");
				}	// end essentiality K1
				
				// check Essentiality K2
				EvolutionNProteins NE2 = new EvolutionNProteins();
				NE2.initialize();
				double[] PVR_Ess_K2 = new double[NE2.N];
				for (int i = 0; i < NE2.numberIterations; i++) 
				{					
					NE2.Essentiality(1, 1e-5, 1, 1);
					NE2.keff2D[0][0] = NE.computeK(in1, out1); NE2.keff2D[0][1] = NE.computeK(in1, out2); NE2.keff2D[0][2] = NE.computeK(in1, out3); NE2.keff2D[0][3] = NE.computeK(in1, out4);
					NE2.keff2D[1][0] = NE.computeK(in2, out1); NE2.keff2D[1][1] = NE.computeK(in2, out2); NE2.keff2D[1][2] = NE.computeK(in2, out3); NE2.keff2D[1][3] = NE.computeK(in2, out4);
					NE2.keff2D[2][0] = NE.computeK(in3, out1); NE2.keff2D[2][1] = NE.computeK(in3, out2); NE2.keff2D[2][2] = NE.computeK(in3, out3); NE2.keff2D[2][3] = NE.computeK(in3, out4);
					NE2.keff2D[3][0] = NE.computeK(in4, out1); NE2.keff2D[3][1] = NE.computeK(in4, out2); NE2.keff2D[3][2] = NE.computeK(in4, out3); NE2.keff2D[3][3] = NE.computeK(in4, out4);			
					NE2.RKSolveN();	
				}
				for(int i=0;i<NE2.N;i++)
				{
					PVR_Ess_K2[i] = maxValue(NE2.y_solution,NE2.steps,NE2.numberIterations,NE2.N)[i]/minValue(NE2.y_solution,NE2.steps,NE2.numberIterations,NE2.N)[i];
					// if stops oscillating, essentiality = 1
				}
				
				if( (PVR_Ess_K2[0]-1)<1e-2 & (PVR_Ess_K2[1]-1)<1e-2 & (PVR_Ess_K2[2]-1)<1e-2 & (PVR_Ess_K2[3]-1)<1e-2 )
				{
					essentiality_N_T[1][time] = 1;
					System.out.print(" 1 ");
				}
				else
				{
					essentiality_N_T[1][time] = 0;
					System.out.print(" 0 ");
				}	// end essentiality K2
				
				// check Essentiality P1
				EvolutionNProteins NE3 = new EvolutionNProteins();
				NE3.initialize();
				double[] PVR_Ess_P1 = new double[NE3.N];
				for (int i = 0; i < NE3.numberIterations; i++) 
				{					
					NE3.Essentiality(1, 1, 1e-5, 1);
					NE3.keff2D[0][0] = NE.computeK(in1, out1); NE3.keff2D[0][1] = NE.computeK(in1, out2); NE3.keff2D[0][2] = NE.computeK(in1, out3); NE3.keff2D[0][3] = NE.computeK(in1, out4);
					NE3.keff2D[1][0] = NE.computeK(in2, out1); NE3.keff2D[1][1] = NE.computeK(in2, out2); NE3.keff2D[1][2] = NE.computeK(in2, out3); NE3.keff2D[1][3] = NE.computeK(in2, out4);
					NE3.keff2D[2][0] = NE.computeK(in3, out1); NE3.keff2D[2][1] = NE.computeK(in3, out2); NE3.keff2D[2][2] = NE.computeK(in3, out3); NE3.keff2D[2][3] = NE.computeK(in3, out4);
					NE3.keff2D[3][0] = NE.computeK(in4, out1); NE3.keff2D[3][1] = NE.computeK(in4, out2); NE3.keff2D[3][2] = NE.computeK(in4, out3); NE3.keff2D[3][3] = NE.computeK(in4, out4);			
					NE3.RKSolveN();	
				}
				for(int i=0;i<NE3.N;i++)
				{
					PVR_Ess_P1[i] = maxValue(NE3.y_solution,NE3.steps,NE3.numberIterations,NE3.N)[i]/minValue(NE3.y_solution,NE3.steps,NE3.numberIterations,NE3.N)[i];
					// if stops oscillating, essentiality = 1
				}
				
				if( (PVR_Ess_P1[0]-1)<1e-2 & (PVR_Ess_P1[1]-1)<1e-2 & (PVR_Ess_P1[2]-1)<1e-2 & (PVR_Ess_P1[3]-1)<1e-2 )
				{
					essentiality_N_T[2][time] = 1;
					System.out.print(" 1 ");
				}
				else
				{
					essentiality_N_T[2][time] = 0;
					System.out.print(" 0 ");
				}	// end essentiality P1
				// check Essentiality P2
				EvolutionNProteins NE4 = new EvolutionNProteins();
				NE4.initialize();
				double[] PVR_Ess_P2 = new double[NE4.N];
				for (int i = 0; i < NE4.numberIterations; i++) 
				{					
					NE4.Essentiality(1, 1, 1, 1e-5);
					NE4.keff2D[0][0] = NE.computeK(in1, out1); NE4.keff2D[0][1] = NE.computeK(in1, out2); NE4.keff2D[0][2] = NE.computeK(in1, out3); NE4.keff2D[0][3] = NE.computeK(in1, out4);
					NE4.keff2D[1][0] = NE.computeK(in2, out1); NE4.keff2D[1][1] = NE.computeK(in2, out2); NE4.keff2D[1][2] = NE.computeK(in2, out3); NE4.keff2D[1][3] = NE.computeK(in2, out4);
					NE4.keff2D[2][0] = NE.computeK(in3, out1); NE4.keff2D[2][1] = NE.computeK(in3, out2); NE4.keff2D[2][2] = NE.computeK(in3, out3); NE4.keff2D[2][3] = NE.computeK(in3, out4);
					NE4.keff2D[3][0] = NE.computeK(in4, out1); NE4.keff2D[3][1] = NE.computeK(in4, out2); NE4.keff2D[3][2] = NE.computeK(in4, out3); NE4.keff2D[3][3] = NE.computeK(in4, out4);			
					NE4.RKSolveN();	
				}
				for(int i=0;i<NE4.N;i++)
				{
					PVR_Ess_P2[i] = maxValue(NE4.y_solution,NE4.steps,NE4.numberIterations,NE4.N)[i]/minValue(NE4.y_solution,NE4.steps,NE4.numberIterations,NE4.N)[i];
					// if stops oscillating, essentiality = 1
				}
				
				if( (PVR_Ess_P2[0]-1)<1e-2 & (PVR_Ess_P2[1]-1)<1e-2 & (PVR_Ess_P2[2]-1)<1e-2 & (PVR_Ess_P2[3]-1)<1e-2 )
				{
					essentiality_N_T[3][time] = 1;
					System.out.print(" 1 ");
				}
				else
				{
					essentiality_N_T[3][time] = 0;
					System.out.print(" 0 ");
				}	// end essentiality P2
				*/
				
				/*
				// for drawing figure 1 in paper 2s
				if(acceptedMutations%1000==0)
				{
					System.out.println(acceptedMutations);
					
					for(int j=0;j<NE.N;j++)
					{				
						for(int i=0;i<NE.N;i++)
						{
							//System.out.print(" "+NE.keff2D[j][i]+" ");
							System.out.println(NE.keff2D[j][i]);
						}					
					}
					
					
					double scaleFactor = 30.0;
					System.out.println("\\draw[edge,loop,line width ="+ NE.keff2D[0][0]/scaleFactor+"] (1) to (1);");
					System.out.println("\\draw[edge, line width="+ NE.keff2D[0][1]/scaleFactor+"] (2) to (1);");
					System.out.println("\\draw[-|,shorten >=0.07 cm, line width = "+ NE.keff2D[0][2]/scaleFactor+"] (3.80) to (1.280);");
					System.out.println("\\draw[-|, shorten >=0.05cm, line width="+ NE.keff2D[0][3]/scaleFactor+"] (4.145) to (1.305);");
					System.out.println("\\draw[edge, line width = "+ NE.keff2D[1][0]/scaleFactor+"] (1.20) to (2.160);");
					System.out.println("\\draw[edge,loop, line width = "+ NE.keff2D[1][1]/scaleFactor+"] (2) to (2);");
					System.out.println("\\draw[-|, shorten >=0.05cm,line width = "+ NE.keff2D[1][2]/scaleFactor+"] (3.35) to (2.235);");
					System.out.println("\\draw[-|,shorten >=0.07 cm, ,line width = "+ NE.keff2D[1][3]/scaleFactor+"] (4.80) to (2.280);");
					System.out.println("\\draw[edge, line width="+ NE.keff2D[2][0]/scaleFactor+"] (1.260) to (3.100);");
					System.out.println("\\draw[edge,shorten <=0.1cm, line width = "+ NE.keff2D[2][1]/scaleFactor+"] (2.215) to (3.55);");
					System.out.println("\\draw[-|,out=315, in=225, loop,shorten >=0.05 cm, line width ="+ NE.keff2D[2][2]/scaleFactor+"] (3) to (3);");
					System.out.println("\\draw[-|,shorten >=0.07 cm, line width = "+ NE.keff2D[2][3]/scaleFactor+"] (4) to (3);");
					System.out.println("\\draw[edge,shorten <=0.1cm, line width="+ NE.keff2D[3][0]/scaleFactor+"] (1.325) to (4.125);");
					System.out.println("\\draw[edge, line width ="+ NE.keff2D[3][1]/scaleFactor+"] (2.260) to (4.100);");
					System.out.println("\\draw[-|,shorten >=0.07 cm, line width = "+ NE.keff2D[3][2]/scaleFactor+"] (3.20) to (4.160);");
					System.out.println("\\draw[-|,out=315, in=225, loop,shorten >=0.05 cm, line width ="+ NE.keff2D[3][3]/scaleFactor+"] (4) to (4);");
					
				if(acceptedMutations%20000==0)	
					System.exit(0);
					
				}
				*/
				
				for(int i=0;i<NE.M;i++)
				{
					 NE.In1[i] = in1[i];
					 NE.In2[i] = in2[i];
					 NE.In3[i] = in3[i];
					 NE.In4[i] = in4[i];
					 NE.Out1[i] = out1[i];
					 NE.Out2[i] = out2[i];
					 NE.Out3[i] = out3[i];
					 NE.Out4[i] = out4[i];
				}
				
				if(
						!(Math.abs(NE.keff2D[0][2]-NE.keff2D[0][3])<20 & Math.abs(NE.keff2D[0][2]-NE.keff2D[0][3])>1 & Math.abs(NE.keff2D[1][2]-NE.keff2D[1][3])<20 & Math.abs(NE.keff2D[1][2]-NE.keff2D[1][3])>1)
				  )
						 
				{
					
				//List b = Arrays.asList(ArrayUtils.toObject(NE.y_solution[0]));

				for(int j=0;j<NE.N;j++)
				{
					double max = NE.y_solution[j][0];
	
					for (int i = 0; i < NE.y_solution[0].length; i++) {
						if (NE.y_solution[j][i] > max) 
						{
							max = NE.y_solution[j][i];
						}
						
					}	
					System.out.print(max+" ");
				}
				System.out.println();
				rejectedMutations=0;
				//System.out.println(NE.keff2D[0][2]-NE.keff2D[0][3]);
				//System.out.print(time+" ");
				// print k_ij
				/*
				for(int j=0;j<NE.N;j++)
				{				
					for(int i=0;i<NE.N;i++)
					{
						System.out.println("k"+j+""+i+"= "+NE.keff2D[j][i]+"; ");
						//System.out.print(" "+NE.keff2D[j][i]+" ");
					}					
				}
				*/
				//System.out.println();
				}
				//System.out.println(PVR[0]+" "+PVR[1]+" "+PVR[2]+" "+PVR[3]);
				
				/*
				for(int i=0;i<NE.M;i++)
				{
					average_Seq[0][i]+=in1[i];
					average_Seq[1][i]+=out1[i];
					average_Seq[2][i]+=in2[i];
					average_Seq[3][i]+=out2[i];
					average_Seq[4][i]+=in3[i];
					average_Seq[5][i]+=out3[i];
					average_Seq[6][i]+=in4[i];
					average_Seq[7][i]+=out4[i];
				}
				*/
				
			}	// move accepted
			
			// reject move
			else
			{
					time--;
					rejectedMutations++;
			}
				
			
			if(graphicsEnabled && moveAccepted)
			{
				String evolutionTime = Integer.toString(time);
				// Plot properties
				PlottingPanel plot = new PlottingPanel("time", "", evolutionTime);
				DrawingFrame frame = new DrawingFrame(plot); // create frame
				Dataset dataset = new Dataset(); // create dataset
				plot.addDrawable(dataset); // add dataset to plot
			    //dataset.setConnected(true);
				dataset.setMarkerSize(1);
				//frame.setLocation(xPos, yPos);
				frame.setBounds(xPos, yPos, 250,250);
				xPos+=250;
				if(xPos%1750==0){xPos=0;yPos+=250;}
				
				// end plot properties
				
				for (int i = 0; i < NE.numberIterations*NE.steps; i++) {
					for(int j=0;j<NE.N;j++)
					{
						dataset.append(i/100.0, NE.y_solution[j][i]);
					}
				}
				
				frame.setVisible(true); // show the frame
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				
			}	// end graphics
			
		}	// for time
		
		
		/*
		// in 1
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[0][i]/(1.0*acceptedMutations)+", ");}System.out.println();
		// out 1
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[1][i]/(1.0*acceptedMutations)+", ");}System.out.println();
		// in 2
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[2][i]/(1.0*acceptedMutations)+", ");}System.out.println();
		// out 2
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[3][i]/(1.0*acceptedMutations)+", ");}System.out.println();
		// in 3
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[4][i]/(1.0*acceptedMutations)+", ");}System.out.println();
		// out 3
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[5][i]/(1.0*acceptedMutations)+", ");}System.out.println();
		// in 4
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[6][i]/(1.0*acceptedMutations)+", ");}System.out.println();
		// out 4
		for(int i=0;i<NE.M;i++) {System.out.print(average_Seq[7][i]/(1.0*acceptedMutations)+", ");}System.out.println();		
		
		// get raw sequences heatmap
		heatMap hm_2A2D = new heatMap(average_Seq,"raw_seq_2A2D");
		*/
		
		/*
		NE4.keff2D[0][0] = NE.computeK(in1, out1); NE4.keff2D[0][1] = NE.computeK(in1, out2); NE4.keff2D[0][2] = NE.computeK(in1, out3); NE4.keff2D[0][3] = NE.computeK(in1, out4);
		NE4.keff2D[1][0] = NE.computeK(in2, out1); NE4.keff2D[1][1] = NE.computeK(in2, out2); NE4.keff2D[1][2] = NE.computeK(in2, out3); NE4.keff2D[1][3] = NE.computeK(in2, out4);
		NE4.keff2D[2][0] = NE.computeK(in3, out1); NE4.keff2D[2][1] = NE.computeK(in3, out2); NE4.keff2D[2][2] = NE.computeK(in3, out3); NE4.keff2D[2][3] = NE.computeK(in3, out4);
		NE4.keff2D[3][0] = NE.computeK(in4, out1); NE4.keff2D[3][1] = NE.computeK(in4, out2); NE4.keff2D[3][2] = NE.computeK(in4, out3); NE4.keff2D[3][3] = NE.computeK(in4, out4);		
		
		System.out.println();
		
		// k_A1->A1
		k_ij_seq_average[0] = dotProduct(average_Seq[0],average_Seq[1],acceptedMutations);
		// k_A2->A1
		k_ij_seq_average[1] = dotProduct(average_Seq[0],average_Seq[3],acceptedMutations);
		// k_D1->A1
		k_ij_seq_average[2] = dotProduct(average_Seq[0],average_Seq[5],acceptedMutations);
		// k_D1->A1
		k_ij_seq_average[3] = dotProduct(average_Seq[0],average_Seq[7],acceptedMutations);
		// K_A1->A2
		k_ij_seq_average[4] = dotProduct(average_Seq[2],average_Seq[1],acceptedMutations);
		// K_A2->A2
		k_ij_seq_average[5] = dotProduct(average_Seq[2],average_Seq[3],acceptedMutations);
		// K_D1->A2
		k_ij_seq_average[6] = dotProduct(average_Seq[2],average_Seq[5],acceptedMutations);
		// K_D2->A2
		k_ij_seq_average[7] = dotProduct(average_Seq[2],average_Seq[7],acceptedMutations);
		// K_A1->D1
		k_ij_seq_average[8] = dotProduct(average_Seq[4],average_Seq[1],acceptedMutations);
		// K_A2->D1
		k_ij_seq_average[9] = dotProduct(average_Seq[4],average_Seq[3],acceptedMutations);
		// K_D1->D1
		k_ij_seq_average[10] = dotProduct(average_Seq[4],average_Seq[5],acceptedMutations);
		// K_D2->D1
		k_ij_seq_average[11] = dotProduct(average_Seq[4],average_Seq[7],acceptedMutations);
		// K_A1->D2
		k_ij_seq_average[12] = dotProduct(average_Seq[6],average_Seq[1],acceptedMutations);
		// K_A2->D2
		k_ij_seq_average[13] = dotProduct(average_Seq[6],average_Seq[3],acceptedMutations);
		// K_D1->D2
		k_ij_seq_average[14] = dotProduct(average_Seq[6],average_Seq[5],acceptedMutations);
		// K_D2->D2
		k_ij_seq_average[15] = dotProduct(average_Seq[6],average_Seq[7],acceptedMutations);
		
		//for(int i=0;i<NE.M;i++) {System.out.print(k_ij_seq_average[0][i]+" ");}System.out.println();
		heatMap kij_hm_2A2D = new heatMap(k_ij_seq_average,"kij_seq_2A2D");
		System.out.println("Finished");
		*/
		
		/*
		for(int j=0;j<NE.N;j++)
		{				
			for(int i=0;i<NE.N;i++)
			{
				System.out.println(" "+NE.keff2D[j][i]+" ");
			}					
		}
		*/
		
	}	// main
	
	public static double[] dotProduct(double[] In, double[] Out,int acceptedMutations)
	{
		double[] result = new double[In.length];
		for(int i=0;i<In.length;i++)
		{
			result[i]=(In[i]/(1.0*acceptedMutations))*(Out[i]/(1.0*acceptedMutations));
		}
		return result;		
		
	}
	
}	// class


	
//uncomment for rate constants
/*
for(int j=0;j<NE.N;j++)
{				
	for(int i=0;i<NE.N;i++)
	{
		System.out.print(" "+NE.keff2D[j][i]+" ");
	}					
}
*/
