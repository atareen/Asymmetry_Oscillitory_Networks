package Neutral_Evolution_Correct_labels;
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
		
		int[] in_A_average = new int[NE.M];
		int[] out_A_average = new int[NE.M];
		int[] in_D_average = new int[NE.M];
		int[] out_D_average = new int[NE.M];
		
		double fitness = 0.0;

		
		int xPos = 0; int yPos = 0;
		
		// initialize diffeqs
		NE.initialize();
		// form sequences of In and Out Strings
		//NE.formStrings();
		// duplicates sequences for the extra two proteins.
		//NE.duplicateStrings();
		int acceptedMutations = 0;
		
		for(int time=0;time<evolutionTimeTotal;time++)
		//while(1>0)
		{
			
			boolean graphicsEnabled = false;
			boolean moveAccepted=false;
			NE.time=0;
			
			// compute all rate constants
			NE.keff2D[0][0] = NE.computeK(NE.In_A, NE.Out_A); NE.keff2D[0][1] = NE.computeK(NE.In_A, NE.Out_D);
			NE.keff2D[1][0] = NE.computeK(NE.In_D, NE.Out_A); NE.keff2D[1][1] = NE.computeK(NE.In_D, NE.Out_D);


	
			// the temporary array in contains the latest values of the In1 string.
			int[] inA = new int[NE.In_A.length];
			int[] inD = new int[NE.In_D.length];
			int[] outA = new int[NE.Out_A.length];
			int[] outD = new int[NE.Out_D.length];
			
			for(int i=0;i<NE.M;i++)
			{
				inA[i] = NE.In_A[i];
				inD[i] = NE.In_D[i];
				outA[i] = NE.Out_A[i];
				outD[i] = NE.Out_D[i];
			}
			
			int sequence=0;
			if(time>0)
			{
				Random randomString = new Random();
				sequence = randomString.nextInt(4);
				if(sequence==0){mutateString(inA);}
				if(sequence==1){mutateString(inD);}
				if(sequence==2){mutateString(outA);}
				if(sequence==3){mutateString(outD);}
			}
			
			
			// k_A_A
			NE.keff2D[0][0] = NE.computeK(inA, outA);
			// k_D_A
			NE.keff2D[0][1] = NE.computeK(inA, outD);
			// k_A_D
			NE.keff2D[1][0] = NE.computeK(inD, outA);
			// k_D_D
			NE.keff2D[1][1] = NE.computeK(inD, outD);
					
				 			
			//essentiality_N_T[0][time] = 0 or 1;
			for (int i = 0; i < NE.numberIterations; i++) 
			{
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
			
			// if all PVR's > then a threshold, system oscillating, else roll back time.
			// all accepted (oscillator) results
			
			if( ( (PVR[0]-1)<1e-2) && ((PVR[1]-1)<1e-2))
			{
				time--;				
				moveAccepted = false;
				//System.out.println(" Rejecting Move!!! ");
				/*
					String seqA = "";
					String seqD = "";
					
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(inA[i]+"");
						seqA+=Integer.toString(inA[i]);						
					}//System.out.println();
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(outA[i]+"");
						seqA+=Integer.toString(outA[i]);
					}//System.out.println();
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(inD[i]+"");
						seqD+=Integer.toString(inD[i]);
					}//System.out.println();
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(outD[i]+"");
						seqD+=Integer.toString(outD[i]);
					}//System.out.println();
					//System.out.println(seq);
					
					Long decimalValueA = Long.parseLong(seqA, 2);
					Long decimalValueD = Long.parseLong(seqA, 2);
					System.out.println(decimalValueA+" "+decimalValueD);
				*/
				
			}	// move rejected
			
			// else non decaying oscillator
			//else if(pvrTest[0]-pvrTest[1]<0.5 & pvrTest[2]-pvrTest[3]< 0.5)
			else if(pvrTest[0]-pvrTest[1]<0.5 & pvrTest[2]-pvrTest[3]< 0.5)
			{
				// relative PVR
					//System.out.println(" Move Accepted! ");
					moveAccepted=true;
					
					for(int i=0;i<NE.M;i++)
					{
						 NE.In_A[i] = inA[i];
						 NE.In_D[i] = inD[i];
						 NE.Out_A[i] = outA[i];
						 NE.Out_D[i] = outD[i];
					}
					//System.out.print(acceptedMutations+" ");
					acceptedMutations++;
					
					if(acceptedMutations%100==0)
						//System.out.println(acceptedMutations);
					
					for(int i=0;i<NE.M;i++) 
					{
						in_A_average[i] +=inA[i];
						in_D_average[i] +=inD[i];
						out_A_average[i]+=outA[i];
						out_D_average[i]+=outD[i];
						//System.out.print(in1[i]);
					}
					
					
					String seqA = "";
					String seqD = "";
					
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(inA[i]+"");
						seqA+=Integer.toString(inA[i]);						
					}//System.out.println();
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(outA[i]+"");
						seqA+=Integer.toString(outA[i]);
					}//System.out.println();
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(inD[i]+"");
						seqD+=Integer.toString(inD[i]);
					}//System.out.println();
					for(int i=0;i<NE.M;i++) 
					{
						//System.out.print(outD[i]+"");
						seqD+=Integer.toString(outD[i]);
					}//System.out.println();
					//System.out.println(seq);
					
					Long decimalValueA = Long.parseLong(seqA, 2);
					Long decimalValueD = Long.parseLong(seqA, 2);
					System.out.println(decimalValueA+" "+decimalValueD);
					
					
					
					/*
					for(int j=0;j<NE.N;j++)
					{				
						for(int i=0;i<NE.N;i++)
						{
							System.out.print(" "+NE.keff2D[j][i]+" ");
						}					
					}
					*/
					//System.out.println();
					//System.out.println(PVR[0]+" "+PVR[1]+" "+PVR[2]+" "+PVR[3]);

			}	// move accepted: oscillator
			
			
			// reject move
			else
			{
					//System.out.println(" Rejecting Move! ");		
					moveAccepted = false;
					time--;
			}

			
			if(graphicsEnabled && moveAccepted)
			{
				String evolutionTime = Integer.toString(acceptedMutations);
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
		
		for(int i=0;i<NE.M;i++) {System.out.print(in_A_average[i]/(1.0*acceptedMutations)+", ");}System.out.println();
		for(int i=0;i<NE.M;i++) {System.out.print(out_A_average[i]/(1.0*acceptedMutations)+", ");}System.out.println();
		for(int i=0;i<NE.M;i++) {System.out.print(in_D_average[i]/(1.0*acceptedMutations)+", ");}System.out.println();
		for(int i=0;i<NE.M;i++) {System.out.print(out_D_average[i]/(1.0*acceptedMutations)+", ");}System.out.println();
		
	}	// main
}
	
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
