package GreedyEvolution;
import java.util.Random;


import javax.swing.JFrame;
import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.display.DrawingFrame;
import org.opensourcephysics.display.PlottingPanel;

public class NeutralEvolutionApp{

	public static int evolutionTimeTotal = 10000000;
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
		
		double fitness = 0.0;

		
		int xPos = 0; int yPos = 0;
		
		// initialize diffeqs
		NE.initialize();
		// form sequences of In and Out Strings
		//NE.formStrings();
		// duplicates sequences for the extra two proteins.
		//NE.duplicateStrings();
		int acceptedMutations = 0;
		double beta = 0.005;
		
		for(int time=0;time<evolutionTimeTotal;time++)
		//while(1>0)
		{
			
			boolean graphicsEnabled = false;
			boolean moveAccepted=false;
			NE.time=0;
			double newFitness = 0;
			
			// compute all rate constants
			NE.keff2D[0][0] = NE.computeK(NE.In1, NE.Out1); NE.keff2D[0][1] = NE.computeK(NE.In1, NE.Out2);
			NE.keff2D[1][0] = NE.computeK(NE.In2, NE.Out1); NE.keff2D[1][1] = NE.computeK(NE.In2, NE.Out2);


	
			// the temporary array in contains the latest values of the In1 string.
			int[] in1 = new int[NE.In1.length];
			int[] in2 = new int[NE.In2.length];
			int[] out1 = new int[NE.Out1.length];
			int[] out2 = new int[NE.Out2.length];
			
			for(int i=0;i<NE.M;i++)
			{
				in1[i] = NE.In1[i];
				in2[i] = NE.In2[i];
				out1[i] = NE.Out1[i];
				out2[i] = NE.Out2[i];
			}
			
			int sequence=0;
			if(time>0)
			{
				Random randomString = new Random();
				sequence = randomString.nextInt(4);
				if(sequence==0){mutateString(in1);}
				if(sequence==1){mutateString(in2);}
				if(sequence==2){mutateString(out1);}
				if(sequence==3){mutateString(out2);}
			}

			NE.keff2D[0][0] = NE.computeK(in1, out1); NE.keff2D[0][1] = NE.computeK(in1, out2);
			NE.keff2D[1][0] = NE.computeK(in2, out1); NE.keff2D[1][1] = NE.computeK(in2, out2);
					
				 			
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
				//time--;				
				moveAccepted = false;
				//System.out.println(" Rejecting Move!!! ");
			}	// move rejected
			
			// else non decaying oscillator
			//else if(pvrTest[0]-pvrTest[1]<0.5 & pvrTest[2]-pvrTest[3]< 0.5)
			else if(pvrTest[0]-pvrTest[1]<0.5 & pvrTest[2]-pvrTest[3]< 0.5)
			{
				// relative PVR
				newFitness = (Math.abs(PVR[0]-PVR[1]))/(PVR[0]+PVR[1]);
				Random mcRandGen = new Random();
				double mcRandomNumber = mcRandGen.nextDouble();
				if(newFitness>fitness)
				{
					acceptedMutations++;
					//System.out.println(" Move Accepted! ");
					fitness = newFitness;
					moveAccepted=true;
					System.out.print(acceptedMutations+" ");
					
					System.out.print(" "+fitness+" ");
												
					for(int i=0;i<NE.M;i++)
					{
						 NE.In1[i] = in1[i];
						 NE.In2[i] = in2[i];
						 NE.Out1[i] = out1[i];
						 NE.Out2[i] = out2[i];
					}
					
					for(int j=0;j<NE.N;j++)
					{				
						for(int i=0;i<NE.N;i++)
						{
							System.out.print(" "+NE.keff2D[j][i]+" ");
						}					
					}
					System.out.println();
					//System.out.println(PVR[0]+" "+PVR[1]+" "+PVR[2]+" "+PVR[3]);
				
				}	// move accepted: fitness
				/*
					else if(mcRandomNumber > Math.exp(-newFitness*beta))
					{
						acceptedMutations++;
						//System.out.println(" Move Accepted! ");
						fitness = newFitness;
						moveAccepted=true;
						System.out.print(acceptedMutations+" ");
						
						System.out.print(" "+fitness+" ");
													
						for(int i=0;i<NE.M;i++)
						{
							 NE.In1[i] = in1[i];
							 NE.In2[i] = in2[i];
							 NE.Out1[i] = out1[i];
							 NE.Out2[i] = out2[i];
						}
						
						for(int j=0;j<NE.N;j++)
						{				
							for(int i=0;i<NE.N;i++)
							{
								System.out.print(" "+NE.keff2D[j][i]+" ");
							}					
						}
						System.out.println();
						//System.out.println(PVR[0]+" "+PVR[1]+" "+PVR[2]+" "+PVR[3]);
					}
				*/
			}	// move accepted: oscillator
			
			
			// reject move
			else
			{
					//System.out.println(" Rejecting Move! ");		
					moveAccepted = false;
					//time--;
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
