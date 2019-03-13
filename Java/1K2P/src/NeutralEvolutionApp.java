
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;


import javax.swing.JFrame;
import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.display.DrawingFrame;
import org.opensourcephysics.display.PlottingPanel;

public class NeutralEvolutionApp 
{

	public static int evolutionTimeTotal = 10000;
	public static double[] PVR;
	public static double[] PVR_index;
	static int[] intervalArray1 = new int[11];
	
	public static double[][] essentiality_N_T;
	
	// for finding PVR
	private static double[] maxValue(double[][] solution, int steps,int numberIterations,int N, boolean recordPeak) 
	{
		double[] max = new double[N];
		// start the max at 10000th point in the array so the initial condition doesn't get counted as the peak or valley.
		for(int i=0;i<N;i++)
		{
			max[i]=solution[i][10000];
		}
		for(int j=0;j<N;j++)
		{
			for (int i = 10000; i < numberIterations*steps; i++) 
			{
				if (solution[j][i] > max[j]) {
					max[j] = solution[j][i];
					if(recordPeak==true)
					{
						PVR_index[j] = i/100.0;
					}
				}
			}
		}
		return max;
	}
	
	private static double[] minValue(double[][] solution, int steps,int numberIterations,int N) 
	{
		double[] min = new double[N];
		// start the max at 10000th point in the array so the initial condition doesn't get counted as the peak or valley.
		for(int i=0;i<N;i++)
		{
			min[i]=solution[i][10000];
		}
		for(int j=0;j<N;j++)
		{
			for (int i = 10000; i < numberIterations*steps; i++) 
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
	
	static int moveAcceptedINT = 0;
	static int moveRejected = 0;
	static int iterations = 0;
			
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException
	{
		//NewtonRasphonNoPlot NR = new NewtonRasphonNoPlot(1,1,1);		// this object determines stables points and eigenvalues.
		EvolutionNProteins NE = new EvolutionNProteins();		// this object has the differential equation solver.
		PVR = new double[NE.N];	// will contain peak to valley ratios for all proteins.
		PVR_index = new double[NE.N];
		essentiality_N_T = new double[NE.N-1][evolutionTimeTotal];	// contains essentiality for P1 and P2, K1 is always essential
		
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
			NE.keff2D[0][0] = NE.computeK(NE.In1, NE.Out1); NE.keff2D[0][1] = NE.computeK(NE.In1, NE.Out2); NE.keff2D[0][2] = NE.computeK(NE.In1, NE.Out3); 
			NE.keff2D[1][0] = NE.computeK(NE.In2, NE.Out1); NE.keff2D[1][1] = NE.computeK(NE.In2, NE.Out2); NE.keff2D[1][2] = NE.computeK(NE.In2, NE.Out3); 
			NE.keff2D[2][0] = NE.computeK(NE.In3, NE.Out1); NE.keff2D[2][1] = NE.computeK(NE.In3, NE.Out2); NE.keff2D[2][2] = NE.computeK(NE.In3, NE.Out3); 
			
			// the temporary array in/out contains the latest values of the In/Out string.
			int[] in2 = new int[NE.In2.length];
			int[] out2 = new int[NE.Out2.length];

			// generate random strings for in2 and out2
			Random random = new Random();
			for(int i=0;i<25;i++)
			{
				in2[i] = random.nextInt(2);
				out2[i] = random.nextInt(2);
			}

			
			// compute all rate constants
			NE.keff2D[0][0] = NE.computeK(NE.In1, NE.Out1); NE.keff2D[0][1] = NE.computeK(NE.In1, out2); NE.keff2D[0][2] = NE.computeK(NE.In1, NE.Out3); 
			NE.keff2D[1][0] = NE.computeK(in2, NE.Out1); NE.keff2D[1][1] = NE.computeK(in2, out2); NE.keff2D[1][2] = NE.computeK(in2, NE.Out3); 
			NE.keff2D[2][0] = NE.computeK(NE.In3, NE.Out1); NE.keff2D[2][1] = NE.computeK(NE.In3, out2); NE.keff2D[2][2] = NE.computeK(NE.In3, NE.Out3); 

							 			
			//essentiality_N_T[0][time] = 0 or 1;
			for (int i = 0; i < NE.numberIterations; i++) 
			{
				NE.Essentiality(1, 1, 1);
				NE.RKSolveN();			
			}			
			for(int i=0;i<NE.N;i++)
			{
				PVR[i] = maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N, true)[i]/minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N)[i];
				//System.out.println(maxValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N)[i]+" "+minValue(NE.y_solution,NE.steps,NE.numberIterations,NE.N)[i]);
			}
			
			// if all PVR's > then a threshold, system oscillating, else roll back time.
			// all accepted (oscillator) results
			if( (PVR[0]-1)<1e-2 & (PVR[1]-1)<1e-2 & (PVR[2]-1)<1e-2 )
			{
				moveRejected++;
				time--;
				//System.out.println("Move rejected ");
			}	// move rejected
			else
			{			
				moveAcceptedINT++;
				moveAccepted=true;
			//if(time%100==0)
				
				/*
				// print rate constants
				for(int j=0;j<NE.N;j++)
				{				
					for(int i=0;i<NE.N;i++)
					{
						System.out.print(" "+NE.keff2D[j][i]+" ");
					}					
				}
				
				// print pvr
				System.out.println(PVR_index[0]+" "+PVR[0]+" "+PVR_index[1]+" "+PVR[1]+" "+PVR_index[2]+" "+PVR[2]);
				*/
				
				/*
				if(time%1000==0)
				{					
					PrintWriter writer = new PrintWriter("/Users/ammar/Dropbox/Clark_Research/essentiality/1K2P_Solutions/1K2P_solution_T_"+Integer.toString(time)+".txt", "UTF-8");
					writer.println("# "+time);
					for (int t = 0; t < NE.numberIterations*NE.steps; t++) 
					{
						writer.println(t + " " + NE.y_solution[0][t]+ " " + NE.y_solution[1][t]+ " " + NE.y_solution[2][t]);
					}
					writer.close();
				}
				*/
			}	// move accepted
			
			if(moveAcceptedINT+moveRejected!=0 & iterations%100==0)
			{
				//System.out.println(iterations+" "+moveAcceptedINT+" "+moveRejected);
				System.out.println(iterations+" "+(moveAcceptedINT)/(1.0*moveAcceptedINT+1.0*moveRejected));
			}
			iterations++;
			
			// Solution graphics
			if(graphicsEnabled)
			{
				String evolutionTime = Integer.toString(time);
				// Plot properties
				//PlottingPanel plot = new PlottingPanel("time", "", evolutionTime);
				PlottingPanel plot = new PlottingPanel("time", "", Integer.toString(iterations));
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
