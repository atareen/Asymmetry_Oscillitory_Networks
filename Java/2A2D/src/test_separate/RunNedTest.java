package ranjan_test_separate;

public class RunNedTest 
{
	
	public static void main(String[] args)
	{
		int success = 0;
		int failure = 0;
		int numberTrials = 1000;
		for(int trial=1;trial<=numberTrials;trial++)
		{
			
			NeutralEvolutionApp nedTesterObject = new NeutralEvolutionApp(false);
			
			if(nedTesterObject.nedTestResult==true)
			{
				success++;
			}
			else
			{
				failure++;
			}
		
			System.out.println("Trial:  "+trial+", Success: "+success+", Fail: "+failure+" success over total: "+(success/(1.0*trial))); 
	
		}
	}
	
}
