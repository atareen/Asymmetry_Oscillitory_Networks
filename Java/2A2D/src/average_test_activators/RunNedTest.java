package ranjan_average_test_activators;

public class RunNedTest 
{
	
	public static void main(String[] args)
	{
		int success = 461;
		int failure = 244;
		int numberTrials = 1000;
		for(int trial=706;trial<=numberTrials;trial++)
		{
			
			// kda
			// results Trial:  999, Success: 655, Fail: 345 success over total: 0.655
			
			
			//kad
			// Trial:  542, Success: 310, Fail: 232 success over total: 0.5719557195571956
			// Trial:  1000, Success: 593, Fail: 407 success over total: 0.593
			
			// ranjan average test:
			//Trial:  1000, Success: 191, Fail: 809 success over total: 0.191
			
			// ranjan average test activators
			//Trial: 	705,	Success:	461,	Fail:	244	success	over	total:	0.653900709
			
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
