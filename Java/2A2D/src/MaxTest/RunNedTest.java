package NedsMaxTest;

public class RunNedTest 
{
	
	public static void main(String[] args)
	{
		int success = 0;
		int failure = 0;
		for(int trial=0;trial<1000;trial++)
		{
			// true = activator test
			// false = deactivator test
			
			//   activator result: 
			// Trial:  99, Success: 38, Fail: 62 ~ 61%
			// Trial:  999, Success: 399, Fail: 601 rate: 0.6638935108153078
			// 250: 128 success, 123 fail
			
			// deactivator result: 
			// Trial:  99, Success: 33, Fail: 67 ~ 49 %
			// Trial:  250, Success: 60, Fail: 191 rate: 0.31413612565445026
			// Trial:  499, Success: 117, Fail: 383 rate: 0.30548302872062666
			
			
			NeutralEvolutionApp nedTesterObject = new NeutralEvolutionApp(false);
			
			if(nedTesterObject.nedTestResult==true)
			{
				success++;
			}
			else
			{
				failure++;
			}
		
			System.out.println("Trial:  "+trial+", Success: "+success+", Fail: "+failure+" rate: "+(success/(1.0*failure))); 
	
		}
	}
	
}
