package oscillations_test_kda;


import java.io.File;
import java.io.IOException;
import org.tc33.jheatchart.HeatChart;

public class heatMap {

	public heatMap(double[][] data, String filename)
	{
		// Create some dummy data.

		// Step 1: Create our heat map chart using our data.
		HeatChart map = new HeatChart(data);

		// Step 2: Customise the chart.
		map.setTitle("Sequence heat map");
		map.setXAxisLabel("position");
		map.setYAxisLabel("String-Interfaces");

		// Step 3: Output the chart to a file.
		try {
			map.saveToFile(new File(filename+".png"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
