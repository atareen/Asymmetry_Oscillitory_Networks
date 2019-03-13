package neutralEvolution;

import java.io.File;
import java.io.IOException;

import org.tc33.jheatchart.HeatChart;

public class heatMap {

	public static void main(String[] args)
	{
		// Create some dummy data.
		double[][] data = new double[][]{
											{
												0.78862, 0.8275, 0.78479, 0.79232, 0.81397, 0.8145, 0.75571, 0.78584, 0.76881, 0.80664, 0.7626, 0.80713, 0.7856, 0.8025, 0.77151, 0.78122, 0.76035, 0.79176, 0.82021, 0.80847, 0.82217, 0.82857, 0.78757, 0.7894, 0.83018
											},
			                                 {
												0.47672, 0.49648, 0.44484, 0.45277, 0.45272, 0.48155, 0.47006, 0.46686, 0.48973, 0.43371, 0.43749, 0.46503, 0.45344, 0.45316, 0.46786, 0.4888, 0.46453, 0.44741, 0.49801, 0.46358, 0.49614, 0.49039, 0.46329, 0.40117, 0.49318
											},
			                                 {
												0.7037, 0.75337, 0.72578, 0.72899, 0.73356, 0.74687, 0.71851, 0.73432, 0.77145, 0.69522, 0.6302, 0.72785, 0.62414, 0.66499, 0.64089, 0.73124, 0.71517, 0.74208, 0.73907, 0.70304, 0.77339, 0.78258, 0.73604, 0.5534, 0.74108
											},
			                                 {
												0.65543, 0.63272, 0.69781, 0.64424, 0.65904, 0.65807, 0.65612, 0.63276, 0.63677, 0.66835, 0.63048, 0.65922, 0.65229, 0.68662, 0.68953, 0.64539, 0.71426, 0.70225, 0.60748, 0.68616, 0.6869, 0.65761, 0.64906, 0.68238, 0.52407
											}
										};

		// Step 1: Create our heat map chart using our data.
		HeatChart map = new HeatChart(data);

		// Step 2: Customise the chart.
		map.setTitle("Sequence heat map");
		map.setXAxisLabel("X Axis");
		map.setYAxisLabel("Y Axis");

		// Step 3: Output the chart to a file.
		try {
			map.saveToFile(new File("java-heat-chart.png"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
