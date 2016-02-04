package jPraat;

import java.util.Collections;
import java.util.ArrayList;

import jPraat.PointProcess;
import jPraat.RealTierPoint;

public class PitchTier
{
    double xmin;
    double xmax;
    int size;
    ArrayList<RealTierPoint> points;
    
    final double NUMundefined = Double.POSITIVE_INFINITY;
    
    public PitchTier(double tmin, double tmax)
    {
    	xmin = tmin;
    	xmax = tmax;
    	size = 0;
    	points = new ArrayList<RealTierPoint>();
    }
    
    // Corresponds to praat function "Add point..."
 	public void addPoint(double time, double frequency)
 	{
 		RealTierPoint point = new RealTierPoint(time, frequency);
 		size++;
 		points.add(point);
 		Collections.sort(points);
 	}
 	
 	PointProcess PitchTier_to_PointProcess ()
 	{
		PointProcess thee = new PointProcess(xmin, xmax, 1000);
		double area = 0.5;   // imagine an event half a period before the beginning
		if (size == 0) return thee;
		for (long interval = 0; interval <= size; interval ++) 
		{
			double t1 = interval == 0 ? xmin : points.get((int) (interval - 1)).number;
			//Melder_assert (NUMdefined (t1));
			double t2 = interval == size ? xmax : points.get((int) interval).number;
			//Melder_assert (NUMdefined (t2));
			double f1 = points.get((int) (interval == 0 ? 0 : (interval - 1))).value;
			//Melder_assert (NUMdefined (f1));
			double f2 = points.get((int) (interval == size ? size-1 : interval)).value;
			//Melder_assert (NUMdefined (f2));
			area += (t2 - t1) * 0.5 * (f1 + f2);
			while (area >= 1.0) 
			{
				double slope = (f2 - f1) / (t2 - t1), discriminant;
				area -= 1.0;
				discriminant = f2 * f2 - 2.0 * area * slope;
				if (discriminant < 0.0) discriminant = 0.0;   // catch rounding errors
				thee.addPoint (t2 - 2.0 * area / (f2 + Math.sqrt(discriminant)));
			}
		}
		return thee;
 	}
 	
 	double RealTier_getValueAtTime (double t) 
 	{
 		long n = size;
 		if (n == 0) return NUMundefined;
 		RealTierPoint pointRight = points.get(0);
 		if (t <= pointRight.number) return pointRight.value;   // constant extrapolation
 		RealTierPoint pointLeft = points.get((int) (n-1));
 		if (t >= pointLeft.number) return pointLeft.value;   // constant extrapolation
 		//Melder_assert (n >= 2);
 		long ileft = AnyTier_timeToLowIndex(t), iright = ileft + 1;
 		//Melder_assert (ileft >= 1 && iright <= n);
 		pointLeft = points.get((int) ileft);
 		pointRight = points.get((int) iright);
 		double tleft = pointLeft.number, fleft = pointLeft.value;
 		double tright = pointRight.number, fright = pointRight.value;
 		return t == tright ? fright   // be very accurate
 			: tleft == tright ? 0.5 * (fleft + fright)   // unusual, but possible; no preference
 			: fleft + (t - tleft) * (fright - fleft) / (tright - tleft);   // linear interpolation
 	}
 	
 	long AnyTier_timeToLowIndex (double time) 
 	{
 		if (size == 0) return 0;   // undefined
 		long ileft = 0, iright = size-1;
 		double tleft = points.get((int) ileft).number;
 		if (time < tleft) return 0;   // offleft
 		double tright = points.get((int) iright).number;
 		if (time >= tright) return iright;
 		//Melder_assert (time >= tleft && time < tright);
 		//Melder_assert (iright > ileft);
 		while (iright > ileft + 1) 
 		{
 			long imid = (ileft + iright) / 2;
 			double tmid = points.get((int) imid).number;
 			if (time < tmid) 
 			{
 				iright = imid;
 				tright = tmid;
 			} 
 			else 
 			{
 				ileft = imid;
 				tleft = tmid;
 			}
 		}
 		/*Melder_assert (iright == ileft + 1);
 		Melder_assert (ileft >= 1);
 		Melder_assert (iright <= my points -> size);
 		Melder_assert (time >= points [ileft] -> number);
 		Melder_assert (time <= points [iright] -> number);*/
 		return ileft;
 	}
 	
 	void PitchTier_multiplyFrequencies (double tmin, double tmax, double factor) 
 	{
 		//Melder_assert (factor > 0.0);
 		for (long i = 0; i < size; i ++) 
 		{
 			RealTierPoint point = points.get((int) i);
 			if (point.number < tmin || point.number > tmax) continue;
 			point.value *= factor;
 		}
 	}
 	
 	void PitchTier_modifyRange_old (double tmin, double tmax, double factor, double fmid) 
 	{
 		for (long i = 0; i < size; i ++) 
 		{
 			
 			RealTierPoint point = points.get((int) i);
 			double f = point.value;
 			if (point.number < tmin || point.number > tmax) 
 			{
 				continue;
 			}
 			f = fmid + (f - fmid) * factor;
 			point.value = f < 0 ? 0 : f;
 		}
 	}
        
    void printInfo()
    {
        System.out.println();

        log("xmin: " + xmin);
        log("xmax: " + xmax);
        log("size: " + size);

        System.out.println();
    }
	
    void log(Object aMsg)
    {
        System.out.println(String.valueOf(aMsg));
    }
 	
}

