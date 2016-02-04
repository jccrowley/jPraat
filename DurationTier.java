package jPraat;

import java.util.ArrayList;
import java.util.Collections;

// Class for praat object DurationTier.
public class DurationTier 
{
	// All members are named the same as they are in praat.
	double xmin;
	double xmax;
	long size;
	ArrayList<RealTierPoint> points;
	
	final double NUMundefined = Double.POSITIVE_INFINITY;
	
	// DurationTier constructor. Corresponds to praat function "Create DurationTier..."
	public DurationTier(double min, double max)
	{
		xmin = min;
		xmax = max;
		size = 0;
		points = new ArrayList<RealTierPoint>();
	}
	
	// Corresponds to praat function "Add point..."
	public void addPoint(double n, double v)
	{
		size++;
		points.add(new RealTierPoint(n, v));
		Collections.sort(points);
	}

	
	// Produce deep copy of DurationTier object
	public DurationTier DurationTierPoint_copy()
	{
		DurationTier durationTierCopy = new DurationTier(xmin, xmax);
		for(int i = 0; i < size; i++)
			durationTierCopy.points.add(points.get(i).RealTierPoint_copy());
		return durationTierCopy;
	}
	
	// Print contents of DurationTier object
	public void print()
	{
		System.out.println("xmin: " + xmin);
 		System.out.println("xmax: " + xmax);
 		System.out.println("NumPoints: " + size);
 		for(int i = 0; i < size; i++)
 			points.get(i).print();
	}
	
	double RealTier_getArea (double tmin, double tmax) 
	{
		long n = size, imin, imax;
		if (n == 0) return NUMundefined;
		if (n == 1) return (tmax - tmin) * points.get(0).value;
		imin = AnyTier_timeToLowIndex(tmin);
		if (imin == n) return (tmax - tmin) * points.get((int) (n-1)).value;
		imax = AnyTier_timeToHighIndex(tmax);
		if (imax == 1) return (tmax - tmin) * points.get(0).value;
		//Melder_assert (imin < n);
		//Melder_assert (imax > 1);
		/*
		 * Sum the areas between the points.
		 * This works even if imin is 0 (offleft) and/or imax is n + 1 (offright).
		 */
		double area = 0.0;
		for (long i = imin; i < imax; i ++) 
		{
			double tleft, fleft, tright, fright;
			if (i == imin) 
			{
				tleft = tmin;
				fleft = RealTier_getValueAtTime(tmin);
			}
			else 
			{
				tleft = points.get((int) i).number;
				fleft = points.get((int) i).value;
			}
			if (i + 1 == imax) 
			{
				tright = tmax;
				fright = RealTier_getValueAtTime(tmax);
			}
			else 
			{
				tright = points.get((int) (i + 1)).number;
				fright = points.get((int) (i + 1)).value;
			}
			area += 0.5 * (fleft + fright) * (tright - tleft);
		}
		return area;
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
 		long ileft = 0, iright = (size-1);
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
 	
 	long AnyTier_timeToHighIndex (double time) 
 	{
 		if (size == 0) return 0;   // undefined; is this right?
 		long ileft = 0, iright = size-1;
 		double tleft = points.get((int) ileft).number;
 		if (time <= tleft) return 1;
 		double tright = points.get((int) iright).number;
 		if (time > tright) return iright + 1;   // offright
 		//Melder_assert (time > tleft && time <= tright);
 		//Melder_assert (iright > ileft);
 		while (iright > ileft + 1) 
 		{
 			long imid = (ileft + iright) / 2;
 			double tmid = points.get((int) imid).number;
 			if (time <= tmid) 
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
 		return iright;
 	}
	
	
}
