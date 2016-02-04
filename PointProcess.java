package jPraat;

import java.util.ArrayList;
import java.util.Collections;

// Class for praat object PointProcess.
public class PointProcess 
{
	// All members are named the same as they are in praat.
	double xmin;
	double xmax;
	int nt;
	int maxnt;
	ArrayList<Double> t;
	final double NUMpi = 3.1415926535897932384626433832795028841972;
	
	// PointProcess constructor. Corresponds to praat function "Create empty PointProcess..."
	public PointProcess(double min, double max)
	{
		xmin = min;
		xmax = max;
		nt = 0;
		maxnt = 0;
		t = new ArrayList<Double>();
	}
	
	// Alternative PointProcess constructor.
	public PointProcess(double tmin, double tmax, int initialMaxnt)
	{
		xmin = tmin;
		xmax = tmax;
		if (initialMaxnt < 1) initialMaxnt = 1;
		nt = 0;
		maxnt = initialMaxnt;
		t = new ArrayList<Double>();
	}
	
	// Corresponds to praat function "Add point..."
	public void addPoint(double point)
	{
		nt++;
		t.add(point);
		Collections.sort(t);
	}
	
	public void printInfo()
	{
        System.out.println();
		System.out.println("xmin: " + xmin);
 		System.out.println("xmax: " + xmax);
 		System.out.println("NumFrames: " + nt);
 		for(int i = 0; i < 10; i++)
 			System.out.println("Point "+(i+1)+": "+t.get(i));
                System.out.println();
	}
	
	// Corresponds to praat function "To Sound (pulse train)..."
	public Wave toSoundPulseTrain(double samplingFrequency, double adaptFactor, double adaptTime, long interpolationDepth)
	{	
		long sound_nt = (long) (1 + Math.floor((xmax - xmin) * samplingFrequency));   // >= 1
		double dt = 1.0 / samplingFrequency;
		double tmid = (xmin + xmax) / 2;
		double t1 = tmid - 0.5 * (sound_nt - 1) * dt;
		Wave sound = new Wave (1, xmin, xmax, (int)sound_nt, dt, t1); // originally thee
		double[] sounds = sound.samples1;
		for (long it = 0; it < nt; it++) 
		{
			double tt = t.get((int) it), amplitude = 0.9, angle, halfampsinangle;
			long mid = sound.Sampled_xToNearestIndex(tt);
			if (it <= 2 || t.get((int)(it - 2)) < t.get((int)it) - adaptTime) 
			{
				amplitude *= adaptFactor;
				if (it == 0 || t.get((int)(it - 1)) < t.get((int)it) - adaptTime)
					amplitude *= adaptFactor;
			}
			long begin = mid - interpolationDepth, end = mid + interpolationDepth;
			if (begin < 0) begin = 0;
			if (end > sound.NumSamples) end = sound.NumSamples;
			angle = NUMpi * (sound.Sampled_indexToX(begin) - tt) / sound.dx;
			halfampsinangle = 0.5 * amplitude * Math.sin(angle);
			for (long j = begin; j < end; j++) 
			{
				if (Math.abs(angle) < 1e-6)
					sounds[(int) j] += amplitude;
				else if (angle < 0.0)
					sounds[(int) j] += halfampsinangle * (1 + Math.cos(angle / (mid - begin + 1))) / angle;
				else
					sounds[(int) j] += halfampsinangle * (1 + Math.cos(angle / (end - mid + 1))) / angle;
				angle += NUMpi;
				halfampsinangle = - halfampsinangle;
			}
		}
		return sound;
	}
	
	// nt must not equal 0
	public Matrix PointProcess_to_Matrix ()
	{
		Matrix thee = new Matrix (1.0, (double)nt, nt, 1.0, 1.0, 1.0, 1.0, 1, 1.0, 1.0);
		for (long i = 0; i < nt; i ++)
			thee.z [0] [(int) i] = t.get((int) i);
		return thee;
	}
	
	long PointProcess_getNearestIndex (double t) 
	{
		if (nt == 0)
			return 0;
		if (t <= this.t.get(0))
			return 1;
		if (t >= this.t.get(nt-1))
			return nt-1;
		/* Start binary search. */
		long left = 0, right = nt-1;
		while (left < right - 1) 
		{
			long mid = (left + right) / 2;
			if (t >= this.t.get((int) mid)) left = mid; else right = mid;
		}
		//Melder_assert (right == left + 1);
		return t - this.t.get((int) left) < this.t.get((int) right) - t ? left : right;
	}
	
}
