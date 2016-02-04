package jPraat;

// class for praat object LPC
public class LPC
{
	// all members are named the same as they are in praat
	double xmin;
	double xmax;
	int nx;
	double dx;
	double x1;
	double SamplingPeriod;
	int maxCoefficients;
	LPC_Frame[] frames;
	final double NUMpi = 3.1415926535897932384626433832795028841972;
	
	
	public LPC(double xmin0, double xmax0, int numFrames, double dt, double t1, int predictionOrder, double samplingPeriod)
	{
		xmin = xmin0;
		xmax = xmax0;
		maxCoefficients = predictionOrder;
		nx = numFrames;
		dx = dt;
		x1 = t1;
		SamplingPeriod = samplingPeriod;
		frames = new LPC_Frame[numFrames];
	}
	
	public void printInfo()
	{
		System.out.println();
		log("xmin: " + xmin);
	    log("xmax: " + xmax);
	    log("NumSamples: " + nx);
	    log("SamplingRate=dx: " + dx);
	    log("x1: " + x1);
	    log("SamplingPeriod: " + SamplingPeriod);
	    log("maxCoefficiients: " + maxCoefficients);
	    System.out.println();
 		frames[nx-1].printInfo();
 		System.out.println();
 		//frames[1].printInfo();
 		//System.out.println();
 		//frames[2].printInfo();
 		//System.out.println();
 		//frames[3].printInfo();
 		//System.out.println();
 		//frames[4].printInfo();
 		//System.out.println();
 		//frames[5].printInfo();
	}
	
	void log(Object aMsg)
	{
	    System.out.println(String.valueOf(aMsg));
	}
	
	public LPC() {} //empty constructor
	
	// LPC constructor. Corresponds to praat function Sound : "To LPC (burg)..."
	LPC Sound_to_LPC_burg (Wave me, int predictionOrder, double analysisWidth, double dt, double preEmphasisFrequency)
	{
		LPC thee = _Sound_to_LPC (me, predictionOrder, analysisWidth, dt, preEmphasisFrequency, 0, 0);
		return thee;
	}
	
	LPC _Sound_to_LPC (Wave me, int predictionOrder, double analysisWidth, double dt, double preEmphasisFrequency, double tol1, double tol2) 
	{
		double[] t1 = new double[1]; 
		double samplingFrequency = 1.0 / me.dx;
		double windowDuration = 2.0 * analysisWidth; /* gaussian window */
		long[] nFrames = new long[1];
		long frameErrorCount = 0;
		
		if (Math.floor(windowDuration / me.NumSamples) < predictionOrder + 1) {} // should return false
		if (windowDuration > me.NumSamples * me.dx) 
		{
			windowDuration = me.NumSamples * me.dx;
		}
		Sampled_shortTermAnalysis (me, windowDuration, dt, nFrames, t1);
		Wave sound = me.Wave_copy();
		Wave sframe = new Wave(); 
		sframe = sframe.Sound_createSimple (1, windowDuration, samplingFrequency);
		Wave window = new Wave();
		window = window.Sound_createGaussian (windowDuration, samplingFrequency);
		//window.printInfo();
		LPC thee = new LPC(me.xmin, me.xmax, (int)nFrames[0], dt, t1[0], predictionOrder, me.dx);
		//System.out.println("sframe num: "+sframe.NumSamples);
		
		if (preEmphasisFrequency < samplingFrequency / 2) 
		{
			sound.Sound_preEmphasis(preEmphasisFrequency);
		}
		
		for (long i = 0; i < nFrames[0]; i++) 
		{
			//LPC_Frame lpcframe = thee.frames[(int) i];
			double t = thee.Sampled_indexToX(i);
			LPC_Frame lpcframe = new LPC_Frame(predictionOrder);
			thee.frames[(int) i] = lpcframe;
			sound.Sound_into_Sound (sframe, t - windowDuration / 2.0);
			sframe.subtractMean();
			sframe.Sounds_multiply (window);
			if (sframe.Sound_into_LPC_Frame_burg (lpcframe) != 0) 
			{
				frameErrorCount++;
			}
		}
		return thee;
	}
	
	// Praat helper function
	long Sampled_xToNearestIndex (double x) 
	{ 
		return (long) Math.round((x - x1) / dx); 
	}
	
	// Praat helper function
	double Sampled_indexToX (long index) 
	{ 
		return x1 + index * dx;
	}
	
	void Sampled_shortTermAnalysis (Wave me, double windowDuration, double timeStep, long[] numberOfFrames, double[] firstTime) 
	{
		//Melder_assert (windowDuration > 0.0);
		//Melder_assert (timeStep > 0.0);
		double myDuration = me.dx * me.NumSamples;
		//if (windowDuration > myDuration)
			//Melder_throw (me, ": shorter than window length."); 
		numberOfFrames[0] = (long) (Math.floor((myDuration - windowDuration) / timeStep) + 1);
		//Melder_assert (numberOfFrames[0] >= 1);
		double ourMidTime = me.x1 - 0.5 * me.dx + 0.5 * myDuration;
		double thyDuration = numberOfFrames[0] * timeStep;
		firstTime[0] = ourMidTime - 0.5 * thyDuration + 0.5 * timeStep;
	}

	

}