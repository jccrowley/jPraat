
package jPraat;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.math.RoundingMode;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import fftpack.RealDoubleFFT;




/* Author: J. Colin Crowley
 * Representation of WAV data in time domain.
 */

//Class representing praat Sound objects.
public class Wave 
{
	// globals
	final double NUMundefined = Double.POSITIVE_INFINITY;
	final double NUMpi = 3.1415926535897932384626433832795028841972;
	final double NUMlog2e = 1.44269504089;
	final double NUM_goldenSection = 0.6180339887498948482045868343656381177203;
	final double step = 0.000001;
	
	// Members with direct links to members in praat.
    int NumChannels;					// ny in praat
    int NumSamples;						// nx in praat
    double[] samples1;					// z in praat
    double xmin;						// xmin in praat. Start time
    double xmax;						// xmax in praat. End time
    double dx;							// dx in praat
    double x1;							// x1 in praat
    int SamplingFrequency;
    
    // Members used to retrieve data from file.
    ByteBuffer bb;
    byte[] bytes;
    
    // Wave constructor. Construct the object from a .wav file.
	public Wave(String filename) throws IOException 
	{
		bb = ByteBuffer.allocate(2);
		bb.order(ByteOrder.LITTLE_ENDIAN);
		File file = new File(filename);
		int size = (int) file.length();
		bytes = new byte[size];
		BufferedInputStream buf = new BufferedInputStream(new FileInputStream(file));
		buf.read(bytes, 0, bytes.length);
		buf.close();
		/*Path path = Paths.get(filename);
		bytes = Files.readAllBytes(path);*/
		byte[] bytes1 = new byte[bytes.length-44];
		NumSamples = (bytes.length-44)/2;
	    int[] samples = new int[NumSamples];
	    samples1 = new double[NumSamples];
	    for(int x = 44; x < bytes.length; x++)
	    {
	    	bytes1[x-44] = bytes[x];
	    }
	    for(int x = 0; x < samples.length; x++)
	    {
	    	bb.put(bytes[(x*2)+44]);
	    	bb.put(bytes[(x*2)+45]);
	    	bb.flip();
	    	samples[x] = bb.getShort();
	    	samples1[x] = (double)samples[x] / 32768;
	    	bb.flip();
	    	bb.clear();
	    }
	    
	    NumChannels = WaveFormat.NumChannels;
	    bb = ByteBuffer.allocate(4);
	    bb.order(ByteOrder.LITTLE_ENDIAN);
	    bb.put(bytes[24]);
    	bb.put(bytes[25]);
    	bb.put(bytes[26]);
    	bb.put(bytes[27]);
    	bb.flip();
    	SamplingFrequency = bb.getInt();
    	bb.flip();
    	bb.clear();
	    dx = 1.0 / (double)SamplingFrequency;
	 	xmin = 0.0;
	 	xmax = NumSamples * dx;
	 	//x1 = ((xmin + xmax) / 2) - 0.5 * (NumSamples - 1) * dx;
	 	x1 = 0.5 * dx;
	}

    // Wave constructor. Construct the object from a .pcm file.
    public Wave(String filename, int samplingFrequency) throws IOException
    {
        bb = ByteBuffer.allocate(2);
        bb.order(ByteOrder.LITTLE_ENDIAN);
        File file = new File(filename);// can pass absolute filename or directory and filename
        int size = (int) file.length();
        byte[] bytes = new byte[size];
        BufferedInputStream buf = new BufferedInputStream(new FileInputStream(file));
        buf.read(bytes, 0, bytes.length);
        buf.close();
        NumSamples = (bytes.length)/2;
        int[] samples = new int[NumSamples];
        samples1 = new double[NumSamples];
        for(int x = 0; x < samples.length; x++)
        {
            bb.put(bytes[(x*2)]);
            bb.put(bytes[(x*2)+1]);
            bb.flip();
            samples[x] = bb.getShort();
            samples1[x] = (double)samples[x] / 32768;
            bb.flip();
            bb.clear();
        }
        NumChannels = WaveFormat.NumChannels;
        SamplingFrequency = samplingFrequency;
        dx = 1.0 / (double)SamplingFrequency;
        xmin = 0.0;
        xmax = NumSamples * dx;
        x1 = ((xmin + xmax) / 2) - 0.5 * (NumSamples - 1) * dx;
    }

    // Wave constructor. Construct the object from an ArrayList of bytes.
    public Wave(ArrayList<Byte> sample, int samplingFrequency) throws IOException
    {
        bb = ByteBuffer.allocate(2);
        bb.order(ByteOrder.LITTLE_ENDIAN);
        Byte[] bytes = sample.toArray(new Byte[sample.size()]);
        NumSamples = (bytes.length)/2;
        int[] samples = new int[NumSamples];
        samples1 = new double[NumSamples];
        for(int x = 0; x < samples.length; x++)
        {
            bb.put(bytes[(x*2)]);
            bb.put(bytes[(x*2)+1]);
            bb.flip();
            samples[x] = bb.getShort();
            samples1[x] = (double)samples[x] / 32768;
            bb.flip();
            bb.clear();
        }
        NumChannels = WaveFormat.NumChannels;
        SamplingFrequency = samplingFrequency;
        dx = 1.0 / (double)SamplingFrequency;
        xmin = 0.0;
        xmax = NumSamples * dx;
        x1 = ((xmin + xmax) / 2) - 0.5 * (NumSamples - 1) * dx;
    }
	
	// Alternative Wave constructor. Sound_create
	public Wave(int channels, double min, double max, int nx, double dxp, double x1p)
	{
		NumChannels = channels;
		xmin = min;
		xmax = max;
		NumSamples = nx;
		dx = dxp;
		x1 = x1p;
		samples1 = new double[NumSamples];
	}
	
	Wave Sound_createSimple (long numberOfChannels, double duration, double samplingFrequency) 
	{
		return Sound_create (numberOfChannels, 0.0, duration, (int)Math.floor(duration * samplingFrequency + 0.5),
			1 / samplingFrequency, 0.5 / samplingFrequency);
	}
	
	Wave Sound_create (long numberOfChannels, double xmin, double xmax, long nx, double dx, double x1) 
	{
		Wave me = new Wave((int)numberOfChannels, xmin, xmax,(int)nx, dx, x1);
		return me;
	}
	
	Wave Sound_createGaussian (double windowDuration, double samplingFrequency) 
	{
		Wave me = Sound_createSimple (1, windowDuration, samplingFrequency);
		double[] s = me.samples1;
		double imid = 0.5 * (me.NumSamples + 1), edge = Math.exp(-12.0);
		for (long i = 0; i < me.NumSamples; i++) 
		{
			s[(int) i] = (Math.exp(-48.0 * ((i+1) - imid) * ((i+1) - imid) / (me.NumSamples + 1) / (me.NumSamples + 1)) - edge) / (1 - edge);
		}
		return me;
	}
	
	public Wave(double min, double max, double dxp, int formula)
	{
		xmin = min;
		xmax = max;
		dx = dxp;
		NumChannels = 1;
		NumSamples = (int) ((xmax - xmin) / dx);
		samples1 = new double[NumSamples];
		if(formula == 1)
		{
			for(int i = 0; i < NumSamples; i++)
				samples1[i] = 0;
		}
		else if(formula == -1)
		{
			Random randGen = new Random();
			for(int i = 0; i < NumSamples; i++)
				samples1[i] = ((double)randGen.nextInt(65535) - 32768.0) / 32768.0;
		}
		x1 = dx / 2;
	}
	
	public Wave() {} //empty constructor
	
	// Make a copy of a Wave object.
	public Wave Wave_copy()
	{
		Wave newWave = new Wave(NumChannels, xmin, xmax, NumSamples, dx, x1);
		newWave.samples1 = new double[NumSamples];
		for(int i = 0; i < NumSamples; i++)
			newWave.samples1[i] = samples1[i];
		return newWave;
	}
	
	/*void updateBytes()
	{
		bytes1 = new byte[samples1.length*2];
		for(int i = 0; i < NumSamples; i++)
		{
			bb = ByteBuffer.allocate(2);
			short tempShort = (short)(Math.round(samples1[i] * 32768));
			if(tempShort > 32767)
				tempShort = 32767;
			if(tempShort < -32768)
				tempShort = -32768;
			bb.putShort(tempShort);
			byte[] temp = bb.array();
			bytes1[2*i] = temp[1];
			bytes1[2*i + 1] = temp[0];
		}
	}*/
	
	// Output Wave data to file.
	public void exportWave(String aFileName) throws IOException
	{
		//updateBytes();
		byte[] temp;
		//Path path = Paths.get(aFileName);
		FileOutputStream os = new FileOutputStream(aFileName);
		byte[] outputBytes = new byte[(samples1.length*2)+44];
		for(int i = 0; i < 44; i++)
		{
			outputBytes[i] = WaveFormat.metaData[i];
		}
		// figure out chunksize, samplerate, byterate, subchunk2size
		int ChunkSize = 36 + (NumSamples * 2);
		int SampleRate = (int) (1.0/dx);
		int ByteRate = (int) (2.0/dx);
		int Subchunk2Size = NumSamples*2;
		bb = ByteBuffer.allocate(4);
		bb.putInt(ChunkSize);
		temp = bb.array();
		outputBytes[4] = temp[3];
		outputBytes[5] = temp[2];
		outputBytes[6] = temp[1];
		outputBytes[7] = temp[0];
		bb = ByteBuffer.allocate(4);
		bb.putInt(SampleRate);
		temp = bb.array();
		outputBytes[24] = temp[3];
		outputBytes[25] = temp[2];
		outputBytes[26] = temp[1];
		outputBytes[27] = temp[0];
		bb = ByteBuffer.allocate(4);
		bb.putInt(ByteRate);
		temp = bb.array();
		outputBytes[28] = temp[3];
		outputBytes[29] = temp[2];
		outputBytes[30] = temp[1];
		outputBytes[31] = temp[0];
		bb = ByteBuffer.allocate(4);
		bb.putInt(Subchunk2Size);
		temp = bb.array();
		outputBytes[40] = temp[3];
		outputBytes[41] = temp[2];
		outputBytes[42] = temp[1];
		outputBytes[43] = temp[0];
		for(int i = 0; i < NumSamples; i++)
		{
			bb = ByteBuffer.allocate(2);
			short tempShort = (short)(Math.round(samples1[i] * 32768));
			if(tempShort > 32767)
				tempShort = 32767;
			if(tempShort < -32768)
				tempShort = -32768;
			bb.putShort(tempShort);
			temp = bb.array();
			outputBytes[44+(2*i)] = temp[1];
			outputBytes[44+(2*i + 1)] = temp[0];
		}
		os.write(outputBytes);
		os.close();
	    //Files.write(path, outputBytes);
	}
	
	public void printInfo()
	{
		System.out.println();
	    log("xmin: " + xmin);
	    log("xmax: " + xmax);
	    log("NumSamples: " + NumSamples);
	    log("SamplingRate=dx: " + dx);
	    log("x1: " + x1);
	    for(int i = 0; i < 10; i++)
	    {
	    	System.out.println("Sample "+(i+1)+": "+samples1[i]);	
	    }
	    		
		
	}
	
	void log(Object aMsg)
	{
	    System.out.println(String.valueOf(aMsg));
	}
	
	// Corresponds to praat function "Filter (pass Hann band)..." (Calls same function in Spectrum class)
	Wave filterPassHannBand(double fmin, double fmax0, double smooth)
	{
	    //Spectrum spectrum = new Spectrum(samples1);
		Spectrum spectrum = Sound_to_Spectrum(true);
		//spectrum.print();
	    spectrum.filterPassHannBand(fmin, fmax0, smooth);
	    //spectrum.print();
	    return spectrum.Spectrum_to_Sound();
	    //double[] tempSamples = spectrum.getSoundSamples();
	    //System.arraycopy(tempSamples, 0, samples1, 0, NumSamples);
	}
	
	// Corresponds to praat functions "Scale..." and "Scale peak..." (same function)
	void scalePeak (double scale) 
	{
		double extremum = 0.0;
		for (long channel = 1; channel <= NumChannels; channel ++) 
		{
			for (long i = 0; i < NumSamples; i ++) 
			{
				if (Math.abs(samples1[(int)i]) > extremum) 
					extremum = Math.abs(samples1[(int)i]);
			}
		}
		if (extremum != 0.0) 
		{
			multiplyByScalar (scale / extremum);
		}
	}
	
	void scaleIntensity (double newAverageIntensity) 
	{
		double currentIntensity = getIntensity_dB(), factor;
		if (currentIntensity == NUMundefined) return;
		factor = Math.pow(10, (newAverageIntensity - currentIntensity) / 20.0);
		for (long channel = 1; channel <= NumChannels; channel++) 
		{
			for (long i = 0; i < NumSamples; i ++) 
			{
				samples1[(int)i] *= factor;
			}
		}
	}
	
	double getIntensity_dB () 
	{
		long[] n = new long[1];
		double sum2 = getSumOfSquares (0, 0, n);
		if(sum2 != NUMundefined && sum2 != 0.0)
			return 10 * Math.log10(sum2 / (n[0] * NumChannels) / 4.0e-10);
		else
			return NUMundefined;
	}
	
	double getSumOfSquares (double xmin0, double xmax0, long[] n) 
	{
		if (xmax0 <= xmin0) 
		{ 
			xmin0 = xmin; 
			xmax0 = xmax; 
		}
		long[] imin = new long[1];
		long[] imax = new long[1];
		n[0] = Sampled_getWindowSamples (xmin0, xmax0, imin, imax);
		if (n[0] < 1) return NUMundefined;
		double sum2 = 0.0;
		for (long channel = 1; channel <= NumChannels; channel ++) 
		{
			double[] amplitude = samples1;
			for (long i = imin[0]; i <= imax[0]; i ++) 
			{
				double value = amplitude [(int)i];
				sum2 += value * value;
			}
		}
		return sum2;
	}
	
	long Sampled_getWindowSamples (double xmin0, double xmax0, long[] ixmin, long[] ixmax) 
	{
		double rixmin = Math.ceil((xmin0 - x1) / dx);
		double rixmax = Math.floor((xmax0 - x1) / dx);   // could be above 32-bit LONG_MAX
		ixmin[0] = rixmin < 0.0 ? 0 : (long) rixmin;
		ixmax[0] = rixmax > (double) (NumSamples-1) ? (NumSamples-1) : (long) rixmax;
		if (ixmin[0] > ixmax[0]) return 0;
		return ixmax[0] - (ixmin[0] - 1);
	}
	
	// Helper praat function
	void multiplyByScalar (double scalar) 
	{
		for (long channel = 1; channel <= NumChannels; channel ++) 
		{
			for (long i = 0; i < NumSamples; i ++) 
			{
				samples1[(int)i] *= scalar;
			}
		}
	}
	
	// Return samples
	double[] getData()
	{
	    return samples1;
	}
	
	// Return samples as bytes
	/*byte[] getBytes()
	{
	    return bytes1;
	}*/
	
	// Corresponds to praat function "Subtract mean"
	void subtractMean()
	{
		for (long channel = 1; channel <= NumChannels; channel++) 
		{
			double sum = 0.0;
			for (long i = 0; i < NumSamples; i++) 
			{
				sum += samples1[(int) i];
			}
			double mean = sum / NumSamples;
			for (long i = 0; i < NumSamples; i ++) 
			{
				samples1[(int) i] -= mean;
			}
		}
	}
	
	// Corresponds to praat function "Get absolute extremum..."
	double getAbsoluteExtrenum()
	{
		double absExtrenum = 0.0;
		for (long channel = 1; channel <= NumChannels; channel++) 
		{
			for (long i = 0; i < NumSamples; i++) 
			{
				if(Math.abs(samples1[(int) i]) > absExtrenum)
					absExtrenum = Math.abs(samples1[(int) i]);
			}
		}
		return absExtrenum;
	}
	
	double Vector_getAbsoluteExtremum (double xmin, double xmax, int interpolation) 
	{
		double minimum = Math.abs(Vector_getMinimum(xmin, xmax, interpolation));
		double maximum = Math.abs(Vector_getMaximum(xmin, xmax, interpolation));
		return minimum > maximum ? minimum : maximum;
	}
	
	double Vector_getMinimum (double xmin, double xmax, int interpolation) 
	{
		double[] minimum = new double[1];
		Vector_getMinimumAndXAndChannel (xmin, xmax, interpolation, minimum, null, null);
		return minimum[0];
	}
	
	double Vector_getMaximum (double xmin, double xmax, int interpolation) 
	{
		double[] maximum = new double[1];
		Vector_getMaximumAndXAndChannel(xmin, xmax, interpolation, maximum, null, null);
		return maximum[0];
	}
	
	void Vector_getMaximumAndXAndChannel (double xmin, double xmax, int interpolation,
			double[] return_maximum, double[] return_xOfMaximum, long[] return_channelOfMaximum)
	{
		double[] maximum = new double[1], xOfMaximum = new double[1];
		long channelOfMaximum = 1;
		Vector_getMaximumAndX (xmin, xmax, 0, interpolation, maximum, xOfMaximum);
		if (return_maximum != null) return_maximum[0] = maximum[0];
		if (return_xOfMaximum != null) return_xOfMaximum[0] = xOfMaximum[0];
		if (return_channelOfMaximum != null )return_channelOfMaximum[0] = channelOfMaximum;
	}
	
	void Vector_getMinimumAndXAndChannel (double xmin, double xmax, int interpolation,
			double[] return_minimum, double[] return_xOfMinimum, long[] return_channelOfMinimum)
	{
		double[] minimum = new double[1], xOfMinimum = new double[1];
		long channelOfMinimum = 1;
		Vector_getMinimumAndX(xmin, xmax, 0, interpolation, minimum, xOfMinimum);
		if (return_minimum != null) return_minimum[0] = minimum[0];
		if (return_xOfMinimum != null) return_xOfMinimum[0] = xOfMinimum[0];
		if (return_channelOfMinimum != null) return_channelOfMinimum[0] = channelOfMinimum;
	}
	
	//Melder_assert (channel >= 1 && channel <= my ny);
	void Vector_getMaximumAndX (double xmin0, double xmax0, long channel, int interpolation,
			double[] return_maximum, double[] return_xOfMaximum)
	{
		long[] imin = new long[1], imax = new long[1];
		long i, n = NumSamples;
		
		double[] y = samples1;
		double maximum, x;
		if (xmax0 <= xmin0) { xmin0 = xmin; xmax0 = xmax; }
		if (Sampled_getWindowSamples (xmin0, xmax0, imin, imax) == 0) 
		{
			/*
			 * No samples between xmin and xmax.
			 * Try to return the greater of the values at these two points.
			 */
			double yleft = Vector_getValueAtX (xmin0, channel,
				interpolation > 0 ? 1 : 0);
			double yright = Vector_getValueAtX (xmax0, channel,
				interpolation > 0 ? 1 : 0);
			maximum = yleft > yright ? yleft : yright;
			x = yleft == yright ? (xmin0 + xmax0) / 2 : yleft > yright ? xmin0 : xmax0;
		} 
		else 
		{
			maximum = y [(int) imin[0]];
			x = imin[0];
			if (y [(int) imax[0]] > maximum) maximum = y [(int) imax[0]];
			x = imax[0];
			if (imin[0] == 0) imin[0] ++;
			if (imax[0] == (NumSamples-1)) imax[0] --;
			for (i = imin[0]; i <= imax[0]; i ++) 
			{
				if (y [(int) i] > y [(int) (i - 1)] && y [(int) i] >= y [(int) (i + 1)]) 
				{
					double[] i_real = new double[1]; 
					double localMaximum = NUMimproveMaximum (y, n, i, interpolation, i_real);
					if (localMaximum > maximum) maximum = localMaximum;
					x = i_real[0];
				}
			}
			x = x1 + x * dx;   /* Convert sample to x. */
			if (x < xmin0) x = xmin0; else if (x > xmax0) x = xmax0;
		}
		if (return_maximum != null) return_maximum[0] = maximum;
		if (return_xOfMaximum != null) return_xOfMaximum[0] = x;
	}
	
	void Vector_getMinimumAndX (double xmin0, double xmax0, long channel, int interpolation,
			double[] return_minimum, double[] return_xOfMinimum)
	{
		long[] imin = new long[1], imax = new long[1];
		int n = NumSamples;
		double[] y = samples1;
		double minimum, x;
		if (xmax0 <= xmin0) { xmin0 = xmin; xmax0 = xmax; }
		if (Sampled_getWindowSamples (xmin0, xmax0, imin, imax) == 0) 
		{
			/*
			 * No samples between xmin and xmax.
			 * Try to return the lesser of the values at these two points.
			 */
			double yleft = Vector_getValueAtX (xmin0, channel,
				interpolation > 0 ? 1 : 0);
			double yright = Vector_getValueAtX (xmax0, channel,
				interpolation > 0 ? 1 : 0);
			minimum = yleft < yright ? yleft : yright;
			x = yleft == yright ? (xmin0 + xmax0) / 2 : yleft < yright ? xmin0 : xmax0;
		} 
		else 
		{
			minimum = y [(int) imin[0]];
			x = imin[0];
			if (y [(int) imax[0]] < minimum) minimum = y [(int) imax[0]];
			x = imax[0];
			if (imin[0] == 0) imin[0] ++;
			if (imax[0] == (NumSamples-1)) imax[0] --;
			for (long i = imin[0]; i <= imax[0]; i ++) 
			{
				if (y [(int) i] < 
						y [(int) (i - 1)] && 
						y [(int) i] <= 
						y [(int) (i + 1)]) 
				{
					double[] i_real = new double[1]; 
					double localMinimum = NUMimproveMinimum(y, n, i, interpolation, i_real);
					if (localMinimum < minimum) minimum = localMinimum;
					x = i_real[0];
				}
			}
			x = x1 + x * dx;   /* Convert sample to x. */
			if (x < xmin0) x = xmin0; else if (x > xmax0) x = xmax0;
		}
		if (return_minimum != null) return_minimum[0] = minimum;
		if (return_xOfMinimum != null) return_xOfMinimum[0] = x;
	}
	
	//Melder_assert (ilevel <= my ny);
	double Vector_getValueAtX (double x, long ilevel, int interpolation) 
	{
		double leftEdge = x1 - 0.5 * dx, rightEdge = leftEdge + NumSamples * dx;
		if (x <  leftEdge || x > rightEdge) return NUMundefined;
		double sum = 0.0;
		for (long channel = 0; channel < NumChannels; channel ++) 
		{
			sum += NUM_interpolate_sinc (samples1, NumSamples, Sampled_xToIndex(x),
				interpolation == 3 ? 70 :
				interpolation == 4 ? 700 :
				interpolation);
		}
		return sum / NumChannels;
	}
	
	// Corresponds to praat function "Override sampling frequency..."
	void overrideSamplingFrequency (double rate) 
	{
		dx = 1 / rate;
		x1 = xmin + 0.5 * dx;
		xmax = xmin + NumSamples * dx;
	}
	
	// Corresponds to praat function "Set part to zero..."
	void Sound_setZero (double tmin_in, double tmax_in, boolean roundTimesToNearestZeroCrossing) 
	{
		double[] tmin_in0 = new double[1];
		double[] tmax_in0 = new double[1];
		tmin_in0[0] = tmin_in;
		tmax_in0[0] = tmax_in;
		Function_unidirectionalAutowindow (tmin_in0, tmax_in0);
		Function_intersectRangeWithDomain (tmin_in0, tmax_in0);
		for (long channel = 1; channel <= NumChannels; channel ++) 
		{
			double tmin = tmin_in0[0], tmax = tmax_in0[0];
			if (roundTimesToNearestZeroCrossing) 
			{
				if (tmin > xmin) tmin = Sound_getNearestZeroCrossing (tmin_in, channel);
				if (tmax < xmax) tmax = Sound_getNearestZeroCrossing (tmax_in, channel);
			}
			if (tmin == NUMundefined) tmin = xmin;
			if (tmax == NUMundefined) tmax = xmax;
			long[] imin = new long[1];
			long[] imax = new long[1];
			Sampled_getWindowSamples (tmin, tmax, imin, imax);
			for (long i = imin[0]; i <= imax[0]; i ++) 
			{
				samples1[(int)i] = 0.0;
			}
		}
	}
	
	// Precondition: my z [1] [i1] != my z [1] [i1 + 1]; 
	double interpolate (long i1, long channel)
	{
		long i2 = i1 + 1;
		double x1 = Sampled_indexToX(i1), x2 = Sampled_indexToX(i2);
		double y1 = samples1[(int) i1], y2 = samples1[(int) i2];
		return x1 + (x2 - x1) * y1 / (y1 - y2);   /* Linear. */
	}
	
	double Sound_getNearestZeroCrossing (double position, long channel) 
	{
		double[] amplitude = samples1;
		long leftSample = Sampled_xToLowIndex(position);
		long rightSample = leftSample + 1, ileft, iright;
		double leftZero = 0, rightZero = 0;
		/* Are we already at a zero crossing? */
		if (leftSample >= 0 && rightSample < NumSamples &&
			(amplitude [(int) leftSample] >= 0.0) !=
			(amplitude [(int) rightSample] >= 0.0))
		{
			return interpolate(leftSample, channel);
		}
		/* Search to the left. */
		if (leftSample >= NumSamples) return NUMundefined;
		for (ileft = leftSample - 1; ileft >= 0; ileft --)
			if ((amplitude [(int) ileft] >= 0.0) != (amplitude [(int) (ileft + 1)] >= 0.0))
			{
				leftZero = interpolate (ileft, channel);
				break;
			}
		/* Search to the right. */
		if (rightSample < 0) return NUMundefined;
		for (iright = rightSample + 1; iright < NumSamples; iright ++)
			if ((amplitude [(int) iright] >= 0.0) != (amplitude [(int) (iright - 1)] >= 0.0))
			{
				rightZero = interpolate (iright - 1, channel);
				break;
			}
		if (ileft < 1 && iright > NumSamples) return NUMundefined;
		return ileft < 1 ? rightZero : iright > NumSamples ? leftZero :
			position - leftZero < rightZero - position ? leftZero : rightZero;
	}
	
	
	
	void Function_unidirectionalAutowindow (double[] xmin0, double[] xmax0) 
	{
		if (xmin0[0] >= xmax0[0]) 
		{
			xmin0[0] = xmin;
			xmax0[0] = xmax;
		}
	}
	
	boolean Function_intersectRangeWithDomain (double[] x10, double[] x2) 
	{
		if (x10[0] == x2[0]) return false;
		if (x10[0] < x2[0]) 
		{
			if (x10[0] < xmin) x10[0] = xmin;   // intersect requested range with logical domain
			if (x2[0] > xmax) x2[0] = xmax;
			if (x2[0] <= x10[0]) return false;   // requested range and logical domain do not intersect
		} 
		else 
		{
			if (x2[0] < xmin) x2[0] = xmin;   // intersect requested range with logical domain
			if (x10[0] > xmax) x10[0] = xmax;
			if (x10[0] <= x2[0]) return false;   // requested range and logical domain do not intersect
		}
		return true;
	}
	
	
	// Precondition: NumChannels must be equal or 1. Sampling frequencies must be equal.
	// Corresponds to praat function "Convolve"
	Wave Sounds_convolve(Wave thee) 
	{
		long n1 = NumSamples, n2 = thee.NumSamples;
		long n3 = n1 + n2 - 1, nfft = 1;
		while (nfft < n3) nfft *= 2;
		double[] data1 = new double[(int) nfft];
		double[] data2 = new double[(int) nfft];
		long numberOfChannels = NumChannels > thee.NumChannels ? NumChannels : thee.NumChannels;
		Wave him = new Wave((int)numberOfChannels, xmin + thee.xmin, xmax + thee.xmax, (int)n3, dx, x1 + thee.x1);
		for (long channel = 1; channel <= numberOfChannels; channel ++) 
		{
			double[] a = samples1;
			for (long i = (n1-1); i >= 0; i--) data1 [(int) i] = a [(int) i];
			for (long i = n1; i < nfft; i++) data1 [(int) i] = 0.0;
			a = thee.samples1;
			for (long i = (n2-1); i >= 0; i--) data2 [(int) i] = a [(int) i];
			for (long i = n2; i < nfft; i++) data2 [(int) i] = 0.0;
			RealDoubleFFT fft1 = new RealDoubleFFT((int) nfft);
			fft1.ft(data1);
			RealDoubleFFT fft2 = new RealDoubleFFT((int) nfft);
			fft1.ft(data2);
			data2 [0] *= data1 [0];
			data2 [1] *= data1 [1];
			double temp;
			for (long i = 2; i < nfft; i += 2) 
			{
				temp = data1[(int) i] * data2[(int) i] - data1[(int) (i + 1)] * data2[(int) (i + 1)];
				data2[(int) (i + 1)] = data1[(int) i] * data2[(int) (i + 1)] + data1[(int) (i + 1)] * data2[(int) i];
				data2[(int) i] = temp;
			}
			RealDoubleFFT bfft = new RealDoubleFFT((int) nfft);
			fft1.bt(data2);
			a = him.samples1;
			for (long i = 0; i < n3; i++) 
			{
				a [(int) i] = data2[(int) i];
			}
		}
		him.multiplyByScalar (1.0 / (double)nfft);
		return him;
		
	}
	
	/* Preconditions: minimumPitch & timStep defined
 	 * timeStep < 0.0 and dx <= 0.0 and minimumPitch <= 0.0 and windowDuration > 0.0
 	 *  Corresponds to praat function "To Intensity..."
 	 */ 
 	Intensity Sound_to_Intensity_ (double minimumPitch, double timeStep, int subtractMeanPressure)
 	{
		if (timeStep == 0.0) timeStep = 0.8 / minimumPitch;   // default: four times oversampling Hanning-wise
		double windowDuration = 6.4 / minimumPitch;
		double halfWindowDuration = 0.5 * windowDuration;
		long halfWindowSamples = (long) (halfWindowDuration / dx);
		double[] amplitude = new double[(int) (halfWindowSamples*2+1)];
		double[] window = new double[(int) (halfWindowSamples*2+1)];
		for (long i = 0; i <= halfWindowSamples*2; i++) 
		{
			double x = (i-halfWindowSamples) * dx / halfWindowDuration, root = 1 - x * x;
			window[(int) i] = root <= 0.0 ? 0.0 : NUMbessel_i0_f ((2 * NUMpi * NUMpi + 0.5) * Math.sqrt(root));
		}
		long[] numberOfFrames = new long[1]; 
		double[] thyFirstTime = new double[1]; 
		double myDuration = dx * NumSamples;
		numberOfFrames[0] = (long) (Math.floor((myDuration - windowDuration) / timeStep) + 1);
		double ourMidTime = x1 - 0.5 * dx + 0.5 * myDuration;
		double thyDuration = numberOfFrames[0] * timeStep;
		thyFirstTime[0] = ourMidTime - 0.5 * thyDuration + 0.5 * timeStep;
		Sampled_shortTermAnalysis (windowDuration, timeStep, numberOfFrames, thyFirstTime);
		Intensity thee = new Intensity(xmin, xmax, (int)numberOfFrames[0], timeStep, thyFirstTime[0]);
		for (long iframe = 0; iframe < numberOfFrames[0]; iframe++) 
		{
			double midTime = thee.Sampled_indexToX(iframe);
			long midSample = Sampled_xToNearestIndex(midTime);
			long leftSample = midSample - halfWindowSamples, rightSample = midSample + halfWindowSamples;
			double sumxw = 0.0, sumw = 0.0, intensity;
			if (leftSample < 0) leftSample = 0;
			if (rightSample >= NumSamples) rightSample = (NumSamples-1);

			for (long channel = 1; channel <= NumChannels; channel ++) 
			{
				for (long i = leftSample; i <= rightSample; i++) 
				{
					amplitude[(int) (i - midSample + halfWindowSamples)] = samples1[(int) i];
				}
				if (subtractMeanPressure == 1) 
				{
					double sum = 0.0;
					for (long i = leftSample; i <= rightSample; i ++) 
					{
						sum += amplitude[(int) (i - midSample + halfWindowSamples)];
					}
					double mean = sum / (rightSample - leftSample + 1);
					for (long i = leftSample; i <= rightSample; i ++) 
					{
						amplitude[(int) (i - midSample + halfWindowSamples)] -= mean;
					}
				}
				for (long i = leftSample; i <= rightSample; i++) 
				{
					sumxw += amplitude[(int) (i - midSample + halfWindowSamples)] * amplitude[(int) (i - midSample + halfWindowSamples)] * window[(int) (i - midSample + halfWindowSamples)];
					sumw += window[(int) (i - midSample + halfWindowSamples)];
				}
			}
			intensity = sumxw / sumw;
			
			if (intensity != 0.0) 
				intensity /= 4e-10;
			thee.samples1[(int) iframe] = intensity < 1e-30 ? -300 : 10 * Math.log10 (intensity);
		}
		return thee;
 		
 	}

	/* Modified Bessel function I0. Abramowicz & Stegun, p. 378.*/
	double NUMbessel_i0_f (double x) 
	{
		if (x < 0.0) return NUMbessel_i0_f (- x);
		if (x < 3.75) 
		{
			/* Formula 9.8.1. Accuracy 1.6e-7. */
			double t = x / 3.75;
			t *= t;
			return 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492
				+ t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
		}
		/*
		otherwise: x >= 3.75
	 	*/
		/* Formula 9.8.2. Accuracy of the polynomial factor 1.9e-7. */
		double t = 3.75 / x;   /* <= 1.0 */
		return Math.exp(x) / Math.sqrt(x) * (0.39894228 + t * (0.01328592
				+ t * (0.00225319 + t * (-0.00157565 + t * (0.00916281
				+ t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633
				+ t * 0.00392377))))))));
	}
	
	// Corresponds to praat function "Sound & LPC: Filter..."
	Wave LPC_and_Sound_filterWithFilterAtTime (LPC me, int channel, double time)
	{
		Wave him = Wave_copy();
		LPC_and_Sound_filterWithFilterAtTime_inline (me, him, channel, time);
		return him;
	}
	
	void LPC_and_Sound_filterWithFilterAtTime_inline (LPC me, Wave thee, int channel, double time) 
	{
		long frameIndex = me.Sampled_xToNearestIndex(time);
		if (frameIndex < 1) 
		{
			frameIndex = 1;
		}
		if (frameIndex > me.nx) 
		{
			frameIndex = (long) me.nx;
		}
		if (channel > NumChannels) 
		{
			channel = 1;
		}
		if (frameIndex < 1 || frameIndex > me.nx) {} // should return false
		if (channel > 0) 
		{
			LPC_Frame_and_Sound_filter (me.frames[(int) frameIndex], thee, channel);
		} 
		else 
		{
			for (long ichan = 1; ichan <= NumChannels; ichan++) 
			{
				LPC_Frame_and_Sound_filter (me.frames[(int) frameIndex], thee, (int)ichan);
			}
		}
	}
	
	static void LPC_Frame_and_Sound_filter (LPC_Frame me, Wave thee, int channel)
	{
		double[] y = thee.samples1, a = me.a;
		for (long i = 1; i <= thee.NumSamples; i++) 
		{
			long m = i > me.nCoefficients ? me.nCoefficients : i - 1;
			for (long j = 1; j <= m; j++) 
			{
				y[(int) i] -= a[(int) j] * y[(int) (i - j)];
			}
		}
	}
	

	
	void Sound_preEmphasis (double preEmphasisFrequency) 
	{
		double preEmphasis = Math.exp(-2.0 * NUMpi * preEmphasisFrequency * dx);
		for (long channel = 1; channel <= NumChannels; channel ++) 
		{
			double[] s = samples1; 
			for (long i = (NumSamples-1); i >= 1; i --) s [(int)i] -= preEmphasis * s [(int) (i - 1)];
		}
	}
	
	// check this
	void Sound_into_Sound (Wave to, double startTime) 
	{
		long index = Sampled_xToNearestIndex(startTime);
		for (long i = 0; i < to.NumSamples; i++) 
		{
			long j = index + i;
			to.samples1[(int) i] = j < 0 || j >= NumSamples ? 0 : samples1[(int) j];
		}
	}
	
	void Sounds_multiply (Wave thee) 
	{
		long n = NumSamples < thee.NumSamples ? NumSamples : thee.NumSamples;
		double[] s1 = samples1, s2 = thee.samples1;
		for (long i = 0; i < n; i++) 
		{
			s1[(int) i] *= s2[(int) i];
		}
	}
	
	int Sound_into_LPC_Frame_burg (LPC_Frame thee) 
	{
		int status = NUMburg1 (samples1, NumSamples, thee.a, thee.nCoefficients, thee.gain);
		thee.gain[0] *= NumSamples;
		for (long i = 0; i < thee.nCoefficients; i++) 
		{
			thee.a[(int)i] = -thee.a[(int)i];
		}
		return status;
	}
	
	int NUMburg1 (double x[], long n, double a[], int m, double[] xms)  
	{
		for (long j = 1; j <= m; j++) 
		{
			int index = (int) (j - 1);
			a[index] = 0;
		}

		/*autoNUMvector<double> b1 (1, n);
		autoNUMvector<double> b2 (1, n);
		autoNUMvector<double> aa (1, m);*/
		
		double[] b1 = new double[(int) n];
		double[] b2 = new double[(int) n];
		double[] aa = new double[m];
		

		// (3)

		double p = 0.0;
		for (long j = 1; j <= n; j++) 
		{
			int index = (int) (j - 1);
			p += x[index] * x[index];
		}

		//*xms = p / n;
		xms[0] = p / n;
		if (xms[0] <= 0) 
		{
			return 0;    // warning empty
		}

		// (9)

		b1[0] = x[0];
		b2[(int) (n - 2)] = x[(int) (n - 1)];
		for (long j = 2; j <= n - 1; j++) 
		{
			int index = (int) (j - 1);
			b1[index] = b2[index - 1] = x[index];
		}

		for (long i = 1; i <= m; i++) 
		{
			int index = (int) (i - 1);
			// (7)

			double num = 0.0, denum = 0.0;
			for (long j = 1; j <= n - i; j++) 
			{
				int index1 = (int) (j - 1);
				num += b1[index1] * b2[index1];
				denum += b1[index1] * b1[index1] + b2[index1] * b2[index1];
			}

			if (denum <= 0) 
			{
				return 0;    // warning ill-conditioned
			}

			a[index] = 2.0 * num / denum;

			// (10)

			xms[0] *= 1.0 - a[index] * a[index];

			// (5)

			for (long j = 1; j <= i - 1; j++) 
			{
				int index1 = (int) (j - 1);
				a[index1] = aa[index1] - a[index] * aa[(index - index1) - 1];
			}

			if (index < (m-1)) 
			{

				// (8)  Watch out: i -> i+1

				for (long j = 1; j <= i; j++) 
				{
					int index1 = (int) (j - 1);
					aa[index1] = a[index1];
				}
				for (long j = 1; j <= n - i - 1; j++) 
				{
					int index1 = (int) (j - 1);
					b1[index1] -= aa[index] * b2[index1];
					b2[index1] = b2[index1 + 1] - aa[index] * b1[index1 + 1];
				}
			}
		}
		return 1;
	}
	
	Spectrum Sound_to_Spectrum (boolean fast) 
	{
		long numberOfSamples = NumSamples;
		if (fast) 
		{
			numberOfSamples = 2;
			while (numberOfSamples < NumSamples) numberOfSamples *= 2;
		}
		long numberOfFrequencies = numberOfSamples / 2 + 1;   // 4 samples -> cos0 cos1 sin1 cos2; 5 samples -> cos0 cos1 sin1 cos2 sin2
		double[] data = new double[(int) numberOfSamples];
		for (long i = 0; i < NumSamples; i ++)
			data [(int) i] = NumChannels == 1 ? samples1[(int) i] : 0.5 * (samples1[(int) i] /*+ samples1[(int) i]*/); // commented part is for two channels
		RealDoubleFFT fft = new RealDoubleFFT((int) numberOfSamples);
		fft.ft(data);	
		Spectrum thee = new Spectrum(0.5 / dx, (int)numberOfFrequencies);
		thee.originalNumberOfSamples = NumSamples; 
		thee.dx = 1.0 / (dx * numberOfSamples);   // override
		double[] re = thee.rl;
		double[] im = thee.im;
		double scaling = dx;
		re [0] = data [0] * scaling;
		im [0] = 0.0;
		for (long i = 1; i < (numberOfFrequencies-1); i ++) 
		{
			re [(int) i] = data [(int) (i + i - 1)] * scaling;
			im [(int) i] = data [(int) (i + i)] * scaling;
		}
		if ((numberOfSamples & 1) != 0) // when fast, this is always false
		{
			if (numberOfSamples > 1) 
			{
				re [(int) (numberOfFrequencies-1)] = data [(int) (numberOfSamples - 2)] * scaling;
				im [(int) (numberOfFrequencies-1)] = data [(int) (numberOfSamples - 1)] * scaling;
			}
		} 
		else 
		{
			re [(int) (numberOfFrequencies-1)] = data [(int) (numberOfSamples-1)] * scaling;
			im [(int) (numberOfFrequencies-1)] = 0.0;
		}
		return thee;
	}
	
	public Wave raspiness(double raspiness)
	{
		if(raspiness > 100)
			raspiness = 100;
		if(raspiness < 1)
			raspiness = 1;
		double value = raspiness / 100000.0;
		Wave wrk = Wave_copy();
		Pitch pitch = wrk.Sound_to_Pitch(0, 40, 600);
		Manipulation manipulation = wrk.Sound_Pitch_to_Manipulation(pitch);
		DurationTier durationtier = new DurationTier(0, xmax);
		durationtier.addPoint(0, 1);
		manipulation.duration = durationtier;
		PointProcess pointprocess = pitch.Pitch_to_PointProcess();
		Matrix matrix = pointprocess.PointProcess_to_Matrix();
		matrix.addSelfRandom(value);
		PointProcess pointprocess2 = matrix.Matrix_to_PointProcess();
		manipulation.pulses = pointprocess2;
		Wave resynthesis = manipulation.synthesize_overlapAdd();
		return resynthesis;
	}
    
    public Wave pitch(double newPitch)
	{
    	if(newPitch > 100)
    		newPitch = 100;
		if(newPitch < 1)
			newPitch = 1;
		newPitch *= 2;
		newPitch += 50;
		Wave wrk = Wave_copy();
        int[] minF0 = new int[1];
        int[] maxF0 = new int[1];
        minMaxF(minF0, maxF0);
        wrk = wrk.Sound_changeGender_old(minF0[0], maxF0[0], 1, newPitch, 1, 1);
        return wrk;
	}
    
    void minMaxF(int[] minF0, int[] maxF0)
    {
        Pitch pitch = Sound_to_Pitch(0.0, 40.0, 600.0);
        long voicedframes = pitch.Pitch_countVoicedFrames();
        if(voicedframes > 0)
        {
            double q25 = pitch.Pitch_getQuantile(0, 0, 0.25, 0);
            double q75 = pitch.Pitch_getQuantile(0, 0, 0.75, 0);
            minF0[0] = (int) Math.round(q25 * 0.75);
            maxF0[0] = (int) Math.round(q75 * 1.5);
        }
        else
        {
            minF0[0] = 40;
            maxF0[0] = 60;
        }
    }
    
    public Wave resonance(double value)
    {
    	Wave wrk = Wave_copy();
    	if(value > 100.0)
    		value = 100.0;
		if(value < 1.0)
			value = 1.0;
		value /= 200.0;
		value += 0.75;
    	wrk = wrk.Sound_changeGender_old(75, 600, value, 0, 1, 1);
    	return wrk;
    }

	public Wave breathiness(double breathiness)
	{
		double gain1 = 0, gain2 = 0;
		if (breathiness > 100)
			breathiness = 100;
		if (breathiness <= 50)
		{
			gain1 = 1;
			gain2 = (breathiness/50.0);
		}
		if (breathiness > 50)
		{
			gain1 = ((100.0-breathiness)/50.0);
			gain2 = 1;
		}
		Wave wrk = Wave_copy();
		wrk = wrk.fixdc();
		Wave whisper = wrk.whisper();
		Wave result = new Wave(xmin, xmax, dx, 1);
		result.addWaveMultiplyScalar(wrk, gain1);
		result.addWaveMultiplyScalar(whisper, gain2);
		result = result.fixdc();
		return result;
	}

	Wave fixdc()
	{
		subtractMean();
		Wave filtered = filterPassHannBand(60, 0, 20);
		filtered.declip();
		return filtered;
	}

	void declip()
	{
		double clip = getAbsoluteExtrenum();
		DecimalFormat df = new DecimalFormat("#.####");
		df.setRoundingMode(RoundingMode.HALF_UP);
		clip = Double.parseDouble(df.format(clip));
		if(clip >= 1)
			scalePeak(0.9999);
	}

	Wave workpre()
	{
		Wave work = Wave_copy();
		work = work.Sound_extractPart(xmin-0.025, xmin+(xmax - xmin)+0.025, 1, false);
		work = work.fixdc();
		return work;
	}

	Wave whisper()
	{
		Wave tmp = workpre();
		tmp.scalePeak(0.9999);
		Wave tmp2 = tmp.gate(20);
		Wave bgnoise = new Wave(xmin, xmax, dx, -1);
		bgnoise.scaleIntensity(0.01);
		tmp2.addWaveMultiplyScalar(bgnoise, 1);
		int pred_order = (int) (Math.round((1.0/dx)/1000.0) +2);
		LPC lpc = new LPC();
		lpc = lpc.Sound_to_LPC_burg(tmp2, pred_order, 0.025, 0.01, 50.0);
		Wave noise = new Wave(xmin, xmax, dx, -1);
		Wave tmp3 = noise.LPC_and_Sound_filter(lpc, 1);
		Wave tmp4 = tmp3.eq10bands(-24, -24, -24, -24, 12, 24, 24, 12, 12, -6);
		tmp4 = tmp4.fixdc();
		//tmp4.workpost();
		return tmp4;
	}
	
	Wave gate(int threshold)
	{
		if(threshold > 100) threshold = 100;
		Wave result = Wave_copy();
		Intensity intensity = result.Sound_to_Intensity_ (100, 0.01, 0);
		double i = 0, zerostart = 0, zeroend = 0, zeroselection = 0;
		int index;
		double dB;
		do
		{
			index = (int) intensity.Sampled_xToNearestIndex(i);
			if(index < 0)
				index = 0;
			if(index >= intensity.NumFrames)
				index = intensity.NumFrames-1;
			dB = intensity.samples1[(int) index];
			if(dB < threshold)
			{
				if(zeroselection == 0)
				{
					zerostart = i;
					zeroselection = 1;
				}
			}
			else
			{
				if(zeroselection == 1)
				{
					zeroend = i;
					if(zerostart != zeroend)
						result.Sound_setZero(zerostart, zeroend, true);
					zeroselection = 0;
				}
			}
			i += 0.01;
		}
		while(i <= xmax);
		if(zeroselection == 1)
		{
			zeroend = i;
			if(zerostart != zeroend)
				result.Sound_setZero(zerostart, zeroend, true);
			zeroselection = 0;
		}
		return result;
	}
	
	Wave eq10bands(int b1, int b2, int b3, int b4, int b5, int b6, int b7, int b8, int b9, int b10)
	{
		Wave wrk = Wave_copy();
		//wrk.fixdc();
		double hifreq = (1.0/dx)/2;
		PointProcess pointprocess = new PointProcess(0, 0.05);
		pointprocess.addPoint(0.025);
		Wave pulse = pointprocess.toSoundPulseTrain((1.0/dx), 1, 0.05, 2000);
		Spectrum sp_pulse = pulse.Sound_to_Spectrum(false);
		Spectrum buffer = sp_pulse.Spectrum_copy();
		buffer.zeroEverything();
		Spectrum sp_eq = buffer.Spectrum_copy();
		sp_eq.EqBand(0, 44.2, b1, 20, buffer, sp_pulse);
		sp_eq.EqBand(44.2, 88.4, b2, 20, buffer, sp_pulse);
		sp_eq.EqBand(88.4, 177, b3, 40, buffer, sp_pulse);
		sp_eq.EqBand(177, 354, b4, 80, buffer, sp_pulse);
		sp_eq.EqBand(354, 707, b5, 100, buffer, sp_pulse);
		sp_eq.EqBand(707, 1414, b6, 100, buffer, sp_pulse);
		sp_eq.EqBand(1414, 2828, b7, 100, buffer, sp_pulse);
		sp_eq.EqBand(2828, 5657, b8, 100, buffer, sp_pulse);
		sp_eq.EqBand(5657, 11314, b9, 100, buffer, sp_pulse);
		sp_eq.EqBand(11314, hifreq, b10, 100, buffer, sp_pulse);
		sp_eq.filterPassHannBand(80, 0, 20);
		sp_eq.filterPassHannBand(0, 20000, 100);
		Wave pulseeqtmp = sp_eq.Spectrum_to_Sound();
	    //extract part nonsense
	    if(pulseeqtmp.dx != dx)
	    {
	    	pulseeqtmp.overrideSamplingFrequency(dx);
	    }

	    Wave resulttmp = wrk.Sounds_convolve(pulseeqtmp);
	    //more extract part nonsense
	    resulttmp.scalePeak(0.9999);
	    return resulttmp;
	}
	
	void addWaveMultiplyScalar(Wave wave, double scalar)
	{
        int numberOfSamples = Math.min(NumSamples, wave.NumSamples);
		for(int i = 0; i < numberOfSamples; i++)
		{
			samples1[i] = samples1[i] + scalar * wave.samples1[i];
		}
	}
	
	// LPC smapling period and Wave dx should match
	Wave LPC_and_Sound_filter(LPC me, int useGain) 
	{
		double xmin0 = me.xmin > xmin ? me.xmin : xmin;
		double xmax0 = me.xmax < xmax ? me.xmax : xmax;
		Wave him = Wave_copy();
		double[] x = him.samples1;
		long ifirst = Sampled_xToHighIndex(xmin0) + 1;
		long ilast = Sampled_xToLowIndex(xmax0) + 1;
		
		for (long i = ifirst; i <= ilast; i++) 
		{
			double t = him.x1 + (i-1) * him.dx; /* Sampled_indexToX (him, i) */
			long iFrame = (long) Math.floor( (t - me.x1) / me.dx + 1.5); /* Sampled_xToNearestIndex (me, t) */
			if (iFrame < 1) 
			{
				continue;
			}
			if (iFrame > me.nx) 
			{
				break;
			}
			
			double[] a = me.frames[(int) (iFrame-1)].a;
			long m = i > me.frames[(int)(iFrame-1)].nCoefficients ? me.frames[(int) (iFrame-1)].nCoefficients : i-1;
			for (long j = 1; j <= m; j++) 
			{
				
				x[(int) (i-1)] -= a[(int) (j-1)] * x[(int) ((i - j) - 1)];
			}
		}
		ifirst--;
		ilast--;
		// Make samples before first frame and after last frame zero.

		for (long i = 0; i < ifirst; i++) 
		{
			x[(int) i] = 0;
		}

		for (long i = ilast + 1; i < him.NumSamples; i++) 
		{
			x[(int) i] = 0;
		}
		if (useGain == 1) 
		{
			for (long i = ifirst; i <= ilast; i++) 
			{
				double t = him.x1 + i * him.dx; /* Sampled_indexToX (him, i) */
				double riFrame = (t - me.x1) / me.dx; /* Sampled_xToIndex (me, t); */
				long iFrame = (long) Math.floor(riFrame);
				double phase = riFrame - iFrame;
				if (iFrame < 0 || iFrame >= me.nx)
				{
					x[(int) i] = 0;
				} 
				else if (iFrame == 0) 
				{
					x[(int) i] *= Math.sqrt(me.frames[1].gain[0]) * phase;
				} 
				else if (iFrame == me.nx-1) 
				{
					x[(int) i] *= Math.sqrt(me.frames[me.nx-1].gain[0]) * (1 - phase);
				} 
				else 
				{
					x[(int) i] *= phase * Math.sqrt(me.frames[(int) (iFrame + 1)].gain[0]) + (1 - phase) * Math.sqrt(me.frames[(int) iFrame].gain[0]);
				}
			}
		}
		return him;
	}
	
	// Praat helper function
	long Sampled_xToHighIndex(double x) 
	{ 
		return (long) Math.ceil((x - x1) / dx /*+ 1.0*/); 
	}
	
	// Praat helper function
	long Sampled_xToLowIndex(double x) 
	{ 
		return (long)Math.floor((x - x1) / dx /*+ 1.0*/); 
	}
	
	// Praat helper function
	long Sampled_xToNearestIndex (double x) 
	{ 
		return (long) Math.round((x - x1) / dx /*+ 1.0*/); 
	}
	
	// Praat helper function
	double Sampled_indexToX (long index) 
	{ 
		return x1 + index * dx;
	}
	
	double Sampled_xToIndex (double x) 
	{ 
		return (x - x1) / dx /*+ 1.0*/; 
	}
	
	Pitch Sound_to_Pitch (double timeStep, double minimumPitch, double maximumPitch) 
	{
		return Sound_to_Pitch_ac (timeStep, minimumPitch, 3.0, 15, 0, 0.03, 0.45, 0.01, 0.35, 0.14, maximumPitch);
	}
	
	Pitch Sound_to_Pitch_ac (double dt, double minimumPitch, double periodsPerWindow, int maxnCandidates, int accurate, double silenceThreshold, double voicingThreshold, double octaveCost, double octaveJumpCost, double voicedUnvoicedCost, double ceiling)
	{
		return Sound_to_Pitch_any(dt, minimumPitch, periodsPerWindow, maxnCandidates, accurate, silenceThreshold, voicingThreshold, octaveCost, octaveJumpCost, voicedUnvoicedCost, ceiling);
	}

	// maxnCandidates must be >= 2
	Pitch Sound_to_Pitch_any (double dt, double minimumPitch, double periodsPerWindow, int maxnCandidates, int method, double silenceThreshold, double voicingThreshold, double octaveCost, double octaveJumpCost, double voicedUnvoicedCost, double ceiling)
	{
		RealDoubleFFT fft; // originally fftTable
		double duration;
		double[] t1 = new double[1];
		double dt_window;   /* Window length in seconds. */
		long nsamp_window, halfnsamp_window;   /* Number of samples per window. */
		long[] nFrames = new long[1]; 
		long minimumLag, maximumLag;
		long nsampFFT;
		double interpolation_depth;
		long nsamp_period, halfnsamp_period;   /* Number of samples in longest period. */
		long brent_ixmax, brent_depth;
		double globalPeak;

		if (maxnCandidates < ceiling / minimumPitch) maxnCandidates = (int) (ceiling / minimumPitch);
		if (dt <= 0.0) dt = periodsPerWindow / minimumPitch / 4.0;   /* e.g. 3 periods, 75 Hz: 10 milliseconds. */

		brent_depth = 3;
		interpolation_depth = 0.5;
		duration = dx * NumSamples;
		// At this point, minimum pitch must not be less than (periodsPerWindow / duration).

		/*
		 * Determine the number of samples in the longest period.
		 * We need this to compute the local mean of the sound (looking one period in both directions),
		 * and to compute the local peak of the sound (looking half a period in both directions).
		 */
		nsamp_period = (long) Math.floor(1 / dx / minimumPitch);
		halfnsamp_period = nsamp_period / 2 + 1;
		if (ceiling > 0.5 / dx) ceiling = 0.5 / dx;

		// Determine window length in seconds and in samples.
		dt_window = periodsPerWindow / minimumPitch;
		nsamp_window = (long) Math.floor(dt_window / dx);
		halfnsamp_window = nsamp_window / 2 - 1;
		// At this point, halfnsamp must not be less than 2.
		nsamp_window = halfnsamp_window * 2;

		// Determine the minimum and maximum lags.
		minimumLag = (long) Math.floor(1 / dx / ceiling);
		if (minimumLag < 2) minimumLag = 2;
		maximumLag = (long) (Math.floor(nsamp_window / periodsPerWindow) + 2);
		if (maximumLag > nsamp_window) maximumLag = nsamp_window;

		/*
		 * Determine the number of frames.
		 * Fit as many frames as possible symmetrically in the total duration.
		 * We do this even for the forward cross-correlation method,
		 * because that allows us to compare the two methods.
		 */
		Sampled_shortTermAnalysis (method >= 2 ? 1 / minimumPitch + dt_window : dt_window, dt, nFrames, t1);
		
		// Create the resulting pitch contour.
		Pitch thee = new Pitch(xmin, xmax, (int) nFrames[0], dt, t1[0], ceiling, maxnCandidates);

		// Create (too much) space for candidates.
		for (long iframe = 0; iframe < nFrames[0]; iframe ++) 
		{
			thee.frame[(int) iframe] = new PitchFrame(maxnCandidates);
		}

		// Compute the global absolute peak for determination of silence threshold.
		globalPeak = 0.0;
		for (long channel = 1; channel <= NumChannels; channel ++) 
		{
			double mean = 0.0;
			for (long i = 0; i < NumSamples; i ++) 
			{
				mean += samples1[(int) i];
			}
			mean /= NumSamples;
			for (long i = 0; i < NumSamples; i ++) 
			{
				double value = Math.abs(samples1[(int) i] - mean);
				if (value > globalPeak) globalPeak = value;
			}
		}
		if (globalPeak == 0.0) 
		{
			return thee;
		}

		double[] window;
		double[] windowR;
		/* For autocorrelation analysis. */

		/*
		* Compute the number of samples needed for doing FFT.
		* To avoid edge effects, we have to append zeroes to the window.
		* The maximum lag considered for maxima is maximumLag.
		* The maximum lag used in interpolation is nsamp_window * interpolation_depth.
		*/
		nsampFFT = 1; while (nsampFFT < nsamp_window * (1 + interpolation_depth)) nsampFFT *= 2;

		// Create buffers for autocorrelation analysis.
		windowR = new double[(int) nsampFFT];
		window = new double[(int) nsamp_window];
		fft = new RealDoubleFFT((int) nsampFFT);

		// Hanning window is applied against phase effects.
		// The Hanning window is 2 to 5 dB better for 3 periods/window.		
		for (long i = 0; i < nsamp_window; i ++)
			window [(int) i] = 0.5 - 0.5 * Math.cos((i+1) * 2 * NUMpi / (nsamp_window + 1));
		
		// Compute the normalized autocorrelation of the window.
		for (long i = 0; i < nsamp_window; i ++) windowR [(int) i] = window [(int) i];
		fft.ft(windowR);
		windowR [0] *= windowR [0];   // DC component
		for (long i = 1; i < (nsampFFT-1); i += 2) 
		{
			windowR [(int) i] = windowR [(int) i] * windowR [(int) i] + windowR [(int) (i+1)] * windowR [(int) (i+1)];
			windowR [(int) (i + 1)] = 0.0;   // power spectrum: square and zero
		}
		windowR [(int) (nsampFFT-1)] *= windowR [(int) (nsampFFT-1)];   // Nyquist frequency
		fft.bt(windowR); // autocorrelation
		for (long i = 1; i < nsamp_window; i ++) windowR [(int) i] /= windowR [0];   // normalize
		windowR [0] = 1.0;   // normalize

		brent_ixmax = (long) (nsamp_window * interpolation_depth);

		long numberOfFramesPerThread = 20;
		int numberOfThreads = (int) ((nFrames[0] - 1) / numberOfFramesPerThread + 1);
		//const int numberOfProcessors = MelderThread_getNumberOfProcessors ();  // perhaps we can make this work later
		final int numberOfProcessors = 1;
		//trace ("%d processors", numberOfProcessors);
		if (numberOfThreads > numberOfProcessors) numberOfThreads = numberOfProcessors;
		if (numberOfThreads > 16) numberOfThreads = 16;
		if (numberOfThreads < 1) numberOfThreads = 1;
		numberOfFramesPerThread = (nFrames[0] - 1) / numberOfThreads + 1;

		//if (! mutex_inited) { MelderThread_MUTEX_INIT (mutex); mutex_inited = true; }
		SoundIntoPitchArgs[] args = new SoundIntoPitchArgs [16];
		long firstFrame = 1, lastFrame = numberOfFramesPerThread;
		int[] cancelled = {0};
		for (int ithread = 1; ithread <= numberOfThreads; ithread ++) 
		{
			if (ithread == numberOfThreads) lastFrame = nFrames[0];
			args[ithread - 1] = new SoundIntoPitchArgs (this, thee,
				firstFrame, lastFrame, minimumPitch, maxnCandidates, method,
				voicingThreshold, octaveCost,
				dt_window, nsamp_window, halfnsamp_window, maximumLag,
				nsampFFT, nsamp_period, halfnsamp_period, brent_ixmax, brent_depth,
				globalPeak, window, windowR,
				ithread == numberOfThreads, cancelled);
			firstFrame = lastFrame + 1;
			lastFrame += numberOfFramesPerThread;
		}
		//MelderThread_run (Sound_into_Pitch, args, numberOfThreads);
		args[0].Sound_into_Pitch();

		thee.Pitch_pathFinder (silenceThreshold, voicingThreshold,
			octaveCost, octaveJumpCost, voicedUnvoicedCost, ceiling, /*Melder_debug == 31 ? true :*/ false);

		
		return thee;
	} 
	

	// windowduration must be greater than 0
	//timestep must be greater than 0
	// windowduration must not be greater than myduration
	// numberofframes must be greater than or equal to 1
	void Sampled_shortTermAnalysis (double windowDuration, double timeStep, long[] numberOfFrames, double[] firstTime) 
	{
		double myDuration = dx * NumSamples;
		numberOfFrames[0] = (long) (Math.floor((myDuration - windowDuration) / timeStep) + 1);
		double ourMidTime = x1 - 0.5 * dx + 0.5 * myDuration;
		double thyDuration = numberOfFrames[0] * timeStep;
		firstTime[0] = ourMidTime - 0.5 * thyDuration + 0.5 * timeStep;
	}
	
	//Melder_assert (startSample >= 1);
	//Melder_assert (endSample <= my nx);
	//Melder_assert (startSample >= 1);
	//Melder_assert (endSample <= my nx);
	void Sound_into_PitchFrame (PitchFrame pitchFrame, double t,
			double minimumPitch, int maxnCandidates, int method, double voicingThreshold, double octaveCost,
			RealDoubleFFT fft, double dt_window, long nsamp_window, long halfnsamp_window,
			long maximumLag, long nsampFFT, long nsamp_period, long halfnsamp_period,
			long brent_ixmax, long brent_depth, double globalPeak,
			double[][] frame, double[] ac, double[] window, double[] windowR,
			double[] r, long[] imax, double[] localMean, long iframe)
	{
		double localPeak;
		long leftSample = Sampled_xToLowIndex(t), rightSample = leftSample + 1;
		long startSample, endSample;

		for (long channel = 0; channel < NumChannels; channel ++) 
		{
			/*
			 * Compute the local mean; look one longest period to both sides.
			 */
			startSample = rightSample - nsamp_period;
			endSample = leftSample + nsamp_period;
			localMean [(int) channel] = 0.0;
			for (long i = startSample; i <= endSample; i ++) 
			{
				localMean [(int) channel] += samples1[(int) i];
			}
			localMean [(int) channel] /= 2 * nsamp_period;

			/*
			 * Copy a window to a frame and subtract the local mean.
			 * We are going to kill the DC component before windowing.
			 */
			startSample = rightSample - halfnsamp_window;
			endSample = leftSample + halfnsamp_window;
			
			for (long j = 0, i = startSample; j < nsamp_window; j ++)
				frame [(int) channel] [(int) j] = (samples1[(int) i ++] - localMean [(int) channel]) * window [(int) j];
			for (long j = nsamp_window; j < nsampFFT; j ++)
				frame [(int) channel] [(int) j] = 0.0;
		}

		/*
		 * Compute the local peak; look half a longest period to both sides.
		 */
		localPeak = 0.0;
		if ((startSample = halfnsamp_window - halfnsamp_period) < 0) startSample = 0;
		if ((endSample = halfnsamp_window + halfnsamp_period - 1) > (nsamp_window - 1)) endSample = (nsamp_window - 1);
		for (long channel = 0; channel < NumChannels; channel ++) 
		{
			for (long j = startSample; j <= endSample; j ++) 
			{
				double value = Math.abs(frame [(int) channel] [(int) j]);
				if (value > localPeak) localPeak = value;
			}
		}
		pitchFrame.intensity = localPeak > globalPeak ? 1.0 : localPeak / globalPeak;

		/*
		 * The FFT of the autocorrelation is the power spectrum.
		 */
		for (long i = 0; i < nsampFFT; i ++) 
		{
			ac [(int) i] = 0.0;
		}
		for (long channel = 0; channel < NumChannels; channel ++) 
		{
			fft.ft(frame [(int) channel]);   /* Complex spectrum. */
			ac [0] += frame [(int) channel] [0] * frame [(int) channel] [0];   /* DC component. */
			for (long i = 1; i < nsampFFT-1; i += 2) 
			{
				ac [(int) i] += frame [(int) channel] [(int) i] * frame [(int) channel] [(int) i] + frame [(int) channel] [(int) (i+1)] * frame [(int) channel] [(int) (i+1)]; /* Power spectrum. */
			}
			ac [(int) (nsampFFT-1)] += frame [(int) channel] [(int) (nsampFFT-1)] * frame [(int) channel] [(int) (nsampFFT-1)];   /* Nyquist frequency. */
		}
		fft.bt(ac);   /* Autocorrelation. */

		/*
		 * Normalize the autocorrelation to the value with zero lag,
		 * and divide it by the normalized autocorrelation of the window.
		 */
		r [(int) nsamp_window] = 1.0;
                for (long i = 1; i <= brent_ixmax; i ++)
			r [(int) (nsamp_window - i)] = r [(int) (nsamp_window + i)] = ac [(int) i] / (ac [0] * windowR [(int) i]);
		/*
		 * Register the first candidate, which is always present: voicelessness.
		 */
		pitchFrame.nCandidates = 1;
		pitchFrame.candidate[0].frequency = 0.0;   /* Voiceless: always present. */
		pitchFrame.candidate[0].strength = 0.0;

		// Shortcut: absolute silence is always voiceless. We are done for this frame.
		if (localPeak == 0) return;
		
		 // Find the strongest maxima of the correlation of this frame, and register them as candidates.
		imax [0] = 0;
		for (long i = 2; i < maximumLag && i < brent_ixmax; i ++)
        {
            if (r [(int) (i+nsamp_window)] > 0.5 * voicingThreshold  && /* Not too unvoiced? */ r [(int) (i+nsamp_window)] > r [(int) (i-1+nsamp_window)] && r [(int) (i+nsamp_window)] >= r [(int) (i+1+nsamp_window)])   /* Maximum? */
            {
                
                int place = 0;

                /*
                 * Use parabolic interpolation for first estimate of frequency,
                 * and sin(x)/x interpolation to compute the strength of this frequency.
                 */
                double dr = 0.5 * (r [(int) (i+1+nsamp_window)] - r [(int) (i-1+nsamp_window)]), d2r = 2 * r [(int) (i+nsamp_window)] - r [(int) (i-1+nsamp_window)] - r [(int) (i+1+nsamp_window)];
                double frequencyOfMaximum = 1 / dx / (i + dr / d2r);
                long offset = - brent_ixmax - 1;
                double strengthOfMaximum = /* method & 1 ? */
                        NUM_interpolate_sinc (Arrays.copyOfRange(r, (int) (offset-1+nsamp_window), r.length), brent_ixmax - offset, 1 / dx / frequencyOfMaximum - offset, 30)
                        /* : r [i] + 0.5 * dr * dr / d2r */;
                /* High values due to short windows are to be reflected around 1. */
                if (strengthOfMaximum > 1.0) strengthOfMaximum = 1.0 / strengthOfMaximum;

                /*
                 * Find a place for this maximum.
                 */
                if (pitchFrame.nCandidates < maxnCandidates) 
                { /* Is there still a free place? */
                        place = ++ pitchFrame.nCandidates;
                } 
                else 
                {
                        /* Try the place of the weakest candidate so far. */
                        double weakest = 1;
                        for (int iweak = 1; iweak < maxnCandidates; iweak ++) 
                        {
                                /* High frequencies are to be favoured */
                                /* if we want to analyze a perfectly periodic signal correctly. */
                                double localStrength = pitchFrame.candidate[iweak].strength - octaveCost *
                                        (Math.log(minimumPitch / pitchFrame.candidate[iweak].frequency) * NUMlog2e);
                                if (localStrength < weakest) { weakest = localStrength; place = iweak; }
                        }
                        /* If this maximum is weaker than the weakest candidate so far, give it no place. */
                        if (strengthOfMaximum - octaveCost * (Math.log(minimumPitch / frequencyOfMaximum) * NUMlog2e) <= weakest)
                                place = 0;
                }
                if (place != 0) 
                {   /* Have we found a place for this candidate? */
                        pitchFrame.candidate[place-1].frequency = frequencyOfMaximum;
                        pitchFrame.candidate[place-1].strength = strengthOfMaximum;
                        imax [place-1] = i;
                }
            }
        }

		/*
         * Second pass: for extra precision, maximize sin(x)/x interpolation ('sinc').
         */
        for (long i = 0; i < pitchFrame.nCandidates-1; i ++) 
        {
            if (pitchFrame.candidate[(int) i].frequency > 0.0 / dx) 
            {
                    double ymid;
                    double[] xmid = new double[1];
                    long offset = - brent_ixmax - 1;
                    ymid = NUMimproveMaximum (Arrays.copyOfRange(r, (int) (offset+nsamp_window), r.length), brent_ixmax - offset, imax [(int) i] - offset,
                            (int) (pitchFrame.candidate[(int) i].frequency > 0.3 / dx ? 4 : brent_depth), xmid);
                    xmid[0] += offset;
                    pitchFrame.candidate[(int) i].frequency = 1.0 / dx / xmid[0];
                    if (ymid > 1.0) ymid = 1.0 / ymid;
                    pitchFrame.candidate[(int) i].strength = ymid;
            }
        }    
	}
	
	double NUMimproveMaximum (double[] y, long nx, long ixmid, int interpolation, double[] ixmid_real)
	{ 
		return NUMimproveExtremum (y, nx, ixmid, interpolation, ixmid_real, 1); 
	}
	
	double NUMimproveMinimum (double[] y, long nx, long ixmid, int interpolation, double[] ixmid_real)
	{ 
		return NUMimproveExtremum (y, nx, ixmid, interpolation, ixmid_real, 0); 
	}
	
	double NUM_interpolate_sinc (double y [], long nx, double x, long maxDepth) 
	{
		long ix, midleft = (long) Math.floor(x), midright = midleft + 1, left, right;
		double result = 0.0, a, halfsina, aa, daa;
		if (nx < 1) return NUMundefined; 
		if (x > nx-1) return y [NumSamples-1]; 
		if (x < 0) return y [0]; 
		if (x == midleft) return y [(int) midleft]; 
		/* 1 < x < nx && x not integer: interpolate. */ 
		if (maxDepth > midright - 1) maxDepth = midright - 1; 
		if (maxDepth > (nx-1) - midleft) maxDepth = (nx-1) - midleft; 
		if (maxDepth <= 0) return y [(int) Math.floor(x + 0.5)]; 
		if (maxDepth == 1) return y [(int) midleft] + (x - midleft) * (y [(int) midright] - y [(int) midleft]); 
		if (maxDepth == 2) 
		{ 
			double yl = y [(int) midleft], yr = y [(int) midright]; 
			double dyl = 0.5 * (yr - y [(int) (midleft - 1)]), dyr = 0.5 * (y [(int) (midright + 1)] - yl); 
			double fil = x - midleft, fir = midright - x; 
			return yl * fir + yr * fil - fil * fir * (0.5 * (dyr - dyl) + (fil - 0.5) * (dyl + dyr - 2 * (yr - yl))); 
		}
		left = midright - maxDepth;
		right = midleft + maxDepth;
		a = NUMpi * (x - midleft);
		halfsina = 0.5 * Math.sin(a);
		aa = a / (x - left + 1);
		daa = NUMpi / (x - left + 1);
		for (ix = midleft; ix >= left; ix --) 
		{
			double d = halfsina / a * (1.0 + Math.cos(aa));
			result += y[(int) ix] * d;
			a += NUMpi;
			aa += daa;
			halfsina = - halfsina;
		}
		a = NUMpi * (midright - x);
		halfsina = 0.5 * Math.sin(a);
		aa = a / (right - x + 1);
		daa = NUMpi / (right - x + 1); 
		for (ix = midright; ix <= right; ix ++)
		{
			double d = halfsina / a * (1.0 + Math.cos(aa));
			result += y [(int) ix] * d;
			a += NUMpi;
			aa += daa;
			halfsina = - halfsina;
		}
		return result;
	}
	
	Manipulation Sound_Pitch_to_Manipulation (Pitch pitch) 
	{
		Manipulation me = new Manipulation(xmin, xmax);
		me.sound = this; //Sound_convertToMono();   	// already mono
		me.sound.subtractMean ();
		me.pulses = me.sound.Sound_Pitch_to_PointProcess_cc(pitch);
		me.pitch = pitch.Pitch_to_PitchTier();
		return me;
	}
	
	PointProcess Sound_Pitch_to_PointProcess_cc (Pitch pitch) 
	{
		PointProcess point = new PointProcess(xmin, xmax, 10);
		double t = pitch.xmin;
		double addedRight = -1e300;
		double globalPeak = Vector_getAbsoluteExtremum(xmin, xmax, 0);
		double[] peak = new double[1];
		
		/*
		 * Cycle over all voiced intervals.
		 */
		for (;;) 
		{
			double[] tleft = new double[1], tright = new double[1];
			if (pitch.Pitch_getVoicedIntervalAfter (t, tleft, tright) == 0) break;
			/*
			 * Go to the middle of the voice stretch.
			 */
			double tmiddle = (tleft[0] + tright[0]) / 2;
			double f0middle = pitch.Pitch_getValueAtTime (tmiddle, 0, 1);

			// f0middle == NUMundefined needs to be false
			
			/*
			 * Our first point is near this middle.
			 */
			double[] tmax = {Sound_findExtremum (tmiddle - 0.5 / f0middle, tmiddle + 0.5 / f0middle, 1, 1)};
			//Melder_assert (NUMdefined (tmax));
			point.addPoint(tmax[0]);

			double tsave = tmax[0];
			for (;;) 
			{
				double f0 = pitch.Pitch_getValueAtTime (tmax[0], 0, 1), correlation;
				if (f0 == NUMundefined) break;
				correlation = Sound_findMaximumCorrelation (tmax[0], 1.0 / f0, tmax[0] - 1.25 / f0, tmax[0] - 0.8 / f0, tmax, peak);
				if (correlation == -1) /*break*/ tmax[0] -= 1.0 / f0;   /* This one period will drop out. */
				if (tmax[0] < tleft[0]) 
				{
					if (correlation > 0.7 && peak[0] > 0.023333 * globalPeak && tmax[0] - addedRight > 0.8 / f0) 
					{
						point.addPoint(tmax[0]);
					}
					break;
				}
				if (correlation > 0.3 && (peak[0] == 0.0 || peak[0] > 0.01 * globalPeak)) 
				{
					if (tmax[0] - addedRight > 0.8 / f0) {   // do not fill in a short originally unvoiced interval twice
						point.addPoint(tmax[0]);
					}
				}
			}
			tmax[0] = tsave;
			for (;;) 
			{
				double f0 = pitch.Pitch_getValueAtTime (tmax[0], 0, 1), correlation;
				if (f0 == NUMundefined) break;
				correlation = Sound_findMaximumCorrelation (tmax[0], 1.0 / f0, tmax[0] + 0.8 / f0, tmax[0] + 1.25 / f0, tmax, peak);
				if (correlation == -1) /*break*/ tmax[0] += 1.0 / f0;
				if (tmax[0] > tright[0]) 
				{
					if (correlation > 0.7 && peak[0] > 0.023333 * globalPeak) 
					{
						point.addPoint(tmax[0]);
						addedRight = tmax[0];
					}
					break;
				}
				if (correlation > 0.3 && (peak[0] == 0.0 || peak[0] > 0.01 * globalPeak)) 
				{
					point.addPoint(tmax[0]);
					addedRight = tmax[0];
				}
			}
			t = tright[0];
		}
		return point;
	}
	
	//Melder_assert (NUMdefined (tmin));
	//Melder_assert (NUMdefined (tmax));
	double Sound_findExtremum (double tmin, double tmax, int includeMaxima, int includeMinima) 
	{
		long imin = Sampled_xToLowIndex(tmin), imax = Sampled_xToHighIndex (tmax);
		if (imin < 1) imin = 1;
		if (imax > NumSamples) imax = NumSamples;
		double iextremum = findExtremum_3 (samples1, null, imin - 1, imax - imin + 1, includeMaxima, includeMinima);
		if (iextremum != 0)
			return x1 + (imin - 1 + iextremum - 1) * dx;
		else
			return (tmin + tmax) / 2;
	}
	
	static double findExtremum_3 (double[] channel1_base, double[] channel2_base, long d, long n, int includeMaxima, int includeMinima) 
	{
		//double[] channel1 = channel1_base + d; channel1 = channel1_base[d+ ];
		double[] channel2 = null;
		boolean includeAll = includeMaxima == includeMinima;
		long imin = 1, imax = 1, i, iextr;
		double minimum, maximum;
		if (n < 3) 
		{
			if (n <= 0) return 0.0;   /* Outside. */
			else if (n == 1) return 1.0;
			else 
			{   /* n == 2 */
				double x1 = channel1_base [(int) d];
				double x2 = channel1_base [(int) (d+1)];
				double xleft = includeAll ? Math.abs(x1) : includeMaxima != 0  ? x1 : - x1;
				double xright = includeAll ? Math.abs (x2) : includeMaxima != 0 ? x2 : - x2;
				if (xleft > xright) return 1.0;
				else if (xleft < xright) return 2.0;
				else return 1.5;
			}
		}
		minimum = maximum = channel1_base [(int) d];
		for (i = 1; i < n; i ++) 
		{
			double value = channel1_base [(int) (i+d)];
			if (value < minimum) { minimum = value; imin = i; }
			if (value > maximum) { maximum = value; imax = i; }
		}
		if (minimum == maximum) 
		{
			return 0.5 * (n + 1.0);   /* All equal. */
		}
		iextr = includeAll ? ( Math.abs (minimum) > Math.abs (maximum) ? imin : imax ) : includeMaxima != 0 ? imax : imin;
		if (iextr == 1) return 1.0;
		if (iextr == n) return (double) n;
		/* Parabolic interpolation. */
		/* We do NOT need fabs here: we look for a genuine extremum. */
		double valueMid = channel1_base [(int) (d + iextr)];
		double valueLeft = channel1_base [(int) (d + iextr - 1)];
		double valueRight = channel1_base [(int) (d + iextr + 1)];
		return iextr + 0.5 * (valueRight - valueLeft) / (2 * valueMid - valueLeft - valueRight);
	}
	
	double Sound_findMaximumCorrelation (double t1, double windowLength, double tmin2, double tmax2, double[] tout, double[] peak) 
	{
		double maximumCorrelation = -1.0, r1 = 0.0, r2 = 0.0, r3 = 0.0, r1_best, r3_best, ir;
		double halfWindowLength = 0.5 * windowLength;
		long ileft1 = Sampled_xToNearestIndex(t1 - halfWindowLength);
		long iright1 = Sampled_xToNearestIndex (t1 + halfWindowLength);
		long ileft2min = Sampled_xToLowIndex (tmin2 - halfWindowLength);
		long ileft2max = Sampled_xToHighIndex (tmax2 - halfWindowLength);
		peak[0] = 0.0;   /* Default. */
		r1_best = r3_best = ir = 0.0;
		for (long ileft2 = ileft2min; ileft2 <= ileft2max; ileft2 ++) 
		{
			double norm1 = 0.0, norm2 = 0.0, product = 0.0, localPeak = 0.0;
			for (long i1 = ileft1, i2 = ileft2; i1 <= iright1; i1 ++, i2 ++) 
			{
				if (i1 < 0 || i1 >= NumSamples || i2 < 0 || i2 >= NumSamples) continue;
				double amp1 = samples1[(int) i1], amp2 = samples1[(int) i2];
				norm1 += amp1 * amp1;
				norm2 += amp2 * amp2;
				product += amp1 * amp2;
				if (Math.abs (amp2) > localPeak)
					localPeak = Math.abs (amp2);
			}
			r1 = r2;
			r2 = r3;
			r3 = product != 0.0 ? product / (Math.sqrt(norm1 * norm2)) : 0.0;
			if (r2 > maximumCorrelation && r2 >= r1 && r2 >= r3) 
			{
				r1_best = r1;
				maximumCorrelation = r2;
				r3_best = r3;
				ir = ileft2 - 1;
				peak[0] = localPeak;  
			}
		}
		/*
		 * Improve the result by means of parabolic interpolation.
		 */
		if (maximumCorrelation > -1.0) 
		{
			double d2r = 2 * maximumCorrelation - r1_best - r3_best;
			if (d2r != 0.0) {
				double dr = 0.5 * (r3_best - r1_best);
				maximumCorrelation += 0.5 * dr * dr / d2r;
				ir += dr / d2r;
			}
			tout[0] = t1 + (ir - ileft1) * dx;
		}
		return maximumCorrelation;
	}
	
	double NUMimproveExtremum (double[] y, long nx, long ixmid, int interpolation, double[] ixmid_real, int isMaximum) 
	{
		ImproveParams params = new ImproveParams();
		double[] result = new double[1];
		if (ixmid <= 1) { ixmid_real[0] = 1; return y [1]; }
		if (ixmid >= nx) { ixmid_real[0] = nx; return y [NumSamples]; }
		if (interpolation <= 0) { ixmid_real[0] = ixmid; return y [(int) ixmid]; }
		if (interpolation == 1) 
		{
			double dy = 0.5 * (y [(int) (ixmid + 1)] - y [(int) (ixmid - 1)]);
			double d2y = 2 * y [(int) ixmid] - y [(int) (ixmid - 1)] - y [(int) (ixmid + 1)];
			ixmid_real[0] = ixmid + dy / d2y;
			return y [(int) ixmid] + 0.5 * dy * dy / d2y;
		}
		 //Sinc interpolation. 
		params.y = y;
		params.depth = interpolation == 3 ? 70 : 700; // int
		params.ixmax = NumSamples; // int
		params.isMaximum = isMaximum; // int
		ixmid_real[0] = NUMminimize_brent (ixmid - 1, ixmid + 1, params, 1e-10, result);
		return isMaximum == 1 ? - result[0] : result[0];
	}
	
	double improve_evaluate (double x, ImproveParams params) 
	{
		double y = NUM_interpolate_sinc(params.y, (long) params.ixmax, x, (long) params.depth);
		return params.isMaximum == 1 ? - y : y;
	}
	
	double NUMminimize_brent (double a, double b, ImproveParams params, double tol, double[] fx) 
	{
		double x, v, fv, w, fw;
		final double golden = 1 - NUM_goldenSection;
		final double sqrt_epsilon = Math.sqrt(1.11e-16); //1.11e-16, 2.22e-16?
		long itermax = 60;
		
		//Melder_assert (tol > 0 && a < b);
		//First step - golden section 
		
		v = a + golden * (b - a);
		fv = improve_evaluate(v, params);
		x = v;  w = v;
		fx[0] = fv;  fw = fv;
		
		for (long iter = 1; iter <= itermax; iter++) 
		{
			double range = b - a;
			double middle_range = (a + b) / 2;
			double tol_act = sqrt_epsilon * Math.abs (x) + tol / 3;
			double new_step;  //Step at this iteration 
			
			if (Math.abs (x - middle_range) + range / 2 <= 2 * tol_act) 
			{
				return x;
			}
			
			//Obtain the golden section step 
			
			new_step = golden * (x < middle_range ? b - x : a - x);
			
			//Decide if the parabolic interpolation can be tried	
			
			if (Math.abs (x - w) >= tol_act) 
			{
			 
			 	//Interpolation step is calculated as p/q;
			 	//division operation is delayed until last moment.
				
				double p, q, t;
				t = (x - w) * (fx[0] - fv);
				q = (x - v) * (fx[0] - fw);
				p = (x - v) * q - (x - w) * t;
				q = 2 * (q - t);
			
				if (q > 0) 
				{
					p = -p;
				} 
				else 
				{
					q = -q;
				}
				
				/*If x+p/q falls in [a,b], not too close to a and b,
				and isn't too large, it is accepted.
				If p/q is too large then the golden section procedure can
				reduce [a,b] range.
				 */
			
				if (Math.abs (p) < Math.abs (new_step * q) && p > q * (a - x + 2 * tol_act) && p < q * (b - x - 2 * tol_act)) 
				{
					new_step = p / q;
				}
			}
			
			//Adjust the step to be not less than tolerance. 
			
			if (Math.abs (new_step) < tol_act) 
			{
				new_step = new_step > 0 ? tol_act : - tol_act;
			}
			
			 //Obtain the next approximation to min	and reduce the enveloping range 
			
			{
				double t = x + new_step;	 //Tentative point for the min	
				double ft = improve_evaluate(t, params);
				
				/*If t is a better approximation, reduce the range so that
				t would fall within it. If x remains the best, reduce the range
				so that x falls within it.*/
				
				if (ft <= fx[0]) 
				{
					if (t < x) 
					{
						b = x;
					} 
					else 
					{
						a = x;
					}
					
					v = w;  w = x;  x = t;
					fv = fw;  fw = fx[0];  
					fx[0] = ft;
				} 
				else 
				{
					if (t < x) 
					{
						a = t;
					} 
					else 
					{
						b = t;
					}
					
					if (ft <= fw || w == x) 
					{
						v = w; w = t;
						fv = fw; fw = ft;
					} 
					else if (ft <= fv || v == x || v == w) 
					{
						v = t;
						fv = ft;
					}
				}
			}
		}
		//Melder_warning (L"NUMminimize_brent: maximum number of iterations (", Melder_integer (itermax), L") exceeded.");
		return x;
	}
	
	Wave Sound_Point_Pitch_Duration_to_Sound (PointProcess pulses, PitchTier pitch, DurationTier duration, double maxT)
	{
		long ipointleft, ipointright;
		double deltat = 0, handledTime = xmin;
		double startOfSourceNoise, endOfSourceNoise, startOfTargetNoise, endOfTargetNoise;
		double durationOfSourceNoise, durationOfTargetNoise;
		double startOfSourceVoice, endOfSourceVoice, startOfTargetVoice, endOfTargetVoice;
		double durationOfSourceVoice, durationOfTargetVoice;
		double startingPeriod, finishingPeriod, ttarget, voicelessPeriod;
		if (duration.size == 0)
		{
			System.out.println("No duration points.");
		}

		/*
		 * Create a Sound long enough to hold the longest possible duration-manipulated sound.
		 */
		Wave thee = new Wave(1, xmin, xmin + 3 * (xmax - xmin), 3 * NumSamples, dx, x1);

		/*
		 * Below, I'll abbreviate the voiced interval as "voice" and the voiceless interval as "noise".
		 */
		if (pitch != null && pitch.size != 0) for (ipointleft = 0; ipointleft < pulses.nt; ipointleft = ipointright + 1) 
		{
			/*
			 * Find the beginning of the voice.
			 */
			startOfSourceVoice = pulses.t.get((int) ipointleft);   /* The first pulse of the voice. */
			startingPeriod = 1.0 / pitch.RealTier_getValueAtTime(startOfSourceVoice);
			startOfSourceVoice -= 0.5 * startingPeriod;   /* The first pulse is in the middle of a period. */

			/*
			 * Measure one noise.
			 */
			startOfSourceNoise = handledTime;
			endOfSourceNoise = startOfSourceVoice;
			durationOfSourceNoise = endOfSourceNoise - startOfSourceNoise;
			startOfTargetNoise = startOfSourceNoise + deltat;
			endOfTargetNoise = startOfTargetNoise + duration.RealTier_getArea(startOfSourceNoise, endOfSourceNoise);
			durationOfTargetNoise = endOfTargetNoise - startOfTargetNoise;

			/*
			 * Copy the noise.
			 */
			voicelessPeriod = NUMrandomUniform (0.008, 0.012);
			ttarget = startOfTargetNoise + 0.5 * voicelessPeriod;
			while (ttarget < endOfTargetNoise) 
			{
				double tsource;
				double tleft = startOfSourceNoise, tright = endOfSourceNoise;
				int i;
				for (i = 1; i <= 15; i ++) 
				{
					double tsourcemid = 0.5 * (tleft + tright);
					double ttargetmid = startOfTargetNoise + duration.RealTier_getArea(startOfSourceNoise, tsourcemid);
					if (ttargetmid < ttarget) tleft = tsourcemid; else tright = tsourcemid;
				}
				tsource = 0.5 * (tleft + tright);
				copyBell (tsource, voicelessPeriod, voicelessPeriod, thee, ttarget);
				voicelessPeriod = NUMrandomUniform (0.008, 0.012);
				ttarget += voicelessPeriod;
			}
			deltat += durationOfTargetNoise - durationOfSourceNoise;

			/*
			 * Find the end of the voice.
			 */
			for (ipointright = ipointleft + 1; ipointright < pulses.nt; ipointright ++)
				if (pulses.t.get((int) ipointright) - pulses.t.get((int) (ipointright - 1)) > maxT)
					break;
			ipointright --;
			endOfSourceVoice = pulses.t.get((int) ipointright);   /* The last pulse of the voice. */
			finishingPeriod = 1.0 / pitch.RealTier_getValueAtTime(endOfSourceVoice);
			endOfSourceVoice += 0.5 * finishingPeriod;   /* The last pulse is in the middle of a period. */
			/*
			 * Measure one voice.
			 */
			durationOfSourceVoice = endOfSourceVoice - startOfSourceVoice;

			/*
			 * This will be copied to an interval with a different location and duration.
			 */
			startOfTargetVoice = startOfSourceVoice + deltat;
			endOfTargetVoice = startOfTargetVoice +
				duration.RealTier_getArea(startOfSourceVoice, endOfSourceVoice);
			durationOfTargetVoice = endOfTargetVoice - startOfTargetVoice;

			/*
			 * Copy the voiced part.
			 */
			ttarget = startOfTargetVoice + 0.5 * startingPeriod;
			while (ttarget < endOfTargetVoice) 
			{
				double tsource, period;
				long isourcepulse;
				double tleft = startOfSourceVoice, tright = endOfSourceVoice;
				int i;
				for (i = 1; i <= 15; i ++) 
				{
					double tsourcemid = 0.5 * (tleft + tright);
					double ttargetmid = startOfTargetVoice + duration.RealTier_getArea(startOfSourceVoice, tsourcemid);
					if (ttargetmid < ttarget) tleft = tsourcemid; else tright = tsourcemid;
				}
				tsource = 0.5 * (tleft + tright);
				period = 1.0 / pitch.RealTier_getValueAtTime(tsource);
				isourcepulse = pulses.PointProcess_getNearestIndex(tsource);
				copyBell2 (pulses, isourcepulse, period, period, thee, ttarget, maxT);
				ttarget += period;
			}
			deltat += durationOfTargetVoice - durationOfSourceVoice;
			handledTime = endOfSourceVoice;
		}

		/*
		 * Copy the remaining unvoiced part, if we are at the end.
		 */
		startOfSourceNoise = handledTime;
		endOfSourceNoise = xmax;
		durationOfSourceNoise = endOfSourceNoise - startOfSourceNoise;
		startOfTargetNoise = startOfSourceNoise + deltat;
		endOfTargetNoise = startOfTargetNoise + duration.RealTier_getArea(startOfSourceNoise, endOfSourceNoise);
		durationOfTargetNoise = endOfTargetNoise - startOfTargetNoise;
		voicelessPeriod = NUMrandomUniform (0.008, 0.012);
		ttarget = startOfTargetNoise + 0.5 * voicelessPeriod;
		while (ttarget < endOfTargetNoise) 
		{
			double tsource;
			double tleft = startOfSourceNoise, tright = endOfSourceNoise;
			for (int i = 1; i <= 15; i ++) 
			{
				double tsourcemid = 0.5 * (tleft + tright);
				double ttargetmid = startOfTargetNoise + duration.RealTier_getArea(startOfSourceNoise, tsourcemid);
				if (ttargetmid < ttarget) tleft = tsourcemid; else tright = tsourcemid;
			}
			tsource = 0.5 * (tleft + tright);
			copyBell (tsource, voicelessPeriod, voicelessPeriod, thee, ttarget);
			voicelessPeriod = NUMrandomUniform (0.008, 0.012);
			ttarget += voicelessPeriod;
		}

		/*
		 * Find the number of trailing zeroes and hack the sound's time domain.
		 */
		thee.xmax = thee.xmin + duration.RealTier_getArea(xmin, xmax);
		if (Math.abs (thee.xmax - xmax) < 1e-12) thee.xmax = xmax;   /* Common situation. */
		thee.NumSamples = (int) thee.Sampled_xToLowIndex(thee.xmax);
		if (thee.NumSamples > 3 * NumSamples) thee.NumSamples = 3 * NumSamples;

		return thee;
	}
	
	double NUMrandomUniform(double lo, double hi)
	{
		int num = (int) ((hi - lo) / step);
		Random rand = new Random();
		int random = rand.nextInt(num);
		double scaled = ((double)random) * step;
		double shifted = scaled + lo;
		return shifted;
	}
	
	void copyBell2 (PointProcess source, long isource, double leftWidth, double rightWidth,
			Wave thee, double tmidTarget, double maxT)
	{
		/*
		 * Replace 'leftWidth' and 'rightWidth' by the lengths of the intervals in the source (instead of target),
		 * if these are shorter.
		 */
		double tmid = source.t.get((int) isource);
		if (isource > 0 && tmid - source.t.get((int) (isource - 1)) <= maxT) 
		{
			double sourceLeftWidth = tmid - source.t.get((int) (isource - 1));
			if (sourceLeftWidth < leftWidth) leftWidth = sourceLeftWidth;
		}
		if (isource < source.nt-1 && source.t.get((int) (isource + 1)) - tmid <= maxT) 
		{
			double sourceRightWidth = source.t.get((int) (isource + 1)) - tmid;
			if (sourceRightWidth < rightWidth) rightWidth = sourceRightWidth;
		}
		copyBell (tmid, leftWidth, rightWidth, thee, tmidTarget);
	}
	
	void copyBell (double tmid, double leftWidth, double rightWidth, Wave thee, double tmidTarget) 
	{
		copyRise (tmid - leftWidth, tmid, thee, tmidTarget);
		copyFall(tmid, tmid + rightWidth, thee, tmidTarget);
	}

	void copyRise (double tmin, double tmax, Wave thee, double tmaxTarget) 
	{
		long imin, imax, imaxTarget, distance, i;
		double dphase;
		imin = Sampled_xToHighIndex(tmin);
		if (imin < 0) imin = 0;
		imax = Sampled_xToHighIndex(tmax) - 1;   /* Not xToLowIndex: ensure separation of subsequent calls. */
		if (imax >= NumSamples) imax = NumSamples-1;
		if (imax < imin) return;
		imaxTarget = thee.Sampled_xToHighIndex(tmaxTarget) - 1;
		distance = imaxTarget - imax;
		dphase = NUMpi / (imax - imin + 1);
		for (i = imin; i <= imax; i ++) 
		{
			long iTarget = i + distance;
			if (iTarget >= 0 && iTarget < thee.NumSamples)
				thee.samples1[(int) iTarget] += samples1[(int) i] * 0.5 * (1 - Math.cos(dphase * (i - imin + 0.5)));
		}
	}
	
	void copyFall (double tmin, double tmax, Wave thee, double tminTarget) 
	{
		long imin, imax, iminTarget, distance, i;
		double dphase;
		imin = Sampled_xToHighIndex(tmin);
		if (imin < 0) imin = 0;
		imax = Sampled_xToHighIndex(tmax) - 1;   /* Not xToLowIndex: ensure separation of subsequent calls. */
		if (imax >= NumSamples) imax = NumSamples-1;
		if (imax < imin) return;
		iminTarget = thee.Sampled_xToHighIndex(tminTarget);
		distance = iminTarget - imin;
		dphase = NUMpi / (imax - imin + 1);
		for (i = imin; i <= imax; i ++) 
		{
			long iTarget = i + distance;
			if (iTarget >= 0 && iTarget < thee.NumSamples)
				thee.samples1[(int) iTarget] += samples1[(int) i] * 0.5 * (1 + Math.cos (dphase * (i - imin + 0.5)));
		}
	}
	
	Wave Sound_changeGender_old (double fmin, double fmax, double formantRatio, double new_pitch, double pitchRangeFactor, double durationFactor) 
	{
		Pitch pitch = Sound_to_Pitch (0.8 / fmin, fmin, fmax);
		Wave thee = Sound_and_Pitch_changeGender_old (pitch, formantRatio, new_pitch, pitchRangeFactor, durationFactor);
		return thee;
	}
	
	Wave Sound_and_Pitch_changeGender_old (Pitch him, double formantRatio, double new_pitch, double pitchRangeFactor, double durationFactor) 
	{
		double samplingFrequency_old = 1 / dx;
		
		if (xmin != him.xmin || xmax != him.xmax) 
		{
			System.out.println("The Pitch and the Sound object must have the same starting times and finishing times.");
		}
		if (new_pitch < 0) 
		{
			System.out.println("The new pitch median must not be negative.");
		}

		Wave sound = Wave_copy();
		sound.subtractMean();

		if (formantRatio != 1) 
		{
			// Shift all frequencies (inclusive pitch!)
			sound.overrideSamplingFrequency (samplingFrequency_old * formantRatio);
		}

		Pitch pitch = him.Pitch_scaleTime_old (1 / formantRatio);
		PointProcess pulses = sound.Sound_Pitch_to_PointProcess_cc (pitch);
		PitchTier pitchTier = pitch.Pitch_to_PitchTier();

		double median = pitch.Pitch_getQuantile(0, 0, 0.5, 0);
		if (median != 0 && median != NUMundefined) 
		{
			// Incorporate pitch shift from overriding the sampling frequency
			if (new_pitch == 0) 
			{
				new_pitch = median / formantRatio;
			}
			double factor = new_pitch / median;
			pitchTier.PitchTier_multiplyFrequencies (sound.xmin, sound.xmax, factor);
			pitchTier.PitchTier_modifyRange_old (sound.xmin, sound.xmax, pitchRangeFactor, new_pitch);
		} 
		else 
		{
			System.out.println("There were no voiced segments found.");
		}
		DurationTier duration = new DurationTier(xmin, xmax);
		duration.addPoint ((xmin + xmax) / 2, formantRatio * durationFactor);

		Wave thee = sound.Sound_Point_Pitch_Duration_to_Sound (pulses, pitchTier, duration, 1.25 / pitch.Pitch_getMinimum (0.0, 0.0, 0, 0));

		// Resample to the original sampling frequency

		if (formantRatio != 1) 
		{
			thee = thee.Sound_resample(samplingFrequency_old, 10);
		}
		return thee;
		
	}


	Wave Sound_resample (double samplingFrequency, long precision) 
	{
		double upfactor = samplingFrequency * dx;
		if (Math.abs (upfactor - 2) < 1e-6) return Sound_upsample();
		if (Math.abs (upfactor - 1) < 1e-6) return Wave_copy();
		long numberOfSamples = (long) Math.floor((xmax - xmin) * samplingFrequency + 0.5);
		if (numberOfSamples < 1)
			System.out.println("The resampled Sound would have no samples.");
		Wave filtered = null;
		if (upfactor < 1.0) /* Need anti-aliasing filter? */
		{   
			long nfft = 1, antiTurnAround = 1000;
			while (nfft < NumSamples + antiTurnAround * 2) nfft *= 2;
			double[] data = new double[(int) nfft];
			//autoNUMvector <double> data (1, nfft);
			filtered = new Wave(NumChannels, xmin, xmax, NumSamples, dx, x1);
			for (long channel = 1; channel <= NumChannels; channel ++) 
			{
				for (long i = 0; i < nfft; i ++) 
				{
					data [(int) i] = 0;
				}
				//NUMvector_copyElements (my z [channel], & data [antiTurnAround], 1, my nx);
				System.arraycopy(samples1, 0, data, (int) antiTurnAround, NumSamples);
				RealDoubleFFT fft = new RealDoubleFFT((int) nfft);
				fft.ft(data);
				//NUMrealft (data.peek(), nfft, 1);   // go to the frequency domain
				for (long i = (long) Math.floor(upfactor * nfft); i < nfft; i ++) 
				{
					data [(int) i] = 0;   /* Filter away high frequencies. */
				}
				data [1] = 0.0;
				fft.bt(data);
				//NUMrealft (data.peek(), nfft, -1);   // return to the time domain
				double factor = 1.0 / nfft;
				double[] to = filtered.samples1;
				for (long i = 0; i < NumSamples; i ++) 
				{
					to [(int) i] = data [(int) (i + antiTurnAround)] * factor;
				}
			}
			//me = filtered;   // reference copy; remove at end
		}
		else
			filtered = this;
		Wave thee = new Wave(filtered.NumChannels, filtered.xmin, filtered.xmax, (int) numberOfSamples, 1.0 / samplingFrequency,
			0.5 * (filtered.xmin + filtered.xmax - (numberOfSamples - 1) / samplingFrequency));
		for (long channel = 1; channel <= filtered.NumChannels; channel ++) 
		{
			double[] from = filtered.samples1;
			double[] to = thee.samples1;
			if (precision <= 1) 
			{
				for (long i = 0; i < numberOfSamples; i ++) 
				{
					double x = thee.Sampled_indexToX (i);
					double index = filtered.Sampled_xToIndex (x);
					long leftSample = (long) Math.floor(index);
					double fraction = index - leftSample;
					to [(int) i] = leftSample < 0 || leftSample >= filtered.NumSamples-1 ? 0.0 :
						(1 - fraction) * from [(int) leftSample] + fraction * from [(int) (leftSample + 1)];
				}
			} 
			else 
			{
				for (long i = 0; i < numberOfSamples; i ++) 
				{
					double x = thee.Sampled_indexToX(i);
					double index = filtered.Sampled_xToIndex(x);
					to [(int) i] = filtered.NUM_interpolate_sinc (filtered.samples1, filtered.NumSamples, index, precision);
				}
			}
		}
		return thee;
	}
	
	Wave Sound_upsample () 
	{
		long nfft = 1;
		while (nfft < NumSamples + 2000) nfft *= 2;
		Wave thee = new Wave(NumChannels, xmin, xmax, NumSamples * 2, dx / 2, x1 - dx / 4);
		double[] data = new double[(int) (nfft*2)];
		//autoNUMvector <double> data (1, 2 * nfft);
		for (long channel = 1; channel <= NumChannels; channel ++) 
		{
			System.arraycopy(samples1, 0, data, 1000, NumSamples);
			//NUMvector_copyElements (my z [channel], & data [1000], 1, my nx);
			RealDoubleFFT fft = new RealDoubleFFT((int) nfft*2);
			fft.ft(data);
			//NUMrealft (data.peek(), nfft, 1);
			long imin = (long) (nfft * 0.95);
			for (long i = imin; i < nfft; i ++) 
			{
				data [(int) i] *= ((double) (nfft - i)) / (nfft - imin);
			}
			data [1] = 0.0;
			//NUMrealft (data.peek(), 2 * nfft, -1);
			fft.bt(data);
			double factor = 1.0 / nfft;
			for (long i = 0; i < thee.NumSamples; i ++) 
			{
				thee.samples1[(int) i] = data [(int) (i + 2000)] * factor;
			}
		}
		return thee;
	}

	public Wave Sound_extractPart (double t1, double t2, double relativeWidth, boolean preserveTimes)
	{
		/*
		 * We do not clip to the Sound's time domain.
		 * Any samples outside it are taken to be zero.
		 */

		/*
		 * Autowindow.
		 */
		if (t1 == t2) { t1 = xmin; t2 = xmax; };
		/*
		 * Allow window tails outside specified domain.
		 */
		if (relativeWidth != 1.0)
		{
			double margin = 0.5 * (relativeWidth - 1) * (t2 - t1);
			t1 -= margin;
			t2 += margin;
		}
		/*
		 * Determine index range. We use all the real or virtual samples that fit within [t1..t2].
		 */
		long ix1 = 1 + (long) Math.ceil ((t1 - x1) / dx);
		long ix2 = 1 + (long) Math.floor ((t2 - x1) / dx);
		//if (ix2 < ix1) Melder_throw (U"Extracted Sound would contain no samples.");
		/*
		 * Create sound, optionally shifted to [0..t2-t1].
		 */
		Wave thee = new Wave (NumChannels, t1, t2, (int) (ix2 - ix1 + 1), dx, x1 + (ix1 - 1) * dx);
		if (! preserveTimes) { thee.xmin = 0.0; thee.xmax -= t1; thee.x1 -= t1; }
		/*
		 * Copy only *real* samples into the new sound.
		 * The *virtual* samples will remain at zero.
		 */
		int a = (int) ( ix1 < 1 ? 1 : ix1 );
		int b = (int) ( ix2 > NumSamples ? NumSamples : ix2 );
		int amount = b - a + 1;
		System.arraycopy(samples1, 0, thee.samples1, (int) (1 - ix1), amount);

		//NUMvector_copyElements (samples1, thee.samples1 + 1 - ix1, ( ix1 < 1 ? 1 : ix1 ), ( ix2 > NumSamples ? NumSamples : ix2 ));
		/*
		 * Multiply by a window that extends throughout the target domain.
		 */
		//thee.Sound_multiplyByWindow ();
		return thee;
	}

}
