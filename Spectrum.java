package jPraat;

import fftpack.RealDoubleFFT;

/* Author: J. Colin Crowley
 * Representation of WAV data in frequency domain.
 */

// Class representing praat Spectrum objects.
public class Spectrum 
{
	// Members with direct links to praat members.
	int NumFrequencies;			// nx in praat
	double xmin;				// xmin in praat
	double xmax;				// xmax in praat
	double x1;					// x1 in praat
	double dx;					// dx in praat
	double[] rl;				// z[1] in praat
	double[] im;				// z[2] in praat
	
	// Members used for fourier transform
	double[] fft;

	int originalNumberOfSamples;

	final double NUMpi = 3.1415926535897932384626433832795028841972;
	
	Spectrum Spectrum_copy()
	{
		Spectrum newSpectrum = new Spectrum(xmax, NumFrequencies);
		newSpectrum.xmin = xmin;
		newSpectrum.dx = dx;
		newSpectrum.x1 = x1;
		newSpectrum.originalNumberOfSamples = originalNumberOfSamples;
		for(int i = 0; i < NumFrequencies; i++)
		{
			newSpectrum.rl[i] = rl[i];
			newSpectrum.im[i] = im[i];
		}
		return newSpectrum;
	}
	
	// Alternative Spectrum constructor in which the fourier data is provided.
	public Spectrum(double[] r, double[] i)
	{
		NumFrequencies = r.length;
		rl = new double[NumFrequencies];
		im = new double[i.length];
		System.arraycopy(r, 0, rl, 0, r.length);
		System.arraycopy(i, 0, im, 0, i.length);
		xmin = 0;
	    xmax = 22050;
	    x1 = 0;
	    //dx = (xmax - xmin) / (double)581184;
	    dx = 0.042057037353515625;
	}
	
	public Spectrum(double fmax, int numF)
	{
		xmin = 0;
		xmax = fmax;
		NumFrequencies = numF;
		dx = fmax / (numF-1);
		x1 = 0;
		rl = new double[NumFrequencies];
		im = new double[NumFrequencies];
	}
	
	void check()
	{
		for(int i = 0; i < im.length; i++)
		{
			//im[i] = 0;
			rl[i] = 0;
		}
	}
	
	// Get the real component of the fourier data.
	double[] getReal()
	{
		return rl;
	}
	
	// Get the imaginary component of the fourier data.
	double[] getImaginary()
	{
		return im;
	}
	
	// Get dx.
	double getDx()
	{
		return dx;
	}
	
	void log(Object aMsg)
	{
	    System.out.println(String.valueOf(aMsg));
	}
	
	void print()
	{
		System.out.println();
		log("SamplingRate=dx: " + dx);
	    /*log("ByteRate: " + ByteRate);
	    log("BlockAlign: " + BlockAlign);
	    log("BitsPerSample: " + BitsPerSample);
	    log("Subchunk2ID: " + Subchunk2ID);
	    log("Subchunk2Size: " + Subchunk2Size);*/
	    log("NumSamples: " + NumFrequencies);
	    //log("Size of file: " + bytes.length);
	    log("xmin: " + xmin);
	    log("xmax: " + xmax);
	    log("x1: " + x1);
		System.out.println("First ten real:");
		for(int i = 0; i < 10; i++)
			System.out.println(i+": "+rl[i]);
		System.out.println("");
		System.out.println("First ten imag:");
		for(int i = 0; i < 10; i++)
			System.out.println(i+": "+im[i]);
		System.out.println();
	}
	
	// Corresponds to praat function "Filter (pass Hann band)..."
	void filterPassHannBand (double fmin, double fmax0, double smooth) 
	{
		double fmax = fmax0 == 0.0 ? xmax : fmax0;
		double f1 = fmin - smooth, f2 = fmin + smooth, f3 = fmax - smooth, f4 = fmax + smooth;
		double halfpibysmooth = smooth != 0.0 ? NUMpi / (2 * smooth) : 0.0;
		for (int i = 0; i < NumFrequencies; i ++) 
		{
			double frequency = x1 + i * dx;
			if (frequency < f1 || frequency > f4) 
			{
				rl[i] = 0.0;
				im[i] = 0.0;
			}
			if (frequency < f2 && fmin > 0.0) 
			{
				double factor = 0.5 - 0.5 * Math.cos(halfpibysmooth * (frequency - f1));
				rl [i] *= factor;
				im [i] *= factor;
			} 
			else if (frequency > f3 && fmax < xmax) 
			{
				double factor = 0.5 + 0.5 * Math.cos(halfpibysmooth * (frequency - f3));
				rl [i] *= factor;
				im [i] *= factor;
			}
		}
	}
	
	void zeroEverything()
	{
		for(int i = 0; i < NumFrequencies; i++)
		{
			rl[i] = 0.0;
			im[i] = 0.0;
		}
	}
	
	void EqBand(double bnd1, double bnd2, double dB, double smooth, Spectrum buffer, Spectrum sp_pulse)
	{
		double factor = Math.pow(10,(dB/20));
		//System.out.println("factor" +factor);
		buffer.copySpectrum(sp_pulse);
		buffer.filterPassHannBand(bnd1, bnd2, smooth);
		buffer.multiplyScalar(factor);
		addSpectrumMultiplyScalar(buffer, 1.0);
	}
	
	void multiplyScalar(double scalar)
	{
		for(int i = 0; i < NumFrequencies; i++)
		{
			rl[i] *= scalar;
			im[i] *= scalar;
		}
	}
	
	void addSpectrumMultiplyScalar(Spectrum spectrum, double scalar)
	{
		for(int i = 0; i < NumFrequencies; i++)
		{
			rl[i] = rl[i] + scalar * spectrum.rl[i];
			im[i] = im[i] + scalar * spectrum.im[i];
		}
	}
	
	void copySpectrum(Spectrum spectrum)
	{
		for(int i = 0; i < NumFrequencies; i++)
		{
			rl[i] = spectrum.rl[i];
			im[i] = spectrum.im[i];
		}
	}
	
	//A Fourier-transformable Spectrum must have a first frequency of 0 Hz
	Wave Spectrum_to_Sound ()
	{
		double lastFrequency = x1 + (NumFrequencies - 1) * dx;
		boolean originalNumberOfSamplesProbablyOdd = im[NumFrequencies-1] != 0.0 || xmax - lastFrequency > 0.25 * dx;
		long numberOfSamples = 2 * NumFrequencies - ( originalNumberOfSamplesProbablyOdd ? 1 : 2 );
		Wave thee = new Wave(); 
		thee = thee.Sound_createSimple (1, 1.0 / dx, numberOfSamples * dx);
		double[] amp = thee.samples1;
		double scaling = dx;
		amp[0] = rl[0] * scaling;
		for (long i = 1; i < (NumFrequencies-1); i++) 
		{
			amp [(int) (i + i - 1)] = rl[(int) i] * scaling;
			amp [(int) (i + i)] = im[(int) i] * scaling;
		}
		if (originalNumberOfSamplesProbablyOdd) 
		{
			amp [(int) (numberOfSamples-1)] = rl[NumFrequencies-1] * scaling;
			if (numberOfSamples > 1) 
				amp [1] = im[NumFrequencies-1] * scaling;
		} 
		else 
		{
			amp [1] = rl[NumFrequencies-1] * scaling;
		}
		
		RealDoubleFFT fft = new RealDoubleFFT((int) numberOfSamples);	
		fft.bt(amp);
		
		for(int i = originalNumberOfSamples; i < numberOfSamples; i++)
			amp[i] = 0.0;
		
		return thee;
	}
	
	
	
	

}
