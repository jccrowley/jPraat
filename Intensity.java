package jPraat;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import java.util.ArrayList;



/* Author: J. Colin Crowley
 * Representation of WAV data in intensity-time domain.
 */

//Class representing praat Sound objects.
public class Intensity
{
	// globals
	final double NUMundefined = Double.POSITIVE_INFINITY;
	final double NUMpi = 3.1415926535897932384626433832795028841972;
	
	// Members with direct links to members in praat.
    int NumChannels;					// ny in praat
    int NumFrames;						// nx in praat
    double[] samples1;					// z in praat
    double xmin;						// xmin in praat. Start time
    double xmax;						// xmax in praat. End time
    double dx;							// dx in praat
    double x1;							// x1 in praat

    // Intensity object constructor.
 	public Intensity(double xmin0, double xmax0, int numberOfFrames, double timeStep, double thyFirstTime)
 	{
 		xmin = xmin0;
 		xmax = xmax0;
 		NumFrames = numberOfFrames;
 		dx = timeStep;
 		x1 = thyFirstTime;
 		samples1 = new double[NumFrames];
 	}
 	
 	// Praat helper function
 	double Sampled_indexToX (long index) 
 	{ 
 		return x1 + index * dx;
 	}
 	
 	// Praat helper function
 	long Sampled_xToNearestIndex (double x) 
 	{ 
 		return (long) Math.round((x - x1) / dx /*+ 1.0*/); 
 	}
 	
 	public void printInfo()
 	{
 		System.out.println("xmin: " + xmin);
 		System.out.println("xmax: " + xmax);
 		System.out.println("NumFrames: " + NumFrames);
 		System.out.println("dx: " + dx);
 		System.out.println("x1: " + x1);
 		for(int i = 0; i < NumFrames; i++)
 			System.out.println("Frame "+(i+1)+": "+samples1[i]);
 	}

}
