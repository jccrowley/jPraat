package jPraat;

import java.util.Random;

import jPraat.PointProcess;

//Class representing praat Matrix objects.
public class Matrix
{	
	final double step = 0.000001;
	
	// All members are named the same as they are in praat.
    int ny;
    int nx;
    double[][] z;
    double xmin;
    double xmax;
    double dx;
    double x1;
    
    double ymin;
    double ymax;
    double dy;
    double y1;
    
    public Matrix(double xmin0, double xmax0, int nx0, double dx0, double x10, double ymin0, double ymax0, int ny0, double dy0, double y10)
    {
    	ny = ny0;
    	nx = nx0;
    	xmin = xmin0;
    	xmax = xmax0;
    	ymin = ymin0;
    	ymax = ymax0;
    	dy = dy0;
    	dx = dx0;
    	x1 = x10;
    	y1 = y10;
    	z = new double[ny][nx];
    	
    }
    
    PointProcess Matrix_to_PointProcess ()
    {
		PointProcess thee = new PointProcess(z [0] [0], z [0] [nx-1], nx);
		for (long i = 0; i < nx; i ++) 
		{
			thee.addPoint (z[0] [(int) i]);
		}
		return thee;
    }
    
    double NUMrandomUniform(double lo, double hi)
	{
		int num = (int) ((hi - lo) / step);
		//System.out.println("num: "+num);
		Random rand = new Random();
		int random = rand.nextInt(num);
		//System.out.println("random: "+random);
		double scaled = ((double)random) * step;
		//System.out.println("scaled: "+scaled);
		double shifted = scaled + lo;
		return shifted;
	}
    
    void addSelfRandom(double value)
    {
    	double rand;
    	for(int i = 0; i < nx; i++)
    	{
    		for(int j = 0; j < ny; j++)
    		{
    			rand = NUMrandomUniform(0, value);
    			//System.out.println(rand);
    			z[j][i] += rand;
    		}
    	}
    }
    
    public void print()
    {
    	System.out.println("xmin: " + xmin);
 		System.out.println("xmax: " + xmax);
 		System.out.println("nx: " + nx);
 		System.out.println("dx: " + dx);
 		System.out.println("x1: " + x1);
 		System.out.println("ymin: " + ymin);
 		System.out.println("ymax: " + ymax);
 		System.out.println("ny: " + ny);
 		System.out.println("dy: " + dy);
 		System.out.println("y1: " + y1);
 		
 		
 		for(int i = 0; i < ny; i++)
 			for(int j = 0; j < nx; j++)
 				System.out.println(i+", "+j+": "+z[i][j]);
    }

}
