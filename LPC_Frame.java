package jPraat;

public class LPC_Frame 
{
	int nCoefficients;
	double[] a;
	double[] gain;

	
	public LPC_Frame(int nCoeff)
	{
		nCoefficients = nCoeff;
		a = new double[nCoeff];
		gain = new double[1];
	}
	
	public void printInfo()
	{
		System.out.println("nCoefficients: "+nCoefficients);
		System.out.println("gain: "+gain[0]);
		for(int i = 0; i < nCoefficients; i++)
			System.out.println("Thingy "+i+": "+a[i]);
	}

}
