package jPraat;

public class RealTierPoint implements Comparable<RealTierPoint> 
{
	double number;
	double value;
	
	// RealTierPoint constructor
	public RealTierPoint(double t, double f)
	{
		number = t;
		value = f;
	}
	
	public int compareTo(RealTierPoint other) 
	{
        return (int) (number - other.number);
    }
	
	// Produce deep copy of DurationTierPoint object
	public RealTierPoint RealTierPoint_copy()
	{
		return new RealTierPoint(number, value);
	}
	
	// Print contents of DurationTierPoint
	public void print()
	{
		System.out.println("number: " + number);
 		System.out.println("value: " + value);
	}
}
