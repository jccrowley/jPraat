package jPraat;

public class PitchCandidate 
{
    double frequency;
    double strength;
    
    // PitchCandidate constructor
    public PitchCandidate(double f, double s)
    {
    	frequency = f;
    	strength = s;
    }
    
    // Empty constructor
    public PitchCandidate() {}
    
    // Produce a deep copy of a PitchCandidate object
    public PitchCandidate copy()
    {
    	return new PitchCandidate(frequency, strength);
    }
    
    // Replace this PitchCandidate with another (maintain existing pointers)
    public void set(PitchCandidate pc)
    {
    	frequency = pc.frequency;
    	strength = pc.strength;
    }
    
    public void printInfo()
    {
        log("frequency: " + frequency);
        log("strength: " + strength);
    }
    
    void log(Object aMsg)
    {
        System.out.println(String.valueOf(aMsg));
    }
}
