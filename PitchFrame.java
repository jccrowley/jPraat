package jPraat;


public class PitchFrame
{
    int nCandidates;
    double intensity;
    PitchCandidate[] candidate;
    //int maxnCandidates;
    
    public PitchFrame(int n)
    {
    	nCandidates = n;
    	candidate = new PitchCandidate[n];
    	for(int i = 0; i < nCandidates; i++)
    		candidate[i] = new PitchCandidate();
    }
    
    public void printInfo()
    {
        log("intensity: " + intensity);
        log("nCandidates: " + nCandidates);
        for(int i = 0; i < nCandidates; i++)
            candidate[i].printInfo();
    }
    
    void log(Object aMsg)
    {
        System.out.println(String.valueOf(aMsg));
    }
      
    /*public void addCandidate(double f, double s)
    {
        PitchCandidate c = new PitchCandidate(f, s);
        nCandidates++;
        candidate.add(c);
    }*/
    
}

