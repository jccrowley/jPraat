package jPraat;


public class Manipulation
{
	double xmin, xmax;
	DurationTier duration;
	Wave sound;
	PitchTier pitch;
	PointProcess pulses;
	
	public Manipulation(double tmin, double tmax)
	{
		xmin = tmin;
		xmax = tmax;
		duration = new DurationTier(tmin, tmax);
	}
	
	Wave synthesize_overlapAdd () 
	{
		Wave thee = sound.Sound_Point_Pitch_Duration_to_Sound (pulses, pitch, duration, 0.02000000001);
		return thee;
	}
}
