package jPraat;


public class Pitch
{
    double ceiling; //max pitch
    int maxnCandidates;
    long NumFrames; // Number frames
    double xmin;
    double xmax;
    double dx;
    double x1;
    PitchFrame[] frame;
    
    final double NUMlog2e = 1.44269504089;
    final double NUMundefined = Double.POSITIVE_INFINITY;
    
    // Pitch constructor
    public Pitch(double tmin, double tmax, int nt, double dt, double t1, double ceiling0, int maxnCandidates0)
    {   
        xmin = tmin;
        xmax = tmax;
        NumFrames = nt;
        dx = dt;
        x1 = t1;
		ceiling = ceiling0;
		maxnCandidates = maxnCandidates0;
		frame = new PitchFrame[(int) NumFrames];
                
		 //Put one candidate in every frame (unvoiced, silent). 
		for (int i = 0; i < nt; i++) 
		{
            frame[i] = new PitchFrame(1);
		}
    }
    
    int Pitch_getMaxnCandidates () 
    {
    	int result = 0;
    	for (long i = 0; i < NumFrames; i ++) 
    	{
    		int nCandidates = frame[(int) i].nCandidates;
    		if (nCandidates > result) result = nCandidates;
    	}
    	return result;
    }
    
    @SuppressWarnings("unused")
	void Pitch_pathFinder (double silenceThreshold, double voicingThreshold,
    		double octaveCost, double octaveJumpCost, double voicedUnvoicedCost,
    		double ceiling, boolean pullFormants)
	{
		long maxnCandidates = Pitch_getMaxnCandidates();
		long place;
		double maximum, value;
		double ceiling2 = pullFormants ? 2 * ceiling : ceiling;
		/* Next three lines 20011015 */
		double timeStepCorrection = 0.01 / dx;
		octaveJumpCost *= timeStepCorrection;
		voicedUnvoicedCost *= timeStepCorrection;

		this.ceiling = ceiling;
		double[][] delta = new double[(int) NumFrames][(int) maxnCandidates];
		long[][] psi = new long[(int) NumFrames][(int) maxnCandidates];

		for (long iframe = 0; iframe < NumFrames; iframe ++) 
		{
			PitchFrame frame1 = frame[(int) iframe];
			double unvoicedStrength = silenceThreshold <= 0 ? 0 : 2 - frame1.intensity / (silenceThreshold / (1 + voicingThreshold));
			unvoicedStrength = voicingThreshold + (unvoicedStrength > 0 ? unvoicedStrength : 0);
			for (long icand = 0; icand < frame1.nCandidates; icand ++) 
			{
				PitchCandidate candidate = frame1.candidate[(int) icand];
				boolean voiceless = candidate.frequency == 0 || candidate.frequency > ceiling2;
				delta [(int) iframe] [(int) icand] = voiceless ? unvoicedStrength : candidate.strength - octaveCost * Math.log(ceiling / candidate.frequency) * NUMlog2e; // log (x) * NUMlog2e
			}
		}

		/* Look for the most probable path through the maxima. */
		/* There is a cost for the voiced/unvoiced transition, */
		/* and a cost for a frequency jump. */

		for (long iframe = 1; iframe < NumFrames; iframe ++) 
		{
			PitchFrame prevFrame = frame[(int) (iframe - 1)], curFrame = frame[(int) iframe];
			double[] prevDelta = delta [(int) (iframe - 1)], curDelta = delta [(int) iframe];
			long[] curPsi = psi[(int) iframe];
			for (long icand2 = 0; icand2 < curFrame.nCandidates; icand2 ++) 
			{
				double f2 = curFrame.candidate[(int) icand2].frequency;
				maximum = -1e30;
				place = 0;
				for (long icand1 = 0; icand1 < prevFrame.nCandidates; icand1 ++) 
				{
					double f1 = prevFrame.candidate[(int) icand1].frequency;
					double transitionCost;
					boolean previousVoiceless = f1 <= 0 || f1 >= ceiling2;
					boolean currentVoiceless = f2 <= 0 || f2 >= ceiling2;
					if (currentVoiceless) 
					{
						if (previousVoiceless) 
						{
							transitionCost = 0;   // both voiceless
						}
						else 
						{
							transitionCost = voicedUnvoicedCost;   // voiced-to-unvoiced transition
						}
					}
					else 
					{
						if (previousVoiceless) 
						{
							transitionCost = voicedUnvoicedCost;   // unvoiced-to-voiced transition
							if (/*Melder_debug == 30*/ false) 
							{
								/*
								 * Try to take into account a frequency jump across a voiceless stretch.
								 */
								long place1 = icand1;
								for (long jframe = iframe - 1; jframe >= 0; jframe --) 
								{
									place1 = psi [(int) (jframe + 1)] [(int) place1];
									f1 = frame[(int) jframe].candidate[(int) place1].frequency;
									if (f1 > 0 && f1 < ceiling) 
									{
										transitionCost += octaveJumpCost * Math.abs(NUMlog2e * (f1 / f2)) / (iframe - jframe);
										break;
									}
								}
							}
						}
						else 
						{
							transitionCost = octaveJumpCost * Math.abs(NUMlog2e * (f1 / f2));   // both voiced
						}
					}
					
					
					value = prevDelta [(int) icand1] - transitionCost + curDelta [(int) icand2];
					//if (Melder_debug == 33) Melder_casual ("Frame %ld, current candidate %ld (delta %g), previous candidate %ld (delta %g), "
					//	"transition cost %g, value %g, maximum %g", iframe, icand2, curDelta [icand2], icand1, prevDelta [icand1], transitionCost, value, maximum);
					if (value > maximum) 
					{
						maximum = value;
						place = icand1;
					} 
					else if (value == maximum) 
					{
						//if (Melder_debug == 33) Melder_casual ("A tie in frame %ld, current candidate %ld, previous candidate %ld", iframe, icand2, icand1);
					}
				}
				curDelta [(int) icand2] = maximum;
				curPsi [(int) icand2] = place;
			}
		}

		/* Find the end of the most probable path. */

		place = 0;
		maximum = delta [(int) (NumFrames - 1)] [(int) place];
		for (long icand = 1; icand < frame[(int) NumFrames-1].nCandidates; icand++) 
		{
			if (delta[(int) (NumFrames - 1)][(int) icand] > maximum) 
			{
				place = icand;
				maximum = delta [(int) (NumFrames - 1)][(int) place];
			}
		}

		/* Backtracking: follow the path backwards. */

		for (long iframe = NumFrames - 1; iframe >= 0; iframe --) 
		{
			//if (Melder_debug == 33) Melder_casual ("Frame %ld: swapping candidates 1 and %ld", iframe, place);
			PitchFrame frame1 = frame[(int) iframe];
			PitchCandidate help = frame1.candidate[0].copy();
			frame1.candidate[0] = frame1.candidate[(int) place];
			frame1.candidate[(int) place] = help;
			place = psi [(int) iframe][(int) place];   // This assignment is challenging to CodeWarrior 11.
		}

		/* Pull formants: devoice frames with frequencies between ceiling and ceiling2. */

		if (ceiling2 > ceiling) 
		{
			//if (Melder_debug == 33) Melder_casual ("Pulling formants...");
			for (long iframe = NumFrames - 1; iframe >= 0; iframe --) 
			{
				PitchFrame frame1 = frame[(int) iframe];
				PitchCandidate winner = frame1.candidate[0];
				double f = winner.frequency;
				if (f > ceiling && f <= ceiling2) 
				{
					for (long icand = 1; icand < frame1.nCandidates; icand ++) 
					{
						PitchCandidate loser = frame1.candidate[(int) icand];
						if (loser.frequency == 0.0) 
						{
							PitchCandidate help = winner.copy();
							winner.set(loser);
							loser.set(help);
							break;
						}
					}
				}
			}
		}
	}
    
    double Sampled_indexToX (long index) 
    { 
            return x1 + index * dx;
    }
    
    double Sampled_xToIndex (double x) 
    { 
            return (x - x1) / dx /*+ 1.0*/; 
    }
    
    long Sampled_xToHighIndex(double x) 
    { 
            return (long) Math.ceil((x - x1) / dx /*+ 1.0*/); 
    }
    
    double Sampled_getValueAtSample (long isamp, long ilevel, int unit) 
    {
    	if (isamp < 0 || isamp >= NumFrames) return NUMundefined;
    	return v_getValueAtSample(isamp, ilevel, unit);
    }
    
    double Pitch_getValueAtTime (double time, int unit, int interpolate) 
    {
    	return Sampled_getValueAtX(time, 1, unit, interpolate);
    }
    
    // Praat helper function
 	long Sampled_xToNearestIndex (double x) 
 	{ 
 		return (long) Math.round((x - x1) / dx /*+ 1.0*/); 
 	}
    
    double Sampled_getValueAtX (double x, long ilevel, int unit, int interpolate) 
    {
    	if (x < xmin || x > xmax) return NUMundefined;
    	if (interpolate == 1) 
    	{
    		double ireal = Sampled_xToIndex (x);
    		long ileft = (long) Math.floor(ireal), inear, ifar;
    		double phase = ireal - ileft;
    		if (phase < 0.5) 
    		{
    			inear = ileft;
    			ifar = ileft + 1;
    		} 
    		else 
    		{
    			ifar = ileft;
    			inear = ileft + 1;
    			phase = 1.0 - phase;
    		}
    		if (inear < 0 || inear >= NumFrames) return NUMundefined;   // x out of range?
    		double fnear = v_getValueAtSample (inear, ilevel, unit);
    		if (fnear == NUMundefined) return NUMundefined;   // function value not defined?
    		if (ifar < 0 || ifar >= NumFrames) return fnear;   // at edge? Extrapolate
    		double ffar = v_getValueAtSample (ifar, ilevel, unit);
    		if (ffar == NUMundefined) return fnear;   // neighbour undefined? Extrapolate
    		return fnear + phase * (ffar - fnear);   // interpolate
    	}
    	return Sampled_getValueAtSample (Sampled_xToNearestIndex(x), ilevel, unit);
    }
    
    double v_getValueAtSample (long iframe, long ilevel, int unit) 
    {
    	double f = frame [(int) iframe].candidate[0].frequency;
    	if (f <= 0.0 || f >= ceiling) return NUMundefined;   // frequency out of range (or NUMundefined)? Voiceless
    	return ilevel == 1 ? f : frame [(int) iframe].candidate [0]. strength;
    }
    
    boolean Pitch_isVoiced_i (long iframe) 
    {
    	if(Sampled_getValueAtSample (iframe, 1, 0) != NUMundefined)
    		return true;
    	else
    		return false;
    }
    
    boolean Pitch_isVoiced_t (double time) 
    {
    	if(Sampled_getValueAtX (time, 1, 0, 0) != NUMundefined)
    		return true;
    	else
    		return false;
    }
    
    int Pitch_getVoicedIntervalAfter (double after, double[] tleft, double[] tright) 
    {
    	long ileft = Sampled_xToHighIndex (after), iright;
    	if (ileft >= NumFrames) return 0;   /* Offright. */
    	if (ileft < 0) ileft = 0;   /* Offleft. */

    	/* Search for first voiced frame. */
    	for (; ileft < NumFrames; ileft ++)
    		if (Pitch_isVoiced_i (ileft)) break;
    	if (ileft > NumFrames) return 0;   /* Offright. */

    	/* Search for last voiced frame. */
    	for (iright = ileft; iright < NumFrames; iright ++)
    		if (! Pitch_isVoiced_i (iright)) break;
    	iright --;

    	tleft[0] = Sampled_indexToX (ileft) - 0.5 * dx;   // the whole frame is considered voiced
    	tright[0] = Sampled_indexToX (iright) + 0.5 * dx;
    	if (tleft[0] >= xmax - 0.5 * dx) return 0;
    	if (tleft[0] < xmin) tleft[0] = xmin;
    	if (tright[0] > xmax) tright[0] = xmax;
    	return 1;
    }
    
    PitchTier Pitch_to_PitchTier () 
    {
		PitchTier thee = new PitchTier (xmin, xmax);
		for (long i = 0; i < NumFrames; i++) 
		{
			double frequency = frame[(int) i].candidate[0].frequency;
			/*
			 * Count only voiced frames.
			 */
			if (frequency > 0.0 && frequency < ceiling) 
			{
				double time = Sampled_indexToX(i);
				thee.addPoint(time, frequency);
			}
		}
		return thee;
    }
    
    PointProcess Pitch_to_PointProcess ()
    {
		PitchTier pitchTier = Pitch_to_PitchTier();
		PointProcess point = PitchTier_Pitch_to_PointProcess(pitchTier);
		return point;
    }
    
    PointProcess PitchTier_Pitch_to_PointProcess (PitchTier me) 
    {
		PointProcess fullPoint = me.PitchTier_to_PointProcess();
		PointProcess thee = new PointProcess(xmin, xmax, fullPoint.maxnt);
		/*
		 * Copy only voiced parts to result.
		 */
		for (long i = 0; i < fullPoint.nt; i ++) 
		{
			double t = fullPoint.t.get((int) i);
			if (Pitch_isVoiced_t (t)) 
			{
				thee.addPoint(t);
			}
		}
		return thee;
    }
    
    Pitch Pitch_scaleTime_old (double scaleFactor) 
    {
		double dx0 = dx, x10 = x1, xmax0 = xmax;
		if (scaleFactor != 1) 
		{
			dx0 = dx * scaleFactor;
			x10 = xmin + 0.5 * dx;
			xmax0 = xmin + NumFrames * dx;
		}
		Pitch thee = new Pitch(xmin, xmax0, (int) NumFrames, dx0, x10, ceiling, 2);

		for (long i = 0; i < NumFrames; i++) 
		{
			double f = frame[(int) i].candidate[0].frequency;
			thee.frame[(int) i].candidate[0].strength = frame[(int) i].candidate[0].strength;
			f /= scaleFactor;
			if (f < ceiling) 
			{
				thee.frame[(int) i].candidate[0].frequency = f;
			}
		}
		return thee;	
    }
    
    double Pitch_getQuantile (double tmin, double tmax, double quantile, int unit) 
    {
    	double value = Sampled_getQuantile (tmin, tmax, quantile, 1, unit);
    	if (value <= 0.0) 
    	{
    		value = NUMundefined;
    	}
    	return value;
    }
    
    double Sampled_getQuantile (double xminP, double xmaxP, double quantile, long ilevel, int unit) 
    {
		double[] values = new double[(int) NumFrames];
		double[] xmin0 = {xminP};
		double[] xmax0 = {xmaxP};
		Function_unidirectionalAutowindow (xmin0, xmax0);
		if (! Function_intersectRangeWithDomain (xmin0, xmax0))
        {
            return NUMundefined;
        }
		long[] imin = new long[1], imax = new long[1]; 
		long numberOfDefinedSamples = 0;
		Sampled_getWindowSamples(xmin, xmax, imin, imax);
		for (long i = imin[0]; i <= imax[0]; i ++) 
		{
			double value = v_getValueAtSample (i, ilevel, unit);
			if (value != NUMundefined) 
			{
				values [(int) ++ numberOfDefinedSamples] = value;
			}
		}
		double result = NUMundefined;
		if (numberOfDefinedSamples >= 1) 
		{
			NUMsort_d (numberOfDefinedSamples, values);
			result = NUMquantile (numberOfDefinedSamples, values, quantile);
		}
		return result;
    }
    
    double Pitch_getMinimum (double tmin, double tmax, int unit, int interpolate) 
    {
    	double[] minimum = new double[1];
    	Pitch_getMinimumAndTime (tmin, tmax, unit, interpolate, minimum, null);
    	return minimum[0];
    }
    
    void Pitch_getMinimumAndTime (double tmin, double tmax, int unit, int interpolate,
    		double[] return_minimum, double[] return_timeOfMinimum)
	{
		Sampled_getMinimumAndX (tmin, tmax, 1, unit, interpolate, return_minimum, return_timeOfMinimum);
		if (return_minimum != null && return_minimum[0] <= 0.0)
		{
			return_minimum[0] = NUMundefined;   /* Not so unlikely. */
		}
	}
    
    void Sampled_getMinimumAndX (double xmin, double xmax, long ilevel, int unit, int interpolate,
    		double[] return_minimum, double[] return_xOfMinimum)
	{
    	boolean end = false;
		long[] imin = new long[1], imax = new long[1];
		long i;
		double[] xmin0 = {xmin};
		double[] xmax0 = {xmax};
		double minimum = 1e301, xOfMinimum = 0.0;
		if (xmin == NUMundefined || xmax == NUMundefined) 
		{
			minimum = xOfMinimum = NUMundefined;
			end = true;
		}
		if(!end)
			Function_unidirectionalAutowindow (xmin0, xmax0);
		if (! Function_intersectRangeWithDomain (xmin0, xmax0) && !end) 
		{
			minimum = xOfMinimum = NUMundefined;   // requested range and logical domain do not intersect
			end = true;
		}
		if(!end)
		{
			if (0 == Sampled_getWindowSamples (xmin0[0], xmax0[0], imin, imax)) 
			{
				/*
				 * No sample centres between xmin and xmax.
				 * Try to return the lesser of the values at these two points.
				 */
				double fleft = Sampled_getValueAtX (xmin, ilevel, unit, interpolate);
				double fright = Sampled_getValueAtX (xmax, ilevel, unit, interpolate);
				if (NUMundefined != fleft && fleft < minimum) 
				{
					minimum = fleft;
					xOfMinimum = xmin;
				}
				if (NUMundefined != fright && fright < minimum)
				{
					minimum = fright;
					xOfMinimum = xmax;
				}
			} 
			else
			{
				for (i = imin[0]; i <= imax[0]; i ++) 
				{
					double fmid = v_getValueAtSample (i, ilevel, unit);
					if (fmid == NUMundefined) continue;
					if (interpolate == 0) 
					{
						if (fmid < minimum) 
						{
							minimum = fmid;
							xOfMinimum = i;
						}
					} 
					else 
					{
						/*
						 * Try an interpolation, possibly even taking into account a sample just outside the selection.
						 */
						double fleft = i < 1 ? NUMundefined : v_getValueAtSample (i - 1, ilevel, unit);
						double fright = i >= NumFrames - 1 ? NUMundefined : v_getValueAtSample (i + 1, ilevel, unit);
						if (fleft == NUMundefined || fright == NUMundefined) 
						{
							if (fmid < minimum) 
							{
								minimum = fmid;
								xOfMinimum = i;
							}
						} 
						else if (fmid < fleft && fmid <= fright) 
						{
							double[] y = new double[4], i_real = new double[1];
							double localMinimum;
							y [1] = fleft;
							y [2] = fmid;
							y [3] = fright;
							Wave wave = new Wave();
							localMinimum = wave.NUMimproveMinimum (y, 3, 2, 1, i_real);
							if (localMinimum < minimum)
							{
								minimum = localMinimum;
								xOfMinimum = i_real[0] + i - 2;
							}
						}
					}
				}
				xOfMinimum = this.x1 + (xOfMinimum - 1) * this.dx;   /* From index plus phase to time. */
				/* Check boundary values. */
				if (interpolate != 0) 
				{
					double fleft = Sampled_getValueAtX (xmin, ilevel, unit, 1);
					double fright = Sampled_getValueAtX (xmax, ilevel, unit, 1);
					if (NUMundefined != fleft && fleft < minimum) 
					{
						minimum = fleft;
						xOfMinimum = xmin;
					}
					if (NUMundefined != fright && fright < minimum)
					{
						minimum = fright;
						xOfMinimum = xmax;
					}
				}
				if (xOfMinimum < xmin) xOfMinimum = xmin;
				if (xOfMinimum > xmax) xOfMinimum = xmax;
			}
		}
		if (minimum == 1e301 && !end) minimum = xOfMinimum = NUMundefined;
	
		if (return_minimum != null) return_minimum[0] = minimum;
		if (return_xOfMinimum != null) return_xOfMinimum[0] = xOfMinimum;
	}

	double NUMquantile (long n, double a [], double factor) 
	{
		double place = factor * n + 0.5;
		long left = (long) Math.floor(place);
		if (n < 1) return 0.0;
		if (n == 1) return a [1];
		if (left < 1) left = 1;
		if (left >= n) left = n - 1;
		if (a [(int) (left + 1)] == a [(int) left]) return a [(int) left];
		return a [(int) left] + (place - left) * (a [(int) (left + 1)] - a [(int) left]);
	}
    
    void NUMsort_d (long n, double a [])
    {
     
    	long l, r, j, i; 
    	double k; 
    	if (n < 2) return;   /* Already sorted. */ 
    	/* This n<2 step is absent from Press et al.'s implementation, */ 
    	/* which will therefore not terminate on if(--ir==1). */ 
    	/* Knuth's initial assumption is now fulfilled: n >= 2. */ 
    	if (n == 2) 
    	{ 
    		if (a [1] > a [2]) { double min = a [2]; a [2] = a [1]; a [1] = min; } 
    			return; 
    	}	 
    	if (n <= 12) 
    	{ 
    		for (i = 1; i < n; i ++) 
    		{ 
    			double min = a [(int) i]; 
    			long pmin = i; 
    			for (j = i + 1; j <= n; j ++) if (a [(int) j] < min) 
    			{ 
    				min = a [(int) j]; 
    				pmin = j; 
    			} 
    			a [(int) pmin] = a [(int) i]; 
    			a [(int) i] = min; 
    		} 
    	return; 
    	} 
    	l = (n >> 1) + 1; 
    	r = n; 
    	for (;;) 
    	{ 
    		if (l > 1) 
    		{ 
    			l --; 
    			k = a [(int) l]; 
    		} 
    		else /* l == 1 */ 
    		{ 
    			k = a [(int) r]; 
    			a [(int) r] = a [1]; 
    			r --; 
    			if (r == 1) { a [1] = k; return; } 
    		} 
    		j = l; 
    		for (;;) 
    		{ 
    			i = j; 
    			j = j << 1; 
    			if (j > r) break; 
    			if (j < r && a [(int) j] < a [(int) (j + 1)]) j ++; 
    			if (k >= a [(int) j]) break; 
    			a [(int) i] = a [(int) j]; 
    		} 
    		a [(int) i] = k; 
    	} 
    }
    
    long Sampled_getWindowSamples (double xmin0, double xmax0, long[] ixmin, long[] ixmax) 
	{
		double rixmin = Math.ceil((xmin0 - x1) / dx);
		double rixmax = Math.floor((xmax0 - x1) / dx);   // could be above 32-bit LONG_MAX
		ixmin[0] = rixmin < 0.0 ? 0 : (long) rixmin;
		ixmax[0] = rixmax > (double) (NumFrames-1) ? (NumFrames-1) : (long) rixmax;
		if (ixmin[0] > ixmax[0]) return 0;
		return ixmax[0] - ixmin[0] + 1;
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
    
    long Pitch_countVoicedFrames () 
    {
    	return Sampled_countDefinedSamples(1, 0);
    }
    
    long Sampled_countDefinedSamples (long ilevel, int unit) 
    {
    	long numberOfDefinedSamples = 0;
    	for (long isamp = 0; isamp < NumFrames; isamp ++) 
    	{
    		double value = v_getValueAtSample (isamp, ilevel, unit);
    		if (value == NUMundefined) continue;
    		numberOfDefinedSamples += 1;
    	}
    	return numberOfDefinedSamples;
    }
    
    void printInfo()
    {
        System.out.println();

        log("xmin: " + xmin);
        log("xmax: " + xmax);
        log("NumFrames: " + NumFrames);
        log("dx: " + dx);
        log("x1: " + x1);
        log("ceiling: " + ceiling);
        log("maxnCandidates: " + maxnCandidates);
        frame[0].printInfo();
        System.out.println();
    }
	
    void log(Object aMsg)
    {
        System.out.println(String.valueOf(aMsg));
    }
    
}