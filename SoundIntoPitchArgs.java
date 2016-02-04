package jPraat;

import fftpack.RealDoubleFFT;

public class SoundIntoPitchArgs 
{
	Wave sound;
	Pitch pitch;
	long firstFrame, lastFrame;
	double minimumPitch;
	int maxnCandidates, method;
	double voicingThreshold, octaveCost, dt_window;
	long nsamp_window, halfnsamp_window, maximumLag, nsampFFT, nsamp_period, halfnsamp_period, brent_ixmax, brent_depth;
	double globalPeak; 
	double[] window, windowR;
	boolean isMainThread;
	int[] cancelled;
	
	public SoundIntoPitchArgs (Wave sound, Pitch pitch,
			long firstFrame, long lastFrame, double minimumPitch, int maxnCandidates, int method,
			double voicingThreshold, double octaveCost,
			double dt_window, long nsamp_window, long halfnsamp_window, long maximumLag, long nsampFFT,
			long nsamp_period, long halfnsamp_period, long brent_ixmax, long brent_depth,
			double globalPeak, double[] window, double[] windowR,
			boolean isMainThread, int[] cancelled)
		{
			this.sound = sound;
			this.pitch = pitch;
			this.firstFrame = firstFrame;
			this.lastFrame = lastFrame;
			this.minimumPitch = minimumPitch;
			this.maxnCandidates = maxnCandidates;
			this.method = method;
			this.voicingThreshold = voicingThreshold;
			this.octaveCost = octaveCost;
			this.dt_window = dt_window;
			this.nsamp_window = nsamp_window;
			this.halfnsamp_window = halfnsamp_window;
			this.maximumLag = maximumLag;
			this.nsampFFT = nsampFFT;
			this.nsamp_period = nsamp_period;
			this.halfnsamp_period = halfnsamp_period;
			this.brent_ixmax = brent_ixmax;
			this.brent_depth = brent_depth;
			this.globalPeak = globalPeak;
			this.window = window;
			this.windowR = windowR;
			this.isMainThread = isMainThread;
			this.cancelled = cancelled;
		}
	
	public void Sound_into_Pitch()
	{
		RealDoubleFFT fft;
		double[][] frame;
		double[] ac, r, localMean;
		long[] imax;
		
		{// scope
			//MelderThread_LOCK (mutex);
			
			  // autocorrelation
			fft = new RealDoubleFFT((int) nsampFFT);
			frame = new double[sound.NumChannels][(int) nsampFFT];
			ac = new double[(int) nsampFFT];
			
			r = new double[(int) ((nsamp_window * 2) + 1)];
			//r.reset (- my nsamp_window, my nsamp_window);
			imax = new long[maxnCandidates];
			localMean = new double[sound.NumChannels];
			//MelderThread_UNLOCK (mutex);
		}
                //System.out.println("right here");
		for (long iframe = firstFrame-1; iframe < lastFrame; iframe ++) 
		{
			PitchFrame pitchFrame = pitch.frame[(int) iframe];
			double t = pitch.Sampled_indexToX(iframe);
			
			sound.Sound_into_PitchFrame (pitchFrame, t,
				minimumPitch, maxnCandidates, method, voicingThreshold, octaveCost,
				fft, dt_window, nsamp_window, halfnsamp_window,
				maximumLag, nsampFFT, nsamp_period, halfnsamp_period,
				brent_ixmax, brent_depth, globalPeak,
				frame, ac, window, windowR,
				r, imax, localMean, iframe);
		}
	}

}
