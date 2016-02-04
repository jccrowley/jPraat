package jPraat;


public final class WaveFormat  
{

  public static final String ChunkID = "RIFF";
  public static final String Format = "WAVE";
  public static final String Subchunk1ID = "fmt ";
  public static final String Subchunk2ID = "data";
  public static final int Subchunk1Size = 16;
  public static final int AudioFormat = 1;
  public static final int NumChannels = 1;
  public static final int BlockAlign = 2;
  public static final int BitsPerSample = 16;
  
  public static final byte[] metaData = { 82, 73, 70, 70,	 // ChunkID = "RIFF"
	  									  0, 0, 0, 0,		 // ChunkSize depends (36 + NumSamples*2)
	  									  87, 65, 86, 69,	 // Format = "WAVE"
	  									  102, 109, 116, 32, // Subchunk1ID = "fmt "
	  									  16, 0, 0, 0,		 // Subchunk1Size = 16
	  									  1, 0,				 // AudioFormat = 1
	  									  1, 0,				 // NumChannels = 1
	  									  0, 0, 0, 0,		 // SampleRate depends (1/dx)
	  									  0, 0, 0, 0,		 // ByteRate depends (2/dx)
	  									  2, 0,				 // BlockAlign = 2
	  									  16, 0,			 // BitsPerSample = 16
	  									  100, 97, 116, 97,  // Subchunk2ID = "data"
	  									  0, 0, 0, 0 };		 // Subchunk2Size depends (NumSamples*2)
  


  private WaveFormat()
  {
    throw new AssertionError();
  }
  
}
