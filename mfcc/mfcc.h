#include <algorithm>
#include <numeric>
#include <complex>
#include <vector>
#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdint>

// vector of doubles (useful for single MFCC frame)
typedef std::vector<double> v_d;
// Complex double (useful for FFT, frequency domain results) 
typedef std::complex<double> c_d;
// matrix of doubles (usefult for MFCC output, filter bank matrix)
typedef std::vector<v_d> m_d;
// vector of complex doubles (useful for complex FFT output from a single frame)
typedef std::vector<c_d> v_c_d;
// 2D sparse map of complex numbers (useful for sparse spectrogram, storing selected FFT bins per frame)
typedef std::map<int,std::map<int,c_d> > twmap;

struct wavHeader {
    /* RIFF Chunk Descriptor */
    uint8_t         RIFF[4];        // RIFF Header Magic header
    uint32_t        ChunkSize;      // RIFF Chunk Size
    uint8_t         WAVE[4];        // WAVE Header
    /* "fmt" sub-chunk */
    uint8_t         fmt[4];         // FMT header
    uint32_t        Subchunk1Size;  // Size of the fmt chunk
    uint16_t        AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw,257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
    uint16_t        NumOfChan;      // Number of channels 1=Mono 2=Stereo
    uint32_t        SamplesPerSec;  // Sampling Frequency in Hz
    uint32_t        bytesPerSec;    // bytes per second
    uint16_t        blockAlign;     // 2=16-bit mono, 4=16-bit stereo
    uint16_t        bitsPerSample;  // Number of bits per sample
    /* "data" sub-chunk */
    uint8_t         Subchunk2ID[4]; // "data"  string
    uint32_t        Subchunk2Size;  // Sampled data length
};

class MFCC
{
  public:

    // constructor
    MFCC(int sampFreq=48000, int nCep=12, int winLength=25, int framShift=10, int numFilt=40, double lf=50, double hf=6500);

    // Process each frame and extract MFCC
    std::string processFrame(int16_t* samples, size_t N);

    // Read input file stream, extract MFCC's and write to output file stream
    int process(std::ifstream &wavFp, std::ofstream &mfcFp);
    
  private:

  // Private Functions
    
    // Hertz to Mel Conversion
    inline double hz2mel (double f)
    {
        return 2595*std::log10 (1+f/700);   
    }

    // Mel to Hertz Conversion
    inline double mel2hz (double m)
    {
        return 700*(std::pow(10,m/2595)-1);   
    }

    // Twiddle factor computation
    void compTwiddle(void);    

    // Cooley-Tukey DIT_FFT recursive function
    v_c_d fft(v_c_d x);

    /////////////// Frame Processing routines ////////////////////////
    
    // Pre-emphasis and hamming window
    void preEmphHam(void);

    // Power spectrum computation
    void computePowerSpec(void);
    
    // Apply Log Mel Filterbank (LMFB)
    void applyLMFB(void);

    // Computing discrete cosine transform (DCT)
    void applyDCT(void);

    /////////////// Initialization Routines //////////////////////////

    // Precomputing Hamm Window and dct matrix
    void initHamDCT(void);
    
    // Precompute filterbank
    void initFilterbank();
    
    // Convert vector of double to string (for writing MFCC file output)
    std::string v_d_to_string(v_d vec);


  // Private Constants

    const double PI = 4*atan(1.0);

  // Private Variables

    // sampling frequency
    int fs; 
    // 
    twmap twiddle;

    // length of the analysis window in seconds
    size_t winlength;
    // How much to shift the frame between frames
    size_t frameShiftSamples;
    // size of the FFT
    size_t numFFT;
    // number of filters
    size_t numFilters;
    // number of Cepstra
    size_t numCepstra;
    // number of FFT bins
    size_t numFFTBins;
    
    // Pre-emphasis ceofficient. 0 is no filter
    double preEmph;
    // lowest band edge of mel filters in Hz
    double lowFreq;
    // highest band edge of mel filters in Hz
    double highFreq;

    // a single frame
    v_d frame;
    // Coefficient for the Power Spectra
    v_d powerSpectralCoef;
    // Log Mel Filterbank Coefficient
    v_d lmfbCoef;
    // vector for hamming function
    v_d hamming;
    // output for Mel Frequency Cepstral Coefficient
    v_d mfcc;
    // previous samples
    v_d prevsamples;
    
    // matrix for filterbank
    m_d fbank;
    // matrix for the Discrete Cosine Transform
    m_d dct;

   
};