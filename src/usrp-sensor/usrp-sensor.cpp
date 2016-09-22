/*Copyright 2013,2012 Joseph "Mitch" Davis mitchd@vt.edu
 *
 *This file is part of usrp-utils.
 *
 *   usrp-sensor is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   usrp-sensor is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with usrp-utils.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *This program is designed to leverage single-precision fftw3 libraries and
 *pthreads to compute spectral periodigrams using an applied window function
 *and the magnitude of a 1D DFT.  A USRP serves as the input data stream for
 *this program.  As of July 2012, this has only been built and tested with an
 *N210 and SBX RX/TX daughtercard.
 *
 *The commandline options are:
 *
 * -o [file]       Output File  -The output file contains raw float data
 *                               representing the computed spectral periodigram
 *                               of the recorded signal.  The peridigram is
 *                               computed with fftw3f's 1D DFT.
 *
 * -s [size]       FFT Size     -Size of the FFT to compute
 *
 * -l [overlap]    Overlap      -How many FFTs overlap in the space of a
 *                               transform.  That is, a transform of FFT Size is
 *                               computed every (FFT Size)/Overlap samples.
 *                               Sanity checking is performed to ensure that the
 *                               speficied overlap is appropriate
 *
 * -w [file]       Window File  -The window file is expected to contain raw
 *                               floats representing the window function to
 *                               apply to each FFT.  It is acceptable for the
 *                               window function to be smaller than the FFT.
 *                               If this is the case, then the program
 *                               left-justifies the window with respect to the
 *                               input array to fftw3f, and zero-pads the
 *                               remaining values.
 *
 * -c [children]   Child Threads-This specifies the number of child worker
 *                               threads to spawn for .  If your
 *                               processor has N cores (including hyperthreading
 *                               cores), it is recommended to specify N-1
 *                               children.  This allows for the parent thread
 *                               to process the input data for a total of N
 *                               threads.
 *
 * -a [args]       USRP Args    -Specify the address for the input USRP. See
 *               http://files.ettus.com/uhd_docs/manual/html/identification.html
 *                               for information about identifying USRPs.
 *
 * -f [frequency]  Center Freq  -The center frequency for the FFT process.
 *                               Future plans including adding sweep
 *                               functionality to the program.
 *
 * -r [rate]       Sample Rate  -The sample rate of the USRP.  Because of the
 *                               I-Q data from the USRP, this correlates to the
 *                               bandwidth of the FFT calculation.
 *
 * -t [time]       Runtime      -The runtime in seconds for the sensing process.
 *
 * -g [gain]       RX Gain      -Gain in DB of the rx chain
 *
 *Description of error messages:
 *
 *Need at least one child thread
 *  Worker threads spawn from the parent process.  You must specify at least
 *  one child thread to do the FFT calculations.
 *
 *Cannot open window file
 *  There was a problem opening the window file.
 *
 *Incompatible FFT Size and Overlap factor
 *  The overlap must be and integral divisor of the FFT Size.  The ambiguity
 *  shown by the error output is the modulus of the parameters.
 *
 *Window is too large!
 *  The window function size must be less-than-or-equal-to the FFT size.
 *
 *Window is smaller than FFT Size, assuming zero-padding.
 *  The program is assuming you intend to use a smaller window function for a
 *  large FFT.  This output is provided in case you accidentally provided the
 *  wrong window function.
 *
 *Error performing calculations.
 *  There was a gross error with the calculation routine.  This could be because
 *  the program ran out of memory, or there was some failure with the
 *  multithreading system
 *
 *ERROR; Cannot create message queue
 *  Is your system POSIX compliant?  This program requires the use of pthreads
 *  and POSIX message queues.  There was a problem creating the message queue.
 *  Investigate /dev/mqueue
 *
 *ERROR; return code from pthread_creat() is [xx]
 *  There was a serious problem creating child threads.  The program could have
 *  run out of system resources, or any number of reasons.
 *
 *Input data terminated with unaligned data
 *  It is entirely possible to have an input data file that is not some integral
 *  mutiple of (FFT Size)/Overlap.  If this is the case, the program runs until
 *  it doesn't have enough data to complete an FFT.  After encountering the tail
 *  of the file, the program finishes remaining calculations and exits
 *  gracefully.  This message is provided to indicate that data was lost in the
 *  FFT process (The bytes at the end of the file that didn't fit into the last
 *  FFT.
 *
 *
 *BENCHMARK MODE
 *
 *Complide with -DBENCHMARK to get Elapsed Time (ET) output of the realtime
 *used to calculate all the FFTs
 *
 *
 *
 *
 * Changelog
 *
 * 0.1 - Initial release
 * 0.2 - Changed USRP streaming mode to continuous and changed the counter datatype
 *        to accomodate large sample periods (unsigned long long int = 2^128) 2012
 * 0.3 - Changed to cmake and upload to github 20130701
 */

//Define some of the values we use to setup the USRP and FFT process
#define __USRP_CPU_FMT  "fc32"

#ifdef WIRE_SC8
  #define __USRP_WIRE_FMT "sc8"
#else
  #define __USRP_WIRE_FMT "sc16"
#endif

#define __USRP_CLK_SRC  "internal"

//USRP Headers
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <boost/thread.hpp>
#include <csignal>
#include <complex.h>

#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include <mqueue.h>
#include <unistd.h>

#include "fft_thread.h"


#ifdef BENCHMARK
  #include <ctime>
  #include <unistd.h>
  #include <sys/time.h>
#endif




using namespace std;

/*useage()
 *
 *Display program useage information
 *
 */
void useage();

/*openFiles( char*, FILE*&, char*, FILE*& )
 *
 *Open the input and output files, and perform error handling if there is a
 *catastrophic failure
 */
int openFiles( const char*  outputFileName,
               FILE*&       outputFile );

/*createFFTPlans( fftwf_plan*&, _Complex float**&, _Complex float**&,
                  int, int)
 *
 *Initialize the FFT plans.  This creates max_children plans for FFTSize FFTs
 */
void createFFTPlans( fftwf_plan*&     plans,
                     _Complex float**& inputData,
                     _Complex float**& outputData,
                     int              FFTSize,
                     int              max_children );

/*initializeThreads(...)
 *
 *Initialize and spawn all the child threads.  Returns 0 upon thread or mqueue
 *failure.
 */
int initializeThreads( pthread_mutex_t&         mutex,
                       pthread_mutexattr_t&     attr,
                       pthread_t*&              fft_children,
                       mqd_t*&                  fft_mq,
                       struct mq_attr*&         ma,
                       int                      max_children,
                       struct fft_thread_data*& fft_child_args,
                       FILE*&                   outputFile,
                       fftwf_plan*&             plans,
                       _Complex float**&         inputData,
                       _Complex float**&         outputData,
                       float*&                  window,
                       int                      FFTSize,
                       int*&                    thread_control );

/*calculateTask(...)
 *
 *This is the main work of the program, performing the specified overlapped
 *FFT transforms
 */
int calculateTask(  const char*                   outputFileName,
                    const int                     FFTSize,
                    const int                     FFTOverlap,
                    const int                     max_children,
                    float*                        window,
                    const unsigned long long int  maximum_samples,
                    uhd::usrp::multi_usrp::sptr&  usrp );

/*setupUSRP(...)
 *
 *Setup the USRP for receiving at the specified freq and rate
 */
int setupUSRP(  uhd::usrp::multi_usrp::sptr&  usrp,
                const float                   center_freq,
                const float                   sample_rate,
                const int                     rx_gain,
                const char*                   dev_addr);







/*******************************************************************************


*******************************************************************************/
int main( int argc, char* argv[])
{
  //First things first, try to set realtime priority for the parent thread
  uhd::set_thread_priority_safe();

  const int FLOAT_SIZE = sizeof(float); //Size of single-precision float
  //in bytes

  //Ensure the correct number of arguments were passed
  if( argc != 19)
  {
    cout << "Only " << argc << " parameters entered" << endl;
    useage();
    return -1;
  }

  char  *outputFileName = NULL;
  char  *windowFileName = NULL;
  char  *usrpArgs       = NULL;
  int   usrpGain        = 0;
  int   FFTSize         = 0;
  int   FFTOverlap      = 0;
  int   arg             = 0;
  int   max_children    = 0;
  float usrpCenterFreq  = 0.0f;
  float usrpSampleRate  = 0.0f;
  float usrpRecordTime  = 0.0f;

  //argument parsing
  while( (arg = getopt( argc, argv, "o:s:l:c:w:a:f:r:t:g:")) != -1 )
  {
    switch (arg)
    {
      case 'o':
        outputFileName = new char[strlen(optarg)+1];
        strcpy(outputFileName,optarg);
        break;
      case 's':
        FFTSize = atoi(optarg);
        break;
      case 'l':
        FFTOverlap = atoi(optarg);
        break;
      case 'c':
        max_children = atoi(optarg);
        break;
      case 'w':
        windowFileName = new char[strlen(optarg)+1];
        strcpy(windowFileName,optarg);
        break;
      case 'a':
        usrpArgs = new char[strlen(optarg)+1];
        strcpy(usrpArgs,optarg);
        break;
      case 'f':
        usrpCenterFreq = atof(optarg);
        break;
      case 'r':
        usrpSampleRate = atof(optarg);
        break;
      case 't':
        usrpRecordTime = atof(optarg);
        break;
      case 'g':
        usrpGain = atoi(optarg);
        break;
      case '?':
        useage();
        if( outputFileName )
          delete [] outputFileName;
        if( usrpArgs )
          delete [] usrpArgs;
        if( windowFileName )
          delete [] windowFileName;
        return -1;
      }
  }
  //Check FFT and Overlap compatibility
  if( FFTSize % FFTOverlap )
  {
    cout  << "Incompatible FFT Size and Overlap factor " << endl
          << "FFT Size: " << FFTSize << endl
          << "Overlap: " << FFTOverlap << endl
          << "Modulus: " << FFTSize % FFTOverlap
          << endl;
    return -1;
  }

  //Check multithreading options
  if( max_children < 1 )
  {
    cout  << "Need at least one child thread" << endl;
    return -1;
  }

  float* window = new float[FFTSize];
  FILE* window_file;
  window_file = fopen( windowFileName, "r" );
  if( !window_file )
  {
    cout << "Cannot open window file" << endl
         << "Assuming uniform window" << endl;
    for( int i = 0; i < FFTSize; i++ )
      window[i] = 1.0f;
  }
  else
  {
    //This will get the number of window elements in the specified file.
    //We want this info because its possible that our window function is shorter
    //than the FFTSize.  In which case, we need to pad the remainder of the 
    //window array with zeros
    fseek( window_file, 0L, SEEK_END );
    int window_file_size = ftell( window_file );
    fseek( window_file, 0L, SEEK_SET );

    //We only support windows that are >= FFTSize
    int window_size = window_file_size / FLOAT_SIZE;
    if( window_size > FFTSize )
    {
      cout  << "Window is too large!" << endl
            << "FFT Size: " << FFTSize << endl
            << "Window Size: " << window_size << endl;
      delete [] outputFileName;
      delete [] windowFileName;
      delete [] window;
      delete [] usrpArgs;
      fclose( window_file );
      return 1;
    }

    //Finally read the window
    fread( window, FLOAT_SIZE , window_size, window_file );
    fclose(window_file);

    //Accomadate zero-padding
    if (window_size < FFTSize )
    {
      cout  << "Window is smaller than FFT Size, assuming zero-padding."
            << endl << "FFT Size: " << FFTSize << endl
                    << "Window Size: " << window_size << endl;
      //Zero-pad the rest of the window
      for(int i = window_size; i < FFTSize; i++ )
        window[i] = 0.0f;
    }
  }
  cout << "Initializing USRP device" << endl;
  //Initialize the USRP hardware
  uhd::usrp::multi_usrp::sptr the_usrp;
  if( !setupUSRP( the_usrp,
                  usrpCenterFreq,
                  usrpSampleRate,
                  usrpGain,
                  usrpArgs ))
  {
    cout << "Error initializing the USRP device." << endl;
    delete [] outputFileName;
    delete [] windowFileName;
    delete [] window;
    delete [] usrpArgs;
    return 1;
  }
  //Perform the actual work
  if( !calculateTask( outputFileName,
                      FFTSize,
                      FFTOverlap,
                      max_children,
                      window,
                      static_cast<unsigned long long int>(usrpSampleRate*usrpRecordTime),
                      the_usrp ) )
  {
    cout << "Error performing calculations" << endl;
    delete [] outputFileName;
    delete [] windowFileName;
    delete [] window;
    delete [] usrpArgs;
    return 1;
  }

  uhd::stream_cmd_t       usrp_stream_stop(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
  usrp_stream_stop.stream_now = false;
  usrp_stream_stop.time_spec = uhd::time_spec_t();
  the_usrp->issue_stream_cmd( usrp_stream_stop );
  delete [] outputFileName;
  delete [] windowFileName;
  delete [] window;
  delete [] usrpArgs;
  return 0;
}








/*******************************************************************************


*******************************************************************************/
int setupUSRP(  uhd::usrp::multi_usrp::sptr&  usrp,
                const float                   center_freq,
                const float                   sample_rate,
                const int                     rx_gain,
                const char*                   dev_addr)
{
  //Initialize the USRP to the specified address
  usrp = uhd::usrp::multi_usrp::make(string(dev_addr));
  //Define the clock reference
  usrp->set_clock_source(__USRP_CLK_SRC);

  //Output some useful information
  cout  << "Using the following USRP device: " << endl
        << usrp->get_pp_string() << endl;

  //Try setting the sample rate.  If the rate we get is not the same as the
  //requested rate, we will return with a warning to ensure the user is aware
  //of the actual sample rate
  usrp->set_rx_rate( sample_rate );

  if( usrp->get_rx_rate() != sample_rate )
  {
    ios_base::fmtflags originalFlags = cout.flags();
    cout.setf(ios_base::left,ios_base::floatfield);
    cout.precision(15);
    cout  << "WARNING! Requested rate = " << sample_rate << endl
          << "WARNING! Actual rate = " << usrp->get_rx_rate() << endl;
    cout.setf(originalFlags);
  }

  //Try setting the center frequency.  Like above, if we get a different
  //frequency than the one we're requesting, we will spit out a warning for the
  //user
  usrp->set_rx_freq( center_freq );

  if( usrp->get_rx_freq() != center_freq )
  {
    ios_base::fmtflags originalFlags = cout.flags();
    cout.setf(ios_base::left,ios_base::floatfield);
    cout.precision(15);
    cout  << "WARNING! Requested frequency = " << center_freq << endl
          << "WARNING! Actual frequency = " << usrp->get_rx_freq() << endl;
    cout.setf(originalFlags);
  }

  //Set the RX gain.  There really shouldn't be any problems here, but the user
  //might request something silly like 50dB of gain when the module can't
  //accomodate.  So we'll perform a similar check here.
  usrp->set_rx_gain( rx_gain );

  if( usrp->get_rx_gain() != rx_gain )
  {
    cout  << "WARNING! Requested gain = " << rx_gain << endl
          << "WARNING! Actual gain = " << usrp->get_rx_gain() << endl;
  }

  //Ensure the LO locked
  vector<string> sensor_names;
  sensor_names = usrp->get_rx_sensor_names(0);
  if( find(sensor_names.begin(), sensor_names.end(), "lo_locked")
      != sensor_names.end() )
  {
    uhd::sensor_value_t lo_locked = usrp->get_rx_sensor("lo_locked",0);
    cout  << "Checking RX: " << endl
          << lo_locked.to_pp_string() << endl;
    UHD_ASSERT_THROW(lo_locked.to_bool());      //We should probably catch this
  }

  return 1;
}








/*******************************************************************************


*******************************************************************************/
void useage( )
{
  cout  << "Useage:\t USRP-Sensor [args]" << endl
        << "-o <file>\t Output File" << endl
        << "-s <size>\t FFT Size" << endl
        << "-l <number>\t FFT Overlap" << endl
        << "-c <number>\t Number of Child Processes" << endl
        << "-w <file>\t Window File" << endl
        << "-a <args>\t USRP Address" << endl
        << "-f <freq>\t USRP Center Frequency" << endl
        << "-r <rate>\t USRP Sample Rate" << endl
        << "-g <gain>\t USRP RX Gain" << endl
        << "-t <time>\t Time to record" << endl;
}







/*******************************************************************************


*******************************************************************************/
int openFiles( const char*  outputFileName,
               FILE*&       outputFile )
{
  outputFile = fopen( outputFileName, "w" );

  /*Make sure the files actually opened.
   *
   *If not... gracefully exit.
   *We should really be throwing exceptions...
   */
  if( !outputFile)
    return 0;
  return 1;
}









/*******************************************************************************


*******************************************************************************/
void createFFTPlans( fftwf_plan*&     plans,
                     _Complex float**& inputData,
                     _Complex float**& outputData,
                     int              FFTSize,
                     int              max_children )
{
  plans       = new fftwf_plan[ max_children ];
  inputData   = new _Complex float*[ max_children ];
  outputData  = new _Complex float*[ max_children ];

  //Setup the FFT plans -- create one for each child thread
  for(int i = 0; i < max_children; i++ )
  {
    inputData[i]  = new _Complex float[FFTSize];
    outputData[i] = new _Complex float[FFTSize];
    plans[i]      = fftwf_plan_dft_1d( FFTSize, (fftwf_complex*)inputData[i], 
                                       (fftwf_complex*)outputData[i],
                                       FFTW_FORWARD, FFTW_EXHAUSTIVE );
  }

  return;
}










/*******************************************************************************


*******************************************************************************/
int initializeThreads( pthread_mutex_t&         mutex,
                       pthread_mutexattr_t&     attr,
                       pthread_t*&              fft_children,
                       mqd_t*&                  fft_mq,
                       struct mq_attr*&         ma,
                       int                      max_children,
                       struct fft_thread_data*& fft_child_args,
                       FILE*&                   outputFile,
                       fftwf_plan*&             plans,
                       _Complex float**&         inputData,
                       _Complex float**&         outputData,
                       float*&                  window,
                       int                      FFTSize,
                       int*&                    thread_control )
{
  //mutex attribute = PTHREAD_MUTEX_NORMAL
  pthread_mutexattr_init( &attr );
  pthread_mutexattr_settype( &attr, PTHREAD_MUTEX_NORMAL );

  //Initialize the mutex
  pthread_mutex_init( &mutex, &attr );

  //Setup the maximum number of child threads
  fft_children    = new pthread_t[max_children];
  fft_mq          = new mqd_t[max_children];
  ma              = new mq_attr[max_children];
  thread_control  = new int;
  //Setup the arguments for the children
  fft_child_args  = new fft_thread_data[max_children];

  //Ensure that we got memory
  if( !fft_children || !fft_mq || !ma || !thread_control || !fft_child_args )
    return 0;

  for(int i = 0; i < max_children; i++ )
  {
    fft_child_args[i].outputFile      = outputFile;
    fft_child_args[i].output_mutex    = &mutex;
    fft_child_args[i].plan            = plans[i];
    fft_child_args[i].fft_size        = FFTSize;
    fft_child_args[i].inputData       = inputData[i];
    fft_child_args[i].outputData      = outputData[i];
    fft_child_args[i].is_running      = false;
    fft_child_args[i].mq_name         = new char[255];
    //Ensure that we got the memory
    if( !fft_child_args[i].mq_name )
      return 0;
    fft_child_args[i].my_id           = i;
    fft_child_args[i].window          = window;
    fft_child_args[i].next_thread     = thread_control;
    fft_child_args[i].max_children    = max_children;
    sprintf(fft_child_args[i].mq_name,"/fft_thread_%i",i);

    ma[i].mq_flags    = 0;
    ma[i].mq_maxmsg   = 1;
    ma[i].mq_msgsize  = sizeof(char);
    ma[i].mq_curmsgs  = 0;
  }

  //Initialize the child threads
  int rc = 0;
  for( int i = 0; i < max_children; i++ )
  {
    //Open the message queue
    fft_mq[i] = mq_open( fft_child_args[i].mq_name, O_WRONLY | O_CREAT,
                         0700, &(ma[i]));
    if( fft_mq[i] == -1 )
    {
      cout << "ERROR; Cannot create message queue\n" << endl;
      return 0;
    }
    rc = pthread_create( &fft_children[i], NULL, fft_thread_start,
                         reinterpret_cast<void *>(&(fft_child_args[i])) );

    if( rc )
    {
      cout << "ERROR; return code from pthread_create() is " << rc << endl;
      return 0;
    }
  }
  return 1;
}










/*******************************************************************************


*******************************************************************************/
int calculateTask(  const char*                   outputFileName,
                    const int                     FFTSize,
                    const int                     FFTOverlap,
                    const int                     max_children,
                    float*                        window,
                    const unsigned long long	  maximum_samples,
                    uhd::usrp::multi_usrp::sptr&  usrp )
{
  ///////////////////////////////////////////////////////////
  //
  //Initialization Section
  ///////////////////////////////////////////////////////////

  //Number of bytes in a single-precision float
  const int   FLOAT_COMPLEX_SIZE  = sizeof(_Complex float);
  const char  msg_thread_start    = static_cast<char>(__FFT_THREAD_START);
  const char  msg_thread_kill     = static_cast<char>(__FFT_THREAD_KILL);

  //Initialize and open the input/output files
  FILE *outputFile;

  if(!openFiles( outputFileName, outputFile ))
    return 0;

  //Create some FFT Plans
  fftwf_plan     *plans       = NULL;
  _Complex float **inputData   = NULL;
  _Complex float **outputData  = NULL;

  createFFTPlans( plans, inputData, outputData, FFTSize, max_children );

  //Setup the mutex for the output file... we don't want a race condition
  pthread_mutex_t     output_mutex;
  pthread_mutexattr_t output_mutex_attr;

  pthread_t               *fft_children   = NULL;
  mqd_t                   *fft_mq         = NULL;
  struct mq_attr          *ma             = NULL;
  struct fft_thread_data  *fft_child_args = NULL;
  int                     *thread_control = NULL;

  if( !initializeThreads( output_mutex,
                          output_mutex_attr,
                          fft_children,
                          fft_mq,
                          ma,
                          max_children,
                          fft_child_args,
                          outputFile,
                          plans,
                          inputData,
                          outputData,
                          window,
                          FFTSize,
                          thread_control ))
  {
    //Cleanup for a graceful exit
    //We have to check to see if things exist before deleting them because
    //initialization failed.  So its possible that some things exist, and others
    //do not.
    if( plans )
    {
      for(int i = 0; i < max_children; i++)
        fftwf_destroy_plan(plans[i]);
      delete [] plans;
    }
    if( inputData )
    {
      for(int i = 0; i < max_children; i++)
	if( inputData[i] )
          delete [] inputData[i];
      delete [] inputData;
    }
    if( outputData )
    {
      for(int i = 0; i < max_children; i++)
  	if( outputData[i] )
	  delete [] outputData[i];
      delete [] outputData;
    }
    if( fft_children )
    {
      for(int i = 0; i < max_children; i++)
      {
        mq_send( fft_mq[i], &msg_thread_kill,
                 __FFT_THREAD_MSG_LENGTH, __FFT_THREAD_MSG_PRIO );
        pthread_join( fft_children[i], NULL );
      }
      delete [] fft_children;
    }
    if( fft_mq )
    {
      for(int i = 0; i < max_children; i++)
      {
        mq_close( fft_mq[i] );
        mq_unlink( fft_child_args[i].mq_name );
      }
      delete [] fft_mq;
    }
    if( ma )
      delete [] ma;
    if( fft_child_args )
    {
      for(int i = 0; i < max_children; i++)
        if( fft_child_args[i].mq_name )
          delete [] fft_child_args[i].mq_name;
      delete [] fft_child_args;
    }
    return 0;
  }



  //Setup the input buffer and tracking variables
  int                   fft_interval_size = FFTSize / FFTOverlap;
  _Complex float*        input_buffer      = new _Complex float[FFTSize];
  int                   head              = 0;
  int                   child_tracker     = 0;
  bool                  isFull            = false;
  int                   return_code       = 1;
  //Check memory allocation
  if( !input_buffer )
    return_code = 0;

  //Setup the USRP for streaming
  vector<_Complex float>   usrpBuffer( fft_interval_size );
  uhd::stream_args_t      stream_args(__USRP_CPU_FMT, __USRP_WIRE_FMT );
  uhd::rx_streamer::sptr  usrp_rx_stream = usrp->get_rx_stream(stream_args);
  uhd::rx_metadata_t      rx_md;
  unsigned long long int  samples_recorded = 0;
  unsigned long long int  buffer_samples_recorded = 0;
  uhd::stream_cmd_t       usrp_stream_command(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);

  usrp_stream_command.stream_now  = true;
  usrp_stream_command.time_spec   = uhd::time_spec_t();

  ///////////////////////////////////////////////////////////
  //
  //Work Section
  ///////////////////////////////////////////////////////////
#ifdef BENCHMARK
  timeval a;
  timeval b;
#endif

  //Start at the first child thread
  (*thread_control) = 0;
  cout << "Begin Data Collection" << endl;
  //Start streaming!
  usrp->issue_stream_cmd( usrp_stream_command );

#ifdef BENCHMARK
  gettimeofday(&a, 0);
#endif

  while( (samples_recorded < maximum_samples) and return_code )
  {
    //Read in the I-Q of fft_interval_size samples...
    buffer_samples_recorded = usrp_rx_stream->recv( &usrpBuffer.front(),
                                                    usrpBuffer.size(),
                                                    rx_md );

    //Check the USRP for errors (including Overflow indication)
    if( rx_md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE )
    {
      //There was a USRP-related problem
      switch( rx_md.error_code ){
        case uhd::rx_metadata_t::ERROR_CODE_OVERFLOW:
          cout << "O";
          break;
        case uhd::rx_metadata_t::ERROR_CODE_TIMEOUT:
          cout << "USRP Timeout" << endl;
          return_code = 0;
          break;
        default:
          cout << "Unexpected USRP Error: " << rx_md.error_code;
          return_code = 0;
      }
    }

    //Ensure we grabbed the correct number of samples
    if( buffer_samples_recorded == static_cast<unsigned int>(fft_interval_size))
    {
      samples_recorded += buffer_samples_recorded;

      //copy into the fftsize-long input buffer at the proper spot
      memmove( input_buffer+head, &usrpBuffer.front(),
                FLOAT_COMPLEX_SIZE * fft_interval_size );

      //Time to take an FFT yet?
      if( !isFull )
      {
        //If we're not full, then we're still filling the input buffer for the
        //very first FFT
        if( head == (FFTSize - fft_interval_size) )
        {
          //Now we are full and ready to take an FFT
          isFull = true;
        }
      }
      //Increment the head pointer before we take the FFT so the head is pointing
      //to the oldest data in the buffer
      head += fft_interval_size;
      if( head == FFTSize )
        head = 0;
      if( isFull )
      {
        //Check to see if the thread is still active before copying data
        //This ideally shouldn't happen... but if it does, we can potentially lose
        //data from the USRP.
        if( fft_child_args[child_tracker].is_running )
        {
          while( fft_child_args[child_tracker].is_running )
          {
            ;
              //We have to block until that thread finishes
              //We can't pthread_join because the thread doesn't
              //terminate until we generate a _THREAD_KILL command
              //So instead we wait for the running flag to go false
          }
        }
        //Copy the buffer into the FFT input data using a 2-part memmove
        memmove( fft_child_args[child_tracker].inputData,
                 input_buffer+head, FLOAT_COMPLEX_SIZE * (FFTSize-head) );
        memmove( fft_child_args[child_tracker].inputData+FFTSize-head,
                 input_buffer, FLOAT_COMPLEX_SIZE * head );

        //Send a message to the worker to compute the FFT
        mq_send( fft_mq[child_tracker], &msg_thread_start,
                 __FFT_THREAD_MSG_LENGTH, __FFT_THREAD_MSG_PRIO );
        //Move the child_tracker to point to the next child thread in the pool
        child_tracker = child_tracker == max_children - 1 ? 0 : child_tracker + 1;
      }
    }
  }
#ifdef BENCHMARK
  gettimeofday(&b, 0);

  float seconds = b.tv_sec - a.tv_sec;
  seconds += b.tv_usec/1000000.0f - a.tv_usec/1000000.0f;
  printf("ET: %.6f s\n",seconds);
#endif
  cout << "End data collection" << endl;

  usrp->issue_stream_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);

  ///////////////////////////////////////////////////////////
  //
  //Cleanup Section
  ///////////////////////////////////////////////////////////

  //Toss out any leftovers and cleanup
  fclose(outputFile);

  //Destroy plans
  for(int i = 0; i < max_children; i++)
  {
    //Send a message for the child to kill itself
    mq_send( fft_mq[i], &msg_thread_kill,
             __FFT_THREAD_MSG_LENGTH, __FFT_THREAD_MSG_PRIO );
    fftwf_destroy_plan(plans[i]);
    delete [] inputData[i];
    delete [] outputData[i];
    mq_close( fft_mq[i] );
    mq_unlink( fft_child_args[i].mq_name );
    //Wait for the child to kill itself
    pthread_join( fft_children[i], NULL );
  }

  //cleanup mutex
  pthread_mutexattr_destroy( &output_mutex_attr );
  pthread_mutex_destroy( &output_mutex );

  //free more memory
  for( int i = 0; i < max_children; i++ )
    delete [] fft_child_args[i].mq_name;

  delete [] plans;
  delete [] inputData;
  delete [] outputData;
  delete [] ma;
  delete [] fft_mq;
  delete [] fft_child_args;
  delete [] input_buffer;

  return 1;
}

