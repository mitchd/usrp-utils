/*Copyright 2013,2012 Joseph "Mitch" Davis mitchd@vt.edu
 *
 *This file is part of usrp-utils.
 *
 *   usrp-utils is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   usrp-energy is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with usrp-utils.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *This program is designed to record the energy of samples streamed from
 *a USRP device.
 *
 *The commandline options are:
 *
 * -o [file]       Output File  -The output file contains raw float data
 *                               representing the computed energy in each bin.
 *
 * -b [size]       Bin Size     -Energy bin size in samples
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
 *
 * -g [gain]       RX Gain      -Gain in DB of the rx chain
 *
 * Changelog
 *
 * 0.1 - Initial release 2012
 * 0.3 - Change to cmake and upload to github 20130701
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
#include <unistd.h>

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


/*setupUSRP(...)
 *
 *Setup the USRP for receiving at the specified freq and rate
 */
int setupUSRP(  uhd::usrp::multi_usrp::sptr&  usrp,
                const float                   center_freq,
                const float                   sample_rate,
                const int                     rx_gain,
                const char*                   dev_addr);




int calculateTask(  const char*                   outputFileName,
                    const int                     binSize,
                    const unsigned long long	    maximum_samples,
                    uhd::usrp::multi_usrp::sptr&  usrp );



/*******************************************************************************


*******************************************************************************/
int main( int argc, char* argv[])
{
  //First things first, try to set realtime priority for the parent thread
  uhd::set_thread_priority_safe();

  //Ensure the correct number of arguments were passed
  if( argc != 15)
  {
    cout << "Only " << argc << " parameters entered" << endl;
    useage();
    return -1;
  }

  char  *outputFileName = NULL;
  char  *usrpArgs       = NULL;
  int   usrpGain        = 0;
  int   arg             = 0;
  int   binSize         = 0;
  float usrpCenterFreq  = 0.0f;
  float usrpSampleRate  = 0.0f;
  float usrpRecordTime  = 0.0f;

  //argument parsing
  while( (arg = getopt( argc, argv, ":g:o:a:f:r:t:b:")) != -1 )
  {
    switch (arg)
    {
      case 'b':
        binSize = atoi(optarg);
        break;
      case 'o':
        outputFileName = new char[strlen(optarg)+1];
        strcpy(outputFileName,optarg);
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
        return -1;
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
    delete [] usrpArgs;
    return 1;
  }

  //Perform the actual work
  if( !calculateTask( outputFileName,
                      binSize,
                      static_cast<unsigned long long int>(usrpSampleRate*usrpRecordTime),
                      the_usrp ) )
  {
    cout << "Error performing calculations" << endl;
    delete [] outputFileName;
    delete [] usrpArgs;
    return 1;
  }

  uhd::stream_cmd_t       usrp_stream_stop(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
  usrp_stream_stop.stream_now = false;
  usrp_stream_stop.time_spec = uhd::time_spec_t();
  the_usrp->issue_stream_cmd( usrp_stream_stop );
  delete [] outputFileName;
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
        << "-b <size>\t Size (in samples) of energy bins" << endl
        << "-a <args>\t USRP Address" << endl
        << "-f <freq>\t USRP Center Frequency" << endl
        << "-r <rate>\t USRP Sample Rate" << endl
        << "-g <gain>\t USRP Rx Gain" << endl
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
int calculateTask(  const char*                   outputFileName,
                    const int                     binSize,
                    const unsigned long long	    maximum_samples,
                    uhd::usrp::multi_usrp::sptr&  usrp )
{
  ///////////////////////////////////////////////////////////
  //
  //Initialization Section
  ///////////////////////////////////////////////////////////

  //Number of bytes in a single-precision float
  const int   FLOAT_SIZE = sizeof(float);
  //Initialize and open the input/output files
  FILE *outputFile;

  if(!openFiles( outputFileName, outputFile ))
    return 0;

  //Setup the input buffer and tracking variables
  int                   return_code       = 1;

  //Setup the USRP for streaming
  vector<_Complex float >   usrpBuffer( binSize );
  float                   energy;
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
  cout << "Begin Data Collection" << endl;
  //Start streaming!
  usrp->issue_stream_cmd( usrp_stream_command );

  while( (samples_recorded < maximum_samples) and return_code )
  {
    energy = 0.0f;
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
    samples_recorded += buffer_samples_recorded;
    //Compute magnitude (we don't want to store phase information)
    for(int i = 0; i < binSize; i++ )
      energy += pow(cabsf( usrpBuffer[i] ), 2);

    //Write results to the output file
    fwrite(&energy, FLOAT_SIZE, 1, outputFile );
  }

  ///////////////////////////////////////////////////////////
  //
  //Cleanup Section
  ///////////////////////////////////////////////////////////

  //Toss out any leftovers and cleanup
  fclose(outputFile);

  return 1;
}

