/*Copyright 2012 Joseph "Mitch" Davis mitchd@vt.edu
 *Copyright 2012 Virginia Tech
 *
 *This file is part of fftcompute.
 *   fftcompute is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   fftcompute is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with fftcompute.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *This program is designed to leverage single-precision fftw3 libraries and
 *pthreads to compute spectral periodigrams using an applied window function
 *and the magnitude of a 1D DFT.
 *
 *The commandline options are:
 *
 * -i [file]       Input File   -This input file is expected to contain
 *                               raw complex float data representing recorded
 *                               samples (from MATLAB, GNU Radio, etc.)
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
 *                               threads to spawn for FFTCompute.  If your
 *                               processor has N cores (including hyperthreading
 *                               cores), it is recommended to specify N-1
 *                               children.  This allows for the parent thread
 *                               to process the input data for a total of N
 *                               threads.
 *
 *
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
 *On an Intel(R) Xeon(R) W3530 processor with 6GB of RAM, this program will
 *compute about 120,000 FFTs of size 1024 per second with 7 worker threads.
 *
 *This program is the post-processing equivalent of usrp-sensor.  usrp-record
 *data may be used if configured to record fc32 data
 *
 */

#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include <mqueue.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>

#include "fft_thread.h"

#ifdef BENCHMARK
#include <ctime>
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
int openFiles( char* inputFileName, FILE*& inputFile,
               char* outputFileName, FILE*& outputFile );

/*createFFTPlans( fftwf_plan*&, float complex**&, float complex**&,
                  int, int)
 *
 *Initialize the FFT plans.  This creates max_children plans for FFTSize FFTs
 */
void createFFTPlans( fftwf_plan*& plans, float complex**& inputData,
                     float complex**& outputData, int FFTSize,
                     int max_children );

/*initializeThreads(...)
 *
 *Initialize and spawn all the child threads.  Returns 0 upon thread or mqueue
 *failure.
 */
int initializeThreads( pthread_mutex_t& mutex, pthread_mutex_t& attr,
                       pthread_t*& fft_children, mqd_t*& fft_mq,
                       struct mq_attr*& ma, int max_children,
                       struct fft_thread_data*& fft_child_args,
                       FILE*& outputFile, fftwf_plan*& plans,
                       float complex**& inputData, float complex**& outputData,
                       float*& window, int FFTSize, int*& thread_control );

/*calculateTask(...)
 *
 *This is the main work of the program, performing the specified overlapped
 *FFT transforms
 */
int calculateTask(  char* inputFileName, char* outputFileName,
                    int FFTSize, int FFTOverlap, int max_children,
                    float* window );










/*******************************************************************************


*******************************************************************************/
int main( int argc, char* argv[])
{
    const int FLOAT_SIZE = sizeof(float); //Size of single-precision float
    //in bytes

    //Ensure the correct number of arguments were passed
    if( argc != 11)
    {
        cout << "Only " << argc << " parameters entered" << endl;
        useage();
        return -1;
    }

    char  *inputFileName  = NULL;
    char  *outputFileName = NULL;
    char  *windowFileName = NULL;
    int   FFTSize         = 0;
    int   FFTOverlap      = 0;
    int   arg             = 0;
    int   max_children    = 0;

    //argument parsing
    while( (arg = getopt( argc, argv, "i:o:s:l:c:w:")) != -1 )
    {
        switch (arg)
        {

        case 'i':
            inputFileName = new char[strlen(optarg)];
            strcpy(inputFileName,optarg);
            break;

        case 'o':
            outputFileName = new char[strlen(optarg)];
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
            windowFileName = new char[strlen(optarg)];
            strcpy(windowFileName,optarg);
            break;

        case '?':
            useage();
            if( inputFileName )
                delete [] inputFileName;
            if( outputFileName )
                delete [] outputFileName;
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
             << "Assuming Uniform Window" << endl;
        for(int i = 0; i < FFTSize; i++)
          window[i] = 1.0f;
    }
    else
    {
      //This will get the number of window elements in the specified file.
      //We want this info because its possible that our window function is
      //shorter than the FFTSize.  In which case, we need to pad the remainder
      //of the window array with zeros
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
          delete [] inputFileName;
          delete [] outputFileName;
          delete [] windowFileName;
          delete [] window;
          fclose( window_file );
          return 1;
      }
      //Finally read the window
      fread( window, FLOAT_SIZE , window_size, window_file );
      fclose(window_file);

      //Accomadate zero-padding
      if (window_size < FFTSize )
      {
          cout  << "Window is smaller than FFT Size, assuming zero-padding. " 
                << endl << "FFT Size: " << FFTSize << endl
                        << "Window Size: " << window_size << endl;
          //Zero-pad the rest of the window
          for(int i = window_size; i < FFTSize; i++ )
              window[i] = 0.0f;
      }  
    }



    if( !calculateTask( inputFileName, outputFileName, FFTSize, FFTOverlap,
                        max_children, window ) )
    {
        cout << "Error performing calculations" << endl;
        delete [] inputFileName;
        delete [] outputFileName;
        delete [] windowFileName;
        delete [] window;
        return 1;
    }

    delete [] inputFileName;
    delete [] outputFileName;
    delete [] windowFileName;
    delete [] window;
    return 0;
}









/*******************************************************************************


*******************************************************************************/
void useage( )
{
    cout  << "Useage:\t FFTCompute [args]" << endl
          << "-i <file>\t Input File" << endl
          << "-o <file>\t Output File" << endl
          << "-s <size>\t FFT Size" << endl
          << "-l <number>\t FFT Overlap" << endl
          << "-c <number>\t Number of Child Processes" << endl
          << "-w <file>\t Window File" << endl;
}







/*******************************************************************************


*******************************************************************************/
int openFiles( char* inputFileName, FILE*& inputFile,
               char* outputFileName, FILE*& outputFile )
{
    inputFile = fopen( inputFileName, "r");
    outputFile = fopen( outputFileName, "w" );

    /*Make sure the files actually opened.
     *
     *If not... gracefully exit.
     *We should really be throwing exceptions...
     */
    if( !inputFile || !outputFile)
    {
        if( inputFile )
            fclose( inputFile );
        if( outputFile )
            fclose( outputFile );
        return 0;
    }
    return 1;
}









/*******************************************************************************


*******************************************************************************/
void createFFTPlans( fftwf_plan*& plans, float complex**& inputData,
                     float complex**& outputData, int FFTSize,
                     int max_children )
{
    plans       = new fftwf_plan[ max_children ];
    inputData   = new float complex*[ max_children ];
    outputData  = new float complex*[ max_children ];

    //Setup the FFT plans -- create one for each child thread
    for(int i = 0; i < max_children; i++ )
    {
        inputData[i]  = new float complex[FFTSize];
        outputData[i] = new float complex[FFTSize];
        plans[i]      = fftwf_plan_dft_1d( FFTSize, inputData[i], outputData[i],
                                           FFTW_FORWARD, FFTW_EXHAUSTIVE );
    }

    return;
}










/*******************************************************************************


*******************************************************************************/
int initializeThreads( pthread_mutex_t& mutex, pthread_mutexattr_t& attr,
                       pthread_t*& fft_children, mqd_t*& fft_mq,
                       struct mq_attr*& ma, int max_children,
                       struct fft_thread_data*& fft_child_args,
                       FILE*& outputFile, fftwf_plan*& plans,
                       float complex**& inputData, float complex**& outputData,
                       float*& window, int FFTSize, int*& thread_control )
{
    //mutex attribute -> PTHREAD_MUTEX_NORMAL
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
        fft_child_args[i].outputFile    = outputFile;
        fft_child_args[i].output_mutex  = &mutex;
        fft_child_args[i].plan          = plans[i];
        fft_child_args[i].fft_size      = FFTSize;
        fft_child_args[i].inputData     = inputData[i];
        fft_child_args[i].outputData    = outputData[i];
        fft_child_args[i].is_running    = false;
        fft_child_args[i].mq_name       = new char[255];
        //Ensure that we got the memory
        if( !fft_child_args[i].mq_name )
            return 0;
        fft_child_args[i].my_id         = i;
        fft_child_args[i].window        = window;
        fft_child_args[i].next_thread   = thread_control;
        fft_child_args[i].max_children  = max_children;
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
int calculateTask(  char* inputFileName, char* outputFileName, int FFTSize,
                    int FFTOverlap, int max_children,
                    float* window )
{
    ///////////////////////////////////////////////////////////
    //
    //Initialization Section
    ///////////////////////////////////////////////////////////

    //Number of bytes in a single-precision float
    const int   FLOAT_COMPLEX_SIZE  = sizeof(float complex);
    const char  msg_thread_start    = static_cast<char>(__FFT_THREAD_START);
    const char  msg_thread_kill     = static_cast<char>(__FFT_THREAD_KILL);

    //Initialize and open the input/output files
    FILE *inputFile;
    FILE *outputFile;

    if(!openFiles( inputFileName, inputFile, outputFileName, outputFile ))
        return 0;

    //Create some FFT Plans
    fftwf_plan     *plans       = NULL;
    float complex **inputData   = NULL;
    float complex **outputData  = NULL;

    createFFTPlans( plans, inputData, outputData, FFTSize, max_children );

    //Setup the mutex for the output file... we don't want a race condition
    pthread_mutex_t     output_mutex;
    pthread_mutexattr_t output_mutex_attr;

    pthread_t               *fft_children   = NULL;
    mqd_t                   *fft_mq         = NULL;
    struct mq_attr          *ma             = NULL;
    struct fft_thread_data  *fft_child_args = NULL;
    int                     *thread_control = NULL;

    if( !initializeThreads( output_mutex, output_mutex_attr, fft_children, fft_mq,
                            ma, max_children, fft_child_args, outputFile, plans,
                            inputData, outputData, window, FFTSize, thread_control ))
    {
        //Cleanup for a graceful exit
        //We have to check to see if things exist before deleting them because
        //initialization failed.  So its possible that some things exist, and others
        //do not.
        if( plans )
        {
            for(int i = 0; i < max_children; i++)
            {
                fftwf_destroy_plan(plans[i]);
            }
            delete [] plans;
        }
        if( inputData )
        {
            for(int i = 0; i < max_children; i++)
                delete [] inputData[i];
            delete [] inputData;
        }
        if( outputData )
        {
            for(int i = 0; i < max_children; i++)
                delete [] outputData[i];
            delete [] outputData;
        }
        if( fft_children )
        {
            for(int i = 0; i < max_children; i++)
                mq_send( fft_mq[i], &msg_thread_kill,
                         __FFT_THREAD_MSG_LENGTH, __FFT_THREAD_MSG_PRIO );
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
            {
                if( fft_child_args[i].mq_name )
                    delete [] fft_child_args[i].mq_name;
            }
            delete [] fft_child_args;
        }
        return 0;
    }



    //Setup the input buffer and tracking variables
    int             fft_interval_size = FFTSize / FFTOverlap;
    float complex*  input_buffer      = new float complex[FFTSize];
    int             head              = 0;
    int             child_tracker     = 0;
    bool            isFirst           = true;
    bool            isFull            = false;
    int             return_code       = 1;
    size_t          bytes_read        = 0;
    //Check memory allocation
    if( !input_buffer )
        return_code = 0;
    ///////////////////////////////////////////////////////////
    //
    //Work Section
    ///////////////////////////////////////////////////////////
#ifdef BENCHMARK
    timeval a;
    timeval b;

    gettimeofday(&a, 0);
#endif
    //Start at the first child thread
    (*thread_control) = 0;

    while( !feof(inputFile) && return_code )
    {
        //Read in the I-Q of fft_interval_size samples...
        //FLOAT_COMPLEX_SIZE bytes per float
        bytes_read = fread( &input_buffer[head], FLOAT_COMPLEX_SIZE,
                            fft_interval_size, inputFile );

        if( bytes_read != static_cast<unsigned int>(fft_interval_size))
        {
            //encountered EOF before reaching the end of the FFT interval
            return_code = 0;
            break;
        }
        //Time to take an FFT yet?
        if( !isFull )
        {
            if( head == (FFTSize - fft_interval_size) )
                isFull = true;
            else
                head += fft_interval_size;
        }
        else
        {
            //If this isn't the very first FFT
            if( !isFirst )
            {
                //rotate the buffer
                head += fft_interval_size;
                if( head == FFTSize )
                    head = 0;
            }
            else
            {
                isFirst = false;
            }
            //Check to see if the thread is still active before copying data
            if( fft_child_args[child_tracker].is_running )        //This is very bayd
            {
                while( fft_child_args[child_tracker].is_running )
                    ;  //We have to block until that thread finishes
                //We can't pthread_join because the thread doesn't
                //terminate until we generate a _THREAD_KILL command
                //So instead we wait for the running flag to go false
            }
            //Copy the buffer into the FFT input data
            memmove( fft_child_args[child_tracker].inputData,
                     input_buffer+head, FLOAT_COMPLEX_SIZE * (FFTSize-head) );

            memmove( fft_child_args[child_tracker].inputData+FFTSize-head,
                     input_buffer, FLOAT_COMPLEX_SIZE * head );
            //Send a message to the worker to compute the FFT
            mq_send( fft_mq[child_tracker], &msg_thread_start,
                     __FFT_THREAD_MSG_LENGTH, __FFT_THREAD_MSG_PRIO );

            child_tracker = child_tracker == max_children - 1 ? 0 : child_tracker + 1;
        }
    }
#ifdef BENCHMARK
    gettimeofday(&b, 0);

    float seconds = b.tv_sec - a.tv_sec;
    seconds += b.tv_usec/1000000.0f - a.tv_usec/1000000.0f;
    printf("ET: %.6f s\n",seconds);
#endif
    ///////////////////////////////////////////////////////////
    //
    //Cleanup Section
    ///////////////////////////////////////////////////////////

    //If we exited the work section due to premature EOF, wait for worker
    //threads to complete before destroying things
    if( !return_code )
    {
        cout << "Input data terminated with unaligned data" << endl;
        for(int i = 0; i < max_children; i++ )
        {
            while( fft_child_args[i].is_running )
                ;
        }
    }

    //Toss out any leftovers and cleanup
    fclose(inputFile);
    fclose(outputFile);

    //Destroy plans
    for(int i = 0; i < max_children; i++)
    {
        mq_send( fft_mq[i], &msg_thread_kill,
                 __FFT_THREAD_MSG_LENGTH, __FFT_THREAD_MSG_PRIO );
        fftwf_destroy_plan(plans[i]);
        delete [] inputData[i];
        delete [] outputData[i];
        mq_close( fft_mq[i] );
        mq_unlink( fft_child_args[i].mq_name );
    }

    //cleanup fftw residuals
    fftwf_cleanup();

    //cleanup mutex
    pthread_mutexattr_destroy( &output_mutex_attr );
    pthread_mutex_destroy( &output_mutex );

    //free more memory
    for( int i = 0; i < max_children; i++ )
    {
        delete [] fft_child_args[i].mq_name;
    }

    delete [] plans;
    delete [] inputData;
    delete [] outputData;
    delete [] ma;
    delete [] fft_mq;
    delete [] fft_child_args;
    delete [] input_buffer;

    return 1;
}

