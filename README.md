usrp-utils
==========

Some UHD level utility programs for use with USRP devices:


This software includes the following utilities:

	usrp-energy
	usrp-record
	usrp-sensor
        energycalculator
	fftcompute



Documentation for usrp-record:
 *This program is designed to record the samples streamed from
 *a USRP device.  This records the raw 16-bit I/Q samples from the USRP
 *
 *The commandline options are:
 *
 * -o [file]       Output File  -The output file contains raw float data
 *                               representing the computed spectral periodigram
 *                               of the recorded signal.  The peridigram is
 *                               computed with fftw3f's 1D DFT.
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



Documentation for usrp-energy:

 *This program is designed to record the magnitude of samples streamed from
 *a USRP device.
 *
 *The commandline options are:
 *
 * -o [file]       Output File  -The output file contains raw float data
 *                               representing the computed spectral periodigram
 *                               of the recorded signal.  The peridigram is
 *                               computed with fftw3f's 1D DFT.
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



Documentation for usrp-sensor:

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



Documentation for energycalculator:

 *This program takes I-Q data from a specified file, and computes a stream of energy bins.
 *
 *The I-Q Data is assumed to be stored in 'binary' complex<float> format.  (2 floats per 
 *sample). The energy bin size is defined as the number of samples summed to compute energy.
 *
 *This is the post-processing equivalent of usrp-energy.  This program can be run with data 
 *from usrp-record if usrp-record is configured to store fc32 data.
 *
 *The commandline options are:
 *
 * -o [file]       Output File  -The output file contains raw float data
 *                               representing the computed energy in each bin
 *
 * -b [size]       Bin Size     -Energy bin size in samples
 *
 * -i [file]	   Input File	-The input file contains raw float data


Documentation for fftcompute:
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




BUILDING:


To build this suite:

$mkdir build
$cd build
$cmake ..
$make

To build with debugging symbols and additional output:

$mkdir build
$cd build
$cmake -DDEBUG=1 ..
$make

You can change the wire format and host format for the various programs with:

To enable sc8 wire format (default sc16) (Affects usrp-energy usrp-recorder and
usrp-sensor)
-DWIRE_SC8=1

To enable sc16 host format (default fc32) (Only affects usrp-recorder)
-DHOST_SC16=1

