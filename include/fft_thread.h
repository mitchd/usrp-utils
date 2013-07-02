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
 *
 *This is the worker thread definition.  Used in the usrp-sensor program
 */
#ifndef FFT_THREAD_H_INCLUDED
#define FFT_THREAD_H_INCLUDED

#include <pthread.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#include <complex.h>
#include <iostream>
#include <mqueue.h>


#define __FFT_THREAD_START       1
#define __FFT_THREAD_KILL        2
#define __FFT_THREAD_MSG_LENGTH  1
#define __FFT_THREAD_MSG_PRIO    0

struct fft_thread_data
{
    FILE*             outputFile;     //File to output the FFT results

    fftwf_plan        plan;           //fftw3 fft plan

    float complex*    outputData;     //fft output data (only valid after
    //fftw_execute)

    float complex*    inputData;      //fft input data (depending on the fft plan,
    //these data may get destroyed)

    float*            window;         //window function

    pthread_mutex_t*  output_mutex;   //Output file mutex

    int               fft_size;       //FFT Size
    volatile bool     is_running;     //Signal that the thread is busy
    char*             mq_name;        //mqueue name

    int               my_id;          //child thread id number
    volatile int*     next_thread;    //next thread id to execute
    int               max_children;   //maximum number of child processes
};

/*fft_thread_start
 *
 *pthread starting function for the fftw calculation thread
 *see the above struct for key parameters
 */
void* fft_thread_start( void* fft_thread_arg );


#endif // FFT_THREAD_H_INCLUDED
