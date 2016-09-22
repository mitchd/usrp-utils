/*Copyright 2012-2016 Joseph "Mitch" Davis mitchd@vt.edu
 *
 *This file is part of usrp-utils.
 *
 *   usrp-utils is free software: you can redistribute it and/or modify
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
 *
 *This is the worker thread implementation.  Used in the FFTCompute program
 */

#include "fft_thread.h"


void* fft_thread_start( void* fft_thread_arg )
{
    const int FLOAT_SIZE = sizeof(float);


    //Sort out the input data
    struct fft_thread_data* my_thread_data;

    my_thread_data = reinterpret_cast<fft_thread_data*>(fft_thread_arg);

    float *magnitude  = new float[my_thread_data->fft_size];

    //increment is used for adavancing the thread_control to the next thread
    //in line
    int   increment   = 1;
    if( my_thread_data->my_id == my_thread_data->max_children - 1 )
    {
        increment = -1*my_thread_data->my_id;
    }

    //Open up this child's message queue
    mqd_t my_queue = mq_open( my_thread_data->mq_name, O_RDONLY );

    //Variables used for interfacing with the mqueue
    char*     my_msg        = new char[255];
    const int my_msg_length = __FFT_THREAD_MSG_LENGTH;
    unsigned* my_prio       = NULL;


    while( true )
    {
        //mq_receive is blocking, so we'll sit here and idle until we get a
        //message
        mq_receive( my_queue, my_msg, my_msg_length, my_prio );

        //Check to see if we got the message to kill this thread
        if( my_msg[0] == static_cast<char>(__FFT_THREAD_KILL) )
            break;
        else if( my_msg[0] == static_cast<char>(__FFT_THREAD_START) )
        {
            //Got a message to signal an FFT computation

            //Signal that the thread is running
            my_thread_data->is_running = true;

            //Apply the window function
            for(int i = 0; i < my_thread_data->fft_size; i++ )
            {
                my_thread_data->inputData[i] *= my_thread_data->window[i];
            }

            //Compute fft
            fftwf_execute( my_thread_data->plan );

            //Compute magnitude (we don't want to store phase information)
            for(int i = 0; i < my_thread_data->fft_size; i++ )
            {
                magnitude[i] = cabsf( my_thread_data->outputData[i] );
            }

            //Ensure that we're writing data to the output file in order
            while( *(my_thread_data->next_thread) != my_thread_data->my_id )
                ;

            //Write results to the output file

            //Lock mutex in preparation to write
            //this is a normal mutex, so this will block until we can grab it
            pthread_mutex_lock( my_thread_data->output_mutex );

            //Write to file

            //Negative freqs first
            fwrite(magnitude+(my_thread_data->fft_size/2), FLOAT_SIZE,
                   my_thread_data->fft_size/2, my_thread_data->outputFile);
            //Positive freqs next
            fwrite(magnitude, FLOAT_SIZE,
                   my_thread_data->fft_size/2, my_thread_data->outputFile);

            //Advance the next_thread
            *(my_thread_data->next_thread) += increment;

            //Unlock mutex
            pthread_mutex_unlock( my_thread_data->output_mutex );

            //Signal that we're done running and going back to the idle state
            my_thread_data->is_running = false;
        }
    }
    //Kill signal received

    //Close the mqueue
    mq_close( my_queue );
    //Free memory
    delete [] my_msg;
    delete [] magnitude;
    //Kill thread
    pthread_exit(NULL);
}
