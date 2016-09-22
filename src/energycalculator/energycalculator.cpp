/*Copyright 2012-2016 Joseph "Mitch" Davis mitchd@vt.edu
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
 *
 * Changelog
 *
 * 0.1 - Initial release 2012
 * 0.3 - Change to cmake and upload to github 20130701
 *
 *
 */

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>
//Uncomment this to get gratuitous debug information
//#define DEBUG 1
using namespace std;

/*useage()
 *
 *Display program useage information
 *
 */
void useage();

/*calculateTask
 *
 *This is the main computation of the program.  inputFileName must be a valid
 *data file, and we must have write permission for outputFileName.
 *
 *energyBinSize defines the number of samples summed for computing the energy
 */
int calculateTask( char* inputFileName, char* outputFileName, int energyBinSize );


/*wrideData
 *
 *Writes 4 bytes to the specified ofstream object
 *
 */
int writeData( FILE* outputFile, double data );



int main( int argc, char* argv[])
{
  //Ensure the correct number of arguments were passed
  if( argc != 7){
    cout << "Only " << argc << " parameters entered" << endl;
    useage();
    return -1;
  }

  char* inputFileName = NULL;
  char* outputFileName = NULL;
  int   energyBinSize = 0;
  int   arg = 0;

  //argument parsing
  while( (arg = getopt( argc, argv, "i:o:s:")) != -1 ){
#ifdef DEBUG
    cout << "Arg: " << optarg << endl;
#endif
    switch (arg){

    case 'i':
      inputFileName = new char[sizeof(optarg)];
      strcpy(inputFileName,optarg);
#ifdef DEBUG
      cout << "Input Filename: " << inputFileName << endl;
#endif
      break;

    case 'o':
      outputFileName = new char[sizeof(optarg)];
      strcpy(outputFileName,optarg);
#ifdef DEBUG
      cout << "Output Filename: " << outputFileName << endl;
#endif
      break;

    case 's':
      energyBinSize = atoi(optarg);
#ifdef DEBUG
      cout << "Energy Bin Size: " << energyBinSize << endl;
#endif
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

  if( !calculateTask( inputFileName, outputFileName, energyBinSize ) ){
    cout << "Error performing calculations" << endl;
    delete [] inputFileName;
    delete [] outputFileName;
    return 1;
  }

  delete [] inputFileName;
  delete [] outputFileName;
  return 0;
}


void useage( )
{
  cout  << "Useage:\t " << endl
        << "-i <file>\t Input File" << endl
        << "-o <file>\t Output File" << endl
        << "-s <size>\t Energy Bin Size" << endl;
}


int calculateTask( char* inputFileName, char* outputFileName, int energyBinSize )
{
  const int FLOAT_SIZE = 4; //Number of bytes in a single-precision float

  FILE* inputFile;
  FILE* outputFile;

  inputFile = fopen( inputFileName, "r");
#ifdef DEBUG
  cout << "Input file open attempt." << endl;
#endif
  outputFile = fopen( outputFileName, "w" );
#ifdef DEBUG
  cout << "Output file open attempt." << endl;
#endif

  /*Make sure the files actually opened.
   *
   *If not... gracefully exit.
   *We should really be throwing exceptions... but I'm feeling lazy
   */
  if( !inputFile || !outputFile){
#ifdef DEBUG
    cout << "Input File " << inputFileName << " Open? " << (inputFile ? "True" : "False") << endl;
    cout << "Output File " << outputFileName << " Open? " << (outputFile ? "True" : "False") << endl;
#endif
    if( inputFile )
      fclose( inputFile );
    if( outputFile )
      fclose( outputFile );
    return 0;
  }
#ifdef DEBUG
  cout << "Input/Output files opened successfully." << endl;
#endif

  float i = 0;
  float q = 0;

  double mag = 0;
  int    iterator = 0;

  while( !feof(inputFile) ){
      //Read in the I-Q of one sample... FLOAT_SIZE bytes per float
      fread( &i, FLOAT_SIZE, 1, inputFile );
      fread( &q, FLOAT_SIZE, 1, inputFile );

      //Ensure we're calling the double version of pow, otherwise data is lost
      mag += pow(static_cast<double>(i),2.) + pow(static_cast<double>(q),2.);

      //Reached the end of the bin?
      if( ++iterator == energyBinSize){
        iterator = 0;
        //If so, write 4 bytes
        if( !writeData( outputFile, mag ) ){
          fclose(inputFile);
          fclose(outputFile);
          return 0;
        }
        mag = 0;
      }
  }
  //Toss out any leftovers (incomplete energy bin)
  fclose(inputFile);
  fclose(outputFile);
#ifdef DEBUG
  cout << endl << "Files Closed." << endl;
#endif
  return 1;
}

int writeData( FILE* outputFile, double data )
{
  const int DOUBLE_SIZE = 8; //Number of bytes in a single-precision float
  fwrite(&data, DOUBLE_SIZE, 1, outputFile);
#ifdef DEBUG
  cout << ".";
#endif
  return 1;
}
