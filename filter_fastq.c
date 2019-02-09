#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#define MAXLINE 300
#define WinLen 20
#define WinThres 19
#define AllLines 400

int main(int argc,char **argv) {

  // Open the file given from CLI for input
  FILE * Fin = fopen(argv[1], "r");
  // Open the file given from CLI for output
  FILE * Fout = fopen(argv[2], "w");

  // Loop through the file and
  // Do the same thing for all 4-line-groups
  int line_count;
  for(line_count = 0; line_count < AllLines ; line_count+=4) {

    int i;
    int Line;

    // Malloc for a 2-dimensional array of strings with
    // 4 lines and MAXLINE of characters per line
    char ** buffer;
    buffer=(char**)malloc(sizeof(char*)*4);
    for(i=0;i<4;i++) {
      buffer[i]=(char*)malloc(sizeof(char)*MAXLINE);
    }

    size_t len = 0;

    // read line-by-line the first 4 lines of the file
    // and store each in the array named buffer
    for(Line=0;Line<4;Line++) {
      getline(&buffer[Line], &len, Fin);
    }

    // The number of nucleotides in the second line
    // or equally in the last line
    int MaxLen=strlen(buffer[1])-1;

    // length of line[1] and line[3] MUST be equally
    if (strlen(buffer[3])!=strlen(buffer[1])) {
      printf("ERROR Lines 2 and 4 have different length\n");
      exit(-1);
    }

    float Qual=0;

    // start and end position of the sliding window
    int start=0;
    int end=start+WinLen;

    // Add 1 to WinThres for the first iteration
    Qual=WinThres+1;
    // slide the window while:
    // the end position has not reached the end of the line
    // and the mean quality score is above the minimum threshold
    while ((end<=MaxLen)&&Qual>WinThres) {

      // calculate the mean quality score
      Qual=0;
      for (int k=start;k<end;k++) {
        Qual+=buffer[3][k]-33;
      }
      //printf("Accumulated Qual=%f \t Mean Qual=%f \n",Qual,Qual/WinLen);
      Qual/=WinLen;

      //slide the window by one position to the right
      start++;
      end=start+WinLen;
    }
    start--;

    // if first window found low quality score
    if(start == 0){
      // free the allocated memory
      for (i=0;i<4;i++) {
        free(buffer[i]);
      }
      free(buffer);
      // discard the whole 4-line-group and move on to the next group
      continue;
    }

    // trim out the filter positions from
    // the second and the last lines up to
    strncpy(buffer[1],buffer[1],start);
    buffer[1][start]='\0';
    strncpy(buffer[3],buffer[3],start);
    buffer[3][start]='\0';

    //write the filtered fastq to the output file
    fprintf(Fout,"%s",buffer[0] );
    fprintf(Fout,"%s\n",buffer[1] );
    fprintf(Fout,"%s",buffer[2] );
    fprintf(Fout,"%s\n",buffer[3] );

    // free the allocated memory
    for (i=0;i<4;i++) {
      free(buffer[i]);
    }
    free(buffer);
  }
  // close the files opened
  fclose(Fin);
  fclose(Fout);

  exit(0);
}
