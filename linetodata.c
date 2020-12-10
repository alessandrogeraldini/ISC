// This module takes a line from a data file and stores it in an array of doubles
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "isc.h"
#define LEN 200


double *linetodata(char *line, int *size) {
	int i, imin, lenline, words, started_word, save=0;
	char word[LEN];
	double *line_broken;
	/* Iterate over number of characters in the whole line
	 */
	if (*size<1)
	line_broken = malloc(LEN*sizeof(double));
	else 
	line_broken = malloc((*size)*sizeof(double));
	//printf("START of linetodata\n");
	lenline = strlen(line);
	imin = 0;
	words = 0;
	started_word = 0;
	for (i=0; i<lenline; i++)
	{ 
/* Compare each character to a white space (32) or newline (10) character. When a white space is found (or a 
 newline character), we break the string and extract the number from the 
 previous whitepace to the new one.
 */
		if ( (line[i] != 32) && (line[i] != 10) ) 
		{
			started_word = 1;
			//words += 1;
			if ( (isdigit(line[i]) != 0) && (words == 0) ) save = 1;
			if ( strncmp(line, "param", 5) == 0 ) save = 1;
		}
		else 
		{
			if (started_word == 1)
			{
				started_word = 0;
				words += 1;
			}
			imin = i+1;
		}
	}
	if (save == 1) { // CHANGE HERE
		//printf("saving line\n");
		imin = 0;
		words = 0;
		started_word = 0;
		for (i=0; i<lenline; i++)
		{
			/* Compare each character to a white space. When a white space is found (or a 
			 * newline character), we break the string and extract the number from the 
			 * previous whitepace to the new one.
			 */
			if ( (line[i] != 32) && (line[i] != 10) ) 
			{
				word[i-imin] = line[i];
				if (isdigit(word[i-imin]) != 0) {
					started_word = 1;
				}
			}
			else
			{
				if (started_word == 1)
				{	
					word[i-imin] = '\0';
					started_word = 0;
					sscanf(word, "%lf", (line_broken+words));
					words += 1;
					//line_broken += 1;
				}
				imin = i + 1;
			}
		}
	}
	*size = words;

	//printf("*line_broken   = %Le\n", *(line_broken)) ;
	//printf("*line_broken+1 = %Le\n", *(line_broken+1));
	//printf("*line_broken+2 = %Le\n", *(line_broken+2));
	//printf("*line_broken+3 = %Le\n", *(line_broken+3));
	//printf("%s\n", line);

	return line_broken;
}
