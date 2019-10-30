// This module takes a line from a data file and converts it to an array of numbers (in practice, to a double pointer)
// There's probably a better way to do it, but this seems to work
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "linetodata.h"
#define LEN 200


double *linetodata(char *line, int *size)
{
int i, imin, lenline, words, started_word;
char word[LEN];
double *line_broken = NULL;
/* Iterate over number of characters in the whole line
 */
line_broken = calloc(4,sizeof(double));
//printf("START of linetodata\n");
lenline = strlen(line);
imin = 0;
words = 0;
started_word = 0;
for (i=0; i<lenline; i++)
{ /* Compare each character to a white space (32) or newline (10) character. When a white space is found (or a 
	 * newline character), we break the string and extract the number from the 
	 * previous whitepace to the new one.
	 */
	if ( (line[i] != 32) && (line[i] != 10) ) 
	{
		started_word = 1;
		//words += 1;
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
if (words == 4)
{
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
		started_word = 1;
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

return line_broken;
}
