#include <stdio.h>

int main (int argc, char *argv [])
{ FILE *fbin = fopen (argv[1], "rb");
  int i = 0;
  double tobs = 0, freq = 0;
  float acm[288*289] = {0,};
  // for (i=0; i<3; i++)
  while (!feof (fbin))
  { fread (&tobs, sizeof (double), 1, fbin);
    fread (&freq, sizeof (double), 1, fbin);
    fread (acm, sizeof (float), 288*289, fbin);
    fprintf (stderr, "timeslice %d: %12.3f, Freq: %12.3f\n", i, tobs, freq);
    i++;
  }

  return 0;
} 
