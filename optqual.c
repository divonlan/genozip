#include <stdio.h>
#include <stdlib.h>
#include <string.h>

main ()
{
    char *data = malloc (2700000);
    char *out = malloc (2700000);

    fread (data, 2700000, 1, stdin);

    out[0] = data[0];
    for (unsigned i=1; i< 2700000; i++) {
/*        if (data[i] == data[i-1]) out[i] = 193;
        else if (data[i] == data[i-1] + 1) out[i] = 194;
        else if (data[i] == data[i-1] - 1) out[i] = 195;
        //else if (data[i] == data[i-1] + 2) out[i] = 196;
        //else if (data[i] == data[i-1] - 2) out[i] = 197;*/
        if (data[i] > 72) out[i] = 72;
        else out[i] = data[i];
    }

    fwrite(out, 2700000,1,stdout);
}