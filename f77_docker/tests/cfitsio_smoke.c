#include <fitsio.h>

/*
 * Function: main
 * Purpose: Verify that the CFITSIO header and shared library agree on version.
 * Parameters: None.
 * Returns: Zero for CFITSIO 4.3.x; nonzero otherwise.
 */
int main(void)
{
    float version = 0.0f;
    float expected;
    float difference;

    ffvers(&version);
    expected = (float) CFITSIO_MAJOR
        + 0.01f * (float) CFITSIO_MINOR
        + 0.0001f * (float) CFITSIO_MICRO;
    difference = version - expected;
    if (difference < 0.0f) {
        difference = -difference;
    }

    return difference < 0.00001f ? 0 : 1;
}
