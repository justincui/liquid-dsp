//
// firhilb_decim_example.c
//
// Hilbert transform: 2:1 real-to-complex decimator.  This example
// demonstrates the functionality of firhilb (finite impulse response
// Hilbert transform) decimator which converts a real time series
// into a complex one with half the number of samples.  The input
// is a real-valued sinusoid of N samples. The output is a complex-
// valued sinusoid of N/2 samples.
//
// SEE ALSO: firhilb_interp_example.c
//

#include <stdio.h>
#include <complex.h>
#include <math.h>

#include "liquid.h"

#define OUTPUT_FILENAME "firhilb_decim_example.m"

int main() {
    unsigned int m=5;           // filter semi-length
    float As=60.0f;             // stop-band attenuation [dB]
    float fc=0.37f;             // signal center frequency
    unsigned int N=128;         // number of samples

    // create Hilbert transform object
    firhilb f = firhilb_create(m,As);
    firhilb_print(f);

    // open output file
    FILE*fid = fopen(OUTPUT_FILENAME,"w");
    fprintf(fid,"%% %s : auto-generated file\n", OUTPUT_FILENAME);
    fprintf(fid,"clear all;\n");
    fprintf(fid,"close all;\n");
    fprintf(fid,"h_len=%u;\n", 4*m+1);
    fprintf(fid,"N=%u;\n", N);

    unsigned int i;
    float x[2], theta=0.0f, dtheta=2*M_PI*fc;
    float complex y;
    for (i=0; i<N; i++) {
        // compute real input signal
        x[0] = cosf(theta);
        theta += dtheta;
        x[1] = cosf(theta);
        theta += dtheta;

        // execute transform (decimator) to compute complex signal
        firhilb_decim_execute(f, x, &y);

        // print results
        fprintf(fid,"x(%3u) = %8.4f;\n", 2*i+1, x[0]);
        fprintf(fid,"x(%3u) = %8.4f;\n", 2*i+2, x[1]);
        fprintf(fid,"y(%3u) = %8.4f + j*%8.4f;\n", i+1, crealf(y), cimagf(y));

        printf("y(%3u) = %8.4f + j*%8.4f;\n", i+1, crealf(y), cimagf(y));
    }

    // plot results
    fprintf(fid,"nfft=512;\n");
    fprintf(fid,"%% compute normalized windowing functions\n");
    fprintf(fid,"wx = 1.8534/N*hamming(length(x)).'; %% x window\n");
    fprintf(fid,"wy = 1.8534/N*hamming(length(y)).'; %% y window\n");
    fprintf(fid,"X=20*log10(abs(fftshift(fft(x.*wx,nfft))));\n");
    fprintf(fid,"Y=20*log10(abs(        (fft(y.*wy,nfft))));\n");
    fprintf(fid,"f =[0:(nfft-1)]/nfft-0.5;\n");
    fprintf(fid,"fd=[0:(nfft-1)]/(2*nfft);\n");
    fprintf(fid,"figure; plot(f,X,'Color',[0.5 0.5 0.5],fd,Y,'LineWidth',2);\n");
    fprintf(fid,"grid on;\n");
    fprintf(fid,"axis([-0.5 0.5 -80 20]);\n");
    fprintf(fid,"xlabel('normalized frequency');\n");
    fprintf(fid,"ylabel('PSD [dB]');\n");
    fprintf(fid,"legend('original/real','transformed/decimated',1);");

    fclose(fid);
    printf("results written to %s\n", OUTPUT_FILENAME);

    firhilb_destroy(f);
    printf("done.\n");
    return 0;
}