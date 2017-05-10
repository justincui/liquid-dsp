/*
 Copyright (c) 2003-2004, Mark Borgerding

 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "kiss_fft.h"
#include "kiss_fftr.h"

//assume buf, bufout are both of nfft kiss_fft_cpx elements
void fft_cplx(kiss_fft_cpx* buf, kiss_fft_cpx* bufout, int nfft, int isinverse) {
	kiss_fft_cfg cpx_st = kiss_fft_alloc(nfft, isinverse, 0, 0);
	kiss_fft(cpx_st, buf, bufout);
	free(cpx_st);
}

//assume rbuf is of nfft kiss_fft_scalar elements and cbuf has (nfft/2+1) number of kiss_fft_cpx elements.
void fft_real(kiss_fft_scalar * rbuf, kiss_fft_cpx * cbuf, int nfft, int isinverse) {
	kiss_fftr_cfg real_st = kiss_fftr_alloc(nfft, isinverse, 0, 0);
	if (isinverse == 0) {
		kiss_fftr(real_st, rbuf, cbuf);
	}
	else {
		kiss_fftri(real_st, cbuf, rbuf);
	}
	free(real_st);
}

kiss_fft_scalar real_imag_2_amplitude_phase(kiss_fft_cpx * cplx_array, int num) {
    kiss_fft_scalar max_amplitude=0;
	for(int i=0;i<num;i++){
		kiss_fft_cpx* pItem=&cplx_array[i];
		kiss_fft_scalar amplitude= sqrt(pItem->r * pItem->r + pItem->i * pItem->i);
		max_amplitude = max_amplitude > amplitude ? max_amplitude : amplitude;
		kiss_fft_scalar phase = atan2(pItem->i, pItem->r);
		pItem->r = amplitude;
		pItem->i = phase;
	}
	return max_amplitude;
}

kiss_fft_scalar* alloc_scalar_array(int nfft){
	return (kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar)*nfft);
}

kiss_fft_cpx* alloc_cplx_array(int nfft){
	return (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*nfft);
}

int sizeof_scalar(){
	return sizeof(kiss_fft_scalar);
}

int sizeof_cplx(){
	return sizeof(kiss_fft_cpx);
}

void test_real(){
	kiss_fft_scalar a[10];
	const int nfft=sizeof(a)/sizeof(a[0]);

	printf("real input: ");
	for(int i=0;i<nfft;i++){
		a[i]=i;
		printf("%f, ",a[i]);
	}
	printf("\n");

	kiss_fft_cpx b[nfft/2+1];
	fft_real(a,b, nfft, 0);
	real_imag_2_amplitude_phase(b,sizeof(b)/sizeof(b[0]));

	printf("cplx output:\n");
	for(int i=0;i<sizeof(b)/sizeof(b[0]);i++){
		printf("   [%d]:i=%f  q=%f\n",i, b[i].r, b[i].i);
	}
}

void test_cplx(){
	kiss_fft_cpx a[10];
	const int nfft=sizeof(a)/sizeof(a[0]);

	printf("cplx input: ");
	for(int i=0;i<nfft;i++){
		a[i].r=i;
		a[i].i=0;
		printf("(%f,%f) ",a[i].r, a[i].i);
	}
	printf("\n");

	kiss_fft_cpx b[nfft];
	fft_cplx(a,b, nfft, 0);
	real_imag_2_amplitude_phase(b,nfft);

	printf("cplx output:\n");
	for(int i=0;i<sizeof(b)/sizeof(b[0]);i++){
		printf("   [%d]:i=%f  q=%f\n",i, b[i].r, b[i].i);
	}
}

int main(){
	printf("KISS_FFT: bytes of scalar = %d\n", sizeof_scalar());
	test_real();
	test_cplx();
}
