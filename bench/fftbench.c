/*
 * Copyright (c) 2007, 2009, 2012 Joseph Gaeddert
 * Copyright (c) 2007, 2009, 2012 Virginia Polytechnic
 *                              Institute & State University
 *
 * This file is part of liquid.
 *
 * liquid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * liquid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with liquid.  If not, see <http://www.gnu.org/licenses/>.
 */

//
// fftbench.c : benchmark fft algorithms
//


// default include headers
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <sys/resource.h>

#include "liquid.h"

void usage()
{
    // help
    printf("Usage: benchmark [OPTION]\n");
    printf("Execute benchmark scripts for liquid-dsp library.\n");
    printf("  -h            display this help and exit\n");
    printf("  -v/q          verbose/quiet\n");
    printf("  -t[SECONDS]   set minimum execution time (s)\n");
    printf("  -o[FILENAME]  export output\n");
    printf("  -n[NFFT_MIN]  minimum FFT size (benchmark single FFT)\n");
    printf("  -N[NFFT_MAX]  maximum FFT size\n");
    printf("  -m[MODE]      mode: all, radix2, composite, prime, single\n");
}

// benchmark structure
struct benchmark_s {
    // fft options
    unsigned int nfft;          // FFT size
    int direction;              // FFT direction
    int flags;                  // FFT flags/method

    // benchmark results
    unsigned int num_trials;    // number of trials
    float extime;               // execution time

    // derived values
    float time_per_trial;       // execution time per trial
    float flops;                // computation bandwidth
};

// simulation structure
struct fftbench_s {
    enum {RUN_ALL=0,
          RUN_RADIX2,
          RUN_COMPOSITE,
          RUN_PRIME,
          RUN_SINGLE
    } mode;
    int verbose;
    float runtime;   // minimum run time (s)
    
    // min/max sizes for other modes
    unsigned int nfft_min;  // minimum FFT size (also, size for RUN_SINGLE mode)
    unsigned int nfft_max;  // maximum FFT size

    // output file
    char filename[128];     // output filename
    FILE * fid;             // output file pointer
    int output_to_file;     // output file write flag
};

// helper functions:
char convert_units(float * _s);
double calculate_execution_time(struct rusage, struct rusage);

// run all benchmarks
void fftbench_execute(struct fftbench_s * _fftbench);

// execute single benchmark
void execute_benchmark_fft(struct benchmark_s * _benchmark,
                           float                _runtime);

// main benchmark script
void benchmark_fft(struct rusage *      _start,
                   struct rusage *      _finish,
                   struct benchmark_s * _benchmark);

void benchmark_print_to_file(FILE * _fid,
                              struct benchmark_s * _benchmark);

void benchmark_print(struct benchmark_s * _benchmark);

// main function
int main(int argc, char *argv[])
{
    // options
    struct fftbench_s fftbench;
    fftbench.mode       = RUN_RADIX2;
    fftbench.verbose    = 1;
    fftbench.runtime    = 0.1f;
    fftbench.nfft_min   = 2;
    fftbench.nfft_max   = 1024;
    fftbench.filename[0]= '\0';
    fftbench.fid        = NULL;
    fftbench.output_to_file = 0;

    // get input options
    int d;
    while((d = getopt(argc,argv,"uvqn:N:t:o:m:")) != EOF){
        switch (d) {
        case 'h':   usage();        return 0;
        case 'v':   fftbench.verbose = 1;    break;
        case 'q':   fftbench.verbose = 0;    break;
        case 'n':
            fftbench.nfft_min = atoi(optarg);
            break;
        case 'N':
            fftbench.nfft_max = atoi(optarg);
            break;
        case 't':
            fftbench.runtime = atof(optarg);
            if (fftbench.runtime < 1e-3f)     fftbench.runtime = 1e-3f;
            else if (fftbench.runtime > 10.f) fftbench.runtime = 10.0f;
            printf("minimum runtime: %d ms\n", (int) roundf(fftbench.runtime*1e3));
            break;
        case 'o':
            fftbench.output_to_file = 1;
            strcpy(fftbench.filename, optarg);
            break;
        case 'm':
            if      (strcmp(optarg,"all")==0)       fftbench.mode = RUN_ALL;
            else if (strcmp(optarg,"radix2")==0)    fftbench.mode = RUN_RADIX2;
            else if (strcmp(optarg,"composite")==0) fftbench.mode = RUN_COMPOSITE;
            else if (strcmp(optarg,"prime")==0)     fftbench.mode = RUN_PRIME;
            else if (strcmp(optarg,"single")==0)    fftbench.mode = RUN_SINGLE;
            else {
                fprintf(stderr,"error: %s, unknown mode '%s'\n", argv[0], optarg);
                exit(1);
            }
            break;
        default:
            usage();
            return 0;
        }
    }

    // run empty loop; a bug was found that sometimes the first package run
    // resulted in a longer execution time than what the benchmark really
    // reflected.  This loop prevents that from happening.
    unsigned int i;
    for (i=0; i<1e6; i++) {
        // do nothing
    }
    
    // open output file (if applicable)
    if (fftbench.output_to_file) {
        fftbench.fid = fopen(fftbench.filename,"w");
        if (!fftbench.fid) {
            fprintf(stderr,"error: %s, could not open file '%s' for writing\n", argv[0], fftbench.filename);
            exit(1);
        }
        FILE * fid = fftbench.fid;

        // print header
        fprintf(fid,"# %s : auto-generated file\n", fftbench.filename);
        fprintf(fid,"#\n");
        fprintf(fid,"# invoked as:\n");
        fprintf(fid,"#   ");
        for (i=0; i<argc; i++)
            fprintf(fid," %s", argv[i]);
        fprintf(fid,"\n");
        fprintf(fid,"#\n");
        fprintf(fid,"# properties:\n");
        fprintf(fid,"#  verbose             :   %s\n", fftbench.verbose ? "true" : "false");
        fprintf(fid,"#  runtime             :   %12.8f s\n", fftbench.runtime);
        fprintf(fid,"#  mode                :   \n");
        fprintf(fid,"#\n");
        fprintf(fid,"# %12s %12s %12s %12s %12s\n",
                "nfft", "num trials", "ex. time", "us/trial", "M-flops");
    }

    // run benchmarks
    fftbench_execute(&fftbench);

    if (fftbench.output_to_file) {
        fclose(fftbench.fid);
        printf("results written to %s\n", fftbench.filename);
    }

    return 0;
}

// convert raw value into metric units,
//   example: "0.01397s" -> "13.97 ms"
char convert_units(float * _v)
{
    char unit;
    if (*_v < 1e-9)     {   (*_v) *= 1e12;  unit = 'p';}
    else if (*_v < 1e-6){   (*_v) *= 1e9;   unit = 'n';}
    else if (*_v < 1e-3){   (*_v) *= 1e6;   unit = 'u';}
    else if (*_v < 1e+0){   (*_v) *= 1e3;   unit = 'm';}
    else if (*_v < 1e3) {   (*_v) *= 1e+0;  unit = ' ';}
    else if (*_v < 1e6) {   (*_v) *= 1e-3;  unit = 'k';}
    else if (*_v < 1e9) {   (*_v) *= 1e-6;  unit = 'M';}
    else if (*_v < 1e12){   (*_v) *= 1e-9;  unit = 'G';}
    else                {   (*_v) *= 1e-12; unit = 'T';}

    return unit;
}

double calculate_execution_time(struct rusage _start, struct rusage _finish)
{
    return _finish.ru_utime.tv_sec - _start.ru_utime.tv_sec
        + 1e-6*(_finish.ru_utime.tv_usec - _start.ru_utime.tv_usec)
        + _finish.ru_stime.tv_sec - _start.ru_stime.tv_sec
        + 1e-6*(_finish.ru_stime.tv_usec - _start.ru_stime.tv_usec);
}

// run all benchmarks
void fftbench_execute(struct fftbench_s * _fftbench)
{
    // validate input
    if (_fftbench->nfft_min > _fftbench->nfft_max && _fftbench->mode != RUN_SINGLE) {
        fprintf(stderr,"execute_benchmarks_composite(), nfft_min cannot be greater than nfft_max\n");
        exit(1);
    } else if (_fftbench->runtime <= 0.0f) {
        fprintf(stderr,"execute_benchmarks_composite(), runtime must be greater than zero\n");
        exit(1);
    }

    // create benchmark structure
    struct benchmark_s benchmark;

    if (_fftbench->mode == RUN_SINGLE) {
        // run single benchmark and exit

        // initialize benchmark structure
        benchmark.nfft       = _fftbench->nfft_min;
        benchmark.direction  = FFT_FORWARD;
        benchmark.num_trials = 1;
        benchmark.flags      = 0;
        benchmark.extime     = 0.0f;
        benchmark.flops      = 0.0f;

        // run the benchmark
        execute_benchmark_fft(&benchmark, _fftbench->runtime);
        benchmark_print(&benchmark);
        return;
    }
        
    printf("running ");
    switch (_fftbench->mode) {
    case RUN_ALL:       printf("all");              break;
    case RUN_RADIX2:    printf("all power-of-two"); break;
    case RUN_COMPOSITE: printf("all composite");    break;
    case RUN_PRIME:     printf("all prime");        break;
    case RUN_SINGLE:
    default:;
    }
    printf(" FFTs from %u to %u:\n",
        _fftbench->nfft_min,
        _fftbench->nfft_max);

    unsigned int nfft;
    for (nfft=_fftbench->nfft_min; nfft<=_fftbench->nfft_max; nfft++) {
        int isprime  = liquid_is_prime(nfft);
        int isradix2 = (1 << liquid_nextpow2(nfft))==nfft ? 1 : 0;

        // check run mode
        if (_fftbench->mode == RUN_RADIX2 && !isradix2)
            continue;
        if (_fftbench->mode == RUN_COMPOSITE && (isprime || isradix2))
            continue;
        if (_fftbench->mode == RUN_PRIME && !isprime)
            continue;

        // run the transform
        benchmark.nfft       = nfft;
        benchmark.direction  = FFT_FORWARD;
        benchmark.num_trials = 1;
        benchmark.flags      = 0;
        benchmark.extime     = 0.0f;
        benchmark.flops      = 0.0f;
        execute_benchmark_fft(&benchmark, _fftbench->runtime);

        if (_fftbench->verbose)
            benchmark_print(&benchmark);

        if (_fftbench->output_to_file)
            benchmark_print_to_file(_fftbench->fid, &benchmark);
    }
}

// execute single benchmark
void execute_benchmark_fft(struct benchmark_s * _benchmark,
                           float                _runtime)
{
    unsigned long int n = _benchmark->num_trials;
    struct rusage start, finish;

    unsigned int num_attempts = 0;
    do {
        // increment number of attempts
        num_attempts++;

        // set number of trials and run benchmark
        _benchmark->num_trials = n;
        benchmark_fft(&start, &finish, _benchmark);

        // calculate execution time
        _benchmark->extime = calculate_execution_time(start, finish);

        // check exit criteria
        if (_benchmark->extime >= _runtime) {
            break;
        } else if (num_attempts == 30) {
            fprintf(stderr,"warning: benchmark could not execute over minimum run time\n");
            break;
        } else {
            // increase number of trials
            n *= 2;
        }
    } while (1);

    // compute time per trial
    _benchmark->time_per_trial = _benchmark->extime / (float)_benchmark->num_trials;

    // computational bandwidth (see: http://www.fftw.org/speed/)
    _benchmark->flops = 5.0f * _benchmark->nfft * log2f(_benchmark->nfft) / _benchmark->time_per_trial;
}

// main benchmark script
void benchmark_fft(struct rusage *      _start,
                   struct rusage *      _finish,
                   struct benchmark_s * _benchmark)
{
    // initialize arrays, plan
    float complex * x = (float complex *) malloc((_benchmark->nfft)*sizeof(float complex));
    float complex * y = (float complex *) malloc((_benchmark->nfft)*sizeof(float complex));
    fftplan q = fft_create_plan(_benchmark->nfft,
                                x, y,
                                _benchmark->direction,
                                _benchmark->flags);
    
    unsigned long int i;

    // initialize input with random values
    for (i=0; i<_benchmark->nfft; i++)
        x[i] = randnf() + randnf()*_Complex_I;

    // scale number of iterations to keep execution time
    // relatively linear
    unsigned int num_iterations = _benchmark->num_trials;

    // start trials
    getrusage(RUSAGE_SELF, _start);
    for (i=0; i<num_iterations; i++) {
        fft_execute(q);
        fft_execute(q);
        fft_execute(q);
        fft_execute(q);
    }
    getrusage(RUSAGE_SELF, _finish);

    // set actual number of iterations in result
    _benchmark->num_trials = num_iterations * 4;

    fft_destroy_plan(q);
    free(x);
    free(y);
}

void benchmark_print_to_file(FILE * _fid,
                             struct benchmark_s * _benchmark)
{
    fprintf(_fid,"  %12u %12u %12.4e %12.6f %12.3f\n",
            _benchmark->nfft,
            _benchmark->num_trials,
            _benchmark->extime,
            _benchmark->time_per_trial * 1e6f,
            _benchmark->flops * 1e-6f);
}

void benchmark_print(struct benchmark_s * _benchmark)
{
    // format trials (iterations)
    float trials_format = (float)(_benchmark->num_trials);
    char trials_units = convert_units(&trials_format);

    // format time (seconds)
    float extime_format = _benchmark->extime;
    char extime_units = convert_units(&extime_format);

    // format rate (time/trial)
    float time_format = _benchmark->time_per_trial;
    char time_units = convert_units(&time_format);

    // format computational bandwidth (flops)
    float flops_format = _benchmark->flops;
    char flops_units = convert_units(&flops_format);

    printf("  %-8u: %8.2f %c trials / %6.2f %cs (%6.2f %cs/t) > %6.2f %c flops\n",
            _benchmark->nfft,
            trials_format, trials_units,
            extime_format, extime_units,
            time_format,   time_units,
            flops_format,  flops_units);
}
