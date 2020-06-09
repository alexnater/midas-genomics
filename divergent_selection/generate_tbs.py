#!/usr/bin/env python3

import argparse, sys
import numpy as np


def get_uniform_sample(umin, umax, nsamples, islog = False):
    if umin == 0. and umax == 0.:
        return [0.] * nsamples
    if islog:
        if umin == 0. or umax == 0.:
            raise Exception("Zero value for loguniform distribution specified!")
        return np.power(10, np.random.uniform(np.log10(umin), np.log10(umax), nsamples))
    else:
        return np.random.uniform(umin, umax, nsamples)


def print_samples(args):
    samples = []    
    if args.theta:
        samples.append(get_uniform_sample(args.theta[0], args.theta[1], args.nsamples))
    
    if args.rho:
        rhos = np.random.exponential(args.rho[0], int(1.2*args.nsamples))
        rhos = rhos[rhos <= args.rho[1]]
        while len(rhos) < args.nsamples:
            temp = np.random.exponential(args.rho[0], int(0.5*args.nsamples))
            temp = temp[temp <= args.rho[1]]
            rhos += temp
        rhos = rhos[:args.nsamples]
        samples.append(rhos)

    if args.alpha:
        if len(args.alpha) == args.npops * 2:
            for pop in range(args.npops):
                alphas = get_uniform_sample(args.alpha[2 * pop], args.alpha[2 * pop + 1], args.nsamples)
                alphas_het = alphas * args.h
                samples.extend([alphas, alphas_het])
        else:
            alphas = get_uniform_sample(args.alpha[0], args.alpha[1], args.nsamples)
            alphas_het = alphas * args.h
            samples.extend([alphas, alphas_het] * args.npops)            
    
    if args.tau:
        samples.append(get_uniform_sample(args.tau[0], args.tau[1], args.nsamples))
    
    if args.freq:
        if len(args.freq) == args.npops * 2:
            for pop in range(args.npops):
                freqs = get_uniform_sample(args.freq[2 * pop], args.freq[2 * pop + 1], args.nsamples, args.log)
                samples.append(freqs)
        else:
            freqs = get_uniform_sample(args.freq[0], args.freq[1], args.nsamples, args.log)
            samples.extend([freqs] * args.npops)  
    
    if args.smu:
        samples.append(get_uniform_sample(args.smu[0], args.smu[1], args.nsamples))
    
    if args.pos:
        samples.append(get_uniform_sample(args.pos[0], args.pos[1], args.nsamples))
    
    outfile = '/dev/stdout' if args.outfile == '-' or args.outfile == 'stdout' else args.outfile
    with open(outfile, 'w') as outhandle:
        for values in zip(*samples):
            print("\t".join(map(str,values)), file=outhandle)



def main(args):
    try:
        print_samples(args)
    except Exception as err:
        print("Error message:", err, file=sys.stderr)
        raise err


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("npops", help="Number of populations", type=int)
    parser.add_argument("nsamples", help="Number of samples", type=int)
    parser.add_argument("outfile", help="Output tbs file")
    parser.add_argument("-t", "--theta", dest="theta", nargs=2, help="Min and max of uniform distribution for theta", type=float, default=None)
    parser.add_argument("-r", "--rho", dest="rho", nargs=2, help="Mean and max of exponential distribution for rho", type=float, default=None)
    parser.add_argument("-a", "--alpha", dest="alpha", nargs='*', help="Min and max of uniform distribution for alpha for each population", type=float, default=None)
    parser.add_argument("-d", "--dom", dest="h", nargs=1, help="Dominance coefficient", type=float, default=0.5)
    parser.add_argument("-u", "--tau", dest="tau", nargs=2, help="Min and max of uniform distribution for tau", type=float, default=None)
    parser.add_argument("-s", "--smu", dest="smu", nargs=2, help="Min and max of uniform distribution for mutation rate to selected allele", type=float, default=None)
    parser.add_argument("-f", "--freq", dest="freq", nargs='*', help="Min and max of uniform distribution for initial frequency of selected allele for each population", type=float, default=None)
    parser.add_argument("-p", "--pos", dest="pos", nargs=2, help="Min and max of uniform distribution for position of selected site", type=float, default=None)
    parser.add_argument("--log", dest='log', action='store_true', help="Use loguniform distribution for frequency", default=False)
    main(parser.parse_args())

