#!/usr/bin/env python
from __future__ import print_function, division

import argparse
import os
import sys
from shutil import which

import numpy as np
import sh


def unphase_vcf(in_file, out_file=None):
    sed = sh.Command(which('sed'))
    if out_file is None:
        out = sys.stdout
    else:
        out = open(out_file, "wt")
    it = sed(*[r'/^##/! s/|/\//g', in_file], _iter=True)
    for line in it:
        print(line[:-1], file=out)
    out.close()


def demography_from_params(params):
    demography = []
    ct = 0.0
    z = list(zip(*params))
    for ai, bi, si in z[:-1]:
        beta = (np.log(ai) - np.log(bi)) / si
        demography += ['-eN', ct, ai]
        demography += ['-eG', ct, beta]
        ct += si
    demography += ['-eN', ct, z[-1][0]]
    demography += ['-eG', ct, 0.0]
    return demography


SCRM = os.environ.get('SCRM_PATH', False) or which('scrm')

if __name__ == "__main__":
    if not SCRM:
        sys.exit("Can't find scrm. Please set SCRM_PATH.")
    scrm = sh.Command(SCRM)
    parser = argparse.ArgumentParser()
    parser.add_argument("--contig", default="contig1", help="name of contig in VCF")
    # parser.add_argument("--demography", choices=["human", "sawtooth"])
    parser.add_argument("--demography", choices=["human"])
    parser.add_argument("-o", help="output location (default: stdout)")
    parser.add_argument("-T", help="Print the simulated local genealogies in Newick format.")
    parser.add_argument("n", type=int, help="diploid sample size")
    parser.add_argument("t", type=float, help="mutation rate")
    parser.add_argument("rho", type=float, help="recombination rate")
    parser.add_argument("length", type=int, help="length of chromosome to simulate")

    args, scrm_extra_args = parser.parse_known_args()

    if args.o is None:
        out = sys.stdout
    else:
        out = open(args.o, "wt")

    if args.demography is not None:
        # demo = getattr(smcpp.util, args.demography)
        demo = {
            "a": np.array([10.0, 0.5, 1.0, 4.0]),
            "b": np.array([1.0, 0.5, 1.0, 4.0]),
            "s": np.array([10000., 70000. - 10000., 200000. - 70000., 1.0]) / 20000. / 29.0,
            "N0": 10000.,
        }

        a = demo['a']
        b = demo['b']
        s = demo['s'] * 0.5
        scrm_extra_args += demography_from_params((a, b, s))

    scrm_args = [2 * args.n, 1]
    scrm_args.append("--transpose-segsites")
    scrm_args += ["-SC", "abs", "-p", 14]
    scrm_args += ["-t", args.t]
    scrm_args += ["-r", args.rho, args.length]

    # Create a minimal VCF header
    header = ["##fileformat=VCFv4.0", """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""]
    header.append("##contig=<ID={},length={}>".format(args.contig, args.length))
    h = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()
    h += ["sample%d" % i for i in range(1, args.n + 1)]
    header.append("\t".join(h))
    print("\n".join(header), file=out)

    # Iterate over scrm output
    print("Calling scrm with args: %s" % str(scrm_args), file=sys.stderr)
    it = scrm(*scrm_args, _iter=True)
    line = next(it)
    while not line.startswith("position"):
        line = next(it)
    next(it)
    for line in it:
        ary = line.strip().split()
        pos, time = ary[:2]
        gts = ary[2:]
        cols = [args.contig, str(int(float(pos))), ".", "A", "C", ".", "PASS", ".", "GT"]
        cols += ["|".join(gt) for gt in zip(gts[::2], gts[1::2])]
        print("\t".join(cols), file=out)
    out.close()
