# pbwt-founder-population

## Requirements:

Install requirements with `brew`:
```
brew install samtools
brew install htslib  # includes bgzip
brew install bcftools
brew install scrm
brew install vcftools
```

To build `SHAPEIT4` in OS X, modify the `makefile`:
```
CXX=g++ -std=c++11

#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
HTSLIB_INC=/usr/local/include/htslib
HTSLIB_LIB=/usr/local/lib/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
BOOST_INC=/usr/local/include/boost/
BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a

#COMPILER & LINKER FLAGS

#Best performance is achieved with this. Use it if running on the same platform you're compiling, it's definitely worth it!
#CXXFLAG=-O3 -march=native
#Good performance and portable on most intel CPUs
CXXFLAG=-O3 -mavx2 -mfma 
#Portable version without avx2 (much slower)
#CXXFLAG=-O3

LDFLAG=-O3

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma -lcurl

#SHAPEIT SOURCES & BINARY
BFILE=bin/shapeit4
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#COMPILATION RULES
all: $(BFILE)

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
```

## Simulate a VCF Panel Using `scrm2vcf.py`

Generate a simulated phased VCF human panel with `scrm2vcf.py`:

```
$ python scrm2vcf.py <n> <t> <rho> <length> --demography human -o <output path>
```

Parameters:
- _n_ — The present day size of population _i_, _n_*_N0_.
- _t_ — The mutation rate _t_ = 4*_N0_*_mu_, where _mu_ is the neutral mutation rate per locus.
- _rho_ & _length_ — The recombination rate _R_ & locus length _L_.

Where _N0_ is the effective population size & _mu_ is the per-generation mutation rate of the population.

For example:
```
$ python scrm2vcf.py 100 10 20000 50000000 --demography human -o ../resources/panels/scrm/m200.vcf
> Calling scrm with args: [200, 1, '--transpose-segsites', '-SC', 'abs', '-p', 14, '-t', 10.0, '-r', 20000.0, 50000000, '-eN', 0.0, 10.0, '-eG', 0.0, 267.0998707873093, '-eN', 0.008620689655172414, 0.5, '-eG', 0.008620689655172414, 0.0, '-eN', 0.06034482758620689, 1.0, '-eG', 0.06034482758620689, 0.0, '-eN', 0.1724137931034483, 4.0, '-eG', 0.1724137931034483, 0.0]
```

Bgzip compress a VCF file:
```
$ bgzip -c unphased.vcf > unphased.vcf.gz
```

Index a panel using `bcftools`:
```
$ bcftools index unphased.vcf.gz
```

Phase a file using `SHAPEIT4`:
```
$ ./shapeit4/bin/shapeit4 --input unphased.vcf.gz --region contig1 --output shapeit_phased.vcf.gz
```

Find switch errors in a `SHAPEIT` phased panel:
```
$ vcftools --vcf true_phased.vcf --diff shapeit_phased.vcf --diff-switch-error
```
The outputs are stored in `out.diff.indiv.switch` & `out.diff.switch`.

## Extracting Data From VCF Files in Python using `scikit-allel`:
[http://alimanfoo.github.io/2017/06/14/read-vcf.html][extract-vcf-data-python]

[extract-vcf-data-python]: http://alimanfoo.github.io/2017/06/14/read-vcf.html