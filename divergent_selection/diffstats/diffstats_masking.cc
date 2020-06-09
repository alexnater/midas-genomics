/*
  Calculate normalized summary statistics for two populations from ms-like data.
  g++ -DNDEBUG -I. -I.. -I$HOME/libs/include -O3 -std=c++11 -o diffstats_masking diffstats_masking.cc -L$HOME/libs/lib64 -L$HOME/libs/lib -lsequence -lgsl -lgslcblas -lz
*/

#include <fstream>
#include <iostream>
#include <array>
#include <cassert>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <Sequence/Fasta.hpp>
#include <mask.hpp>
#include <sumstats.hpp>


int main(int argc, char ** argv)
{
	if (argc != 13)
	{
	std::cerr << "Too few arguments.\n"
			  << "mask_file(fasta) region_file(bed) "
			  << "locuslength nwindows "
			  << "samplesize1 samplesize2 "
			  << "minpropvalid maxmissing data_is_phased "
			  << "min_alternative_allel_count hapstat_window_proportion seed\n";
	exit(1);
	}
	unsigned argument = 1;
	const char * maskfile = argv[argument++];
	const char * bedfile = argv[argument++];
 	const unsigned locuslength = atoi(argv[argument++]);
 	const unsigned nwindows = atoi(argv[argument++]);
 	const unsigned samplesize1 = atoi(argv[argument++]);
 	const unsigned samplesize2 = atoi(argv[argument++]);
	const double minpropvalid = atof(argv[argument++]);
	const unsigned maxmissing = atoi(argv[argument++]);
	const bool is_phased = (atoi(argv[argument++]) == 1);
 	const unsigned minaltc = atoi(argv[argument++]);
 	const double hap_prop = atof(argv[argument++]);
 	const unsigned seed = atoi(argv[argument++]);

 	// read mask file:
	std::ifstream maskin(maskfile);
	mask mymask(maskin, locuslength, seed);
	maskin.close();

 	// read bed file:
 	if (std::string(bedfile) != "none")
 	{
		std::ifstream bedin(bedfile);
		unsigned nloci = mymask.read_bed(bedin);
		std::cerr << "Read in " << nloci << " loci from bedfile." << std::endl;
		bedin.close();
	}
	
	// read data in ms format from stdin:
	unsigned countit = 0;
	const std::array<unsigned, 2> samplesizes = {samplesize1, samplesize2};

	do
	{
		++countit;
		std::cerr << "Working on iteration " << countit << " ..." << std::endl;
		std::vector<unsigned> validsites;
		if (std::string(bedfile) != "none")
		{
			validsites = mymask.set_mask(countit - 1, nwindows);
			if (!mymask.subwindows_are_valid(validsites, locuslength / nwindows, minpropvalid, maxmissing))
			{
				validsites.assign(nwindows, 0);
				std::cerr << "Locus " << countit << " has not enough valid sites." << std::endl;
			}
		}
		else validsites = mymask.set_random_mask(locuslength, nwindows, minpropvalid, maxmissing);
		auto vm = Sequence::from_msformat(std::cin);
		std::cerr << "Generated variant matrix with " << vm.nsam() << " samples and " << vm.nsites() << " sites." << std::endl;
		auto filtered_vm = filter_positions(vm, mymask);
		if (is_phased)
			StandardizedWindows(std::cout, filtered_vm, samplesizes, nwindows, locuslength,
								validsites, countit, mymask.get_current_region(), minaltc, hap_prop);
		else
			StandardizedWindowsUnphased(std::cout, filtered_vm, samplesizes, nwindows, locuslength,
										validsites, countit, mymask.get_current_region());			
	}
	while (!std::cin.eof());
}

