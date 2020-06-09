#ifndef MASK_H
#define MASK_H

#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <Sequence/Fasta.hpp>


struct region
{
	std::string id;
	std::size_t start;
	std::size_t end;
	region(){};
	region(const std::string & id_, std::size_t start_, std::size_t end_)
		   : id{id_}, start{start_}, end{end_}{}
};


class mask
{
	private:
		using gsl_rng_ptr_t = std::unique_ptr<gsl_rng, void (*)(gsl_rng *)>;
		using gsl_ran_discrete_t_ptr = std::unique_ptr<gsl_ran_discrete_t, void (*)(gsl_ran_discrete_t *)>;
		std::vector<Sequence::Fasta> seqs;
		std::unordered_map<std::string, std::size_t> names;
		std::vector<region> regions;
		std::vector<std::int8_t> vmask;
		region current_region;
		gsl_rng_ptr_t rng;
		gsl_ran_discrete_t_ptr lookup;
		
	public:
		mask(std::ifstream & infile, const unsigned minlength, const unsigned long seed, unsigned nseqs = 10)
			: rng(gsl_rng_alloc(gsl_rng_mt19937), [](gsl_rng *r){ gsl_rng_free(r); }),
			  lookup(nullptr, [](gsl_ran_discrete_t * l){ gsl_ran_discrete_free(l); })
		{
			gsl_rng_set(rng.get(), seed);
			seqs.reserve(nseqs);
			std::vector<double> lengths;
			lengths.reserve(nseqs);
			Sequence::Fasta seq;
			while (infile >> seq)
			{
				if (seq.length() > minlength)
				{
					std::cerr << "Reading in sequence " << seq.name << ", length " << seq.length() << " ..." << std::endl;
					lengths.push_back(static_cast<double>(seq.length() - minlength));
					names[seq.name] = seqs.size();
					seqs.push_back(std::move(seq));
				}
			}
			lookup.reset(gsl_ran_discrete_preproc(lengths.size(), lengths.data()));
		}

		unsigned read_bed(std::ifstream & bed, unsigned nrecords = 100)
		{
			regions.clear();
			regions.reserve(nrecords);
			std::string id;
			std::size_t start;
			std::size_t end;
			unsigned nloci = 0;
			while (bed >> id >> start >> end)
			{
				regions.emplace_back(id, start, end);
				++nloci;
			}
			return nloci;
		}

		bool subwindows_are_valid(const std::vector<unsigned> & nvalid, const unsigned winlength, const double minprop)
		{
			return std::all_of(nvalid.begin(), nvalid.end(), [winlength, minprop](const unsigned val)
								{ return (val / static_cast<double>(winlength) >= minprop); });
		}

		bool subwindows_are_valid(std::vector<unsigned> & nvalid, const unsigned winlength,
								  const double minprop, const unsigned maxmiss)
		{
			unsigned nmissing = 0;
			bool prevmissing = false;
			for (std::size_t i = 0; i < nvalid.size(); ++i)
			{
				if (nvalid[i] / static_cast<double>(winlength) < minprop)
				{
					++nmissing;
					nvalid[i] = 0;	// set number of valid sites to 0 to mark subwindow for interpolation.
					if (i == 0 || i == nvalid.size() - 1 || prevmissing || nmissing > maxmiss) return false;
					prevmissing = true;
				}
				else prevmissing = false;
			}
			return true;
		}

		region get_current_region()
		{
			return current_region;
		}

		region get_random_start(unsigned length)
		{
			std::size_t chridx = gsl_ran_discrete(rng.get(), lookup.get());
			assert(length < seqs[chridx].length());
			std::size_t start = gsl_rng_uniform_int(rng.get(), seqs[chridx].length() - length);
			return region(seqs[chridx].name, start, start + length);
		}

		std::string get_random_seq(unsigned length)
		{
			auto reg = get_random_start(length);
			return seqs[names[reg.id]].substr(reg.start, length);
		}

		unsigned set_mask(region & reg)
		{
			std::cerr << "Setting mask for region " << reg.id << ":" << reg.start + 1 << "-" << reg.end << "." << std::endl;
			current_region = reg;
			std::string seq = seqs[names[reg.id]].substr(reg.start, reg.end - reg.start);
			unsigned nvalid = 0;
			vmask = std::vector<std::int8_t>(seq.length(), 1);
			for (std::string::size_type i = 0; i < seq.length(); ++i)
			{
				if (std::toupper(seq[i]) == 'N' || seq[i] == '0') vmask[i] = 0;
				else ++nvalid;
			}
			return nvalid;
		}

		// overload for subwindows:
		std::vector<unsigned> set_mask(region & reg, const unsigned nwindows)
		{
			std::cerr << "Setting mask for region " << reg.id << ":" << reg.start << "-" << reg.end - 1 << "." << std::endl;
			current_region = reg;
			std::string seq = seqs[names[reg.id]].substr(reg.start, reg.end - reg.start);
			std::vector<unsigned> nvalid(nwindows, 0);
			vmask = std::vector<std::int8_t>(seq.length(), 1);
			double winsize = seq.length() / static_cast<double>(nwindows);
			double winstart = 0.;
		    for (unsigned win = 0; win < nwindows; ++win)
		    {
				for (unsigned idx = static_cast<unsigned>(winstart);
				              idx < static_cast<unsigned>(winstart + winsize); ++idx)
				{
					if (std::toupper(seq[idx]) == 'N' || seq[idx] == '0') vmask[idx] = 0;
					else ++nvalid[win];
				}
				std::cerr << "Subwindow " << static_cast<unsigned>(winstart)  << "-" << static_cast<unsigned>(winstart + winsize)
						  << ", proportion of valid sites: " << nvalid[win] / winsize << std::endl;
				winstart += winsize;
			}
			return nvalid;
		}
		
		unsigned set_mask(std::size_t index)
		{
			if (index >= regions.size()) throw std::invalid_argument("index out of range");
			return set_mask(regions[index]);
		}

		// overload for subwindows:
		std::vector<unsigned> set_mask(std::size_t index, const unsigned nwindows)
		{
			if (index >= regions.size()) throw std::invalid_argument("index out of range");
			return set_mask(regions[index], nwindows);
		}
		
		unsigned set_random_mask(const unsigned length, const double minprop = 0.)
		{
			unsigned nvalid = 0;
			while (!nvalid || nvalid / static_cast<double>(length) < minprop)
			{
				auto reg = get_random_start(length);
				assert(length == reg.end - reg.start);
				nvalid = set_mask(reg);
			}
			return nvalid;
		}

		// overload for subwindows:
		std::vector<unsigned> set_random_mask(const unsigned length, const unsigned nwindows,
											  const double minprop = 0., const unsigned maxmiss = 0)
		{
			unsigned winlength = length / nwindows;
			std::vector<unsigned> nvalid(nwindows, 0);
			while (!subwindows_are_valid(nvalid, winlength, minprop, maxmiss))
			{
				auto reg = get_random_start(length);
				assert(length == reg.end - reg.start);
				nvalid = set_mask(reg, nwindows);
			}
			return nvalid;
		}
		
		inline bool operator()(const double pos) const
		{
			// assuming pos in range [0,1)
			std::size_t idx(pos * vmask.size());
			if (idx == vmask.size()) --idx;	// msms seems to produce both 0.00 and 1.00 sites.
			assert(idx < vmask.size());
			return static_cast<bool>(vmask[idx]);		
		}
};


#endif

