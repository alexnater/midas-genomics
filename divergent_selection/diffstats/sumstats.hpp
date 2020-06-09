#ifndef SUMSTATS_H
#define SUMSTATS_H

#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/StateCounts.hpp>
#include <Sequence/variant_matrix/windows.hpp>
#include <Sequence/variant_matrix/filtering.hpp>
#include <Sequence/summstats/classics.hpp>
#include <Sequence/summstats/garud.hpp>
#include <Sequence/summstats/lhaf.hpp>
#include <Sequence/FST.hpp>
#include <Sequence/SummStats/nSL.hpp>
#include <Sequence/SummStats/Garud.hpp>
#include <mask.hpp>


using const_data_iterator = Sequence::SimData::const_data_iterator;


template<typename iterator_t>
std::string join(iterator_t begin, iterator_t end, const std::string sep = ", ")
{
  std::ostringstream ss;
  if (begin != end) ss << *begin++;
  while (begin != end)
  {
    ss << sep;
    ss << *begin++;
  }
  return ss.str();
}

template <typename iterator>
inline std::pair<iterator, iterator> mismatch_skip_missing(iterator beg, iterator end, iterator beg2)
{
	auto m = std::mismatch(beg, end, beg2);
	while (m.first < end && (*m.first < 0 || *m.second < 0))
	{
		m = std::mismatch(m.first + 1, end, m.second + 1);
	}
	return m;
}

template <typename iterator>
inline std::int32_t ndiff_skip_missing(iterator beg, iterator end, iterator beg2)
{
	std::int32_t ndiffs = 0;
	auto m = mismatch_skip_missing(beg, end, beg2);
	while (m.first < end)
	{
		++ndiffs;
		m = mismatch_skip_missing(m.first + 1, end, m.second + 1);
	}
	return ndiffs;
}


// Check if this could be solved with libsequence built in masking.
Sequence::VariantMatrix filter_positions(const Sequence::VariantMatrix & vm, const mask & mymask)
{
	std::vector<std::int8_t> data;
	std::vector<double> positions;
	data.reserve(vm.nsites() * vm.nsam());
	positions.reserve(vm.nsites());
	for (auto i = 0; i < vm.nsites(); ++i)
	{
		const double pos = vm.position(i);
		if (mymask(pos))
		{
			positions.push_back(pos);
			auto row = get_ConstRowView(vm, i);
			data.insert(data.end(), row.begin(), row.end());
		}
	}
	std::cerr << "Retained " << positions.size() << " out of " << vm.nsites() << " variant positions." << std::endl;
	return Sequence::VariantMatrix(std::move(data), std::move(positions));
}


// only subset columns of VariantMatrix. only needed if remaining individuals are not in a single range. use make_slice otherwise:
Sequence::VariantMatrix subset_individuals(const Sequence::VariantMatrix& vm, const std::vector<std::pair<std::size_t, std::size_t> > & indices)
{
	std::vector<std::int8_t> data;
	for (auto i = 0; i < vm.nsites(); ++i)
	{
		auto row = get_ConstRowView(vm, i);
		for (auto & popi : indices)
		{
			if (!(popi.second > popi.first))
			{
				throw std::invalid_argument("i must be < j");
			}
			if (popi.first >= vm.nsam() || popi.second > vm.nsam())
			{
				throw std::invalid_argument("slice indexes out of range");
			}
			data.insert(data.end(), row.begin() + popi.first, row.begin() + popi.second);
		}
	}
	return Sequence::VariantMatrix(std::move(data), std::vector<double>(vm.cpbegin(), vm.cpend()));
}


Sequence::SimData vm_to_simdata(Sequence::VariantMatrix & vm)
{
	std::vector<std::string> data;
	data.reserve(vm.nsam());
	for (std::size_t i = 0; i < vm.nsam(); ++i)
	{
		std::string haplotype;
		auto col_view = get_ConstColView(vm, i);
		for (auto allele : col_view) haplotype += std::to_string(allele);
		data.push_back(std::move(haplotype));
	}
	return Sequence::SimData(std::vector<double>(vm.pbegin(), vm.pend()), std::move(data));
}


void print_stats(std::ostream & outstream, const std::vector<std::vector<double> > & sumstats,
				 const std::vector<std::string> & ssnames, const unsigned iteration, const region & reg)
{
	if (iteration == 1)
	{
		// print header:
		outstream << "iteration" << '\t' << "chrom" << '\t' << "start" << '\t' << "end";
		for (unsigned i = 0; i < sumstats.size(); ++i)
		{
			for (unsigned j = 0; j < sumstats[i].size(); ++j) outstream << '\t' << ssnames[i] << "_win" << j;
		}
		outstream << '\n';
	}
	// print stats:
	outstream << iteration << '\t' << reg.id << '\t' << reg.start << '\t' << reg.end;
	for (auto & vstats : sumstats)
	{
		for (double ss : vstats) outstream << '\t' << ss;
	}
	outstream << '\n';
	return;
}


template <typename iterator>
double get_mean(iterator beg, iterator end)
{
	double rt = 0.;
	unsigned nc = 0;
	for (; beg != end; ++beg)
	{
		if (std::isfinite(*beg))
		{
			rt += *beg;
			++nc;
		}
	}
	if (!nc) return std::numeric_limits<double>::quiet_NaN();
	return rt / static_cast<double>(nc);
}


void interpolate_windows(std::vector<double> & stats)
{
	for (std::size_t i = 0; i < stats.size(); ++i)
	{
		if (!std::isfinite(stats[i]))
		{
			if (i == 0 || i == stats.size() - 1 ||
				!std::isfinite(stats[i - 1]) ||
				!std::isfinite(stats[i + 1]))
					throw std::runtime_error("Invalid window values for interpolation in window " + std::to_string(i) + "!");
			stats[i] = (stats[i - 1] + stats[i + 1]) / static_cast<double>(2);
		}
	}
	return;
}


std::vector<double> normalize_windows(std::vector<double> & stats)
{
	double minval = *std::min_element(stats.begin(), stats.end());
	std::vector<double> normstats(stats);
	if (minval < 0.) std::transform(normstats.begin(), normstats.end(), normstats.begin(), [minval](double stat){ return stat - minval; });
	double sum_of_stats = std::accumulate(normstats.begin(), normstats.end(), 0.);
	if (sum_of_stats == 0.)
	{
		normstats = std::vector<double>(normstats.size(), 1.0 / normstats.size());
	}
	else
	{
		for (double & stat : normstats) stat /= sum_of_stats;
	}
	return normstats;
}


double calculate_hapdiff(Sequence::VariantMatrix & vm, const unsigned popsize1, const unsigned popsize2)
{
	if (vm.nsam() != popsize1 + popsize2)
	{
		throw std::runtime_error("ERROR: Number of samples in variant matrix doesn't match sum of population sizes! ("
                                 + std::to_string(vm.nsam()) + " vs. " + std::to_string(popsize1 + popsize2) + ")");
	}
	auto labels = Sequence::label_haplotypes(vm);
	std::unordered_map<std::int32_t, std::pair<std::size_t, std::size_t> > counts;
	std::pair<std::size_t, std::size_t> nmissing(0, 0);
	for (unsigned i = 0; i < popsize1; ++i)
	{
		if (labels[i] < 0) nmissing.first++;
		else if (counts.find(labels[i]) != counts.end()) counts[labels[i]].first++;
		else counts[labels[i]] = {1, 0};
	}
	for (unsigned i = popsize1; i < popsize1 + popsize2; ++i)
	{
		if (labels[i] < 0) nmissing.second++;
		else if (counts.find(labels[i]) != counts.end()) counts[labels[i]].second++;
		else counts[labels[i]] = {0, 1};
	}
	
	double hapdiff = 0.;
	const unsigned n1 = popsize1 - nmissing.first, n2 = popsize2 - nmissing.second;
	for (auto & hc : counts)
	{
		double freq1 = hc.second.first / static_cast<double>(n1);
		double freq2 = hc.second.second / static_cast<double>(n2);
		hapdiff += std::pow(freq1 - freq2, 2);
	}
	return hapdiff;
}


double calculate_gmin(Sequence::VariantMatrix & vm, const unsigned popsize1, const unsigned popsize2)
{
	if (vm.nsam() != popsize1 + popsize2)
	{
		throw std::runtime_error("ERROR: Number of samples in variant matrix doesn't match sum of population sizes! ("
                                 + std::to_string(vm.nsam()) + " vs. " + std::to_string(popsize1 + popsize2) + ")");
	}
	std::vector<unsigned> diffs(popsize1 * popsize2, 0);
	unsigned comp = 0;
	for (unsigned i = 0; i < popsize1; ++i)
	{
		auto hap1 = get_ConstColView(vm, i);
		for (unsigned j = popsize1; j < popsize1 + popsize2; ++j)
		{
			auto hap2 = get_ConstColView(vm, j);
			diffs[comp] = ndiff_skip_missing(hap1.cbegin(), hap1.cend(), hap2.cbegin());
			++comp;
		}
	}
	
	double dxy = std::accumulate(diffs.cbegin(), diffs.cend(), 0) / static_cast<double>(comp);
	auto minitr = std::min_element(diffs.cbegin(), diffs.cend());
	return (*minitr) / dxy;
}


void StandardizedWindows(std::ostream & outstream, Sequence::VariantMatrix & vm, 
						 const std::array<unsigned, 2> & samplesizes, const unsigned nwindows,
						 const unsigned chromlength, std::vector<unsigned> & validsites,
						 const unsigned iteration, const region & reg,
						 const unsigned minaltc = 0u, const double hap_prop = 1.0)
{
  try
  {

    // initialize stats containers:
    std::vector<std::string> ssnames = {"pi_1", "pi_2", "pi_tot",
    									"tajD_1", "tajD_2", "tajD_tot",
    									"faywuH_1", "faywuH_2", "faywuH_tot",
    									"H1_1", "H1_2", "H1_tot",
    									"H2H1_1", "H2H1_2", "H2H1_tot",
    									"H12_1", "H12_2", "H12_tot",
    									"SS-H12", "fst", "dxy", "gmin",
    									"1-HAF_1", "1-HAF_2"};
    std::vector<std::vector<double> > sumstats(ssnames.size(), std::vector<double>(nwindows, std::numeric_limits<double>::quiet_NaN()));
    
    // take shortcut if none of the subwindows has data:
	if (validsites.empty()) validsites.resize(nwindows, chromlength / nwindows);
	else if (validsites.size() != nwindows) throw std::runtime_error("Wrong size of validsites vector!");
    else if (std::all_of(validsites.begin(), validsites.end(), [](const unsigned val){ return (val == 0); }))
    {
    	print_stats(outstream, sumstats, ssnames, iteration, reg);
    	std::cerr << "Skipping locus." << std::endl;
    	return;
    }

	// define sample configurations:
    unsigned config12[2] = {samplesizes[0], samplesizes[1]};

    // loop through the windows:
    double start = 0.;
    double stepsize = 1. / static_cast<double>(nwindows);
    for (unsigned i = 0; i < nwindows; ++i)
	{
		auto win_vm = Sequence::make_window(vm, start, start + stepsize);
	    double hstart = start + ((stepsize * (1 - hap_prop)) / 2.);
	    double hend = hstart + (stepsize * hap_prop);
	    auto hwin_vm = Sequence::make_window(vm, hstart, hend);
		std::cerr << "Window: " << i << " - Number of SNPs: " << win_vm.nsites() << "/" << hwin_vm.nsites() << std::endl;

		if (validsites[i] > 0 && win_vm.nsites() > 0 && hwin_vm.nsites() > 0)
		{
	        // create VariantMatrix objects for all pop comparisons:
	        auto win_vm1 = subset_individuals(win_vm, {std::pair<std::size_t, std::size_t>(0, samplesizes[0])});
	        auto win_vm2 = subset_individuals(win_vm, {std::pair<std::size_t, std::size_t>(samplesizes[0], vm.nsam())});
			
	        auto hwin_vm1 = subset_individuals(hwin_vm, {std::pair<std::size_t, std::size_t>(0, samplesizes[0])});
	        auto hwin_vm2 = subset_individuals(hwin_vm, {std::pair<std::size_t, std::size_t>(samplesizes[0], vm.nsam())});
	        
			if (minaltc)
			{
				Sequence::StateCounts c1;
				auto f1 = Sequence::filter_sites(hwin_vm1, [&c1, minaltc](const Sequence::RowView & r){ c1(r); return c1.counts[1] < minaltc; });
				Sequence::StateCounts c2;
				auto f2 = Sequence::filter_sites(hwin_vm2, [&c2, minaltc](const Sequence::RowView & r){ c2(r); return c2.counts[1] < minaltc; });
				Sequence::StateCounts c12;
				auto f12 = Sequence::filter_sites(hwin_vm, [&c12, minaltc](const Sequence::RowView & r){ c12(r); return c12.counts[1] < minaltc; });
				std::cerr << "Filtered " << f1 << "/" << f2 << "/" << f12 << " sites." << std::endl;
			}
			
			// calculate diversity-based statistics:
			Sequence::AlleleCountMatrix ac1(win_vm1);
			Sequence::AlleleCountMatrix ac2(win_vm2);
			Sequence::AlleleCountMatrix ac12(win_vm);
			sumstats[0][i] = Sequence::thetapi(ac1) / static_cast<double>(validsites[i]);
			sumstats[1][i] = Sequence::thetapi(ac2) / static_cast<double>(validsites[i]);
			sumstats[2][i] = Sequence::thetapi(ac12) / static_cast<double>(validsites[i]);
			sumstats[3][i] = Sequence::tajd(ac1);
			sumstats[4][i] = Sequence::tajd(ac2);
			sumstats[5][i] = Sequence::tajd(ac12);
			sumstats[6][i] = Sequence::faywuh(ac1, 0);
			sumstats[7][i] = Sequence::faywuh(ac2, 0);
			sumstats[8][i] = Sequence::faywuh(ac12, 0);
			
	        // generate Garud statistics:
	        auto garud1 = Sequence::garud_statistics(hwin_vm1);
	        auto garud2 = Sequence::garud_statistics(hwin_vm2);
	        auto garud12 = Sequence::garud_statistics(hwin_vm);
	        double h12_anc = garud12.H12 - calculate_hapdiff(hwin_vm, samplesizes[0], samplesizes[1]);
	        auto cfactor = std::minmax(garud1.H12, garud2.H12);
	        sumstats[9][i] = garud1.H1;
	        sumstats[10][i] = garud2.H1;
	        sumstats[11][i] = garud12.H1;
	        sumstats[12][i] = garud1.H2H1;
	        sumstats[13][i] = garud2.H2H1;
	        sumstats[14][i] = garud12.H2H1;
	        sumstats[15][i] = garud1.H12;
	        sumstats[16][i] = garud2.H12;
	        sumstats[17][i] = garud12.H12;
	        sumstats[18][i] = h12_anc * (cfactor.first / cfactor.second);
			
			// calculate Fst-based statistics:
			auto wpop12 = vm_to_simdata(win_vm);
			Sequence::FST fst12(&wpop12, 2, config12);
	        sumstats[19][i] = fst12.HSM();
			sumstats[20][i] = fst12.piB() / static_cast<double>(validsites[i]);
			sumstats[21][i] = calculate_gmin(win_vm, samplesizes[0], samplesizes[1]);
			
			// calculate mean 1-HAF statistic per population:
			std::vector<double> lhaf = Sequence::lhaf(win_vm, 0, 1);
			sumstats[22][i] = get_mean(lhaf.begin(), lhaf.begin() + samplesizes[0]);
			sumstats[23][i] = get_mean(lhaf.begin() + samplesizes[0], lhaf.end());
		}
		start += stepsize;
	}

	// normalize stats:
	std::vector<std::vector<double> > normstats(ssnames.size(), std::vector<double>(nwindows, std::numeric_limits<double>::quiet_NaN()));
	for (unsigned i = 0; i < sumstats.size(); ++i)
	{
		try
		{
			interpolate_windows(sumstats[i]);
			normstats[i] = normalize_windows(sumstats[i]);
		} catch (std::exception & error)
		{
			std::cerr << error.what() << std::endl;
			std::cerr << "Sumstat " << ssnames[i] << ": " << join(sumstats[i].begin(), sumstats[i].end()) << std::endl;
		}
	}
	
	// print stats:
	print_stats(outstream, normstats, ssnames, iteration, reg);
		
  } catch (std::exception & error){ std::cerr << error.what() <<std::endl; }
  return;
}

void StandardizedWindowsUnphased(std::ostream & outstream, Sequence::VariantMatrix & vm, 
						 		 const std::array<unsigned, 2> & samplesizes, const unsigned nwindows,
						 		 const unsigned chromlength, std::vector<unsigned> & validsites,
						 		 const unsigned iteration, const region & reg)
{
  try
  {

    // initialize stats containers:
    std::vector<std::string> ssnames = {"pi_1", "pi_2", "pi_tot",
    									"tajD_1", "tajD_2", "tajD_tot",
    									"faywuH_1", "faywuH_2", "faywuH_tot",
    									"fst", "dxy"};
    std::vector<std::vector<double> > sumstats(ssnames.size(), std::vector<double>(nwindows, std::numeric_limits<double>::quiet_NaN()));
    
    // take shortcut if none of the subwindows has data:
	if (validsites.empty()) validsites.resize(nwindows, chromlength / nwindows);
	else if (validsites.size() != nwindows) throw std::runtime_error("Wrong size of validsites vector!");
    else if (std::all_of(validsites.begin(), validsites.end(), [](const unsigned val){ return (val == 0); }))
    {
    	print_stats(outstream, sumstats, ssnames, iteration, reg);
    	std::cerr << "Skipping locus." << std::endl;
    	return;
    }

	// define sample configurations:
    unsigned config12[2] = {samplesizes[0], samplesizes[1]};

    // loop through the windows:
    double start = 0.;
    double stepsize = 1. / static_cast<double>(nwindows);
    for (unsigned i = 0; i < nwindows; ++i)
	{
		auto win_vm = Sequence::make_window(vm, start, start + stepsize);
		std::cerr << "Window: " << i << " - Number of SNPs: " << win_vm.nsites() << std::endl;

		if (validsites[i] > 0 && win_vm.nsites() > 0)
		{
	        // create VariantMatrix objects for all pop comparisons:
	        auto win_vm1 = subset_individuals(win_vm, {std::pair<std::size_t, std::size_t>(0, samplesizes[0])});
	        auto win_vm2 = subset_individuals(win_vm, {std::pair<std::size_t, std::size_t>(samplesizes[0], vm.nsam())});
	        
			// calculate diversity-based statistics:
			Sequence::AlleleCountMatrix ac1(win_vm1);
			Sequence::AlleleCountMatrix ac2(win_vm2);
			Sequence::AlleleCountMatrix ac12(win_vm);
			sumstats[0][i] = Sequence::thetapi(ac1) / static_cast<double>(validsites[i]);
			sumstats[1][i] = Sequence::thetapi(ac2) / static_cast<double>(validsites[i]);
			sumstats[2][i] = Sequence::thetapi(ac12) / static_cast<double>(validsites[i]);
			sumstats[3][i] = Sequence::tajd(ac1);
			sumstats[4][i] = Sequence::tajd(ac2);
			sumstats[5][i] = Sequence::tajd(ac12);
			sumstats[6][i] = Sequence::faywuh(ac1, 0);
			sumstats[7][i] = Sequence::faywuh(ac2, 0);
			sumstats[8][i] = Sequence::faywuh(ac12, 0);
			
			// calculate Fst-based statistics:
			auto wpop12 = vm_to_simdata(win_vm);
			Sequence::FST fst12(&wpop12, 2, config12);
	        sumstats[9][i] = fst12.HSM();
			sumstats[10][i] = fst12.piB() / static_cast<double>(validsites[i]);
		}
		start += stepsize;
	}

	// normalize stats:
	std::vector<std::vector<double> > normstats(ssnames.size(), std::vector<double>(nwindows, std::numeric_limits<double>::quiet_NaN()));
	for (unsigned i = 0; i < sumstats.size(); ++i)
	{
		try
		{
			interpolate_windows(sumstats[i]);
			normstats[i] = normalize_windows(sumstats[i]);
		} catch (std::exception & error)
		{
			std::cerr << error.what() << std::endl;
			std::cerr << "Sumstat " << ssnames[i] << ": " << join(sumstats[i].begin(), sumstats[i].end()) << std::endl;
		}
	}
	
	// print stats:
	print_stats(outstream, normstats, ssnames, iteration, reg);
		
  } catch (std::exception & error){ std::cerr << error.what() <<std::endl; }
  return;
}


#endif

