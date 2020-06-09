#!/usr/bin/env python3

"""
Takes a file with sample information (sample id, P/F1/F2, family id) and a VCF file and outputs a filtered VCF or JoinMap file with a/b encoding.

"""

import argparse
import sys
import re
import scipy.stats as stats

encoding = ['a', 'b', '.']
genotype_encoding = {('a','a'):'a', ('a','b'):'h', ('b','a'):'h', ('b','b'):'b', ('.','.'):'-'}


def read_samplelist(samplefile):
    """
    Reads in a text file with sample information (sample id, P/F1/F2, family id).
    
    Takes a string with the full path to the text file.
    
    Returns a dict of dict with list of sample ids by generation and family id.
    """
    with open(samplefile, 'r') as inhandle:
        samples = {}
        for line in inhandle:
            if not line.startswith("#"):
                fields = line.strip().split()
                if fields[1] == 'P':
                    if not len(fields) >= 2:
                        raise Exception("ERROR: Incorrectly formated sample file:\n{}".format(line) )
                    samples.setdefault('P', {})
                    samples['P'].setdefault('NA', []).append(fields[0])
                elif fields[1] == 'F1' or fields[1] == 'F2':
                    if not len(fields) == 3:
                        raise Exception("ERROR: Incorrectly formated sample file:\n{}".format(line) )
                    samples.setdefault(fields[1], {})
                    samples[fields[1]].setdefault(fields[2], []).append(fields[0])
                else:
                    raise Exception("ERROR: Invalid category {} specified for sample {}!".format(fields[1], fields[0]) )
    if len(samples['P']['NA']) != 2:
        raise Exception("ERROR: Need to specify two parental samples!")
    if 'F1' in samples:
        for key,val in samples['F1'].items():
            if len(val) != 2:
                raise Exception("ERROR: Need to specify two samples for F1 family {}!".format(key) )
    else:
        samples['F1'] = {}
    return samples


def get_sampleidx(samples, sampleids, startcol = 9):
    """Takes a samples dict and list of sample ids and returns a dict of dicts with lists of sample indices by generation and family."""
    sample_map = {x[1]: x[0] + startcol for x in enumerate(sampleids)}
    sampleidx = {}
    for gen in samples.keys():
        sampleidx[gen] = {}
        for famidx in samples[gen].keys():
            sampleidx[gen][famidx] = [ sample_map[sampleid] for sampleid in samples[gen][famidx] ]
    return sampleidx


def check_site(fields, args, fidx):
    """Checks if site fulfills filtering criteria."""
    reflen = len(fields[3])
    if args.only_snps and reflen > 1:
        return False
    if args.excl_indels and any(len(alt) != reflen for alt in fields[4].split(',')):
        return False
    if args.min_vq and (fields[5] == '.' or float(fields[5]) < args.min_vq):
        return False
    if args.min_dp:
        mobj = re.search('DP=(\d+)', fields[7])
        if not mobj:
            raise Exception("ERROR: Incorrectly formated line in vcf file:\n{}".format("\t".join(fields)) )
        elif int(mobj.group(1)) < args.min_dp:
            return False
    if args.ind_dp:
        if not ('DP' in fidx):
            return False
    if args.ind_gq:
        if not ('GQ' in fidx):
            return False 
    return True



class Genotypes:
    """Class to represent parsed genotypes from a single generation."""
    def __init__(self, indices, dpfilter = 0, gqfilter = 0, glfilter = 0, separator = '/'):
        self.indices = indices
        self.gts = []
        self.alleles = set()
        self.mindp = dpfilter
        self.mingq = gqfilter
        self.mingl = glfilter
        self.sep = separator
        self.ntot = 0
        self.nvalid = 0
        self.nmissing = 0
        self.nonparental = 0
    
    def parse_genotypes(self, fidx, fields, p_alleles = None):
        self.gts = []
        for idx in self.indices:
            try:
                indstring = fields[idx]
            except IndexError:
                print("Error on line {}:".format(sys.exc_info()[-1].tb_lineno), file=sys.stderr)
                print("WARNING: Index {} does not exist in list of fields!".format(idx), file=sys.stderr)
                continue
            
            indfields = indstring.strip().split(':')
            if len(indfields) < len(fidx):
            	self.nmissing += 1
            elif self.mindp and (indfields[fidx['DP']] == '.' or int(indfields[fidx['DP']]) < self.mindp):
                self.nmissing += 1
            elif self.mingq and (indfields[fidx['GC']] == '.' or int(indfields[fidx['GC']]) < self.mingq):
                self.nmissing += 1
            elif self.mingl and (indfields[fidx['GL']] == '.' or self.get_gldiff(indfields[fidx['GL']]) < self.mingl):
                self.nmissing += 1
            else:
                gt = indfields[0].split(self.sep)
                if '.' in gt:
                    self.nmissing += 1
                else:
                    gt = [ int(al) for al in gt ]
                    if p_alleles and not p_alleles.issuperset(gt):
                        self.nonparental += 1
                        self.nmissing += 1
                    else:
                        self.nvalid += 1
                        self.alleles.update(gt)
                        self.gts.append(gt)
            self.ntot += 1
    
    def get_alleles(self):
        return self.alleles
    
    def get_indices(self):
        return tuple(self.alleles)
    
    def get_nalleles(self):
        return len(self.alleles)
    
    def get_nmissing(self):
        return self.nmissing
    
    def get_alleles_pair(self):
        if len(self.gts) != 2:
            raise Exception("WARNING: Need to have exactly two non-missing genotypes!")
        return set(self.gts[0]), set(self.gts[1])
        
    def get_shared_pair(self):
        if len(self.gts) != 2:
            raise Exception("WARNING: Need to have exactly two non-missing genotypes!")
        return set(self.gts[0]) & set(self.gts[1])
    
    def get_unique_pair(self):
        if len(self.gts) != 2:
            raise Exception("WARNING: Need to have exactly two non-missing genotypes!")
        alleles = self.get_alleles_pair()
        ualleles = alleles[0] ^ alleles[1]
        return alleles[0] & ualleles, alleles[1] & ualleles, ualleles
    
    def get_shared(self, genotypes):
        return self.alleles & genotypes.alleles
    
    def get_unique(self, genotypes):
        ualleles = self.alleles ^ genotypes.alleles
        return self.alleles & ualleles, genotypes.alleles & ualleles, ualleles
    
    def get_missingness(self):
        return self.nmissing / self.ntot

    def get_nonparental(self):
        return self.nonparental / self.ntot
    
    def check_compatibility(self, parents):
        if len(self.gts) != 2 and len(parents.gts) != 2:
            raise Exception("WARNING: Need to have exactly two non-missing genotypes in parental and F1 genotpyes!")
        f1_alleles = self.get_alleles_pair()
        p_alleles = parents.get_alleles_pair()
        return bool(f1_alleles[0] & p_alleles[0]) and bool(f1_alleles[0] & p_alleles[1]) and bool(f1_alleles[1] & p_alleles[0]) and bool(f1_alleles[1] & p_alleles[1])
    
    def count_genotypes(self, trmap):
        gtcounts = [0, 0, 0]
        for gt in self.gts:
            trgt = (trmap[gt[0]], trmap[gt[1]])
            if '.' in trgt:
                continue
            elif trgt[0] == trgt[1]:
                if trgt[0] == encoding[0]:
                    gtcounts[0] += 1
                else:
                    gtcounts[2] += 1
            else:
                gtcounts[1] += 1
        return gtcounts
    
    def chisquare_test(self, trmap):
        observed = self.count_genotypes(trmap)
        if not all(x >= 5 for x in observed):
            raise Exception("WARNING: Chi-squared test with insufficent valid genotypes per category!")
        nval = sum(observed)
        expected = (0.25*nval, 0.5*nval, 0.25*nval)
        return stats.chisquare(observed, expected, ddof=1), observed
    
    def get_gldiff(self, glstring):
        gls = sorted(map(float, glstring.split(',')), reverse=True)
        if len(gls) < 3:
            raise Exception("ERROR: Invalid length of GL field ({})!".format(glstring))
        return gls[0] - gls[1]


def recode_genotypes(fields, fidx, indices, mindp, trans_map, annotate = False, family = "all", separator = '/'):
    """Recodes selected genotypes of line in VCF file with a/b allele encoding."""
    valid = any(encoding != '.' for encoding in trans_map)
    rec_gts = []
    nmissing = 0
    for idx in indices:
        indfields = fields[idx].strip().split(':')
        gt = ('.', '.')
        if valid and len(indfields) == len(fidx) and indfields[fidx['DP']] != '.' and int(indfields[fidx['DP']]) >= mindp:
            lgt = indfields[0].split(separator)
            if not '.' in lgt:
                gt = (trans_map[int(lgt[0])], trans_map[int(lgt[1])])
        if '.' in gt:
            gt = ('.', '.')
            nmissing += 1
        indfields[0] = "/".join(gt)
        fields[idx] = ":".join(indfields)
        rec_gts.append(genotype_encoding[gt])
        
    if annotate:
        used_map = {org : rec for org,rec in enumerate(trans_map) if rec != '.'}
        text_to_add = ','.join("{}->{}".format(org, rec) for org, rec in used_map.items())
        fields[7] += ";RECODING=FamilyID_{}:{}".format(family, text_to_add)
    
    return rec_gts, nmissing


def print_line(outhandle, fields, args, subindices):
    """Prints line to output in VCF format."""
    outfields = fields if args.print_all else [field for idx,field in enumerate(fields) if idx < args.startcol or idx in subindices]
    print("\t".join(outfields), file=outhandle)

def print_joinmap(outhandle, rec_markers, args, samples):
    """Prints markers in JoinMap format."""
    settings = "-i_{}_-m_{}".format(args.ind_dp, args.f2_missing)
    f2_ids = [sampleid for fam in samples['F2'].values() for sampleid in fam]
    print("name = {}\npopt = F2\nnloc = {}\nnind = {}\n".format(settings, len(rec_markers), len(f2_ids)), file=outhandle)
    locus = 1
    for marker in rec_markers:
        infostring = ", full locus name: {}".format(marker[0]) if len(marker[0]) > 20 else ""
        print("{} ({}) ; locus{}{}".format(marker[0][-20:], ",".join(['a','h','b']), locus, infostring), file=outhandle)
        print(" ".join(marker[1:]), file=outhandle)
        locus += 1
    print("\nindividual names:", file=outhandle)
    for idx,indid in enumerate(f2_ids):
        print("{} ; {}".format(indid, idx + 1), file=outhandle)



def filter_recode_vcf(samples, args):
    """
    
    Reads in a VCF file, filters each line and outputs a VCF or JoinMap file recoded with a/b allele encoding.
    
    Takes a dictionary of sample ids returned by read_samplelist and a command line arguments object.
    
    Returns nothing.
    
    """
    with open(args.infile, 'r') as inhandle, open(args.outfile, 'w') as outhandle:
        
        # go through vcf file:
        nvar, filtered, pmissing, notinformative, potinformative, notresolved, recoded = 0, 0, 0, 0, 0, 0, 0
        f1_miss, f1_incomp, f2_miss, f2_newal, f2_failseg, fam_not_res = {}, {}, {}, {}, {}, {}
        rec_markers = []
        for line in inhandle:
            if line.startswith("##"):
                if args.outformat == "vcf":
                    print(line, file=outhandle, end="")
                continue
            fields = line.strip().split()
            if line.startswith("#CHROM"):
                sampleids = line.strip().split()[args.startcol:]
                sampleidx = get_sampleidx(samples, sampleids, args.startcol)
                f2_indices = [idx for fam in sampleidx['F2'].values() for idx in fam]
                if args.outformat == "vcf":
                    print_line(outhandle, fields, args, f2_indices)
                continue
            
            if fields[4] == '.':
                continue
            nvar += 1
            if not nvar % 1000:
                print("Processed {} variants.".format(nvar), file=sys.stderr)

            scaff = fields[0]
            pos = fields[1]            
            bases = [fields[3]]
            bases.extend(fields[4].split(','))
            
            # filter sites not satisfying variant type, quality or depth criteria, or not containing indiviudal depth or genotype quality fields:
            fidx = { x[1]: x[0] for x in enumerate(fields[8].split(':')) }
            if not check_site(fields, args, fidx):
                filtered += 1
                continue
                    
            # obtain parental genotypes and check for missing data:
            p_genotypes = Genotypes(sampleidx['P']['NA'], args.ind_dp, args.ind_gq, args.ind_gl)
            p_genotypes.parse_genotypes(fidx, fields)
            if p_genotypes.get_nmissing():
                pmissing += 1
                continue
            
            # check overlap of parental alleles:
            shared_alleles = p_genotypes.get_shared_pair()
            nshared = len(shared_alleles)
            unique_alleles = p_genotypes.get_unique_pair()
            if nshared == 2:
                notinformative += 1
                continue
            
            # check all the families for validity:
            potinformative += 1
            validfam = {fam: True for fam in sampleidx['F2'].keys()}
            f1_genotypes = {}
            f2_genotypes = {}
            rec_gts = ["{}_{}".format(scaff, pos)]
            for key, famindices in sampleidx['F2'].items():
                f1_avail = False
                # if F1s are present, check if their genotypes are compatible with the parental genotpyes:
                if key in sampleidx['F1']:
                    f1_genotypes[key] = Genotypes(sampleidx['F1'][key], args.ind_dp, args.ind_gq, args.ind_gl)
                    f1_genotypes[key].parse_genotypes(fidx, fields)
                    # each F1 genotype needs to share an allele with each of the parents:
                    if not f1_genotypes[key].get_nmissing():
                        f1_avail = True
                        if not f1_genotypes[key].check_compatibility(p_genotypes):
                            validfam[key] = False
                            f1_incomp[key] = f1_incomp.get(key, 0) + 1
                            if args.verbose:
                                print("WARNING: {}:{} - Family {}: Rejected parentally informative variant due to incompatible F1 genotypes.".format(scaff, pos, key), file=sys.stderr)
                                print("Parental genotypes: {}, {}; F1 genotypes: {}, {}.".format(p_genotypes.gts[0], p_genotypes.gts[1], f1_genotypes[key].gts[0], f1_genotypes[key].gts[1]), file=sys.stderr)
                
                # now decide how to handle incomplete or incompatible F1 genotypes:
                if args.require_f1 and not f1_avail:
                    validfam[key] = False
                    f1_miss[key] = f1_miss.get(key, 0) + 1
                    if args.verbose:
                        print("WARNING: {}:{} - Family {}: Rejected parentally informative variant due to missing F1 genotypes.".format(scaff, pos, key), file=sys.stderr)
                
                # check proportion of valid F2 genotypes:
                f2_mindp = args.f2_dp if args.f2_dp else args.ind_dp
                f2_genotypes[key] = Genotypes(famindices, f2_mindp, args.ind_gq, args.ind_gl)
                f2_genotypes[key].parse_genotypes(fidx, fields)
                if f2_genotypes[key].get_missingness() > args.f2_missing:
                    validfam[key] = False
                    f2_miss[key] = f2_miss.get(key, 0) + 1
                    if args.verbose:
                        print("WARNING: {}:{} - Family {}: Rejected parentally informative variant due to high levels of missing F2 genotypes.".format(scaff, pos, key), file=sys.stderr)
                
                # check if all F2 alleles are present in the parents:
                p_alleles = p_genotypes.get_alleles()
                f2_alleles = f2_genotypes[key].get_alleles()
                non_parental = f2_genotypes[key].get_nonparental()
                if (args.excl_new and non_parental) or (args.max_nonparental and non_parental > args.max_nonparental):
                    validfam[key] = False
                    f2_newal[key] = f2_newal.get(key, 0) + 1
                    if args.verbose:
                        print("WARNING: {}:{} - rejected parentally informative variant due to unexpected alleles in the F2 individuals.".format(scaff, pos), file=sys.stderr)
                        print("Parental alleles: {}; F2 alleles: {}.".format(p_alleles, f2_alleles), file=sys.stderr)
                
                # now start with the recoding procedure:
                transl_map = ['.'] * len(bases)
                rcfields = fields.copy()    # make copy of fields list for recoding of genotypes:
                require_seg = False
                
                if validfam[key]:
                    resolved, f2_seg = False, False
                    
                    # add all unique parental alleles to the translation map:
                    for pidx in range(2):
                        for idx in unique_alleles[pidx]:
                            transl_map[idx] = encoding[pidx]
                    
                    # if there is no allele overlap in the parents, recoding is simple:
                    if nshared == 0:
                        resolved = True
                    
                    # if there is only one shared allele, more recoding is possible:
                    elif nshared == 1 and not args.only_strict:
                        if f1_avail:
                            # recoding possible if unique allele from one parent is present in both F1s:
                            f1_alleles = f1_genotypes[key].get_alleles_pair()
                            for pidx in range(2):
                                if bool(unique_alleles[pidx] & f1_alleles[0]) and bool(unique_alleles[pidx] & f1_alleles[1]):
                                    transl_map[next(iter(shared_alleles))] = encoding[1-pidx]
                                    resolved = True
                                    break
                        
                        # if no F1s available, recoding possible if only a unique and a shared allele segregate in the F2 individuals:
                        elif len(f2_alleles & p_alleles) == 2:
                            require_seg = True   # we need to test F2 segregation proportions in any case.
                            for pidx in range(2):
                                if unique_alleles[pidx]:
                                    transl_map[next(iter(shared_alleles))] = encoding[1-pidx]
                                    resolved = True
                                    break
                    
                    # if variant was recoded, check if F2 genotypes segregate in expected proportions:
                    if resolved:
                        if require_seg or args.require_allseg:
                            pvalue = args.pvalue if require_seg else args.pvalue_all
                            nodata = False 
                            try:
                                chi2test, gt_counts = f2_genotypes[key].chisquare_test(transl_map)
                            except Exception as err:
                                nodata = True
                                if args.verbose:
                                    print(err, file=sys.stderr)
                            
                            if nodata or chi2test[1] < pvalue:
                                # reset the translation map to all missing:
                                transl_map = ['.'] * len(bases)
                                validfam[key] = False
                                f2_failseg[key] = f2_failseg.get(key, 0) + 1
                                if args.verbose and not nodata:
                                    print("WARNING: {}:{} - Family {}: Rejected parentally informative variant due to unexpected genotype proportions in the F2 individuals.".format(scaff, pos, key), file=sys.stderr)
                                    print("Genotype counts: {}.".format(",".join(map(str, gt_counts))), file=sys.stderr)
                                    print("Chi-square statistic: {}; p-value: {}.".format(round(chi2test[0], 4), chi2test[1]), file=sys.stderr)
                    
                    else:
                        # reset the translation map to all missing:
                        transl_map = ['.'] * len(bases)
                        validfam[key] = False
                        fam_not_res[key] = fam_not_res.get(key, 0) + 1
                
                # now recode the genotype fields:
                rec_dp = args.rec_dp if args.rec_dp else f2_mindp
                fam_rec_gts, nnotrec = recode_genotypes(rcfields, fidx, famindices, rec_dp, transl_map, args.annot, key)
                if validfam[key] and nnotrec / len(fam_rec_gts) > args.f2_missing + 0.1:
                   print("WARNING: {}:{} - Suspicious proportion of genotypes could not be recoded.".format(scaff, pos), file=sys.stderr)
                   print(line, file=sys.stderr)
                if validfam[key] and nnotrec / len(fam_rec_gts) > args.rec_missing:
                    validfam[key] = False
                    fam_not_res[key] = fam_not_res.get(key, 0) + 1
                    if args.verbose:
                        print("WARNING: {}:{} - Too many genotypes could not be recoded.".format(scaff, pos), file=sys.stderr)
                        print("Translation map: ", transl_map, sep="", file=sys.stderr)
                rec_gts += fam_rec_gts
            
            if args.require_allfam and not all(validfam.values()):
                notresolved += 1
                continue
            elif not any(validfam.values()):
                notresolved += 1
                continue
            else:
                if args.outformat == "vcf":
                    print_line(outhandle, rcfields, args, f2_indices)
                elif args.outformat == "joinmap":
                    rec_markers.append(rec_gts)
                recoded += 1
        
        # print output in JoinMap format:
        if args.outformat == "joinmap":
            print_joinmap(outhandle, rec_markers, args, samples)
        
        print("Processed a total of {} variants, out of which {} were removed; recoded a total of {} variants.".format(nvar, nvar-recoded, recoded), file=sys.stderr)
        print("{} variants were filtered for the following reasons:".format(nvar-recoded), file=sys.stderr)
        print("Quality-filtered sites: {} ({}%)".format(filtered, round(filtered / nvar * 100, 2)), file=sys.stderr)
        print("Parental genotypes missing: {} ({}%)".format(pmissing, round(pmissing / nvar * 100, 2)), file=sys.stderr)
        print("Parental genotypes not informative: {} ({}%)".format(notinformative, round(notinformative / nvar * 100, 2)), file=sys.stderr)
        print("Of {} potentially informative variants, {} ({}%) variants were filtered for the following reasons:".format(potinformative, notresolved, round(notresolved / potinformative * 100, 2)), file=sys.stderr)
        for key in samples['F2']:
            print("Family {}:".format(key), file=sys.stderr)
            print("F1 missing: {} ({}%)".format(f1_miss.get(key, 0), round(f1_miss.get(key, 0) / potinformative * 100, 2)), file=sys.stderr)
            print("F1 incompatible: {} ({}%)".format(f1_incomp.get(key, 0), round(f1_incomp.get(key, 0) / potinformative * 100, 2)), file=sys.stderr)
            print("F2 excess missing: {} ({}%)".format(f2_miss.get(key, 0), round(f2_miss.get(key, 0) / potinformative * 100, 2)), file=sys.stderr)
            print("F2 new alleles: {} ({}%)".format(f2_newal.get(key, 0), round(f2_newal.get(key, 0) / potinformative * 100, 2)), file=sys.stderr)
            print("F2 failed segregation test: {} ({}%)".format(f2_failseg.get(key, 0), round(f2_failseg.get(key, 0) / potinformative * 100, 2)), file=sys.stderr)
            print("Alleles not resolved: {} ({}%)".format(fam_not_res.get(key, 0), round(fam_not_res.get(key, 0) / potinformative * 100, 2)), file=sys.stderr)



def main(args):
    try:
        samples = read_samplelist(args.samplefile)
        filter_recode_vcf(samples, args)
    except Exception as err:
        print(err, file=sys.stderr)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("samplefile", help="List of samples")
    parser.add_argument("infile", help="VCF file")
    parser.add_argument("outfile", help="Output file")
    parser.add_argument("-p", "--pvalue", dest='pvalue', help="P-value threshold for chi-squared test", type=float, default=0.05)
    parser.add_argument("-v", "--min-qual", dest='min_vq', help="Minimum variant quality of site to be considered", type=float, default=0.)
    parser.add_argument("-d", "--min-depth", dest='min_dp', help="Minimum sequencing depth of site to be considered", type=int, default=0)
    parser.add_argument("-i", "--ind-depth", dest='ind_dp', help="Minimum sequencing depth of individuals to count as valid genotypes", type=int, default=5)
    parser.add_argument("--f2-depth", dest='f2_dp', help="Minimum sequencing depth of F2 individuals to count as valid genotypes", type=int)
    parser.add_argument("--rec-depth", dest='rec_dp', help="Minimum sequencing depth for genotypes to be recoded? Set to missing in output otherwise.", type=int)
    parser.add_argument("-q", "--ind-gq", dest='ind_gq', help="Minimum genotype quality of individual to count as valid genotype", type=float, default=0.)
    parser.add_argument("--ind-gl", dest='ind_gl', help="Minimum difference of genotype likelihoods between most and second likely genotype.", type=float, default=0.)
    parser.add_argument("-m", "--f2-missing", dest='f2_missing', help="Maximum proportion of missing genotypes for each F2 family", type=float, default=0.50)
    parser.add_argument("-r", "--rec-missing", dest='rec_missing', help="Maximum proportion of missing genotypes after recoding for each F2 family", type=float, default=0.50)
    parser.add_argument("--startcol", help="0-based start index of first individual in VCF file", type=int, default=9)
    parser.add_argument("--only-strict", dest='only_strict', action='store_true', help="Only record variants without parental allele sharing", default=False)
    parser.add_argument("--only-snps", dest='only_snps', action='store_true', help="Use only SNPs?", default=False)
    parser.add_argument("--include-indels", dest='excl_indels', action='store_false', help="Also include InDels?", default=True)
    parser.add_argument("--exclude-new", dest='excl_new', action='store_true', help="Exclude sites with new alleles in the F2 individuals?", default=False)
    parser.add_argument("--max-nonparental", dest='max_nonparental', help="Maximum proportion of F2 genotypes per family containing alleles not found in the parents", type=float, default=0.01)
    parser.add_argument("--require-f1", dest='require_f1', action='store_true', help="Require valid and compatible F1 genotypes to include variant?", default=False)
    parser.add_argument("--require-allseg", dest='require_allseg', action='store_true', help="Require proper F2 segregation to include variant?", default=False)
    parser.add_argument("--pvalue_all", dest='pvalue_all', help="P-value threshold for chi-squared test for well-resolved variants", type=float, default=0.0001)
    parser.add_argument("--require-allfam", dest='require_allfam', action='store_true', help="Require that all families are valid for printing variant?", default=False)
    parser.add_argument("--annotate", dest='annot', action='store_true', help="Document recoding in annotation field of output VCF file?", default=False)
    parser.add_argument("--output-format", dest='outformat',choices=['vcf', 'joinmap'], help="Desired output format?", default='vcf')
    parser.add_argument("--print-all", dest='print_all', action='store_true', help="Print all individuals instead of just the F2 individuals?", default=False)
    parser.add_argument("--verbose", dest='verbose', action='store_true', help="Print all warning messages to stderr?", default=False)
    main(parser.parse_args())


