#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import pandas as pd
import json
from multiprocessing.pool import Pool
import numpy as np
from time import time
import platform
import gzip

ACMG_DB = 'clinvar dbsnp fathmm fathmm_mkl genocanyon gerp gnomad3 interpro lrt metalr metasvm ' \
          'mutation_assessor mutationtaster omim polyphen2 provean sift siphy spliceai hpo'

parser = argparse.ArgumentParser(description='DIABLO ANNOTATE')

parser.add_argument('-i', '--input', required=True, help='Input as vcf or annotated tsv file', type=str)
parser.add_argument('-o', '--output', required=True, help='Output in tsv format (has to end with ".tsv")', type=str)
parser.add_argument('-s', '--size', required=False, help='Number of lines in output', type=int)
parser.add_argument('-d', '--data', required=True, help='Folder with databases', type=str)
parser.add_argument('-t', '--threads', required=False, help='Number of threads to use', type=int, default=1)
parser.add_argument('-S', '--splice', required=False,
                    help='Create a separate file with predicted splice variants (SpliceAI)', type=bool, default=False)

args = parser.parse_args()

file_input = os.path.abspath(args.input)
file_output = os.path.abspath(args.output)

if args.input is None:
    raise Exception("Input is missing (-i)")
else:
    if platform.system() == 'Windows' and args.input[-3:] != 'tsv':
        mod_n = file_input.split('\\')[-1].split('.')[0]
        mod_d = file_output.strip(file_output.split('\\')[-1])
        os.system(f"oc run {file_input} -l hg38 -t tsv -a {ACMG_DB} -n {mod_n} -d {mod_d}")
    elif platform.system() != 'Windows' and args.input[-3:] != 'tsv':
        os.system(f"oc run {file_input} -l hg38 -t tsv -a {ACMG_DB} -n {file_input.split('/')[-1].split('.')[0]}"
                  f" -d {file_output.strip(file_output.split('/')[-1])}")

start = time()

print('Loading dataframe...')

if args.input[-3:] == 'tsv':
    df = pd.read_csv(str(os.path.abspath(args.input)), sep='\t', comment='#')
else:
    df = pd.read_csv(
        str(file_output.strip(file_output.split('/')[-1]) + file_input.split('/')[-1].split('.')[0] + '.variant.tsv'),
        sep='\t', comment='#')

print('Dataframe loaded')

print('Loading databases...')

PS1_chrst = []
PS1_stop = []
PS1_refs = []
PS1_alts = []
PS1_refa = []
PS1_alta = []
PS1_full = []

for i in open(os.path.join(args.data, 'PS1.txt'), 'r').readlines():
    PS1_chrst.append(str(i.split('\t')[0]) + ':' + str(i.split('\t')[1]))
    PS1_stop.append(i.split('\t')[2])
    PS1_refs.append(i.split('\t')[3])
    PS1_alts.append(i.split('\t')[4])
    PS1_refa.append(i.split('\t')[5])
    PS1_alta.append(i.split('\t')[6])
    PS1_full.append(i.split('\t')[7])

BS2_hom_het_dict = {"het": {'chr1': {},
                            'chr2': {},
                            'chr3': {},
                            'chr4': {},
                            'chr5': {},
                            'chr6': {},
                            'chr7': {},
                            'chr8': {},
                            'chr9': {},
                            'chr10': {},
                            'chr11': {},
                            'chr12': {},
                            'chr13': {},
                            'chr14': {},
                            'chr15': {},
                            'chr16': {},
                            'chr17': {},
                            'chr18': {},
                            'chr19': {},
                            'chr20': {},
                            'chr21': {},
                            'chr22': {},
                            'chrX': {},
                            'chrY': {},
                            'chrM': {},
                            'chrMT': {}},
                    "hom": {'chr1': {},
                            'chr2': {},
                            'chr3': {},
                            'chr4': {},
                            'chr5': {},
                            'chr6': {},
                            'chr7': {},
                            'chr8': {},
                            'chr9': {},
                            'chr10': {},
                            'chr11': {},
                            'chr12': {},
                            'chr13': {},
                            'chr14': {},
                            'chr15': {},
                            'chr16': {},
                            'chr17': {},
                            'chr18': {},
                            'chr19': {},
                            'chr20': {},
                            'chr21': {},
                            'chr22': {},
                            'chrX': {},
                            'chrY': {},
                            'chrM': {},
                            'chrMT': {}}
                    }

for line in gzip.open(os.path.join(args.data, 'BS2_mod.txt.gz'), 'rt').readlines():
    chr = line.split(' ')[0].strip()
    pos = line.split(' ')[1].strip()
    ref = line.split(' ')[2].strip()
    alt = line.split(' ')[3].strip()
    hom = line.split(' ')[4].strip()
    het = line.split(' ')[5].strip()
    if str(het) == '1':
        BS2_hom_het_dict["het"][chr][pos] = ref + '_' + alt
    if str(hom) == '1':
        BS2_hom_het_dict["hom"][chr][pos] = ref + '_' + alt

bs2_hom_het_ad_df = pd.read_csv(os.path.join(args.data, 'BS2_rec_dom_ad.txt.gz'), compression='gzip', sep='\t')

rec_list = set()
dom_list = set()
adult_list = set()
for gene, rec, dom, adult in zip(bs2_hom_het_ad_df['gene'], bs2_hom_het_ad_df['recessive'],
                                 bs2_hom_het_ad_df['dominant'], bs2_hom_het_ad_df['adult_onset']):
    if rec == 1:
        rec_list.add(gene)
    if dom == 1:
        dom_list.add(gene)
    if adult == 1:
        adult_list.add(gene)

BP1_list = []
for i in open(os.path.join(args.data, 'BP1.txt'), 'r').readlines():
    BP1_list.append(str(i.strip('\n')))

PP2_list = []
for i in open(os.path.join(args.data, 'PP2.txt'), 'r').readlines():
    PP2_list.append(i.strip('\n'))

pli = open(os.path.join(args.data, 'pli_dict.json'), 'r')
pli_dict = json.load(pli)

PM1_str = ''

for i in open(os.path.join(args.data, 'PM1.txt'), 'r').readlines():
    PM1_str += i.strip('\n')

repeat_dict = gzip.open(os.path.join(args.data, 'repeat_dict.hg38.gz'), 'rt')
repeat_reg = json.load(repeat_dict)

print('Databases loaded')

print('Starting ACMG assignment')


def PVS1(subdf):
    PVS1_crit = []
    for so, hpo in zip(subdf['so'].to_numpy(), subdf['hpo.id'].to_numpy()):
        if not pd.isna(so) and not pd.isna(hpo) and (
                str(so).find('frameshift') > -1 or str(so) == 'stop_gained' or
                str(so) == 'splice_site_variant'):
            PVS1_crit.append(1)
        else:
            PVS1_crit.append(0)

    subdf['PVS1'] = PVS1_crit

    return print('PVS1 assigned')


def PS1(subdf):
    PS1_crit = []

    for so, chrom, pos, achange in zip(subdf['so'].to_numpy(), subdf['chrom'].to_numpy(), subdf['pos'].to_numpy(),
                                       subdf['achange'].to_numpy()):
        if str(so) == 'missense_variant':

            if str(chrom) + ':' + str(pos) in PS1_chrst:

                if str(achange.split('p.')[1][0:3]) == str(PS1_refa[PS1_chrst.index(
                        str(str(chrom) + ':' + str(pos)))]) and str(
                    achange.split('p.')[1][-3:]) == \
                        str(PS1_alta[
                                PS1_chrst.index(str(str(chrom) + ':' + str(pos)))]):
                    PS1_crit.append(1)
                else:
                    PS1_crit.append(0)
            else:
                PS1_crit.append(0)
        else:
            PS1_crit.append(0)

    subdf['PS1'] = PS1_crit

    return print('PS1 assigned')


def PS3(subdf):
    PS3_crit = []

    for sig, rev_stat in zip(subdf['clinvar.sig'].to_numpy(), subdf['clinvar.rev_stat'].to_numpy()):
        if (not pd.isna(sig)) and (
                str(sig).find('Pathogenic') > -1) and (
                str(rev_stat).find('reviewed by expert panel') > -1 or
                str(rev_stat).find('practice guideline') > -1 or
                str(rev_stat).find('no conflicts') > -1):
            PS3_crit.append(1)
        else:
            PS3_crit.append(0)

    subdf['PS3'] = PS3_crit

    return print('PS3 assigned')


def PM1(subdf):
    PM1_crit = []
    for so, dom in zip(subdf['so'].to_numpy(), subdf['interpro.domain'].to_numpy()):
        counter = 0
        if not pd.isna(dom):
            for domain in str(dom)[2:-2].replace("'", "").split('|'):
                if domain in PM1_str:
                    counter += 1
        if not pd.isna(dom) and so == 'missense_variant' and counter == 0:
            PM1_crit.append(1)
        else:
            PM1_crit.append(0)

    subdf['PM1'] = PM1_crit

    return print('PM1 assigned')


def PM2(subdf):
    PM2_crit = []

    for hugo, gnom_af in zip(subdf['hugo'].to_numpy(), subdf['gnomad3.af'].to_numpy()):
        if not pd.isna(hugo):
            if not pd.isna(gnom_af):
                try:
                    if (float(gnom_af) <= 0.0001 and float(
                            pli_dict[hugo]['pli']) >= 0.85) or (
                            (float(gnom_af) <= 0.05) and (
                            float(pli_dict[hugo]['pli']) <= 0.85 or
                            pli_dict[hugo]['pli'] == 'nan')):
                        PM2_crit.append(1)
                    else:
                        PM2_crit.append(0)
                except KeyError:
                    PM2_crit.append(0)

            elif pd.isna(gnom_af):
                PM2_crit.append(1)

        else:
            PM2_crit.append(0)
    subdf['PM2'] = PM2_crit

    return print('PM2 assigned')


def PM4_BP3(subdf):
    PM4_crit = []
    BP3_crit = []

    for so, chrom, pos, domain in zip(subdf['so'].to_numpy(), subdf['chrom'].to_numpy(), subdf['pos'].to_numpy(),
                                      subdf['interpro.domain'].to_numpy()):

        if str(so).find('inframe') > -1:
            list_of_all_first_elements = [item[0] for item in repeat_reg[str(chrom)]]
            list_of_all_second_elements = [item[1] for item in repeat_reg[str(chrom)]]
            closest1 = min(list_of_all_first_elements, key=lambda x: abs(x - int(pos)))
            closest2 = min(list_of_all_second_elements, key=lambda x: abs(x - int(pos)))
            index_of_first = list_of_all_first_elements.index(closest1)
            index_of_second = list_of_all_second_elements.index(closest2)
            range_1 = range(repeat_reg[str(chrom)][index_of_first][0],
                            repeat_reg[str(chrom)][index_of_first][1])
            range_2 = range(repeat_reg[str(chrom)][index_of_second][0],
                            repeat_reg[str(chrom)][index_of_second][1])

            if int(pos) in range_1 or int(pos) in range_2:
                PM4_crit.append(0)
                if pd.isna(domain):
                    BP3_crit.append(1)
                else:
                    BP3_crit.append(0)
            else:
                PM4_crit.append(1)
                BP3_crit.append(0)

        elif str(so).find('stop_lost') > -1:
            list_of_all_first_elements = [item[0] for item in repeat_reg[str(chrom)]]
            list_of_all_second_elements = [item[1] for item in repeat_reg[str(chrom)]]
            closest1 = min(list_of_all_first_elements, key=lambda x: abs(x - int(pos)))
            closest2 = min(list_of_all_second_elements, key=lambda x: abs(x - int(pos)))
            index_of_first = list_of_all_first_elements.index(closest1)
            index_of_second = list_of_all_second_elements.index(closest2)
            range_1 = range(repeat_reg[str(chrom)][index_of_first][0],
                            repeat_reg[str(chrom)][index_of_first][1])
            range_2 = range(repeat_reg[str(chrom)][index_of_second][0],
                            repeat_reg[str(chrom)][index_of_second][1])

            if int(pos) in range_1 or int(pos) in range_2:
                PM4_crit.append(0)
                BP3_crit.append(0)
            else:
                PM4_crit.append(1)
                BP3_crit.append(0)
        else:
            PM4_crit.append(0)
            BP3_crit.append(0)

    subdf['PM4'] = PM4_crit
    subdf['BP3'] = BP3_crit

    return print('PM4 and BP3 assigned')


def PM5(subdf):
    PM5_crit = []

    for so, chrom, pos, achange in zip(subdf['so'].to_numpy(), subdf['chrom'].to_numpy(), subdf['pos'].to_numpy(),
                                       subdf['achange'].to_numpy()):
        if so == 'missense_variant':
            if str(str(chrom) + ':' + str(pos)) in PS1_chrst:
                index_PS1 = PS1_chrst.index(str(str(chrom) + ':' + str(pos)))

                if str(achange.split('p.')[1][0:3]) == str(PS1_refa[index_PS1]) and str(
                        achange.split('p.')[1][-3:]) != str(PS1_alta[index_PS1]):
                    PM5_crit.append(1)
                else:
                    PM5_crit.append(0)
            else:
                PM5_crit.append(0)
        else:
            PM5_crit.append(0)

    subdf['PM5'] = PM5_crit

    return print('PM5 assigned')


def PP2(subdf):
    PP2_crit = []

    for so, hugo in zip(subdf['so'].to_numpy(), subdf['hugo'].to_numpy()):
        if str(so) == 'missense_variant':
            if str(hugo) in PP2_list:
                PP2_crit.append(1)
            else:
                PP2_crit.append(0)
        else:
            PP2_crit.append(0)

    subdf['PP2'] = PP2_crit

    return print('PP2 assigned')


def PP3_BP4(subdf):
    PP3_crit = []
    BP4_crit = []

    for sift_score, lrt_score, \
        taste, assess, fatmm, pro, \
        mvm, lr, mkl, geno, grp in zip(subdf['sift.score'].to_numpy(),
                                       subdf['lrt.lrt_pred'].to_numpy(),
                                       subdf['mutationtaster.prediction'].to_numpy(),
                                       subdf['mutation_assessor.impact'].to_numpy(),
                                       subdf['fathmm.fathmm_pred'].to_numpy(),
                                       subdf['provean.score'].to_numpy(),
                                       subdf['metasvm.pred'].to_numpy(),
                                       subdf['metalr.pred'].to_numpy(),
                                       subdf['fathmm_mkl.fathmm_mkl_coding_pred'].to_numpy(),
                                       subdf['genocanyon.score'].to_numpy(),
                                       subdf['gerp.gerp_rs'].to_numpy()):
        if not pd.isna(sift_score):
            if float(sift_score) <= 0.05:
                sift = 1
            elif float(sift_score) > 0.05:
                sift = -1
            else:
                sift = 0
        else:
            sift = 0

        if not pd.isna(lrt_score):
            if str(lrt_score)[0] == 'D':
                lrt = 1
            elif str(lrt_score)[0] == 'N':
                lrt = -1
            else:
                lrt = 0
        else:
            lrt = 0

        if not pd.isna(taste):
            if str(taste).find('Disease') > -1 or str(
                    taste).find('Damaging') > -1:
                muttaste = 1

            elif str(taste).find('Polymorphism') > -1:
                muttaste = -1

            else:
                muttaste = 0
        else:
            muttaste = 0

        if not pd.isna(assess):
            if str(assess)[0] == 'h' or \
                    str(assess)[0] == 'm':
                mutassessor = 1
            elif str(assess)[0] == 'l' or \
                    str(assess)[0] == 'n':
                mutassessor = -1
            else:
                mutassessor = 0
        else:
            mutassessor = 0

        if not pd.isna(fatmm):
            if str(fatmm)[0] == 'D':
                fathmm = 1
            elif str(fatmm)[0] == 'T':
                fathmm = -1
            else:
                fathmm = 0
        else:
            fathmm = 0

        if not pd.isna(pro):
            if float(pro) <= -2.5:
                provean = 1
            elif float(pro) > -2.5:
                provean = -1
            else:
                provean = 0
        else:
            provean = 0

        if not pd.isna(mvm):
            if str(mvm) == 'Damaging':
                metasvm = 1
            elif str(mvm) == 'Tolerated':
                metasvm = -1
            else:
                metasvm = 0
        else:
            metasvm = 0

        if not pd.isna(lr):
            if str(lr) == 'Damaging':
                metalr = 1
            elif str(lr) == 'Tolerated':
                metalr = -1
            else:
                metalr = 0

        else:
            metalr = 0

        if not pd.isna(mkl):
            if str(mkl) == 'Damaging':
                fathmm_mkl = 1
            elif str(mkl) == 'Neutral':
                fathmm_mkl = -1
            else:
                fathmm_mkl = 0
        else:
            fathmm_mkl = 0

        if not pd.isna(geno):
            if float(geno) > 0.5:
                genocanyon = 1
            elif float(geno) <= 0.5:
                genocanyon = -1
            else:
                genocanyon = 0
        else:
            genocanyon = 0

        if not pd.isna(grp):
            if float(grp) >= 2.25:
                gerp = 1
            elif float(grp) < 2.25:
                gerp = -1
            else:
                gerp = 0
        else:
            gerp = 0

        total_pred = (
                gerp + genocanyon + fathmm_mkl + metalr + metasvm +
                provean + fathmm + mutassessor + muttaste + lrt + sift
        )

        if total_pred > 2:
            PP3_crit.append(1)
            BP4_crit.append(0)
        elif 0 <= total_pred <= 2:
            PP3_crit.append(0)
            BP4_crit.append(0)
        elif total_pred < 0:
            PP3_crit.append(0)
            BP4_crit.append(1)
        else:
            PP3_crit.append(0)
            BP4_crit.append(0)

    subdf['PP3'] = PP3_crit
    subdf['BP4'] = BP4_crit

    return print('PP3 and BP4 assigned')


def PP5(subdf):
    PP5_crit = []

    for sig, rev in zip(subdf['clinvar.sig'].to_numpy(), subdf['clinvar.rev_stat'].to_numpy()):
        if not pd.isna(sig):
            if (str(sig).find('Pathogenic') > -1 or str(sig).find(
                    'pathogenic') > -1) and (
                    str(sig).find('conflict') < 0 or
                    str(sig).find('Conflict') < 0) and \
                    (str(rev).find('expert') > -1 or
                     str(rev).find('practical') > -1 or
                     str(rev).find('no conflicts') > -1):
                PP5_crit.append(2)
            elif (str(sig).find('Pathogenic') > -1 or
                  str(sig).find('pathogenic') > -1) and \
                    (str(sig).find('conflict') < 0 or
                     str(sig).find('Conflict') < 0):
                PP5_crit.append(1)
            else:
                PP5_crit.append(0)
        else:
            PP5_crit.append(0)

    subdf['PP5'] = PP5_crit

    return print('PP5 assigned')


def BA1(subdf):
    BA1_crit = []

    for af in subdf['gnomad3.af'].to_numpy():
        if not pd.isna(af):
            if float(af) > 0.05:
                BA1_crit.append(1)
            else:
                BA1_crit.append(0)
        else:
            BA1_crit.append(0)

    subdf['BA1'] = BA1_crit

    return print('BA1 assigned')


def BS2(subdf):
    BS2_crit = []

    for hugo, ref, alt, pos, chrom in zip(subdf['hugo'].to_numpy(),
                                          subdf['ref_base'].to_numpy(),
                                          subdf['alt_base'].to_numpy(),
                                          subdf['pos'].to_numpy(),
                                          subdf['chrom'].to_numpy()):
        if not pd.isna(hugo):
            if str(hugo) in adult_list:
                BS2_crit.append(0)
            elif str(hugo) in rec_list:
                key_snv = str(ref) + '_' + str(alt)
                try:
                    if key_snv == BS2_hom_het_dict["hom"][str(chrom)][str(pos)]:
                        BS2_crit.append(1)
                    else:
                        BS2_crit.append(0)
                except KeyError:
                    BS2_crit.append(0)

            elif str(hugo) in dom_list:
                key_snv = str(ref) + '_' + str(alt)
                try:
                    if key_snv == BS2_hom_het_dict["het"][str(chrom)][str(pos)]:
                        BS2_crit.append(1)
                    else:
                        BS2_crit.append(0)
                except KeyError:
                    BS2_crit.append(0)
            else:
                BS2_crit.append(0)
        else:
            BS2_crit.append(0)

    subdf['BS2'] = BS2_crit

    return print('BS2 assigned')


def BS3(subdf):
    BS3_crit = []

    for sig, rev in zip(subdf['clinvar.sig'].to_numpy(), subdf['clinvar.rev_stat'].to_numpy()):
        if not pd.isna(sig):
            if (str(sig).find('Benign') > -1 or
                str(sig).find('benign') > -1) and \
                    (str(rev).find('expert') > -1 or
                     str(rev).find('practical') > -1):
                BS3_crit.append(1)
            else:
                BS3_crit.append(0)
        else:
            BS3_crit.append(0)

    subdf['BS3'] = BS3_crit

    return print('BS3 assigned')


def BP1(subdf):
    BP1_crit = []

    for so, hugo in zip(subdf['so'].to_numpy(), subdf['hugo'].to_numpy()):
        if not pd.isna(so) and not pd.isna(hugo):
            if str(so) == 'missense_variant' and str(hugo) in BP1_list:
                BP1_crit.append(1)
            else:
                BP1_crit.append(0)
        else:
            BP1_crit.append(0)

    subdf['BP1'] = BP1_crit

    return print('BP1 assigned')


def BP6(subdf):
    BP6_crit = []

    for sig in subdf['clinvar.sig'].to_numpy():
        if not pd.isna(sig):
            if str(sig).find('Benign') > -1 or str(sig).find(
                    'benign') > -1:
                BP6_crit.append(1)
            else:
                BP6_crit.append(0)
        else:
            BP6_crit.append(0)

    subdf['BP6'] = BP6_crit

    return print('BP6 assigned')


def BP7(subdf):
    BP7_crit = []

    for so, ag, al, dg, dl in zip(subdf['so'].to_numpy(),
                                  subdf['spliceai.ds_ag'].to_numpy(),
                                  subdf['spliceai.ds_al'].to_numpy(),
                                  subdf['spliceai.ds_dg'].to_numpy(),
                                  subdf['spliceai.ds_dl'].to_numpy()):
        if not pd.isna(so):
            if str(so) == 'synonymous_variant':
                if (not pd.isna(ag) or
                    not pd.isna(al) or
                    not pd.isna(dg) or
                    not pd.isna(dl)) and (float(ag) >= 0.25 or
                                          float(al) >= 0.25 or
                                          float(dg) >= 0.25 or
                                          float(dl) >= 0.25):
                    BP7_crit.append(0)
                else:
                    BP7_crit.append(1)
            else:
                BP7_crit.append(0)
        else:
            BP7_crit.append(0)

    subdf['BP7'] = BP7_crit

    return print('BP7 assigned')


def BS1(subdf):
    BS1_crit = []

    for af in subdf['gnomad3.af'].to_numpy():
        if not pd.isna(af):
            if float(af) > 0.05:
                BS1_crit.append(1)
            else:
                BS1_crit.append(0)
        else:
            BS1_crit.append(0)

    subdf['BS1'] = BS1_crit

    return print('BS1 assigned')


def remove():
    sqlite_file = \
        f"{os.path.join(file_output.strip(file_output.split('/')[-1]), file_input.split('/')[-1].split('.')[0])}.sqlite"
    err_file = \
        f"{os.path.join(file_output.strip(file_output.split('/')[-1]), file_input.split('/')[-1].split('.')[0])}.err"
    log_file = \
        f"{os.path.join(file_output.strip(file_output.split('/')[-1]), file_input.split('/')[-1].split('.')[0])}.log"

    if os.path.isfile(sqlite_file):
        os.remove(sqlite_file)
    if os.path.isfile(err_file):
        os.remove(err_file)
    if os.path.isfile(log_file):
        os.remove(log_file)

    return print('Removing temporary files...')


def pathogenicity_assignment(subdf):
    probabilities = []
    predictions = []

    for pvs1, ps1, ps3, pm1, pm2, pm4, pm5, pp2, pp3, pp5, bs2, bs3, bs1, bp3, bp4, bp1, bp6, bp7, ba1 in \
            zip(subdf['PVS1'].to_numpy(),
                subdf['PS1'].to_numpy(),
                subdf['PS3'].to_numpy(),
                subdf['PM1'].to_numpy(),
                subdf['PM2'].to_numpy(),
                subdf['PM4'].to_numpy(),
                subdf['PM5'].to_numpy(),
                subdf['PP2'].to_numpy(),
                subdf['PP3'].to_numpy(),
                subdf['PP5'].to_numpy(),
                subdf['BS2'].to_numpy(),
                subdf['BS3'].to_numpy(),
                subdf['BS1'].to_numpy(),
                subdf['BP3'].to_numpy(),
                subdf['BP4'].to_numpy(),
                subdf['BP1'].to_numpy(),
                subdf['BP6'].to_numpy(),
                subdf['BP7'].to_numpy(),
                subdf['BA1'].to_numpy()):
        PVS = int(pvs1)
        PS = int(ps1) + int(ps3)
        PM = int(pm1) + int(pm2) + int(pm4) + int(pm5)
        PP = int(pp2) + int(pp3) + int(pp5)
        BS = int(bs2) + int(bs3) + int(bs1)
        BP = int(bp3) + int(bp4) + int(bp1) + int(bp6) + int(bp7)
        BA1 = int(ba1)
        pc = 0.10
        X = 2

        if BA1 != 1:
            odds_path = 350 ** ((PP / (X ** 3)) + (PM / X ** 2) + (PS / X) + (PVS / 1) - (BP / X ** 3) - (BS / X))
            proba = (odds_path * pc) / (((odds_path - 1) * pc) + 1)
            probability = float("%.4f" % proba)

        else:
            probability = 0

        probabilities.append(probability)

        if BA1 == 1:
            prediction = 'Benign auto'
        elif BS > 2:
            prediction = 'Benign'
        elif (BS == 1 and BP >= 1) or (BP >= 2):
            prediction = 'Likely Benign'
        elif (PVS == 1 and (PS >= 1 or PM >= 2 or (PM == 1 and PP == 1) or PP >= 2)) or (PS >= 2) or (
                PS == 1 and (PM >= 3 or (PM == 2 and PP >= 2) or (PM == 1 and PP >= 4))):
            prediction = 'Pathogenic'
        elif (PVS == 1 and PM == 1) or (PS == 1 and PM >= 1) or (PS == 1 and PP >= 2) or (PM >= 3) or (
                PM >= 2 and PP >= 2) or (PM == 1 and PP >= 4):
            prediction = 'Likely Pathogenic'
        else:
            prediction = 'VUS'

        predictions.append(prediction)

    subdf['Score'] = probabilities
    subdf['ACMG'] = predictions

    return print('Pathogenicity assigned')


def splice(subdf):
    splice_df = subdf.query("`spliceai.ds_al` > 0.7 or `spliceai.ds_dg` > 0.7 or "
                            "`spliceai.ds_dl` > 0.7 or `spliceai.ds_ag` > 0.7")
    return splice_df


def main(subdf):
    PVS1(subdf)
    PS1(subdf)
    PS3(subdf)
    PM1(subdf)
    PM2(subdf)
    PM4_BP3(subdf)
    PM5(subdf)
    PP2(subdf)
    PP3_BP4(subdf)
    PP5(subdf)
    BA1(subdf)
    BS2(subdf)
    BS3(subdf)
    BP1(subdf)
    BP6(subdf)
    BP7(subdf)
    BS1(subdf)
    pathogenicity_assignment(subdf)

    return subdf


if args.threads >= 2:
    df = np.array_split(df, args.threads)
    m_list = []
    for part in df:
        tupl = pd.DataFrame(part)
        m_list.append(tupl.reset_index())

    pool = Pool(processes=args.threads)
    result_proc = pool.map(main, m_list)
    final_pool = pd.concat(result_proc, ignore_index=True, sort=False)

    pool.close()
    pool.join()

else:
    df = main(df)

print('Sorting...')

df['ACMG'] = pd.Categorical(df['ACMG'], ['Pathogenic', 'Likely Pathogenic', 'VUS',
                                         'Likely Benign', 'Benign', 'Benign auto'])

df = df.sort_values(by=['ACMG', 'Score'], ascending=[True, False])

if args.size is None:
    df.to_csv(os.path.abspath(args.output), sep='\t', index=False)
else:
    df.head(args.size).to_csv(os.path.abspath(args.output), sep='\t', index=False)

remove()

if args.splice is not False:
    print('Writing file with splice variants...')
    splice(df).to_csv(os.path.abspath(args.output).replace('.tsv', '_splice.tsv'), sep='\t', index=False)
else:
    pass

print('Success!')

stop = time()

print(f'Finished in {str(np.ceil(stop - start))} seconds')
