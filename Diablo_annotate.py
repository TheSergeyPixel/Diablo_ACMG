#!/usr/bin/env python
# coding: utf-8



import os
# from openpyxl import load_workbook
import argparse
import pandas as pd
import json
from multiprocessing.pool import Pool
from numpy import array_split, ceil
from time import time


ACMG_DB = 'clinvar dbsnp fathmm fathmm_mkl genocanyon gerp gnomad gnomad3 interpro lrt metalr metasvm mutation_assessor mutationtaster omim polyphen2 provean sift siphy spliceai hpo' 

parser = argparse.ArgumentParser(description='DIABLO ANNOTATE')

parser.add_argument('-i', '--input', required=True, help='Input a vcf file', type=str )
parser.add_argument('-o', '--output', required=True, help='Output in tsv format', type=str )
parser.add_argument('-s', '--size', required=False, help='Number of lines in output', type=int )
parser.add_argument('-d', '--data', required=True, help='Folder with databases', type=str )
parser.add_argument('-t', '--threads', required=False, help='Number of threads to use', type=int, default=1 )

args = parser.parse_args()


file_input = os.path.abspath(args.input)
file_output = os.path.abspath(args.output)

if args.input != None:
    os.system(f"oc run {file_input} -l hg38 -t tsv -a {ACMG_DB} -n {file_input.split('/')[-1].split('.')[0]} -d {file_output.strip(file_output.split('/')[-1])}")

start = time()

print('Loading dataframe...')

df = pd.read_csv(str(file_output.strip(file_output.split('/')[-1]) +  file_input.split('/')[-1].split('.')[0] + '.variant.tsv') , sep='\t', comment='#', low_memory=False)

print('Dataframe loaded')



#def process(subdf, smth):
#def load_databases():

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


for line in open(os.path.join(args.data, 'BS2_mod.txt'), 'r').readlines():
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

        
bs2_hom_het_ad_df = pd.read_csv(os.path.join(args.data, 'BS2_rec_dom_ad.txt'), sep='\t')

rec_list = []
dom_list = []
adult_list = []
for gene, rec, dom, adult in zip(bs2_hom_het_ad_df['gene'], bs2_hom_het_ad_df['recessive'],
                                     bs2_hom_het_ad_df['dominant'], bs2_hom_het_ad_df['adult_onset']):
    if rec == 1:
        rec_list.append(gene)
    if dom == 1:
        dom_list.append(gene)
    if adult == 1:
        adult_list.append(gene)    

        
BP1_list = []
for i in open(os.path.join(args.data, 'BP1.txt'), 'r').readlines():
    BP1_list.append(str(i.strip('\n')))
        
PP2_list = []
for i in open(os.path.join(args.data, 'PP2.txt'), 'r').readlines():
    PP2_list.append(i.replace('\n', ''))

    
pli = open(os.path.join(args.data, 'pli_dict.json'), 'r')
pli_dict = json.load(pli)
    
PM1_str = ''

for i in open(os.path.join(args.data, 'PM1.txt'), 'r').readlines():
    PM1_str += i.strip('\n')

repeat_dict = open(os.path.join(args.data, 'repeat_dict.hg38'), 'r')
repeat_reg = json.load(repeat_dict)

print('Databases loaded')

def PVS1(subdf):

    PVS1_crit = []
    for i in range(len(subdf)) :
        if pd.isna(subdf.loc[i, 'so']) == False and (str(subdf.loc[i, 'so']).find('frameshift') > -1 or str(subdf.loc[i, 'so']) == 'stop_gained' or str(subdf.loc[i, 'so']) == 'splice_site_variant'):
            PVS1_crit.append(1)
        elif (pd.isna(subdf.loc[i, 'spliceai.ds_ag']) == False or pd.isna(subdf.loc[i, 'spliceai.ds_al']) == False or pd.isna(subdf.loc[i, 'spliceai.ds_dg']) == False or pd.isna(subdf.loc[i, 'spliceai.ds_dl']) == False) and (float(subdf.loc[i, 'spliceai.ds_ag']) >= 0.7 or float(subdf.loc[i, 'spliceai.ds_al']) >= 0.7 or float(subdf.loc[i, 'spliceai.ds_dg']) >= 0.7 or float(subdf.loc[i, 'spliceai.ds_dl']) >= 0.7):
            PVS1_crit.append(1)
        else:
            PVS1_crit.append(0)

    subdf['PVS1'] = PVS1_crit

    return print('PVS1 assigned')

def PS1(subdf):

    PS1_crit = []

    for i in range(len(subdf)):
        if str(subdf.loc[i, 'so']) == 'missense_variant':

            if str(str(subdf.loc[i, 'chrom']) + ':' + str(subdf.loc[i, 'pos'])) in PS1_chrst:

                if str(subdf.loc[i, 'achange'].split('p.')[1][0:3]) == str(PS1_refa[PS1_chrst.index(str(str(subdf.loc[i, 'chrom']) + ':' + str(subdf.loc[i, 'pos'])))]) and  str(subdf.loc[i, 'achange'].split('p.')[1][-3:]) == str(PS1_alta[PS1_chrst.index(str(str(subdf.loc[i, 'chrom']) + ':' + str(subdf.loc[i, 'pos'])))]): 
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

    for i in range(len(subdf)):
        if (pd.isna(subdf.loc[i, 'clinvar.sig']) == False) and (str(subdf.loc[i, 'clinvar.sig']).find('Pathogenic') > -1) and (str(subdf.loc[i, 'clinvar.rev_stat']).find('reviewed by expert panel') > -1 or str(subdf.loc[i, 'clinvar.sig']).find('practice guideline') > -1):
            PS3_crit.append(1) 
        else:
            PS3_crit.append(0) 

    subdf['PS3'] = PS3_crit
        
    return print('PS3 assigned')

def PM1(subdf):

    PM1_crit = []
    for i in range(len(subdf)):    
        counter = 0
        if pd.isna(subdf.loc[i, 'interpro.domain']) == False:
            for domain in str(subdf.loc[i, 'interpro.domain'])[2:-2].replace("'", "").split('|'):
                if domain in PM1_str:
                    counter += 1
        if pd.isna(subdf.loc[i, 'interpro.domain']) == False and subdf.loc[i, 'so'] == 'missense_variant' and counter == 0:
            PM1_crit.append(1)
        else:
            PM1_crit.append(0)
            
    subdf['PM1'] = PM1_crit       

    return print('PM1 assigned')

def PM2(subdf):

    PM2_crit = []

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'hugo']) == False:
            if pd.isna(subdf.loc[i, 'gnomad3.af']) == False: 
                try:
                    if (float(subdf.loc[i, 'gnomad3.af']) <= 0.0001 and float(pli_dict[subdf.loc[i, 'hugo']]['pli']) >= 0.85) or ((float(subdf.loc[i, 'gnomad3.af']) <= 0.05) and (float(pli_dict[subdf.loc[i, 'hugo']]['pli']) <= 0.85 or pli_dict[subdf.loc[i, 'hugo']]['pli'] == 'nan')):
                        PM2_crit.append(1)
                    else:
                        PM2_crit.append(0)
                except KeyError:
                    PM2_crit.append(0)

            elif pd.isna(subdf.loc[i, 'gnomad3.af']) == True:
                PM2_crit.append(1)
            
        else:
            PM2_crit.append(0)
    subdf['PM2'] = PM2_crit

    return print('PM2 assigned')

def PM4_BP3(subdf):

    PM4_crit = []
    BP3_crit = []

    for so, chrom, pos, domain in zip(subdf['so'], subdf['chrom'], subdf['pos'], subdf['interpro.domain']):

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
                if pd.isna(domain) == True:
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


# def PM4_BP3(subdf):

#     PM4_crit = []
#     BP3_crit = []

#     for i in range(len(subdf)):
#         if str(subdf.loc[i, 'so']).find('inframe') > -1:
#             list_of_all_first_elements = [item[0] for item in repeat_reg[str(subdf.loc[i, 'chrom'])]]  
#             list_of_all_second_elements = [item[1] for item in repeat_reg[str(subdf.loc[i, 'chrom'])]]
#             closest1 = min(list_of_all_first_elements, key=lambda x: abs(x - int(subdf.loc[i, 'pos'])))
#             closest2 = min(list_of_all_second_elements, key=lambda x: abs(x - int(subdf.loc[i, 'pos'])))
#             index_of_first = list_of_all_first_elements.index(closest1)
#             index_of_second = list_of_all_second_elements.index(closest2)
#             range_1 = range(repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_first][0],
#                                 repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_first][1])
#             range_2 = range(repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_second][0],
#                                 repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_second][1])
            
#             if int(subdf.loc[i, 'pos']) in range_1 or int(subdf.loc[i, 'pos']) in range_2:
#                 PM4_crit.append(0)
#                 if pd.isna(subdf.loc[i, 'interpro.domain']) == True:
#                     BP3_crit.append(1)
#                 else:
#                     BP3_crit.append(0)
#             else:
#                 PM4_crit.append(1)
#                 BP3_crit.append(0)
#         elif str(subdf.loc[i, 'so']).find('stop_lost') > -1:
#             list_of_all_first_elements = [item[0] for item in repeat_reg[str(subdf.loc[i, 'chrom'])]]  
#             list_of_all_second_elements = [item[1] for item in repeat_reg[str(subdf.loc[i, 'chrom'])]]
#             closest1 = min(list_of_all_first_elements, key=lambda x: abs(x - int(subdf.loc[i, 'pos'])))
#             closest2 = min(list_of_all_second_elements, key=lambda x: abs(x - int(subdf.loc[i, 'pos'])))
#             index_of_first = list_of_all_first_elements.index(closest1)
#             index_of_second = list_of_all_second_elements.index(closest2)
#             range_1 = range(repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_first][0],
#                                 repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_first][1])
#             range_2 = range(repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_second][0],
#                                 repeat_reg[str(subdf.loc[i, 'chrom'])][index_of_second][1])
            
#             if int(subdf.loc[i, 'pos']) in range_1 or int(subdf.loc[i, 'pos']) in range_2:
#                 PM4_crit.append(0) 
#                 BP3_crit.append(0)
#             else:
#                 PM4_crit.append(1) 
#                 BP3_crit.append(0)
#         else:
#             PM4_crit.append(0) 
#             BP3_crit.append(0)

#     subdf['PM4'] = PM4_crit
#     subdf['BP3'] = BP3_crit

#     return print('PM4 and BP3 assigned')

    
def PM5(subdf):

    PM5_crit = []

    for i in range(len(subdf)):
        if subdf.loc[i, 'so'] == 'missense_variant':
            if  str(str(subdf.loc[i, 'chrom']) + ':' + str(subdf.loc[i, 'pos'])) in PS1_chrst:
                index_PS1 = PS1_chrst.index(str(str(subdf.loc[i, 'chrom']) + ':' + str(subdf.loc[i, 'pos'])))

                if str(subdf.loc[i, 'achange'].split('p.')[1][0:3]) == str(PS1_refa[index_PS1]) and  str(subdf.loc[i, 'achange'].split('p.')[1][-3:]) != str(PS1_alta[index_PS1]):
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

    for i in range(len(subdf)):
        if str(subdf.loc[i, 'so']) == 'missense_variant':
            if str(subdf.loc[i, 'hugo']) in PP2_list:
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

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'sift.score']) == False:
            if float(subdf.loc[i, 'sift.score']) <= 0.05:
                sift = 1
            elif float(subdf.loc[i, 'sift.score']) > 0.05:
                sift = -1
            else:
                sift = 0
        else:
            sift = 0 
        
        if pd.isna(subdf.loc[i, 'lrt.lrt_pred']) == False:
            if str(subdf.loc[i, 'lrt.lrt_pred'])[0] == 'D':
                lrt = 1
            elif str(subdf.loc[i, 'lrt.lrt_pred'])[0] == 'N':
                lrt = -1
            else:
                lrt = 0 
        else:
            lrt = 0 
        
        if pd.isna(subdf.loc[i, 'mutationtaster.prediction']) == False:
            if str(subdf.loc[i, 'mutationtaster.prediction']).find('Disease') > -1 or str(subdf.loc[i, 'mutationtaster.prediction']).find('Damaging') > -1:
                muttaste = 1
                
            elif  str(subdf.loc[i, 'mutationtaster.prediction']).find('Polymorphism') > -1:
                muttaste = -1
                
            else:
                muttaste = 0
        else:
            muttaste = 0 
            
        if pd.isna(subdf.loc[i, 'mutation_assessor.impact']) == False:
            if str(subdf.loc[i, 'mutation_assessor.impact'])[0] == 'h' or str(subdf.loc[i, 'mutation_assessor.impact'])[0] == 'm':
                mutassessor = 1 
            elif str(subdf.loc[i, 'mutation_assessor.impact'])[0] == 'l' or str(subdf.loc[i, 'mutation_assessor.impact'])[0] == 'n':
                mutassessor = -1 
            else:
                mutassessor = 0      
        else:
            mutassessor = 0
            
        if pd.isna(subdf.loc[i, 'fathmm.fathmm_pred']) == False:
            if str(subdf.loc[i, 'fathmm.fathmm_pred'])[0] == 'D':
                fathmm = 1
            elif str(subdf.loc[i, 'fathmm.fathmm_pred'])[0] == 'T':
                fathmm = -1
            else:
                fathmm = 0
        else:
            fathmm = 0
            
        if pd.isna(subdf.loc[i, 'provean.score']) == False:
            if float(subdf.loc[i, 'provean.score']) <= -2.5:
                provean = 1
            elif float(subdf.loc[i, 'provean.score']) > -2.5:
                provean = -1
            else:
                provean = 0
        else:
            provean = 0 
            
        if pd.isna(subdf.loc[i, 'metasvm.pred']) == False:
            if str(subdf.loc[i, 'metasvm.pred']) == 'Damaging':
                metasvm = 1
            elif str(subdf.loc[i, 'metasvm.pred']) == 'Tolerated':
                metasvm = -1 
            else:
                metasvm = 0 
        else:
            metasvm = 0 
            
        if pd.isna(subdf.loc[i, 'metalr.pred']) == False:
            if str(subdf.loc[i, 'metalr.pred']) == 'Damaging':
                metalr = 1
            elif str(subdf.loc[i, 'metalr.pred']) == 'Tolerated':
                metalr = -1 
            else:
                metalr = 0
                
        else:
            metalr = 0
                
        if pd.isna(subdf.loc[i, 'fathmm_mkl.fathmm_mkl_coding_pred']) == False:
            if str(subdf.loc[i, 'fathmm_mkl.fathmm_mkl_coding_pred']) == 'Damaging':
                fathmm_mkl = 1 
            elif str(subdf.loc[i, 'fathmm_mkl.fathmm_mkl_coding_pred']) == 'Neutral':
                fathmm_mkl = -1
            else:
                fathmm_mkl = 0
        else:
            fathmm_mkl = 0
            
        if pd.isna(subdf.loc[i, 'genocanyon.score']) == False:
            if float(subdf.loc[i, 'genocanyon.score']) > 0.5:
                genocanyon = 1
            elif float(subdf.loc[i, 'genocanyon.score']) <= 0.5:
                genocanyon = -1
            else:
                genocanyon = 0
        else:
            genocanyon = 0
            
        if pd.isna(subdf.loc[i, 'gerp.gerp_rs']) == False:
            if float(subdf.loc[i, 'gerp.gerp_rs']) >= 2.25:
                gerp = 1
            elif float(subdf.loc[i, 'gerp.gerp_rs']) < 2.25:
                gerp = -1
            else:
                gerp = 0
        else:
            gerp = 0
        
        total_pred = (gerp + genocanyon + fathmm_mkl + metalr + metasvm + provean + fathmm + mutassessor + muttaste + lrt + sift)
        
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
    #print('BP4 assigned')

def PP5(subdf):


    PP5_crit = []
        
    for i in range(len(subdf)):    
        if pd.isna(subdf.loc[i, 'clinvar.sig']) == False:
            if (str(subdf.loc[i, 'clinvar.sig']).find('Pathogenic') > -1 or str(subdf.loc[i, 'clinvar.sig']).find('pathogenic') > -1) and (str(subdf.loc[i, 'clinvar.sig']).find('conflict') < 0 or str(subdf.loc[i, 'clinvar.sig']).find('Conflict') < 0) and (str(subdf.loc[i, 'clinvar.rev_stat']).find('expert') > -1 or str(subdf.loc[i, 'clinvar.rev_stat']).find('practical') > -1):
                PP5_crit.append(2)
            elif (str(subdf.loc[i, 'clinvar.sig']).find('Pathogenic') > -1 or str(subdf.loc[i, 'clinvar.sig']).find('pathogenic') > -1) and (str(subdf.loc[i, 'clinvar.sig']).find('conflict') < 0 or str(subdf.loc[i, 'clinvar.sig']).find('Conflict') < 0):
                PP5_crit.append(1) 
            else:
                PP5_crit.append(0)
        else:
            PP5_crit.append(0)
              
    subdf['PP5'] = PP5_crit 

    return print('PP5 assigned')


def BA1(subdf):

    BA1_crit = []

    for i in range(len(subdf)): 
        if pd.isna(subdf.loc[i, 'gnomad3.af']) == False:
            if  float(subdf.loc[i, 'gnomad3.af']) > 0.05:
                BA1_crit.append(1)
            else:
                BA1_crit.append(0)
        else:
            BA1_crit.append(0) 
            
    subdf['BA1'] = BA1_crit

    return print('BA1 assigned')

def BS2(subdf):

    BS2_crit = []

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'hugo']) == False:
            if str(subdf.loc[i, 'hugo']) in adult_list:  
                BS2_crit.append(0)
            elif str(subdf.loc[i, 'hugo']) in rec_list:
                key_snv = str(subdf.loc[i, 'ref_base']) + '_' + str(subdf.loc[i, 'alt_base'])
                try:
                    if key_snv == BS2_hom_het_dict["hom"][str(subdf.loc[i, 'chrom'])][str(subdf.loc[i, 'pos'])]:
                        BS2_crit.append(1)
                    else:
                        BS2_crit.append(0)
                except KeyError:
                    BS2_crit.append(0)
                    
            elif str(subdf.loc[i, 'hugo']) in dom_list:
                key_snv = str(subdf.loc[i, 'ref_base']) + '_' + str(subdf.loc[i, 'alt_base'])
                try:
                    if key_snv == BS2_hom_het_dict["het"][str(subdf.loc[i, 'chrom'])][str(subdf.loc[i, 'pos'])]:
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

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'clinvar.sig']) == False:
            if (str(subdf.loc[i, 'clinvar.sig']).find('Benign') > -1 or str(subdf.loc[i, 'clinvar.sig']).find('benign') > -1) and (str(subdf.loc[i, 'clinvar.rev_stat']).find('expert') > -1 or str(subdf.loc[i, 'clinvar.rev_stat']).find('practical') > -1):
                BS3_crit.append(1)
            else:
                BS3_crit.append(0)
        else:
            BS3_crit.append(0)
            

    subdf['BS3'] = BS3_crit
        
    return print('BS3 assigned')

def BP1(subdf):

    BP1_crit = []

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'so']) == False and pd.isna(subdf.loc[i, 'hugo']) == False:
            if str(subdf.loc[i, 'so']) == 'missense_variant' and str(subdf.loc[i, 'hugo']) in BP1_list:
                BP1_crit.append(1)
            else:
                BP1_crit.append(0)
        else:
            BP1_crit.append(0)
        
    subdf['BP1'] = BP1_crit
        
    return print('BP1 assigned')

def BP6(subdf):

    BP6_crit = []

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'clinvar.sig']) == False:
            if str(subdf.loc[i, 'clinvar.sig']).find('Benign') > -1 or str(subdf.loc[i, 'clinvar.sig']).find('benign') > -1:
                BP6_crit.append(1)
            else:
                BP6_crit.append(0)
        else:
            BP6_crit.append(0)
        
    subdf['BP6'] = BP6_crit
            
            
    return print('BP6 assigned')

def BP7(subdf):

    BP7_crit = []

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'so']) == False:
            if str(subdf.loc[i, 'so']) == 'synonymous_variant':
                if (pd.isna(subdf.loc[i, 'spliceai.ds_ag']) == False or pd.isna(subdf.loc[i, 'spliceai.ds_al']) == False or pd.isna(subdf.loc[i, 'spliceai.ds_dg']) == False or pd.isna(subdf.loc[i, 'spliceai.ds_dl']) == False) and (float(subdf.loc[i, 'spliceai.ds_ag']) >= 0.25 or float(subdf.loc[i, 'spliceai.ds_al']) >= 0.25 or float(subdf.loc[i, 'spliceai.ds_dg']) >= 0.25 or float(subdf.loc[i, 'spliceai.ds_dl']) >= 0.25):
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

    for i in range(len(subdf)):
        if pd.isna(subdf.loc[i, 'gnomad3.af']) == False:
            if float(subdf.loc[i, 'gnomad3.af']) > 0.005:
                BS1_crit.append(1)
            else:
                BS1_crit.append(0)
        else:
            BS1_crit.append(0)

    subdf['BS1'] = BS1_crit

    return print('BS1 assigned')    

def remove():

    sqlite_file = f"{os.path.join(file_output.strip(file_output.split('/')[-1]), file_input.split('/')[-1].split('.')[0])}.sqlite"
    err_file = f"{os.path.join(file_output.strip(file_output.split('/')[-1]), file_input.split('/')[-1].split('.')[0])}.err"
    log_file = f"{os.path.join(file_output.strip(file_output.split('/')[-1]), file_input.split('/')[-1].split('.')[0])}.log"
    
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


    for i in range(len(subdf)):
        PVS = int(subdf.loc[i, 'PVS1'])
        PS = int(subdf.loc[i, 'PS1']) + int(subdf.loc[i, 'PS3'])
        PM = int(subdf.loc[i, 'PM1']) + int(subdf.loc[i, 'PM2']) + int(subdf.loc[i, 'PM4']) + int(subdf.loc[i, 'PM5'])
        PP = int(subdf.loc[i, 'PP2']) + int(subdf.loc[i, 'PP3']) + int(subdf.loc[i, 'PP5']) 
        BS = int(subdf.loc[i, 'BS2']) + int(subdf.loc[i, 'BS3']) + int(subdf.loc[i, 'BS1'])
        BP = int(subdf.loc[i, 'BP3']) + int(subdf.loc[i, 'BP4']) + int(subdf.loc[i, 'BP1']) +  int(subdf.loc[i, 'BP6']) + int(subdf.loc[i, 'BP7'])
        BA1 = int(subdf.loc[i, 'BA1'])
        pc = 0.10
        X = 2
        if BA1 != 1:
            odds_path = 350**((PP/(X**3)) + (PM/X**2) + (PS/X) + (PVS/1) - (BP/X**2) - (BS/X))
            proba = (odds_path * pc) / (((odds_path - 1) * pc) +1)
            probability = float("%.4f" % proba)
            
        elif BA1 == 1:
            probability = 0
            
        probabilities.append(probability)
        
        if BA1 == 1:
            prediction = 'Benign auto'
        elif BS > 2:
            prediction = 'Benign'
        elif (BS == 1 and BP >= 1) or (BP >= 2):
            prediction = 'Likely Benign'
        elif (PVS == 1 and (PS >= 1 or PM >= 2 or (PM == 1 and PP == 1) or PP >= 2)) or (PS >= 2) or (PS == 1 and (PM >= 3 or (PM == 2 and PP >=2) or (PM == 1 and PP >= 4))): 
            prediction = 'Pathogenic'
        elif (PVS == 1 and PM == 1) or (PS == 1 and PM >= 1) or (PS == 1 and PP >= 2) or (PM >= 3) or (PM >= 2 and PP >= 2) or (PM == 1 and PP >= 4):
            prediction = 'Likely Pathogenic'
        else:
            prediction = 'VUS'

        predictions.append(prediction)
        
    subdf['Score'] = probabilities
    subdf['ACMG'] = predictions

    
    return print('Pathogenicity assigned')   


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
    split_df = array_split(df, args.threads)
    m_list = []
    for part in split_df:
        tupl = pd.DataFrame(part)
        m_list.append(tupl.reset_index())

    pool = Pool(processes=args.threads)
    result_proc = pool.map(main, m_list)
    final_pool = pd.concat(result_proc, ignore_index=True, sort=False)

    pool.close()
    pool.join()

else:
    final_pool = main(df)



# for i in range(len(final_pool)):
#     PVS = int(final_pool.loc[i, 'PVS1'])
#     PS = int(final_pool.loc[i, 'PS1']) + int(final_pool.loc[i, 'PS3'])
#     PM = int(final_pool.loc[i, 'PM1']) + int(final_pool.loc[i, 'PM2']) + int(final_pool.loc[i, 'PM4']) + int(final_pool.loc[i, 'PM5'])
#     PP = int(final_pool.loc[i, 'PP2']) + int(final_pool.loc[i, 'PP3']) + int(final_pool.loc[i, 'PP5']) 
#     BS = int(final_pool.loc[i, 'BS2']) + int(final_pool.loc[i, 'BS3']) + int(final_pool.loc[i, 'BS1'])
#     BP = int(final_pool.loc[i, 'BP3']) + int(final_pool.loc[i, 'BP4']) + int(final_pool.loc[i, 'BP1']) +  int(final_pool.loc[i, 'BP6']) + int(final_pool.loc[i, 'BP7'])
#     BA1 = int(final_pool.loc[i, 'BA1'])
#     pc = 0.10
#     X = 2
#     if BA1 != 1:
#         odds_path = 350**((PP/(X**3)) + (PM/X**2) + (PS/X) + (PVS/1) - (BP/X**2) - (BS/X))
#         proba = (odds_path * pc) / (((odds_path - 1) * pc) +1)
#         probability = float("%.4f" % proba)
        
#     elif BA1 == 1:
#         probability = 0
        
#     probabilities.append(probability)
    
#     if BA1 == 1:
#         prediction = 'Benign auto'
#     elif BS > 2:
#         prediction = 'Benign'
#     elif (BS == 1 and BP >= 1) or (BP >= 2):
#         prediction = 'Likely Benign'
#     elif (PVS == 1 and (PS >= 1 or PM >= 2 or (PM == 1 and PP == 1) or PP >= 2)) or (PS >= 2) or (PS == 1 and (PM >= 3 or (PM == 2 and PP >=2) or (PM == 1 and PP >= 4))): 
#         prediction = 'Pathogenic'
#     elif (PVS == 1 and PM == 1) or (PS == 1 and PM >= 1) or (PS == 1 and PP >= 2) or (PM >= 3) or (PM >= 2 and PP >= 2) or (PM == 1 and PP >= 4):
#         prediction = 'Likely Pathogenic'
#     else:
#         prediction = 'VUS'

#     predictions.append(prediction)
    
# final_pool['Score'] = probabilities
# final_pool['ACMG'] = predictions
# print('Pathogenicity assigned')


print('Sorting...')


pat_final_pool = final_pool.loc[final_pool['ACMG'] == 'Pathogenic'].sort_values(by=['Score'], ascending=False)

like_pat_final_pool = final_pool.loc[final_pool['ACMG'] == 'Likely Pathogenic'].sort_values(by=['Score'], ascending=False)

VUS_final_pool = final_pool.loc[final_pool['ACMG'] == 'VUS'].sort_values(by=['Score'], ascending=False)

like_ben_final_pool = final_pool.loc[final_pool['ACMG'] == 'Likely Benign'].sort_values(by=['Score'], ascending=False)

ben_final_pool = final_pool.loc[final_pool['ACMG'].isin(['Benign', 'Benign auto'])].sort_values( by=['Score'], ascending=False)

final_final_pool = pd.concat([pat_final_pool, like_pat_final_pool, VUS_final_pool, like_ben_final_pool, ben_final_pool])

if args.size != None:
    final_final_pool.head(args.size).to_csv(os.path.abspath(args.output), sep = '\t', index=None)
else:
    final_final_pool.to_csv(os.path.abspath(args.output), sep = '\t', index=None)

remove()

print('Success!')

stop = time()

print(f'Finished in {str(ceil(stop - start))} seconds')
