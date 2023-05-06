#!/usr/bin/python
#-*- coding: utf-8 _*_

"""
@ Description: Parsing the eggnog-mapper annotation file and kegg json file and buidding txt files for making OrgDb by AnnotationForge::makeOrgPackage.
@ Url: https://www.kegg.jp/kegg-bin/get_htext?ko00001
@ Url: http://www.geneontology.org/page/guide-go-evidence-codes # go-evidence-codes
@ Url: https://www.genome.jp/kegg/pathway.htmls # layer description
@ Usage: python preprosessing_OrgDb.py --input_path
@ Author: Yilei
@ Date:2020-02-06
"""

import sys
import argparse
from pathlib import Path
import json

def preprosessing_OrgDb(args):
        p = Path(args.input_path)

        gene_info_list = list()
        GO_list = list()
        KEGG_list = list()

        K_pathway = list()
        pathway_Description = list()
        K_term = list()
        sep = '\n'

        for child in p.iterdir():
                if child.name.endswith('annotations'):
                        with open(child, 'r') as tsv:
                                for line in tsv:
                                        if line.startswith("#"):
                                                continue
                                        else:
                                                records = line.split("\t")
                                                geneid = records[0]

                                                # gene info Preferred_name and description
                                                gene_name = records[8]
                                                description = records[7]
                                                gene_info_list.append(geneid + "\t" + gene_name + "\t" + description)

                                                # the keys of GO
                                                GO = records[9]
                                                GO_split = GO.split(',')
                                                for i in GO_split:
                                                        if i == '-':
                                                                continue
                                                        else:
                                                        # go-evidence-codes
                                                                combined = geneid + "\t" + i + "\t" + "IEA"
                                                                GO_list.append(combined)
                                                # the keys of K number
                                                KEGG = records[11]
                                                KEGG_split = KEGG.split(',')
                                                for i in KEGG_split:
                                                        if i == '-':
                                                                continue
                                                        else:
                                                                i = i.replace("ko:", "")
                                                                combined = geneid + "\t" + i
                                                                KEGG_list.append(combined)
                elif child.name.endswith('json'):
                        with open(child, 'r') as js:
                                # convert json to dict
                                map_dict = json.load(js)
                                maps = map_dict['children']
                                for layer_1 in maps:
                                        # this layer contain the maps
                                        Pathway_Maps = layer_1['name'][6:]

                                        for layer_2 in layer_1['children']:
                                                # Pathway Identifiers
                                                Pathway_Identifiers = layer_2['name'][6:]

                                                for layer_3 in layer_2['children']:
                                                        # this layer contain PATH information, [PATH:ko00010]
                                                        tmp = layer_3['name'][6:].split('[')
                                                        ko_Description = tmp[0]
                                                        try:
                                                                tmp[1]
                                                        except IndexError:
                                                                continue
                                                        else:
                                                                ko_num =  tmp[1].strip(']').replace('PATH:', '').replace('BR:', '')
                                                                combine_tmp = ko_num + "\t" + ko_Description
                                                                pathway_Description.append(combine_tmp)

                                                        try:
                                                                for layer_4 in layer_3['children']:
                                                                        # this is the last layer of kegg json, contain paired of gene and KO
                                                                        tmp = layer_4['name'].split('  ')
                                                                        K_num = tmp[0]
                                                                        K_Description =  tmp[1]
                                                                        # KEGG Universal enrichment analysis
                                                                        K_term_tmp = K_num + "\t" + K_Description
                                                                        K_term.append(K_term_tmp)
                                                                        K_ko = K_num + "\t" + ko_num
                                                                        K_pathway.append(K_ko)

                                                        except KeyError:
                                                                continue

        # return the resluts
        dup_gene_info_list = list(set(gene_info_list))
        dup_gene_info_list.insert(0, "GID\tGENENAME\tDESCRIPTION")

        gene_info = p.joinpath('gene_info.txt')
        with open(gene_info, 'w+') as gene_info:
                #term_file.write(reslut)
                gene_info.write(sep.join(dup_gene_info_list))


        dup_GO_list = list(set(GO_list))
        dup_GO_list.insert(0, "GID\tGO\tEVIDENCE")

        GID_GO = p.joinpath('GID_GO.txt')
        with open(GID_GO, 'w+') as GID2GO:
                GID2GO.write(sep.join(dup_GO_list))

        dup_KEGG_list = list(set(KEGG_list))
        dup_KEGG_list.insert(0, "GID\tK")

        GID_K = p.joinpath('GID_K.txt')
        with open(GID_K, 'w+') as GID2K:
                GID2K.write(sep.join(dup_KEGG_list))

        dup_K_pathway = list(set(K_pathway))
        dup_K_pathway.insert(0, "K\tPATHWAY")

        K_pathway = p.joinpath('K_pathway.txt')
        with open(K_pathway, 'w+') as K2pathway:
                K2pathway.write(sep.join(dup_K_pathway))

        dup_pathway_Description = list(set(pathway_Description))
        dup_pathway_Description.insert(0, "PATHWAY\tdescription")

        pathway_Description = p.joinpath('pathway_Description.txt')
        with open(pathway_Description, 'w+') as pathway2Description:
                pathway2Description.write(sep.join(dup_pathway_Description))

        dup_K_term = list(set(K_term))
        dup_K_term.insert(0, "K_id\tTERM")

        K_term = p.joinpath('K_term.txt')
        with open(K_term, 'w+') as K2term:
                K2term.write(sep.join(dup_K_term))

if __name__=='__main__':
        parser = argparse.ArgumentParser(description = 'parsing the annotations file from eggnog-mapper and buidding a file')
        parser.add_argument('--input_path',
                                                dest = 'input_path',
                                                help = 'The absolute path of the annotations file from eggnog-mapper and kegg josn, the resluts will alse be saved here.')


        if len(sys.argv) <= 1:
                parser.print_help()
                sys.exit()
        args = parser.parse_args()

        preprosessing_OrgDb(args)

