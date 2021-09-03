#!/usr/bin/python
#-*- coding: utf-8 _*_

"""
Description: this a simple statistics-script for the number of different type RNA in the output of INFERNAL
Usage: python infernal_counter.py --infernal_out cmscan.tblout --Rfam Rfam_anno.tsv --genome_size size --output out.tsv
Date: 2020/12/7
"""
import pandas as pd
import argparse
import sys
def infernal_counter(args):
        with open(args.Rfam) as lines:
                next(lines)
                Rfam = dict()
                tmp1 = []
                tmp2 = []
                for line in lines:
                        line_spl = line.split('\t')
                        tmp1.append(line_spl[0])
                        tmp2.append(line_spl[2])
                Rfam['accession'] = tmp1
                Rfam['type'] = tmp2
        Rfamdf = pd.DataFrame(Rfam)
        # pre-prosessing cmsan
        cmscan_tbl = pd.read_table(args.infernal_out, comment='#', sep = '\t')
        cmscan_tbl.columns = ["accession", "chromosome", "starts", "end"]
        # combine data frames
        combined_df = pd.merge(cmscan_tbl, Rfamdf, how = 'left', on = 'accession')
        combined_df['length'] = abs(combined_df['starts'] - combined_df['end'])
        # fix my data frame
        genome_size = int(args.genome_size)
        Copy = combined_df.groupby("type").agg({'accession':['count']})
        Total_length = combined_df.groupby("type").agg({"length":['sum']})
        Percent = (Total_length/genome_size)*100

        mydf = pd.concat([Copy, Total_length, Percent], axis = 1, ignore_index = True)
        #mydf.rename(columns={0:'Copy', 1:'Total_length(bp)',  2:'Percent of genome'})
        mydf.to_csv(args.output, sep='|', index=True, header = True)
        print("{0:'Copy', 1:'Total_length(bp)',  2:'Percent of genome'}")


if __name__=='__main__':
        parser = argparse.ArgumentParser(description = 'parsing the annotations file from eggnog-mapper and buidding a file')
        parser.add_argument('--Rfam', 
                            dest = 'Rfam', 
                            help = 'total annotation information of Rfam, links: http://rfam.xfam.org/search#tabview=tab5')
        parser.add_argument('--infernal_out', 
                            dest = 'infernal_out', 
                            help = 'parseable table of cmscan hits')
        parser.add_argument('--genome_size', 
                            dest = 'genome_size', 
                            help = 'size of target genome')
        parser.add_argument('--output', 
                            dest = "output", 
                            help = 'The result file')

        if len(sys.argv) <= 1:
                parser.print_help()
                sys.exit()
        args = parser.parse_args()

        infernal_counter(args)
