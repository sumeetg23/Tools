# Created by Sumeet Gupta to modify GTF to be used for building 10x reference
# cellranger multi is sensitive to changes in the attribute formating
# This was done for Arabidopsis gtf file

#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import re
import csv

INPUT_FILE = '/lab/htdata/10xgenomes/thaliana_TAIR_10/carly/genes.gtf'
OUTPUT_FILE = '/lab/htdata/10xgenomes/thaliana_TAIR_10/carly/genesmodifiedbysg.gtf'

def add_transcript_feature(input_df):
    transcripts = {}

    R_SEMICOLON = re.compile(r'\s*;\s+')

    for index, row in input_df.iterrows():
        attributes = row['attributes']
        start = row['start']
        end = row['end']
        
        my_dictionary = dict()
        
        # Parse attributes to get gene_id and transcript_id
        attribute_dict = dict(item.split(' "') for item in re.split(R_SEMICOLON, attributes) if item.strip())
        # print(attribute_dict)
        gene_id = attribute_dict.get('gene_id').strip('"')
        my_dictionary['gene_id'] = gene_id
        if "transcript_id" in attribute_dict:
            transcript_id = attribute_dict.get('transcript_id').replace('"','')

            # Add or update transcript_id in the attributes
            my_dictionary['transcript_id'] = transcript_id
            
            # Keep track of unique transcript_ids
            transcripts[transcript_id] = gene_id
        
        updated_attributes = '; '.join([f'{key} \"{value}\"' for key, value in my_dictionary.items()])

        #print(updated_attributes)
        
        # Update the DataFrame with the new attributes
        input_df.at[index, 'attributes'] = updated_attributes + ";"

    print("Transcript features added successfully.")

    return input_df


gftfile = pd.read_csv(
    filepath_or_buffer=INPUT_FILE, 
    sep='\t', 
    header=None,
    names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
    skiprows=[i for i in range(5)])

df_with_transcripts = add_transcript_feature(gftfile)
 
df_with_transcripts.to_csv(OUTPUT_FILE, sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

