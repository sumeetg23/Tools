#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import re
import csv

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
    filepath_or_buffer='/lab/htdata/10xgenomes/thaliana_TAIR_10/carly/genes.gtf', 
    sep='\t', 
    header=None,
    names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
    skiprows=[i for i in range(5)])

df_with_transcripts = add_transcript_feature(gftfile)

output_gtf_file = '/lab/htdata/10xgenomes/thaliana_TAIR_10/carly/genesmodifiedbysg.gtf'
df_with_transcripts.to_csv(output_gtf_file, sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

