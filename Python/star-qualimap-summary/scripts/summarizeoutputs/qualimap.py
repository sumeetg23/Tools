################################################################################
# qualimap.py
# Sumeet Gupta
# Created on Jan 2020
# class to handle files generated from Qualimap for RNA-Seq QC
################################################################################


import os
import re
from basicio import basicio

class qualimap():

    def __init__(self, logdir):
        self.qualimap_data = dict()
        print(logdir)
        for f in basicio().getallfiles(logdir, "results.txt"):
            content = basicio().readfile(f)
            parsed_data = self.parse_report(content)
            if parsed_data is not None:
                filepath, filename = os.path.split(f)
                self.qualimap_data[filepath] = parsed_data

    # reformat for the numpy object
    def reformatfornp(self):
        return basicio().reformatinput(self.qualimap_data)

    def parse_report(self, raw_data):

        regexes = {
            'reads_aligned': r"read(?:s| pairs) aligned\s*=\s*([\d,]+)",
            'total_alignments': r"total alignments\s*=\s*([\d,]+)",
            'non_unique_alignments': r"non-unique alignments\s*=\s*([\d,]+)",
            'reads_aligned_genes': r"aligned to genes\s*=\s*([\d,]+)",
            'ambiguous_alignments': r"ambiguous alignments\s*=\s*([\d,]+)",
            'not_aligned': r"not aligned\s*=\s*([\d,]+)",
            '5_3_bias': r"5'-3' bias\s*=\s*([\d,\.]+)$",
            'reads_aligned_exonic': r"exonic\s*=\s*([\d,]+)",
            'reads_aligned_intronic': r"intronic\s*=\s*([\d,]+)",
            'reads_aligned_intergenic': r"intergenic\s*=\s*([\d,]+)",
            'reads_aligned_overlapping_exon': r"overlapping exon\s*=\s*([\d,]+)",
        }

        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1).replace(',', ''))

        if len(parsed_data) == 0: return None
        return parsed_data

    