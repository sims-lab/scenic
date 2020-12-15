"""
Author: Lucy Garner

This script works with outputs from pipeline_pyscenic.py:
    1. Generate regulons from reg.csv file (output from pyscenic ctx)
    2. Output regulons as csv file
"""

import argparse
import pandas as pd
import pickle
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
import logging

parser = argparse.ArgumentParser()
parser.add_argument('--sample', default = 'merged-all',
                    help = 'sample name')
parser.add_argument('--ctx_output', default = 'pyscenic_results.dir/normalised.dir/merged-all_reg.csv',
                    help = 'Output from pyscenic ctx, dataframe of enriched features')
args = parser.parse_args()

logging.info(args)

sample = args.sample

if "raw" in args.ctx_output:
    datatype = "raw"
elif "normalised" in args.ctx_output:
    datatype = "normalised"

# Transform reg.csv output from pipeline to regulons (in reg.csv, each TF can
# be listed multiple times)
# Create regulons from a dataframe of enriched features
df_motifs = load_motifs(args.ctx_output)
regulons = df2regulons(df_motifs)

# Pickle these regulons
with open("pyscenic_results.dir/" + datatype + ".dir/" + sample + "_regulons.P", 'wb') as f:
    pickle.dump(regulons, f)

# Output regulons as a csv file
regulon_df = pd.DataFrame(columns = ["regulon_name", "transcription_factor", "genes",
                                     "weights", "score", "context"])
for i in range(len(regulons)):
    regulon = regulons[i]
    regulon_dict = dict({"regulon_name": regulon.name,
                        "transcription_factor": regulon.transcription_factor,
                        "genes": list(regulon.genes),
                        "weights": list(regulon.weights),
                        "score": regulon.score,
                        "context": list(regulon.context)})
    regulon_df = regulon_df.append(regulon_dict, ignore_index = True)

regulon_df.genes = regulon_df.genes.apply(lambda x: ", ".join(x))
regulon_df.weights = regulon_df.weights.apply(lambda x: str(x))
regulon_df.weights = regulon_df.weights.apply(lambda x: x.replace('[', ''))
regulon_df.context = regulon_df.context.apply(lambda x: ", ".join(x))

regulon_df.to_csv("pyscenic_results.dir/" + datatype + ".dir/" + sample + "_regulons.csv", index = False)

# # Get motif logo for each regulon
# if not os.path.exists('plots.dir/motifs.dir'):
#     os.makedirs('plots.dir/motifs.dir')

# base_url = "http://motifcollections.aertslab.org/v9/logos/"
# motif_id = regulon_df.context
# for value in motif_id:
#     motif_id = value.split(",")[1].strip()
#     url = base_url + motif_id
#     r = requests.get(url, allow_redirects=True)
#     with open('plots.dir/motifs.dir/' + motif_id, 'wb') as f:
#         f.write(r.content)
