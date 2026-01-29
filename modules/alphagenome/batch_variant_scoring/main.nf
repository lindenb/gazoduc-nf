/*

Copyright (c) 2026 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include {verify   } from '../../../modules/utils/functions'

process ALPHAGENOME_BATCH_VARIANT_SCORING {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
secret 'ALPHAGENOME_API_KEY'
maxForks 1
conda "${moduleDir}/../../../conda/alphagenome.yml"
input:
    tuple val(meta ),path(variants) /* tsv file variant_id\tCHROM\POS\REF\ALT */
output:
    tuple val(meta),path("*.tsv.gz"),emit:tsv
script:
    def organism  = task.ext.organism?:"human"
    def sequence_length = task.ext.sequence_length ?:"1MB" 
	def prefix = task.ext.prefix?:"${meta.id}.variant_scores"
    // https://www.alphagenomedocs.com/api/models.html#variant-scorers
    def scorers = task.ext.scorers?:"variant_scorers.RECOMMENDED_VARIANT_SCORERS"
"""
hostname 1>&2

ulimit -c unlimited



cat << __EOF__ > script.py
from io import StringIO
#from alphagenome import colab_utils
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
#from google.colab import data_table, files
import pandas as pd
from tqdm import tqdm
import sys

#data_table.enable_dataframe_formatter()

# Load the model.
dna_model = dna_client.create("\${ALPHAGENOME_API_KEY}")

# Load VCF file containing variants.
vcf_file = 'placeholder'  # @param

vcf = pd.read_csv("${variants}", sep=${variants.name.endsWith("csv")?"','":"'\t'"})

required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
for column in required_columns:
  if column not in vcf.columns:
    raise ValueError(f'VCF file is missing required column: {column}.')

organism = '${organism}'  # @param ["human", "mouse"] {type:"string"}

# @markdown Specify length of sequence around variants to predict:
sequence_length = '${sequence_length}'  # @param ["16KB", "100KB", "500KB", "1MB"] { type:"string" }
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
    f'SEQUENCE_LENGTH_{sequence_length}'
]

# @markdown Specify which scorers to use to score your variants:
score_rna_seq = False  # @param { type: "boolean"}
score_cage = False  # @param { type: "boolean" }
score_procap = False  # @param { type: "boolean" }
score_atac = False  # @param { type: "boolean" }
score_dnase = False  # @param { type: "boolean" }
score_chip_histone = False  # @param { type: "boolean" }
score_chip_tf = False  # @param { type: "boolean" }
score_polyadenylation = False  # @param { type: "boolean" }
score_splice_sites = False  # @param { type: "boolean" }
score_splice_site_usage = False  # @param { type: "boolean" }
score_splice_junctions = False  # @param { type: "boolean" }

# @markdown Other settings:
download_predictions = True  # @param { type: "boolean" }

# Parse organism specification.
organism_map = {
    'human': dna_client.Organism.HOMO_SAPIENS,
    'mouse': dna_client.Organism.MUS_MUSCULUS,
}
organism = organism_map[organism]

# Parse scorer specification.
scorer_selections = {
    'rna_seq': score_rna_seq,
    'cage': score_cage,
    'procap': score_procap,
    'atac': score_atac,
    'dnase': score_dnase,
    'chip_histone': score_chip_histone,
    'chip_tf': score_chip_tf,
    'polyadenylation': score_polyadenylation,
    'splice_sites': score_splice_sites,
    'splice_site_usage': score_splice_site_usage,
    'splice_junctions': score_splice_junctions,
}

all_scorers = ${scorers}
selected_scorers = []
for key in all_scorers:
    if  scorer_selections.get(key.lower(), False):
        print(f"OK Scorer '{key}' selected.", file=sys.stderr)
        selected_scorers.append(all_scorers[key])
    else:
        print(f"Scorer '{key}' not selected. Removing from list.", file=sys.stderr)

# Remove any scorers or output types that are not supported for the chosen organism.
unsupported_scorers = [
    scorer
    for scorer in selected_scorers
    if (
        organism.value
        not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
    )
    | (
        (scorer.requested_output == dna_client.OutputType.PROCAP)
        & (organism == dna_client.Organism.MUS_MUSCULUS)
    )
]
if len(unsupported_scorers) > 0:
  print(
      f'Excluding {unsupported_scorers} scorers as they are not supported for'
      f' {organism}.' , 
  )
  for unsupported_scorer in unsupported_scorers:
    selected_scorers.remove(unsupported_scorer)


# Score variants in the VCF file.
results = []

for i, vcf_row in tqdm(vcf.iterrows(), total=len(vcf)):
  variant = genome.Variant(
      chromosome=str(vcf_row.CHROM),
      position=int(vcf_row.POS),
      reference_bases=vcf_row.REF,
      alternate_bases=vcf_row.ALT,
      name=vcf_row.variant_id,
  )
  interval = variant.reference_interval.resize(sequence_length)

  variant_scores = dna_model.score_variant(
      interval=interval,
      variant=variant,
      variant_scorers=selected_scorers,
      organism=organism,
  )
  results.append(variant_scores)

df_scores = variant_scorers.tidy_scores(results)

if download_predictions:
  df_scores.to_csv('variant_scores.tsv', index=False, sep='\t')
  #files.download('variant_scores.csv')

df_scores
__EOF__



python3 script.py

gzip --best variant_scores.tsv
mv variant_scores.tsv.gz ${prefix}.tsv.gz

cat << END_VERSIONS > versions.yml
"${task.process}":
	alphagenome: ""
END_VERSIONS
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}.variant_scores"
"""
touch versions.yml ${prefix}.tsv
gzip ${prefix}.tsv
"""
}
