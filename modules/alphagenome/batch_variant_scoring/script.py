from io import StringIO
#from alphagenome import colab_utils
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
#from google.colab import data_table, files
import pandas as pd
from tqdm import tqdm
import sys

#data_table.enable_dataframe_formatter()
def run_batch_variant_scoring(api_key, variant_filename, organism, sequence_length, scorers):
    # Load the model.
    dna_model = dna_client.create(api_key)

    # Load VCF file containing variants.
    vcf_file = 'placeholder'  # @param

    delim = "\t"

    if variant_filename.endswith(".csv"):
        delim = ","

    vcf = pd.read_csv(variant_filename, sep=delim)

    required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
    for column in required_columns:
        if column not in vcf.columns:
            raise ValueError(f'VCF file is missing required column: {column}.')

      # @param ["human", "mouse"] {type:"string"}

    # @markdown Specify length of sequence around variants to predict:
    # sequence_length : @param ["16KB", "100KB", "500KB", "1MB"] { type:"string" }
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

    all_scorers = scorers
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


if __name__ == '__main__':
    try:
        run_batch_variant_scoring(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4],variant_scorers.RECOMMENDED_VARIANT_SCORERS)
    except Exception:
        print(traceback.format_exc(), file=sys.stderr)
