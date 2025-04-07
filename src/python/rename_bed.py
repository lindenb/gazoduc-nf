import sys

def build_chrom_alias_dict(fai_path):
    alias_dict = {}
    with open(fai_path, 'r') as f:
        for line in f:
            chrom = line.strip().split('\t')[0]
            aliases = set()
            aliases.add(chrom)  # original name

            if chrom.startswith('chr'):
                aliases.add(chrom[3:])   # stripped 'chr' prefix
            else:
                aliases.add('chr' + chrom) # add 'chr' version

            # Handle mitochondrial chromosome aliases
            if chrom in ('M', 'MT','chrM','chrMT'):
                aliases.update(['chrM', 'MT', 'chrM', 'chrMT'])

            for alias in aliases:
                alias_dict[alias] = chrom  # map all aliases to canonical name

    return alias_dict

def process_bed_input(chrom_dict):
    for line in sys.stdin:
        if line.startswith('#'):
            continue
        line = line.strip()
        if line == "":
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        if chrom in chrom_dict:
            fields[0] = chrom_dict[chrom]
            print('\t'.join(fields))
        # else: skip the line silently

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python2 script.py reference.fasta.fai < input.bed\n")
        sys.exit(1)

    fai_file = sys.argv[1]
    chrom_aliases = build_chrom_alias_dict(fai_file)
    process_bed_input(chrom_aliases)
