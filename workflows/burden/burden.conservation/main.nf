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
nextflow.enable.dsl=2



include {runOnComplete;moduleLoad;dumpParams;slurpJsonFile} from '../../../modules/utils/functions.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}




workflow 	{
		BURDEN_CONSERVATION(
			params.genomeId,
			file(params.vcf),
			file(params.cases),
			file(params.controls),
			file(params.beds),
			params.bigwig,
			file(params.json),
			file(params.blacklist)
			)
		}


runOnComplete(workflow)

workflow BURDEN_CONSERVATION {
	take:
		genomeId
		vcf
		cases
		controls
		beds
		bigwig
		json_f
		blacklist
	main:
		json = slurpJsonFile(json_f)
		conditions_ch = Channel.of(json.conditions).flatten()
	
		if(beds.name.equals("NO_FILE")) {
			beds_ch = VCF2BEDS(vcf).output
			}
		else
			{
			beds_ch = Channel.fromPath(beds)
			}
		each_bed_ch = beds_ch.splitText().map{it.trim()}

		if(file(bigwig).name.equals("NO_FILE")) {
			bigwig_in = DOWNLOAD_BIGWIG(genomeId).output
			}
		else
			{
			bigwig_in = bigwig;
			}

		genes_ch  = GENES_TO_BED(genomeId)

		per_bed_ch = PER_BED(
			genomeId,
			vcf,
			cases,
			controls,
			bigwig_in,
			blacklist,
			conditions_ch.combine(each_bed_ch)
			)

		merge1_ch=MERGE_CONDITION_TSV(genomeId,genes_ch.output,per_bed_ch.map{[it[0],it[1]]}.groupTuple())

		MERGE_ALL_TSVS(genomeId, merge1_ch.output.collect())

		MERGE_VCFS(per_bed_ch.map{it[2]}.collect())

	}

process VCF2BEDS {
afterScript "rm -rf TMP"
input:
	val(vcf)
output:
	path("beds.list"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP BEDS

if ${file(vcf).name.endsWith(".list")} ; then

	grep -v '^#' '${vcf}' | while read F; do bcftools index -s "\${F}" | awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' ; done |\\
	sort -T TMP | \\
	uniq |\\
	split --lines=1 --additional-suffix=.bed - BEDS/line.
	
else

	bcftools index	-s "${vcf}" | awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' |\\
	split --lines=1 --additional-suffix=.bed - BEDS/line.

fi

find \${PWD}/BEDS/ -type f -name "line*.bed" | sort > beds.list
test -s beds.list
"""
}

process DOWNLOAD_BIGWIG {
input:
	val(genomeId)
script:
	def genome = params.genomeId
	def name= genome.ucsc_name
"""
mkdir -p TMP

"""
}


process PER_BED {
tag "${bed.name}"
memory "5g"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	val(vcf)
	path(cases)
	path(controls)
	val(bigwig)
	val(blacklist)
	tuple val(condition),path(bed)
output:
	tuple val(condition),path("burden.tsv"),path("burden.vcf.gz"),emit:output
script:
	def genome = params.genomes[genomeId]
	def fasta = genome.fasta
	def tag = "CONS"
	def minDP=10
	def maxDP=300
	def lowGQ=60
	def minGQsingleton=90
	def minRatioSingleton=0.2
	def exclude_bed_arg = (file(blacklist).name.equals("NO_FILE")?"":" --targets-file ^${blacklist}")
	if(!condition.containsKey("max_gnomad_AF")) throw new IllegalArgumentException("json.condition.max_gnomad_AF missing");
	if(!condition.containsKey("max_internal_AF")) throw new IllegalArgumentException("json.condition.max_internal_AF missing");
	if(!condition.containsKey("window_size")) throw new IllegalArgumentException("json.condition.window_size missing");
	if(!condition.containsKey("window_shift")) throw new IllegalArgumentException("json.condition.window_shift missing");
	if(!condition.containsKey("min_score")) throw new IllegalArgumentException("json.condition.min_score missing");
"""
hostname 1>&2
${moduleLoad("bedtools bcftools jvarkit")}
mkdir -p TMP BEDS

if ${file(vcf).name.endsWith(".list")} ; then

	grep -v '^#' '${vcf}' | while read F; do bcftools index -s "\${F}" | awk -F '\t' -vF=\${F} '{printf("%s\t0\t%s\t%s\\n",\$1,\$2,F);}' ; done |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n | \\
	uniq > TMP/jeter.vcfs.bed
	
else

	bcftools index	-s "${vcf}" | awk -F '\t' '{printf("%s\t0\t%s\t${vcf}\\n",\$1,\$2);}' |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n | \\
	uniq > TMP/jeter.vcfs.bed

fi

test -s TMP/jeter.vcfs.bed

bedtools intersect -a TMP/jeter.vcfs.bed  -b "${bed}" | cut -f 4 | sort | uniq > TMP/jeter.vcfs.list

test -s TMP/jeter.vcfs.list

bcftools query -l  `head -n1 TMP/jeter.vcfs.list`  | sort | uniq > TMP/jeter.a
cat "${cases}" "${controls}"  | sort | uniq > TMP/jeter.b
comm -12 TMP/jeter.a TMP/jeter.b > TMP/jeter.samples.txt
test -s TMP/jeter.samples.txt


cat << __EOF__ > TMP/jeter.code

if(variant.getNAlleles()<2 || variant.getNAlleles()>3) {
	return false;
	}

/** missing */
final double  missing = variant.getGenotypes().stream().
	filter(G->G.isNoCall()).
	count();

if(missing/variant.getNSamples() > 0.05) return false;


/** filters 
if(variant.isFiltered()) {
	for(final String f: variant.getFilters()) {
		if(f.equals("HENG_LI_MASK")) return false;
		if(f.equals("DUKEMAP20_LT_1")) return false;
		if(f.startsWith("VQSRTranche")) return false;
	}
}
*/

/** low DP */
final double dp= variant.getGenotypes().stream().
	filter(G->G.isCalled() && G.hasDP()).
	mapToInt(G->G.getDP()).average().orElse(${minDP}); 

if(dp <${minDP} || dp > ${maxDP}) return false;

/** low GQ */

final long count_alt = variant.getGenotypes().stream().
        filter(g->g.isCalled() && !g.isHomRef()).
	count();

if(count_alt==0L) return false;

final double count_low_gq = variant.getGenotypes().stream().
	filter(g->g.isCalled() && !g.isHomRef() && g.hasGQ()).
	filter(g->g.getGQ()<${lowGQ}).
	count();

if(count_low_gq/count_alt >= 0.25) {
	return false;
	}


Genotype singleton=null;
for(final Genotype g: variant.getGenotypes()) {
	if(g.isCalled() && !g.isHomRef()) {
		if(singleton!=null) return true;
		singleton=g;
		}
	}
//if(singleton!=null && singleton.isFiltered()) return false;

if(singleton!=null && singleton.isHet() && singleton.hasGQ() && singleton.getGQ()<${minGQsingleton}) {
	return false; 
	}
if(singleton !=null && singleton.hasAD() && singleton.isHet() && singleton.getAD().length==2) 
	{
	int array[]=singleton.getAD();double r= array[1]/(double)(array[0]+array[1]);
	if(r< ${minRatioSingleton} || r>(1.0 - ${minRatioSingleton})) return false;
	}
return true;
__EOF__

set -o pipefail

bcftools concat  --allow-overlaps --regions-file "${bed}" --file-list TMP/jeter.vcfs.list -O u |\\
	bcftools view ${exclude_bed_arg} --samples-file TMP/jeter.samples.txt -O v |\\
	java -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterjdk --nocode  -f TMP/jeter.code |\\
	java -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfbigwig -B "${bigwig}" --tag ${tag}  |\\
	bcftools view -i '${tag} >=${condition.min_score}' -O u |\\
	bcftools norm -f "${fasta}" --multiallelics -any -O u |\\
	bcftools annotate -x 'FILTER,INFO' -O v |\\
	java -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad  --bufferSize 10000 --gnomad "${genome.gnomad_genome}" --fields "AF_POPMAX"  --max-af '${condition.max_gnomad_AF}' |\\
	bcftools view  -e  'ALT=="*" || AC==0 ||  AF> ${condition.max_internal_AF} || FILTER ~ "GNOMAD_GENOME_BAD_AF" || FILTER ~ "GNOMAD_GENOME_InbreedingCoeff" || FILTER ~ "GNOMAD_GENOME_RF" || FILTER ~ "GNOMAD_GENOME_InbreedingCoeff"'  |\\
	java -Djava.io.tmpdir=TMP -jar ${JVARKIT_DIST}/jvarkit.jar vcfburdenslidingwindow --cases  ${cases} --controls ${controls}  -t '${condition.treshold?:1E-5}' -w ${condition.window_size} -s ${condition.window_shift} --save-vcf TMP/burden.vcf.gz > TMP/burden.tsv

set +o pipefail

mv TMP/burden.tsv ./
mv TMP/burden.vcf.gz ./
mv TMP/burden.vcf.gz.tbi ./
"""
}


process GENES_TO_BED {
afterScript "rm -rf TMP"
input:
	val(genomeId)
output:
	path("genes.bed"),emit:output
script:
	def genome = params.genomes[genomeId]
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP BEDS

wget -O - https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz | gunzip -c | cut -f3,9 |\
	LC_ALL=C sort -T TMP -k1,1 --unique > TMP/ncbi.tsv

gunzip -c "${genome.gtf}" |\\
	awk -F '\t' '(\$3=="gene")' |\\
	java  -Djava.io.tmpdir=.  -jar \${JVARKIT_DIST}/jvarkit.jar gtf2bed -c 'gene_name,gene_id,gene_biotype' |\
	LC_ALL=C sort -T TMP -k4,4 > TMP/gtf.txt

LC_ALL=C join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6,2.2' -a 1 -e '.' TMP/gtf.txt TMP/ncbi.tsv |\
	LC_ALL=C sort -T TMP. -t \$'\\t' -k1,1 -k2,2n |\\
	uniq > genes.bed

test -s genes.bed

"""
}

process MERGE_CONDITION_TSV {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(genes)
	tuple val(condition),val(L)
output:
	path("*.closest.bed"),emit:output
script:
	def col = condition.entrySet().collect{it.getKey()+":"+it.getValue()}.join(";")
	def md5 = col.md5() 
"""
hostname 1>&2
${moduleLoad("bedtools")}
set -o pipefail
mkdir -p TMP

head -n 1 "${L[0]}" | tr "\\n" "\t" > TMP/jeter.bed
echo "gene.chrom;gene.start;gene.end;gene_name;gene_id;gene_biotype;gene_desc;distance;analysis_name;analysis_params" | tr ";" "\t" >> TMP/jeter.bed


cat ${L.join(" ")} |\\
	grep -v '^#' |\\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools closest -a - -b '${genes}' -d |\
	awk -F '\t' '{printf("%s\t${md5}\t${col}\\n",\$0);}' |\
	LC_ALL=C sort -t '\t' -k6,6g >> TMP/jeter.bed

mv -v TMP/jeter.bed "${params.prefix?:""}${md5}.closest.bed"
"""
}


process MERGE_ALL_TSVS {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
        val(genomeId)
        val(L)
output:
        path("${params.prefix?:""}closest.bed"),emit:output
        path("${params.prefix?:""}closest.md")
script:
"""
set -o pipefail
mkdir -p TMP
head -n 1 "${L[0]}"  > TMP/jeter.bed

cat ${L.join(" ")} |\\
        grep -v '^#' |\\
        LC_ALL=C sort -T . -t '\t' -k6,6g >> TMP/jeter.bed


cat << 'EOF' > TMP/jeter.awk
        {
        printf("| ");
        for(i=1;i<=NF;i++) printf("%s%s",(i==1?"":" | "),\$i);
        printf(" |\\n");
        if(NR==1) {
                printf("| ");
                for(i=1;i<=NF;i++) printf("%s:----------",(i==1?"":" | "));
                printf(" |\\n");
                }
        }
EOF

cut -f 4-11,15-19,21 TMP/jeter.bed | awk -F '\t' -f TMP/jeter.awk >  "${params.prefix?:""}closest.md"

mv -v TMP/jeter.bed "${params.prefix?:""}closest.bed"




"""
}

process MERGE_VCFS {
tag "N=${L.size()}"
input:
	val(L)
output:
	path("${params.prefix?:""}concat.vcf.gz"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

bcftools concat -a --remove-duplicates  -O u   ${L.join(" ")}  |\
	bcftools sort -T ./tmp -O z -o "${params.prefix?:""}concat.vcf.gz"
"""
}
