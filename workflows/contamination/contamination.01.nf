/*

Copyright (c) 2022 Pierre Lindenbaum

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



/** path to indexed fasta reference */
params.reference = ""
params.references = "NO_FILE"
params.prefix = ""
params.publishDir = ""
/* give bams AND/OR fastqs */
params.bams = "NO_FILE"
params.fastqs = "NO_FILE"
params.virus_acn = "NO_FILE"
params.help=false
params.min_gc = 0.4
params.max_gc = 0.6
params.min_qual = 30
params.min_read_len = 50
params.max_reads= 1_000_000
params.max_repeat = 6
params.max_N = 2
/* use spades ? */
params.with_spades=false

include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {getVersionCmd;moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {NCBI_DOWNLOAD_TAXDUMP_01} from '../../modules/ncbi/ncbi.download.taxdump.01.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'


workflow {
	CONTAMINATION(params,params.reference,file(params.references),file(params.virus_acn),file(params.bams),file(params.fastqs))
	}


runOnComplete(workflow);

workflow CONTAMINATION {
	take:
		meta
		reference
		references
		virus_acn
		bams
		fastqs
	main:
		version_ch = Channel.empty()

		compile_ch = COMPILE(meta)
		version_ch= version_ch.mix(compile_ch.version)
	
		taxdump_ch = NCBI_DOWNLOAD_TAXDUMP_01(meta)
		version_ch= version_ch.mix(taxdump_ch.version)


		to_zip = Channel.empty()	
		fastq_ch = Channel.empty()

		if(!bams.name.equals("NO_FILE")) {
			all_samples_ch = SAMTOOLS_SAMPLES_01(meta.plus("with_header":true), reference, references, bams)
			version_ch= version_ch.mix(all_samples_ch.version)
			sn_bam_ch = all_samples_ch.output.splitCsv(header:true,sep:'\t')

			unmapped_ch = EXTRACT_BAM_UNMAPPED(meta, compile_ch.executable, sn_bam_ch)
			version_ch= version_ch.mix(unmapped_ch.version)
		
			fastq_ch = fastq_ch.mix(unmapped_ch.output.map{T->[sample:T[0],R1:T[1]]})
			}


		if(!fastqs.name.equals("NO_FILE")) {
			sn2fastq_ch = Channel.fromPath(fastqs).splitCsv(header:true,sep:'\t').
				map{T->{
					if(T.containsKey("R2") && !T.R2.isEmpty()) {
						return [ [T.sample,T.R1], [T.sample,T.R2]  ];
						}
					else {
						return [ [T.sample,T.R1] ];			
					}
				}}.flatMap().groupTuple()

			funmapped_ch = EXTRACT_FASTQ_UNMAPPED(meta, reference, compile_ch.executable, sn2fastq_ch)
			version_ch= version_ch.mix(funmapped_ch.version)
		
			fastq_ch = fastq_ch.mix(funmapped_ch.output.map{T->[sample:T[0],R1:T[1]]})
			}

		vir_ch = BUILD_VIRUS_REF(meta, virus_acn)
		version_ch= version_ch.mix(vir_ch.version)

		compile2_ch = COMPILE_TAX(
			meta,
			vir_ch.acn2taxon, 
			taxdump_ch.nodes,
			taxdump_ch.names
			)
		version_ch= version_ch.mix(compile2_ch.version)

		bwa_vir_ch = BWA_VIRUS(
			meta,
			vir_ch.reference,
			compile2_ch.jar,
			fastq_ch
			)
		to_zip = to_zip.mix(bwa_vir_ch.idxstats)
		to_zip = to_zip.mix(bwa_vir_ch.svg)
		to_zip = to_zip.mix(bwa_vir_ch.freq)
		to_zip = to_zip.mix(bwa_vir_ch.spades)
		version_ch= version_ch.mix(bwa_vir_ch.version)

		dig_ch = DIGEST_IDX(meta, compile2_ch.jar, bwa_vir_ch.idxstats.collect())
		version_ch= version_ch.mix(dig_ch.version)
		to_zip = to_zip.mix(dig_ch.idxstats)

		version_ch = MERGE_VERSION(meta, "Contamination", "Contamination", version_ch.collect())
		to_zip = to_zip.mix(version_ch)

		html =  VERSION_TO_HTML(params,version_ch)
		to_zip = to_zip.mix(html.html)

		zip_ch = SIMPLE_ZIP_01(meta,to_zip.collect())
	emit:
		version = version_ch
		zip = zip_ch.zip
	}


process COMPILE {
input:
	val(meta)
output:
	path("a.out"),emit:executable
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools/1.15.1 gcc")} 

echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}" 1>&2

cat << EOF > jeter.cpp
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cerrno>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <htslib/sam.h>
//#include "htslib/kseq.h"
#include <zlib.h>

//KSEQ_INIT(gzFile, gzread)

using namespace std;

int main(int argc,char** argv) {
	int i=0;
	long n_reads = 0L;
	vector<bool> virus_tids;
	vector<long> counts;
	int ret;
	samFile *fp = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *b = NULL;
	gzFile out = NULL;

	fp = hts_open(argv[1], "r");

	if(fp==NULL) {
	    fprintf(stderr,"Cannot open %s (%s).\\n",argv[1],strerror(errno));
            exit(EXIT_FAILURE);
	    }

	hdr = sam_hdr_read(fp);
	if (hdr == NULL) {
	    fprintf(stderr,"Cannot read header for %s.\\n",argv[1]);
            exit(EXIT_FAILURE);
            }

	hts_set_opt(fp,CRAM_OPT_REFERENCE,argv[2]);

       virus_tids.reserve(hdr->n_targets);

       for(i=0;i< hdr->n_targets; i++) {
		counts.push_back(0L);
                if(strcmp("NC_007605",hdr->target_name[i])==0 ||
		   strcmp("EBV",hdr->target_name[i])==0
		) {
			virus_tids.push_back(true);
			} else {
			virus_tids.push_back(false);
			}
                }
	counts.push_back(0L);/* last is for unmapped */


        out = gzdopen(fileno(stdout), "wb9");
        if (out==NULL) {
	    fprintf(stderr,"Cannot open FASTQ for writing.\\n");
            exit(EXIT_FAILURE);
            }

	b = bam_init1();

	if(b==NULL) {
	    	fprintf(stderr,"out of memory\\n");
		exit(EXIT_FAILURE);		
		}

#define HAS_FLAG(FLAG) (b->core.flag & (FLAG))
	string seq;

	while ((ret = sam_read1(fp, hdr,b)) >= 0) {
		int tid = b->core.tid;
		if(tid < 0 ) {
			counts[hdr->n_targets]++;/* unmapped */
			} else {
			counts[tid]++;
			}

		double gc=0;
		uint8_t *s = bam_get_seq(b);
		uint8_t *q = bam_get_qual(b);
		int len = b->core.l_qseq;
		int i = 0;
		if(${params.max_reads}>0 && n_reads >= ${params.max_reads}) break;
		if (HAS_FLAG(BAM_FSECONDARY)) continue;
		if (HAS_FLAG(BAM_FSUPPLEMENTARY)) continue;
		if (len < ${params.min_read_len} ) continue;
		if (!(HAS_FLAG(BAM_FUNMAP) || virus_tids[tid])) continue;

		seq.resize(len, ' ');
		int N=0;
		for (i = 0; i < len; i++) {
			seq[i] =seq_nt16_str[bam_seqi(s, i)];
			switch(seq[i]) {
				case 'C': case 'c': case 'G': case 'g': gc++; break;
				case 'N': case 'n': N++; break;
				default: break;				
				}
			}
		if (N > ${params.max_N} ) continue;
		double gc_percent = gc/seq.size();
		if ( gc_percent < ${params.min_gc} || gc_percent > ${params.max_gc} ) continue;

		unsigned int j=0;
		int rep=0;
		while(j<seq.size() && rep < ${params.max_repeat}) {
			rep=0;
			char c = seq[j];
			while(j < seq.size() && c == seq[j]) {
				rep++;
				++j;
				}
			}
		if (rep >= ${params.max_repeat}) continue;

		double mean_qual=0.0;
		if (q[0] != 0xff) {
            		for (i = 0; i < len; ++i) {
	                	mean_qual += (int)(q[i]);/* don't add 33 */
				}
			mean_qual /= len;
			if ( mean_qual < ${params.min_qual} ) continue;
			}


		gzprintf(out,"@%ld\\n", n_reads);
		gzwrite(out,(void*)(seq.c_str()),seq.size());
		gzputc(out,'\\n');
		gzputc(out,'+');
		gzputc(out,'\\n');

		// Quality
		if (q[0] == 0xff) {
			for (i = 0; i < len; i++) {
				gzputc(out, 'N');
			}
		} else {
            		for (i = 0; i < len; ++i) {
	                	gzputc(out, q[i]+33);
				}
			}

		if(gzputc(out,'\\n')==-1) {
			fprintf(stderr,"I/O error\\n");
			return EXIT_FAILURE;
			}

		n_reads++;
		}
	/* save counts */
	FILE* xo = fopen(argv[3],"w");
	if(xo==NULL) {
		fprintf(stderr,"I/O error cannot open %s\\n",argv[3]);
                return EXIT_FAILURE;
		}
	for(i=0;i< (int)counts.size();i++) {
		if(i<  hdr->n_targets) {
			fprintf(xo,"%s\t%ld",hdr->target_name[i],(long)hdr->target_len[i]);
			} 
		else
			{
			fprintf(xo,"*\t0");//unmappd
			}
		
		fprintf(xo,"\t%ld\t0\\n",counts[i]);
		}
	fflush(xo);
	fclose(xo);



	sam_hdr_destroy(hdr);
	hts_close(fp);
	bam_destroy1(b);
	gzflush(out,Z_FINISH);
        gzclose(out);


	return EXIT_SUCCESS;
	}

EOF

## FIX THIS FULL PATH!!
g++ -o a.out -Wall -L /CONDAS/users/ifb-tools/miniconda3/envs/samtools-1.15.1/lib -O3 jeter.cpp -lz -lhts

rm jeter.cpp

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Compile filter</entry>
	<entry key="version">${getVersionCmd("gcc samtools")}</entry>
	<entry key="min_gc">${params.min_gc}</entry>
	<entry key="max_gc">${params.max_gc}</entry>
	<entry key="max_reads">${params.max_reads}</entry>
	<entry key="max_repeat">${params.max_repeat}</entry>
	<entry key="min_qual">${params.min_qual}</entry>
	<entry key="min_read_len">${params.min_read_len}</entry>
</properties>
EOF
"""
}


process EXTRACT_BAM_UNMAPPED {
	tag "${row.sample} / ${file(row.bam).name}"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		val(fastqfilter)
		val(row)
	output:
		tuple val("${row.new_sample}"),path("${row.new_sample}.fastq.gz"),path("${row.new_sample}.ctg.count"),emit:output
		path("version.xml"),emit:version
	script:
		def sample = row.new_sample
		def bam = row.bam
	"""
	hostname 1>&2
	${moduleLoad("samtools/1.15")}
        mkdir TMP
	echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}"  1>&2

	${fastqfilter} "${bam}" "${row.reference}" "${sample}.ctg.count" > TMP/tmp.fastq.gz

	mv TMP/tmp.fastq.gz "${sample}.fastq.gz"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract unmapped bams from BAM</entry>
	<entry key="sample">${row.new_sample}</entry>
	<entry key="bam">${row.bam}</entry>
</properties>
EOF
	"""
	}



process EXTRACT_FASTQ_UNMAPPED {
	tag "${sample} N=${L.size()}"
	afterScript "rm -rf TMP"
	cpus 4
	memory "10g"
	input:
		val(meta)
		val(reference)
		val(fastqfilter)
		tuple val(sample),val(L)
	output:
		tuple val("${sample}"),path("${sample}.fastq.gz"),path("${sample}.ctg.count"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("bwa samtools/1.15")}
        mkdir TMP
	echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}"  1>&2
	set -o pipefail

	for F in ${L.join(" ")}
	do
		echo "#\${F}" 1>&2
		bwa mem -t 3 -R '@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tCN:Nantes\\tPL:ILLUMINA' "${reference}" "\${F}" |\
		${fastqfilter} "-" "${reference}" "${sample}.ctg.count" >> TMP/tmp.fastq.gz
	done
	mv TMP/tmp.fastq.gz "${sample}.fastq.gz"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract unmapped bams from FASTQ</entry>
	<entry key="reference">${reference}</entry>
	<entry key="sample">${sample}</entry>
	<entry key="fastqs.count">${L.size()}</entry>
	<entry key="version">${getVersionCmd("bwa")}</entry>
</properties>
EOF
	"""
	}

process BUILD_VIRUS_REF {
	input:
		val(meta)
		path(acns)
	output:
		path("${meta.prefix?:""}virusdb.fa"),emit:reference
		path("acn2taxon.tsv"),emit:acn2taxon
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("samtools bwa")}
	set -x

	echo -n 'db=nucleotide&id=' > jeter.http 


if ${acns.name.equals("NO_FILE")} ; then

## via https://viralzone.expasy.org/678
## "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core" 
cat << EOF | grep -v "#" | tr "," "\\n" | sort | uniq | paste -sd, >> jeter.http
# homo sapiens Y
CP086569
AF033819.3
CP003914
AB009864.2
AB095955.1
AF061923.1
AF061931.1
AF129072.1
AF310136.1
AF327711.1
AF327712.1
AY390769.1
AY753578.1
DD038280.1
DQ058733.1
DQ091275.1
DQ133906.1
DQ391279.1
DQ493876.2
DQ493888.2
DQ666274.1
DQ883650.3
DQ883651.3
DQ883652.3
DQ883653.3
DQ883654.3
DQ883656.3
DQ883658.3
DQ883659.3
DQ883661.3
DQ883663.3
DQ883664.3
DQ883669.3
DQ883670.3
DQ883671.3
DQ883673.3
DQ883682.3
DQ883685.3
DQ883689.3
DQ996013.1
EF437940.1
EF694056.1
EU101022.1
EU140755.1
EU182593.1
EU496089.1
EU496100.1
EU496103.1
EU545992.1
FJ160466.1
FJ710109.1
GQ231553.1
GU339046.1
GU593054.1
GU593055.1
HM162834.1
HQ245711.1
HQ418395.1
J01636.1
J01749.1
J02400.1
J02459.1
J02482.1
JF837313.1
JN204881.1
JN204882.1
JN204883.1
JN204884.1
JN204885.1
JN204886.1
JN204887.1
JN900237.1
JQ927446.1
JX025643.1
JX310867.1
JX560321.2
JX560326.2
JX560334.1
JX560336.2
KC702164.1
KC702165.1
KC702166.1
KC702167.1
KC702168.1
KC702169.1
KC702170.1
KC702171.1
KC702172.1
KC702173.1
KC702175.1
KC702176.1
KC702177.1
KC702178.1
KC702179.1
KC702180.1
KC702181.1
KC702182.1
KC702183.1
KC702184.1
KC702185.1
KC702187.1
KC702189.1
KC702190.1
KC702191.1
KC702192.1
KC702193.1
KC702194.1
KC702195.1
KC702196.1
KC702197.1
KC702198.1
KC702200.1
KC702201.1
KC702202.1
KC702203.1
KC702204.1
KC702205.1
KC702206.1
KC702207.1
KC702208.1
KC702209.1
KC702210.1
KC702211.1
KC702212.1
KC702213.1
KC702214.1
KC702215.1
KC702216.1
KC702217.1
KC702218.1
KC702219.1
KC702220.1
KC702222.1
KC702223.1
KC702224.1
KC702225.1
KC702226.1
KC702227.1
KC702228.1
KC702229.1
KC702230.1
KC702231.1
KC702232.1
KC702233.1
KC702234.1
KC702235.1
KC702236.1
KC702237.1
KC702238.1
KF680544.1
KF680545.1
KP057689.1
KY031724.1
KY031725.1
KY031726.1
KY031732.1
KY031733.1
KY031734.1
KY031735.1
KY031736.1
KY031737.1
L05081.1
L05083.1
L08776.1
L08817.1
L08860.1
L08862.1
L08919.1
L19900.1
LC121524.1
M13163.1
M28829.1
M34519.1
M99566.1
NGB00029.1
NGB00030.1
NGB00034.1
NGB00412.1
NGB00413.1
NGB00585.1
NGB00609.1
NGB01027.1
U01668.1
U02436.1
U02447.1
U09128.1
U09365.1
U09476.1
U13189.1
U29899.1
U37573.1
U56996.1
U75325.1
U75992.1
U89960.1
U89961.1
V00604.2
X65279.1
X66730.1
X70276.1
NC_000913
NC_001401
NC_000883
NC_000943
NC_001347
NC_001348
NC_001352
NC_001356
NC_001364
NC_001405
NC_001430
NC_001434
NC_001436
NC_001437
NC_001449
NC_001477
NC_001479
NC_001489
NC_001498
NC_001512
NC_001526
NC_001538
NC_001542
NC_001544
NC_001545
NC_001547
NC_001560
NC_001563
NC_001608
NC_001611
NC_001612
NC_001617
NC_001653
NC_001664
NC_001699
NC_001710
NC_001716
NC_001731
NC_001781
NC_001786
NC_001798
NC_001802
NC_001806
NC_001809
NC_001897
NC_001906
NC_001918
NC_001943
NC_001959
NC_002031
NC_002058
NC_002076
NC_002200
NC_002549
NC_002642
NC_002645
NC_002728
NC_003215
NC_003243
NC_003417
NC_003461
NC_003663
NC_003687
NC_003690
NC_003899
NC_003908
NC_003977
NC_004102
NC_004162
NC_004718
NC_005179
NC_005336
NC_006429
NC_006430
NC_006554
NC_006560
NC_006998
NC_007580
NC_007605
NC_009238
NC_009333
NC_009527
NC_009539
NC_010277
NC_012532
NC_012800
NC_012957
NC_019843
NC_020805
NC_020806
NC_020807
NC_020810
NC_024070
NC_039025
NC_045512
NC_063383
EOF

else
	cat "${acns}" | paste -sd, >> jeter.http
fi

# EPOST ####################################

wget -O jeter.xml --post-file=jeter.http "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi" 

QUERYKEY=`xmllint --xpath '/ePostResult/QueryKey/text()' jeter.xml`
test ! -z "\${QUERYKEY}"
WEBENV=`xmllint --xpath '/ePostResult/WebEnv/text()' jeter.xml`
test ! -z "\${WEBENV}"

# ESUMMARY ##################################

wget -O jeter.xml --post-file=jeter.http "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi" 

cat << EOF > jeter.xsl
<?xml version='1.0' encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform' version='1.0'>
<xsl:output method="text"  encoding="UTF-8"/>
<xsl:template match="/">
<xsl:apply-templates select="/eSummaryResult/DocSum"/>
</xsl:template>

<xsl:template match="DocSum">
<xsl:value-of select="Item[@Name='Caption']/text()"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="Item[@Name='TaxId']/text()"/>
<xsl:text>
</xsl:text>
</xsl:template>

</xsl:stylesheet>
EOF

xsltproc --novalid jeter.xsl jeter.xml > acn2taxon.tsv

rm jeter.xml jeter.http jeter.xsl

# EFETCH ###################################################
wget -O "${meta.prefix?:""}virusdb.fa" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&query_key=\${QUERYKEY}&WebEnv=\${WEBENV}&rettype=fasta&retmode=text"


# BWA INDEX
bwa index "${meta.prefix?:""}virusdb.fa"
samtools faidx "${meta.prefix?:""}virusdb.fa"
samtools dict -A --output "${meta.prefix?:""}virusdb.dict" "${meta.prefix?:""}virusdb.fa"
cut -f 1 "${meta.prefix?:""}virusdb.fa" > virus.chr.txt

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Build virus database</entry>
	<entry key="version">${getVersionCmd("wget samtools bwa")}</entry>
</properties>
EOF
	"""
	}

process COMPILE_TAX {
input:
	val(meta)
	path(acn2taxon)
	path(taxnodes)
	path(taxnames)
output:
	path("minikit.jar"),emit:jar
	path("version.xml"),emit:version	
script:
"""
mkdir -p TMP
cat << EOF > Minikit.java
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.regex.*;

public class Minikit
	{
	
	private final Map<Integer,Node> id2node=new HashMap<>(100_000);
	private final Map<String,Integer> acn2taxid =new HashMap<>();
	private long n_reads=0L;

	private class Node
		{
		final int id;
		int parent_id=-1;
		final Set<Node> children=new HashSet<Node>();;
		String scientific_name=null;
		long count=0L;
		
		private Node(int id)
			{
			this.id=id;
			this.scientific_name = String.valueOf(this.id);
			}	
		@Override
		public int hashCode()
			{
			return id;
			}
		@Override
		public boolean equals(Object obj)
			{
			if(obj==this) return true;
			return Node.class.cast(obj).id==this.id;
			}
		void visit(final long count) {
			if(count<=0L) return;
			System.err.println("add "+count+" to "+this.scientific_name+" "+id);
			Node v=this;
			while(v!=null) {
				v.count+=count;
				v = id2node.get(v.parent_id);
				}
			}
		
		String id() {
			if(id<0) return "unmapped";
			return "n"+id;
			}
		
		void print(PrintWriter w) {
			if(this.count==0L) return;
			int g=100;
			if(n_reads>0) g= 100 - (int)(( (this.count)/(double)(Minikit.this.n_reads) )*50.0);
			w.println(id()+"[style=filled;fillcolor=\\"gray"+g+"\\";label=\\""+this.scientific_name+" ("+this.count+")\\"];");
			for(Node c:this.children) {
				if(c.count==0L) continue;
				c.print(w);
				w.print(c.id());
				w.print(" -> ");
				w.print(this.id());
				w.println(";");
				}
			}
		
		}
	
	private Node getNodeById(int id)
		{
		Node node=id2node.get(id);
		if(node==null)
			{
			node=new Node(id);
			id2node.put(id,node);
			}
		return node;
		}
	
	private void readNames(BufferedReader in) throws IOException {
		final Pattern delim=Pattern.compile("\\t\\\\\\\\|(\\t)?");
		String line;
		while((line=in.readLine())!=null)
			{
			final String tokens[]=delim.split(line);
			final int id= Integer.parseInt(tokens[0].trim());
			final Node node=this.id2node.get(id);
			if(node==null) continue;
			if(tokens[3].equals("scientific name"))
				{
				node.scientific_name=tokens[1].trim();
				}
			}
		}	
	private void readNodes(BufferedReader in) throws IOException
		{
		final Pattern delim=Pattern.compile("\\t\\\\\\\\|(\\t)?");
		String line;
		while((line=in.readLine())!=null)
			{
			final String tokens[]=delim.split(line);
			
			final Node node = getNodeById(Integer.parseInt( tokens[0].trim()));
			final Node parent= getNodeById(Integer.parseInt( tokens[1].trim()));
			if(parent!=node)
				{
				parent.children.add(node);
				node.parent_id=parent.id;
				}
			}
		}	
	
	
	
	public int instanceMainWithExit(String[] args) {
		try
			{
			try(BufferedReader in=Files.newBufferedReader(Paths.get("${taxnodes.toRealPath()}"))) {
				readNodes(in);
				}
			System.err.println("Number of nodes "+ this.id2node.size());
			try(BufferedReader in=Files.newBufferedReader(Paths.get("${taxnames.toRealPath()}"))) {
				readNames(in);
				}
			final int root_id=1;
			final Node root= this.id2node.get(root_id);
			if(root==null)
				{
				System.err.println("Cannot get node id."+root_id);
				return -1;
				}
			
			final int unmapped_id = -1;
			final Node unmapped_node = new Node(unmapped_id);
			unmapped_node.scientific_name="*";
			root.children.add(unmapped_node);
			unmapped_node.parent_id= root.id;
			
			final Pattern tab = Pattern.compile("[\\t]");
			try(BufferedReader in=Files.newBufferedReader(Paths.get("${acn2taxon.toRealPath()}"))) {
				in.lines().
					map(S->tab.split(S)).
					forEach(T->acn2taxid.put(T[0],Integer.parseInt(T[1])));
				}
			
			try(BufferedReader in=Files.newBufferedReader(Paths.get(args[0]))) {
				String line;
				while((line=in.readLine())!=null) {
					final String[] tokens = tab.split(line);
					String contig = tokens [0];
					final long n = Long.parseLong(tokens[2]) + Long.parseLong(tokens[3]);
					n_reads +=n;
					if(contig.equals("*")) {//unmapped
						unmapped_node.visit(n);
						}
					else
						{
						final int dot = contig.indexOf(".");
						if(dot!=-1 && !acn2taxid.containsKey(contig)) {
							contig= contig.substring(0,dot);
							}
						if(!acn2taxid.containsKey(contig)) {
							System.err.println("undefined acn2taxid "+contig);
							continue;
							}
						final int taxid = acn2taxid.get(contig);
						final Node node  = id2node.get(taxid);
						if(node==null) {
							System.err.println("where is "+taxid+"??");
							continue;
							}
						node.visit(n);
						}
					}
				}
			
			try(PrintWriter pw= new PrintWriter(System.out)) {
				pw.println("digraph G {");
				root.print(pw);
				pw.println("}");
				pw.flush();
			}
			return 0;
			}
		catch(Throwable err)
			{
			err.printStackTrace();
			return -1;
			}
		}
	public static void main(String[] args)
		{
		System.exit(new Minikit().instanceMainWithExit(args));
		}
	}
 
EOF

javac -d TMP Minikit.java
jar cvf minikit.jar -C TMP .

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile taxon2dot</entry>
	<entry key="version">${getVersionCmd("javac")}</entry>
</properties>
EOF
"""

}


process BWA_VIRUS {
	tag "${row.sample}"
	afterScript "rm -rf TMP"
	memory "15g"
	cpus ((params.with_spades as boolean)?16:4)
	input:
		val(meta)
		val(reference)
		val(jar)
		val(row)
	output:
		path("${meta.prefix?:""}${row.sample}.contaminations.idx.stats"),emit:idxstats
		path("${meta.prefix?:""}${row.sample}.contaminations.svg"),emit:svg
		path("${meta.prefix?:""}${row.sample}.freq.txt"),emit:freq
		path("${meta.prefix?:""}${row.sample}.spades.fa.gz"),optional:true,emit:spades
		path("version.xml"),emit:version
	script:
		def with_spades = (params.with_spades as boolean)
	"""
	hostname 1>&2
	${moduleLoad("samtools bwa "+(with_spades?"spades/3.15.2":""))}

        mkdir TMP
	set -x
	bwa mem  -R '@RG\\tID:${row.sample}\\tSM:${row.sample}\\tLB:${row.sample}\\tCN:Nantes\\tPL:ILLUMINA' "${reference}" "${row.R1}" |\
	samtools view -O BAM --uncompressed -F 2304 - |\
	samtools sort -o TMP/sorted.bam  -O BAM -T TMP/st.tmp.sort 
	samtools index "TMP/sorted.bam"

	samtools idxstats "TMP/sorted.bam" > "${meta.prefix?:""}${row.sample}.contaminations.idx.stats" 

	java -Xmx${task.memory.giga}G -cp ${jar}  Minikit "${meta.prefix?:""}${row.sample}.contaminations.idx.stats" > TMP/jeter.dot

	samtools view -f 4 TMP/sorted.bam | samtools fastq - | gzip --best > TMP/unmapped.fastq.gz

	gunzip -c TMP/unmapped.fastq.gz |\
		paste - - - - |\
		cut -f 2 |\
		LC_ALL=C sort -T TMP |\
		LC_ALL=C uniq -c |\
		LC_ALL=C sort -T TMP -n |\
		tail > "${meta.prefix?:""}${row.sample}.freq.txt"

	neato -T svg -o TMP/jeter.svg TMP/jeter.dot
	mv TMP/jeter.svg "${meta.prefix?:""}${row.sample}.contaminations.svg"
	

	if ${with_spades} ; then

		export TMPDIR=\${PWD}/TMP

		spades.py -k 21,33,55,77,99,127 --careful --sc \
			-t ${task.memory.giga} \
			-m ${task.cpus} \
			-o TMP/${row.sample}_assembly \
			-s "TMP/unmapped.fastq.gz"

		gzip --best TMP/${row.sample}_assembly/misc/assembled_contigs.fasta 
		mv "TMP/${row.sample}_assembly/misc/assembled_contigs.fasta.gz" "${meta.prefix?:""}${row.sample}.spades.fa.gz"
	fi

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">map contaminantes</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="bam">${row.R1}</entry>
	<entry key="version">${getVersionCmd("samtools bwa java")}</entry>
</properties>
EOF
	"""
	}


process DIGEST_IDX {
        tag "N=${L.size()}"
	memory "5g"
        input:
              	val(meta)
                val(jar)
		val(L)
	output:
		path("${meta.prefix?:""}ALL_SAMPLES.contaminations.svg"),emit:idxstats
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("samtools bwa spades/3.15.2")}

cat ${L.join(" ")} > jeter.txt

java -Xmx${task.memory.giga}G -cp ${jar}  Minikit jeter.txt > jeter.dot
neato -T svg -o jeter.svg jeter.dot
mv jeter.svg "${meta.prefix?:""}ALL_SAMPLES.contaminations.svg"

rm jeter.dot jeter.txt

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">digest samtools idxstats</entry>
	<entry key="samples.count">${L.size()}</entry>
	<entry key="version">${getVersionCmd("java")}</entry>
</properties>
EOF
	"""
	}
