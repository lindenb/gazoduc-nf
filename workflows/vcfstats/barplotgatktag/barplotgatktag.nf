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

include {runOnComplete;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'
include {GS_SIMPLE_01} from '../../../modules/gs/gs.simple.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'

params.reference=""
params.vcf="NO_FILE"
params.prefix=""
params.publishDir=""
params.help=false


def helpMessage() {
  log.info"""
## About

bcftools stats on mulitple files

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf (file) VCF/BCF [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume vcfstats.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--vcfs /path/to/vcfs.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also


"""
}


if( params.help ) {
    helpMessage()
    exit 0
}



workflow {
	ch1 = BARPLOT_VCF_STATS(params, params.reference,file(params.vcf))
	//PUBLISH(ch1.zip)
	}

workflow BARPLOT_VCF_STATS {
	take:
		meta
		reference
		vcf
	main:
		version_ch = Channel.empty()
		ch1 = SCAN_HEADER(meta,reference,vcf)
		version_ch = version_ch.mix(ch1.version)

		rows_ch = ch1.output.splitCsv(header:true,sep:';')
		
		ch2 = VCF_BARPLOT_TAG_01(meta, rows_ch)
		version_ch = version_ch.mix(ch2.version)
	

		pdf_ch = GS_SIMPLE_01(meta,ch2.output.map{T->["gatktag",T[2]]}.groupTuple())
		version_ch = version_ch.mix(pdf_ch.version)
		
		version_ch = MERGE_VERSION(meta, "BarplotGATKTags", "Barplot GATK tags", version_ch.collect())

	}

def mkloop(def a,def b,def dx) {
	StringBuilder sb = new StringBuilder();
	while(a<=b) {
		if(sb.length()>0) sb.append(",");
		sb.append(a);
		a += dx;
		}
	return sb.toString();
	}

/*
*/

process SCAN_HEADER {
	tag "${vcf.name}"
	executor "local"
	afterScript "rm -rf TMP"
	input:
		val(meta)
		val(reference)
		path(vcf)
	output:
		path("tags.tsv"),emit:output
		path("version.xml"),emit:version	
	script:
	"""
	hostname 1>&2
	${moduleLoad("bcftools")}
	mkdir -p TMP
	
	if ${vcf.name.endsWith(".list")} ; then
		bcftools concat -a --files-list "${vcf}" -O u  | bcftools view --header-only > jeter.header
	else
		bcftools view --header-only "${vcf}" > jeter.header
	fi

	awk -F '[<,=]+' '/^##INFO/ {I="";T="";for(i=1;i<=NF;i++) {if(\$i=="Type") T=\$(i+1);if(\$i=="ID") I=\$(i+1);if(T!="" && I!="") break;} if(T!="" && T!="Flag" && I!="") print I;}' < jeter.header | sort > TMP/jeter.a

	cat <<- EOF | sort -T TMP -t ';' -k1,1 | awk '{printf("%s;${vcf.toRealPath()}\\n",\$0);}' > TMP/jeter.b
	QD;QualityByDepth;${mkloop(0.0, 50.0, 1.0)}
	FS;FisherStrand;${mkloop(0.0, 120.0, 5.0)}
	SOR;StrandOddsRatio;${mkloop(0.0 ,10.0 ,0.5)}
	MQ;RMSMappingQuality;${mkloop(0.0, 100.0, 5.0)}
	MQRankSum;MappingQualityRankSumTest;${mkloop(-10.0, 10.0, 0.5)}
	ReadPosRankSum;ReadPosRankSumTest;${mkloop(-10.0, 10.0, 0.5)}
	EOF


	echo "tag;label;range;vcf" > tags.tsv
	join -t ';' -1 1 -2 1 TMP/jeter.a TMP/jeter.b >> tags.tsv

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Scan INFO header name/type</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="version">${getVersionCmd("bcftools awk")}</entry>
	</properties>
	EOF
	"""
	}


process VCF_BARPLOT_TAG_01 {
executor "local"
tag "${row.tag} EXECUTOR LOCAL FIX IT"
maxForks 1
afterScript "rm -rf TMP"
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("${meta.prefix?:""}${row.tag}.tsv"),path("${meta.prefix?:""}${row.tag}.pdf"),emit:output
	path("version.xml"),emit:version
script:
	def datatype = row.datatype?:"BCF_HT_REAL"
	"""
	hostname 1>&2
	set -o pipefail
	${moduleLoad("bcftools/1.15.1 htslib")}
	set -x


	mkdir -p TMP

	# header only
        if ${row.vcf.endsWith(".list")} ; then
                bcftools concat -a --file-list "${row.vcf}" -O u | bcftools view --header-only > TMP/jeter.header
        else
                 bcftools view --header-only "${row.vcf}" > TMP/jeter.header
        fi

	# extract datatype from header
	DATATYPE=`grep -m1 -F '##INFO=<ID=${row.tag},' TMP/jeter.header | tr "," "\\n" | grep -m1 "^Type=" | cut -d '=' -f 2 | sed 's/Float/BCF_HT_REAL/;s/Integer/BCF_HT_INT/'`
	test ! -z "\${DATATYPE}"

cat << EOF > TMP/jeter.cpp
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <getopt.h>
#include "htslib/vcf.h"
#include "htslib/synced_bcf_reader.h"

using namespace std;

#if DATATYPE== BCF_HT_INT
	typedef  int value_type_t;
	#define from_string(A) (int)strtol(A,&p2,10)
	#define extract_info bcf_get_info_int32	
#else
	typedef float value_type_t;
	#define from_string(A) strtof(A,&p2);
	#define extract_info bcf_get_info_float
#endif


struct Range {
	bool first,last;
	value_type_t a,b;
	unsigned long count;
	bool contains(value_type_t v) const {
		if(first) {
			if(last) return true;
			return v < b;
			}
		if(last) {
			return v >= a;
			}
		return a <= v && v<=b;
		}
	void print(std::ostream& out) const {
		out << "[";
		if(!first) out << a;
		if(!first && !last)  out << "/";
		if(!last) out << b;
		out << "[";
		}
	void print() const {
		print(std::cout);
		}
	};

int main(int argc,char** argv) {
	char* p2= NULL;
	vector<Range> ranges;
	vector<value_type_t> values;
	char number[]="${row.range}";
	char* rest = number;
	for(;;) {
		char* token = strtok_r(rest,",", &rest);
		if(token==NULL) break;
		value_type_t v = from_string(token);
		if(*p2!=0) {
			cerr << "bad string " << token << endl;
			exit(EXIT_FAILURE);
			}
		values.push_back(v);
		}
	sort(values.begin(),values.end());

	Range r0;
	r0.first = true;
	r0.last = false;
	r0.count = 0UL;

	Range r1;
	r1.first = false;
	r1.last = true;
	r1.count = 0UL;

	if(values.empty()){
		r0.first = true;
		r0.last = true;
		ranges.push_back(r0);
		}
	else if(values.size()==1UL)
		{
		value_type_t x1=values[0];
		r0.b = x1;
		r1.a = x1;
		ranges.push_back(r0);
		ranges.push_back(r1);
		}
	else
		{
		r0.b = values[0];
		ranges.push_back(r0);
		for(unsigned int i=0;i+1< values.size();i++) {
			Range r;
			r.count=0UL;
			r.first=false;
			r.last=false;
			r.a=(values[i  ]);
			r.b=(values[i+1]);
			ranges.push_back(r);
			}

		r1.a = values[values.size()-1];
		ranges.push_back(r1);
		}
	bcf_srs_t *sr = bcf_sr_init();
        bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, BCF_SR_PAIR_BOTH_REF);

	if(!bcf_sr_add_reader(sr,"-")) {
		cerr << "Cannot add VCF from stdin.." << endl;
	        if ( sr->errnum ) cerr << "Error: " << bcf_sr_strerror(sr->errnum) << "." << endl;
		return EXIT_FAILURE;
		}
	value_type_t* array=NULL;
	while(bcf_sr_next_line (sr)) {
		bcf1_t *line = sr->readers[0].buffer[0];
		int nvalues=0;
		int ret = extract_info(
			sr->readers[0].header,
			line,
			"${row.tag}",
			&(array),
			&nvalues
			);
            	if (ret<0 ) {
			//cerr << "skip "<< ret << " " << nvalues << endl;
			continue;
			}
		
		for(int j=0; j< nvalues;++j) {
			for(unsigned int i=0UL ;i< ranges.size();i++) {
				if(ranges[i].contains(array[j])) {
					ranges[i].count++;
					break;
					}
				}
			}
		}
	free(array);
        if ( sr->errnum ) {
		cerr << "Error: " << bcf_sr_strerror(sr->errnum) << "." << endl;
		return EXIT_FAILURE;
		}
	bcf_sr_destroy(sr);


	ofstream out("${meta.prefix?:""}${row.tag}.tsv");
	for(unsigned int i=0UL ;i< ranges.size();i++) {
                if(i>0UL) out << "\t";
                ranges[i].print(out);
                }
        out << endl;
	for(unsigned int i=0UL ;i< ranges.size();i++) {
                if(i>0UL) out << "\t";
                out <<  ranges[i].count;
                }
	out << endl;
	out.close();


	while(ranges.size()>1 && ranges[0].count==0UL) {
		ranges.erase(ranges.begin());
		}
	while(ranges.size()>1 && ranges[ranges.size()-1].count==0UL) {
		ranges.pop_back();
		}

	cout << "T1 <- c(";
	for(unsigned int i=0UL ;i< ranges.size();i++) {
		if(i>0UL) cout << ",";
		cout <<  ranges[i].count;
		}	
	cout << ")" << endl;
	cout << "N1 <- c(";
	for(unsigned int i=0UL ;i< ranges.size();i++) {
		if(i>0UL) cout << ",";
		cout << "\\"";
		ranges[i].print();
		cout << "\\"";
		}	
	cout << ")" << endl;
	cout << "pdf(\\"${meta.prefix?:""}${row.tag}.pdf\\")" << endl;
	cout << "barplot(T1,main=\\"INFO/${row.tag} ${row.label}\\",sub=\\"${file(row.vcf).name}\\",names.arg=c(N1),las=2,xlab=\\"${row.tag}\\",ylab=\\"Count\\")" << endl;
	cout << "dev.off()" << endl;
	return EXIT_SUCCESS;
	}
EOF

	g++ -Wall -DDATATYPE=\${DATATYPE} -O3 -std=c++11 -o TMP/a.out TMP/jeter.cpp -lhts

	if ${row.vcf.endsWith(".list")} ; then
		bcftools concat -a --file-list "${row.vcf}" -O u | TMP/a.out > TMP/jeter.R
	else
		./TMP/a.out < "${row.vcf}" > TMP/jeter.R
	fi

	${moduleLoad("r")}
	R --vanilla < TMP/jeter.R
	mv TMP/jeter.R "${meta.prefix?:""}${row.tag}.R"

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">barplot for INFO/${row.tag}</entry>
		<entry key="vcf">${row.vcf}</entry>
		<entry key="version">${getVersionCmd("R bcftools awk r gcc")}</entry>
	</properties>
	EOF
	"""
	}


process PUBLISH {
tag "${zip.name}"
input:
	path(zip)
when:
	!isBlank(params.getOrDefault("publishDir",""))
script:
"""
mkdir -p "${params.publishDir}"
cp -v "${zip}" "${params.publishDir}/"
"""
}

runOnComplete(workflow);
