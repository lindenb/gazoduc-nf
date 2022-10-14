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

include {isHg19;runOnComplete;moduleLoad;getKeyValue;hasFeature;getVersionCmd} from '../../../modules/utils/functions.nf'
include {BURDEN_SAMPLE_WGSELECT_PART_01}  from '../../../subworkflows/burden/burden.samples.wgselect.part.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {RVTESTS_REHEADER_01} from '../../../modules/rvtests/rvtests.reheader.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'

params.reference=""
params.pedigree=""
params.vcf=""
params.help=false
params.bed_cluster_method = " --size 1mb "
params.bed=""
params.levels=3

if(params.help) {
  log.info"""
## About

Burden for 1st intron.

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --pedigree <file> jvarkit formatted pedigree. phenotype MUST be case|control. Sex MUST be male|female|unknown
  * --bed <file> CHROM/START/END/GENE_NAME.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume worfklow.nf \\
        --publishDir output \\
        --prefix "analysis." \\
        --reference /path/to/reference.fasta \\
        --vcf /path/to/my.vcf.gz \\
        --pedigree /path/to/input.ped \
```

## Workflow

![workflow](./workflow.svg)

"""
exit 0
}

workflow {
		burden_ch = BURDEN_PAIRS(params, params.reference, params.vcf, file(params.pedigree),file(params.bed))
		ZIPIT(params,burden_ch.zip.collect())
		}

workflow BURDEN_PAIRS {
	take:
		meta
		reference
		vcf
		pedigree
		bed
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()

		pairs_ch = DIGEST_PAIRS(meta, reference, bed)
		version_ch = version_ch.mix(pairs_ch.version)
	
		ch1_ch = BURDEN_SAMPLE_WGSELECT_PART_01(meta,reference,vcf, pedigree, pairs_ch.merged_bed)
		version_ch = version_ch.mix(ch1_ch.version)
		

		header_ch = RVTESTS_REHEADER_01(meta, reference)
		version_ch = version_ch.mix(header_ch.version)

		merged_ch = MERGE_VCF_CONTIGS(meta, ch1_ch.contig_vcfs)
		version_ch = version_ch.mix(merged_ch.version)


		assoc_ch = RVTEST_BY_PAIR(meta,
				reference,
				merged_ch.vcf,
				ch1_ch.rvtest_pedigree,
				header_ch.output, 
				pairs_ch.setFiles.splitText().map{it.trim()}
				)
		version_ch = version_ch.mix(assoc_ch.version)
		
		concat_ch = CONCAT_FILES_01(meta,assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference, file(vcf).name ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "burden 1st intron", "Burden 1st intron ${vcf}", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
	emit:
		version = version_ch
		zip = to_zip
	}

process DIGEST_PAIRS {
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(bed)
output:
	path("merged.bed"),emit:merged_bed
	path("setfile.list"),emit:setFiles
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bedtools")}
set -o pipefail


mkdir -p TMP BEDS

${bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${bed}" | cut -f1-4 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/digest.bed

# test no empty name
test `awk -F '\t' '(\$4=="")' TMP/digest.bed |wc -l` -eq 0

cat << EOF > TMP/Minikit.java

import java.util.*;
import java.io.*;
import java.util.function.Function;
import java.util.stream.Collectors;


public class Minikit {

private final Map<String,Gene> name2gene = new HashMap<>();
private PrintWriter writer = null;
private int num_chuncks = 0;
private long chunk_size =0L;
private final static long MAX_CHUNK_SIZE = 5_000_000L;
private final static int MAX_DEPTH= ${params.levels} ;

private class Interval implements Comparable<Interval> {
	final String contig;
	final int start;
	final int end;
	Interval(String contig,int start,int end) {
		this.contig = contig;
		this.start =start;
		this.end=end;
		}
	@Override
	public int compareTo(final Interval o) {
		int i=this.contig.compareTo(o.contig);
		if(i!=0) return i;
		i= Integer.compare(start,o.start);
		if(i!=0) return i;
		return Integer.compare(end,o.end);
		}
	boolean intersects(final Interval o) {
		return this.contig.equals(o.contig) && !(this.end <= o.start || this.start >= o.end);
		}
	int size() {
		return 1+(end - start);
		}
	@Override
	public String toString() {
		return this.contig+":"+start+"-"+end;
		}
	}

private static class Gene {
	final String name;
	final List<Interval> intervals = new ArrayList<>();
	Gene(final String name) {
		this.name = name;
		}
	}

private void merge(final List<Interval> intervals) {
		Collections.sort(intervals);
		int i=0;
		while(i+1<intervals.size()) {
			if(intervals.get(i).intersects(intervals.get(i+1))) {
				intervals.set(i,new Interval(
					intervals.get(i).contig,
					Math.min(intervals.get(i).start, intervals.get(i+1).start),
					Math.max(intervals.get(i).end, intervals.get(i+1).end)
					));
				intervals.remove(i+1);
				}
			else
				{
				i++;
				}
			}
		}


private void recursive(
	final List<Gene> genes,
	final int[] indexes,
	int index,
	int depth
	) throws IOException
	{
	if(depth>=indexes.length)
		{
		final List<Interval> intervals = new ArrayList<>();
		final Set<String> names = new TreeSet<>( );
		
		//for(int i=0;i< indexes.length;i++) System.err.print(genes.get(indexes[i]).name+" ");System.err.println();
		
		for(int i=0;i< indexes.length;i++) {
			final Gene g = genes.get(indexes[i]);
			names.add(g.name);
			intervals.addAll(g.intervals);
			}
		merge(intervals);
		
		if(this.writer==null) {
			this.writer = new PrintWriter("BEDS/chunk."+(this.num_chuncks++)+".setfile");
			this.chunk_size = 0;
			}
		this.writer.print(String.join("_",names));
		this.writer.print("\t");
		this.writer.println(intervals.stream().map(T->T.toString()).collect(Collectors.joining(",")));
		
		this.chunk_size+= intervals.stream().mapToInt(R->R.size()).sum();
		
		if(this.chunk_size > MAX_CHUNK_SIZE && this.writer!=null) {
			//System.err.println("Closing "+this.chunk_size);
			this.writer.flush();
			this.writer.close();
			this.writer = null;
			this.chunk_size=0L;
			}
		return;
		}
	final int[] indexes2 = Arrays.copyOf(indexes,indexes.length);
	for(int i=0;i< genes.size();i++)
		{
		if(depth>0 && indexes2[depth-1]>i) continue;
		indexes2[depth]=i;
		recursive(genes,indexes2,0,depth+1);
		}
	
	}

private int instanceMain(final String[] args) {
	try {
		try(BufferedReader br=new BufferedReader(new InputStreamReader(System.in))) {
			String line;
			while((line=br.readLine())!=null) {
				final String[] tokens = line.split("[\\t]");
				final String name = tokens[3];
				Gene g = name2gene.get(name);
				if(g==null) {
					g=new Gene(name);
					name2gene.put(name,g);
					}
				final int start = Integer.parseInt(tokens[1])+1;
				final int end = Integer.parseInt(tokens[2]);
				g.intervals.add(new Interval(tokens[0],start,end));
				}
			}
		final List<Gene> genes = new ArrayList<>(name2gene.values());
		
		for(Gene g: genes ) {
			merge(g.intervals);
			}
		
		this.num_chuncks=0;
		this.writer = null;
		final int[] indexes = new int[Math.min(MAX_DEPTH,genes.size())];
		Arrays.fill(indexes,0);
		recursive(genes,indexes,0,0);
		if(this.writer!=null) {
			this.writer.flush();
			this.writer.close();
			}	
		
		return 0;
		}
	catch(Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}

private void instanceMainWithExit(final String[] args) {
	System.exit(instanceMain(args));
	}

public static void main(String[] args) {
	new Minikit().instanceMainWithExit(args);
	}
}
EOF

javac -d TMP -sourcepath TMP TMP/Minikit.java
java -cp "TMP" Minikit < TMP/digest.bed


## find all setfiles
find \${PWD}/BEDS -type f -name "*.setfile" > setfile.list
test -s setfile.list

# merge all regions for wgselect
cut -f1,2,3 TMP/digest.bed |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > merged.bed



###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">digest input bed, create gene list</entry>
	<entry key="versions">${getVersionCmd("bedtools jvarkit/bedrenamechr")}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
}

process MERGE_VCF_CONTIGS {
tag "${vcfs}"
cpus 10
input:
	val(meta)
	path(vcfs)
output:
	path("merged.vcf.gz"),emit:vcf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}

bcftools concat --threads ${task.cpus} -a --remove-duplicates --file-list "${vcfs}" -O z -o "merged.vcf.gz"
bcftools index  --threads ${task.cpus} --tbi merged.vcf.gz


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge vcf</entry>
	<entry key="versions">${getVersionCmd("bcftools")}</entry>
	<entry key="vcfs">${vcfs}</entry>
</properties>
EOF
"""
}

process RVTEST_BY_PAIR {
tag "${file(setfile).name}"
input:
	val(meta)
	val(reference)
	path(vcf)
	path(pedigree)
	path(reheader)
	val(setfile)
output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
when:
	true
script:
	def  rvtest_params = "--burden 'cmc,exactCMC,zeggini' --kernel 'skato'"

	// --burden cmc,zeggini,mb,fp,exactCMC,cmcWald,rarecover,cmat --vt price,analytic --kernel 'skat[nPerm=1000],kbac,skato'

"""
hostname 1>&2
${moduleLoad("rvtests")}
set -o pipefail
set -x
mkdir -p ASSOC


rvtest  --noweb \
        --inVcf ${vcf.toRealPath()} \
	--setFile "${setfile}" \
	--pheno "${pedigree}" \
	--out "ASSOC/part" \
	${rvtest_params} 1>&2 2> last.rvtest.log

find \${PWD}/ASSOC -type f -name "part.*assoc" > assoc.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">rvtest for set of genes</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("rvtest")}</entry>
	<entry key="setfile">${setfile}</entry>
</properties>
EOF
"""
}


process ZIPIT {
tag "N=${files.size()}"
publishDir "${meta.publishDir}" , mode: 'copy', overwrite: true
input:
	val(meta)
        val(files)
output:
        path("${meta.prefix?:""}archive.zip")
when:
        !meta.getOrDefault("publishDir","").trim().isEmpty()
script:
        prefix = meta.getOrDefault("prefix","")
"""

mkdir "${prefix}archive"

cat << EOF | while read F ; do ln -s "\${F}" "./${prefix}archive/" ; done
${files.join("\n")}
EOF

zip -r "${prefix}archive.zip" "${prefix}archive"
"""
}


runOnComplete(workflow);
