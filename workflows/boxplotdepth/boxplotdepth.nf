/*

Copyright (c) 2023 Pierre Lindenbaum

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
params.mapq = 1
params.bams = ""
params.bed = ""
params.help = false
params.publishDir = ""
params.prefix = ""
params.sample2group="NO_FILE"
params.proportional_width=true

include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;runOnComplete;parseBoolean} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'

def helpMessage() {
  log.info"""
## About

apply mosdepth to a set of bams.

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --bed (file) required bed CHROM(tab)START(tab)END(tab)TRANSCRIPT(tab)(INTERVAL_NAME). default: ""
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume boxplotdepth.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--bed /path/to/file.bed
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
	ch1 = BOXPLOTDEPTH(params,params.reference,file(params.references),file(params.bams),file(params.bed),file(params.sample2group))
	PUBLISH(ch1.zip)
	}

workflow BOXPLOTDEPTH {
	take:
		meta
		reference
		references
		bams
		bed
		sample2group
	main:
		version_ch  = Channel.empty()
		bams_ch = SAMTOOLS_SAMPLES_01(meta.plus(["with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false]),reference,file("NO_FILE"),bams)
		version_ch = version_ch.mix(bams_ch.version)
	
		sn_col_ch = JOIN_GROUPS(meta,bams_ch.output,sample2group)
		version_ch = version_ch.mix(sn_col_ch.version)

		merge_bed_ch = MERGE_BED(meta,reference,bed)
		version_ch = version_ch.mix(merge_bed_ch.version)


		dp_ch = SAMPLEDEPTH(meta,reference, merge_bed_ch.bed, sn_col_ch.output.splitCsv(header:true,sep:'\t'))
		version_ch = version_ch.mix(dp_ch.version)

		ch1_ch = CONCAT_FILES_01([:],dp_ch.output.collect())
		version_ch = version_ch.mix(ch1_ch.version)

		kit_ch = COMPILE_PLOTTER(meta)
		version_ch = version_ch.mix(kit_ch.version)

		each_transcript = merge_bed_ch.transcripts.splitText().map{it.trim()}

		to_zip = Channel.empty()
		plot_ch = PLOT(meta,reference,bed,each_transcript,ch1_ch.output, kit_ch.jar )
		version_ch = version_ch.mix(plot_ch.version)
		to_zip = to_zip.mix(plot_ch.pdf)
		to_zip = to_zip.mix(plot_ch.R)

		version_ch = MERGE_VERSION(meta, "BoxplotDEPTH", "Box Plot DEPTH", version_ch.collect())

		to_zip = to_zip.mix(version_ch.version)
		
		html = VERSION_TO_HTML(meta,version_ch.version)
		to_zip = to_zip.mix(html.html)

		zip_ch = SIMPLE_ZIP_01([:],to_zip.collect())
	emit:
		version= version_ch
		zip = zip_ch.zip
}

process MERGE_BED {
input:
	val(meta)
	val(reference)
	path(bed)
output:
	path("merged.bed"),emit:bed
	path("transcripts.txt"),emit:transcripts
	path("version.xml"),emit:version
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
set -o pipefail

grep -v -E "^(browser|track|#)" "${bed}" |\
	cut -f  1,2,3 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	sort -T . -t '\t' -k1,1 -k2,2n  |\
	bedtools merge > merged.bed

test -s merged.bed

grep -v -E "^(browser|track|#)" "${bed}"  | cut -f 4 | sort | uniq > transcripts.txt
test -s transcripts.txt


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge bed intervals</entry>
	<entry key="bed">${bed}</entry>
	<entry key="bedtools.version">\$(bedtools --version)</entry>
</properties>
EOF
"""
}

                sn_col_ch = 
process JOIN_GROUPS {
executor "local"
input:
	val(meta)
	path(snbams)
	path(sample2group)
output:
	path("samples.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
set -o pipefail
if ${sample2group.name.equals("NO_FILE")} ; then

	awk -F '\t' '{printf("%s\t%s\t%s\tALL\\n",\$1,\$3,\$4);}' "${snbams}" > samples.tsv

else

	join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,2.2' \
		<( cut -f1,3,4 "${snbams}" | sort -T . -t '\t' -k1,1) \
		<( cut -f1,2 "${sample2group}" | sort -T . -t '\t' -k1,1 --unique | uniq) > samples.tsv

fi


test -s samples.tsv

echo -e "sample\tbam\treference\tcollection" > jeter.txt
cat samples.tsv >> jeter.txt
mv jeter.txt samples.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">join samples and groups</entry>
	<entry key="groups">${sample2group}</entry>
</properties>
EOF
"""
}

process SAMPLEDEPTH {
tag "${row.sample} ${row.collection} ${bed.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(bed)
	val(row)
output:
	path("${row.sample}.txt"),emit:output
	path("version.xml"),emit:version
script:
	def mapq=meta.mapq?:10
"""
hostname 1>&2
${moduleLoad("samtools tabix jvarkit")}
set -o pipefail
mkdir TMP

java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${row.reference}" --column 1 --convert SKIP "${bed}" > TMP/jeter.bed

samtools view -F 3844 -q "${mapq}" -O BAM --uncompressed  -T "${row.reference}" -M -L "TMP/jeter.bed" "${row.bam}"  |\
	samtools depth -Q "${mapq}" -a -b "TMP/jeter.bed" - > TMP/jeter.tsv

java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP TMP/jeter.tsv > "TMP/${row.sample}.depth.tsv"


bgzip "TMP/${row.sample}.depth.tsv"
tabix -s 1 -b 2 -e 2 "TMP/${row.sample}.depth.tsv.gz"

mv "TMP/${row.sample}.depth.tsv.gz" ./
mv "TMP/${row.sample}.depth.tsv.gz.tbi" ./

echo "${row.sample}\t\${PWD}/${row.sample}.depth.tsv.gz\t${row.collection}" > "${row.sample}.txt"


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">samtools depth</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="mapq">${mapq}</entry>
	<entry key="samtools.version">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
</properties>
EOF
"""
}


process COMPILE_PLOTTER {
executor "local"
afterScript "rm -rf TMP"
input:
	val(meta)
output:
	path("minikit.jar"),emit:jar
	path("version.xml"),emit:version
script:
"""
${moduleLoad("jvarkit")}

mkdir TMP

cat << "__EOF__" > TMP/Minikit.java

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.function.ToIntFunction;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;

public class Minikit {
private static final Pattern TAB = Pattern.compile("[\t]");	
	
private static class Sample {
	@SuppressWarnings("unused")
	String name;
	String tabix;
	String collection;
	
	public double getMedianDepth(Interval r) {
		final List<Integer> array = new ArrayList<>(r.getLengthOnReference());
		try(TabixReader tbr = new TabixReader(this.tabix)) {
			final TabixReader.Iterator tr = tbr.query(r.getContig(), r.getStart(), r.getEnd());
			for(;;) {
				final String line=tr.next();
				if(line==null) break;
				array.add(Integer.parseInt(TAB.split(line)[2]));
				}
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		if(array.isEmpty()) return 0;
		Collections.sort(array);
		final int mid_x= array.size()/2;
		if(array.size()%2==0)
	        {
			return (array.get(mid_x-1)+ array.get(mid_x))/2.0;
	        }
		else
	        {
	        return array.get(mid_x);
	        }
		}
	}

private String quote(Object o) {
	return "\\""+ o +"\\"";
}

private void instanceMainWithExit(final List<String> args) {
	Set<String> collectionsSet = new TreeSet<>();
	final List<Sample> samples = new ArrayList<>();
	final List<Interval> intervals = new ArrayList<>();
	try {
		try (BufferedReader br=Files.newBufferedReader(Paths.get(args.get(0)))) {
			br.lines().map(S->TAB.split(S)).forEach(T->
				{
				final Sample sn = new Sample();
				sn.name=T[0];
				sn.tabix= T[1];
				sn.collection=T[2];
				samples.add(sn);
				collectionsSet.add(sn.collection);
				});
			}
		final String[] titles= new String[1];
		try (BufferedReader br=Files.newBufferedReader(Paths.get(args.get(1)))) {
			br.lines().map(S->TAB.split(S)).forEach(T->
				{
				final int start = Integer.parseInt(T[1])+1;
				final int end  = Integer.parseInt(T[2]);
				final Interval r = new Interval(
						T[0],
						start,
						end,
						false,
						(T.length >4?T[4]:T[0]+":"+start+"-"+end)
						);
				titles[0] = T[3];
				intervals.add(r);
				});
			}
		Collections.sort(intervals);
		final String title= titles[0];
		final List<String> collections = new ArrayList<>(collectionsSet);
		final int min_length_on_reference  = intervals.stream().mapToInt(R->R.getLengthOnReference()).min().orElse(1);
		final ToIntFunction<Integer> toWidth = L->(int)((L/(double)min_length_on_reference)*50.0);

		final boolean proportional_width=false;//${parseBoolean(meta.proportional_width)};

		try(PrintWriter w = new PrintWriter(System.out)) {
			final List<String> components= new ArrayList<>();
			final List<String> widths= new ArrayList<>();
			final List<String> names= new ArrayList<>();
			w.println("# intervals:");
			w.println("# " + intervals.stream().map(R->R.toString()).collect(Collectors.joining(",")));
			w.println("# groups:");
			w.println("# " + collections.stream().map(R->R.toString()).collect(Collectors.joining(",")));
			for(int ii=0;ii< intervals.size();ii++) {
				for(int g = 0;g< collections.size();++g) {
					final String col = collections.get(g);
					final Interval the_interval = intervals.get(ii);
					final String s = "G"+g+"E"+ii;
					components.add(s);
					widths.add(String.valueOf(toWidth.applyAsInt(the_interval.getLengthOnReference())));
					names.add(quote(g==0?the_interval.getName():""));
					w.println("# "+col+" "+the_interval);
					w.print(s+" <- c(");
					w.print(samples.stream().
							filter(S->S.collection.equals(col)).
							map(R->String.valueOf(R.getMedianDepth(the_interval))).
							collect(Collectors.joining(",")));
					w.println(")");
					}
				}
			w.println("cols<-rainbow("+collections.size()+")");
			w.println("pdf("+quote("TMP/jeter.pdf")+",width=max(20,"+(collections.size()*intervals.size())+"*0.2),height=8)");
			w.print("boxplot(");
			w.print(String.join(",",components));
			if(proportional_width) {
				w.print(",width=c("+String.join(",",widths)+")");
				}
		
			w.print(",col=cols,"
					+ "las=2,ylab="+quote("Median Depth")+","
					+ "names=c("+String.join(",",names)+"),"
					+ "main="+quote(title)
					+ ")\\n");
			if(collections.size()>1) {
				w.print("legend("+quote("topleft")+",legend=c(");
				w.print(collections.stream().map(S->quote(S)).collect(Collectors.joining(",")));
				w.println("),fill=cols)");
				}
			
			w.println("dev.off()");
			w.flush();
			}
		}
	catch(Throwable err) {
		err.printStackTrace();
		System.exit(-1);
		}
	}
	
public static void main(final String args[])  {
	new Minikit().instanceMainWithExit(Arrays.asList(args));
	}
}

__EOF__


cat << EOF > TMP/tmp.mf
Manifest-Version: 1.0
Main-Class: Minikit
EOF


javac -cp \${JVARKIT_DIST}/coverageplotter.jar -d TMP -sourcepath TMP TMP/Minikit.java
jar cfm minikit.jar TMP/tmp.mf -C TMP .

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile minikit</entry>
	<entry key="javac.version">\$(javac -version 2>&1)</entry>
</properties>
EOF
"""
}

process PLOT {
tag "${transcript}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(bed)
	val(transcript)
	path(samples)
	path(minikit)
output:
	path("${meta.prefix?:""}${transcript}.pdf"),emit:pdf
	path("${meta.prefix?:""}${transcript}.R"),emit:R
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit R")}
set -o pipefail
mkdir TMP
awk -F '\t' '(\$4=="${transcript}")' "${bed}" |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	sort -T . -k1,1 -k2,2n > TMP/jeter.bed

test -s "TMP/jeter.bed"

java -cp ${minikit}:\${JVARKIT_DIST}/coverageplotter.jar Minikit "${samples}" TMP/jeter.bed > TMP/jeter.R

R --vanilla < TMP/jeter.R
mv -v TMP/jeter.pdf "${meta.prefix?:""}${transcript}.pdf"
mv -v TMP/jeter.R "${meta.prefix?:""}${transcript}.R"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
</properties>
EOF
"""
}



process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(zip)
output:
	path("*.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
ln -s ${zip} ./
"""
}

runOnComplete(workflow);

