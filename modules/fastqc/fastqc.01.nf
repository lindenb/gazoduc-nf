/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {moduleLoad;isBlank} from '../utils/functions.nf'

String sampleOf(String s) {
	int i= s.lastIndexOf("/");
	if(i!=-1) s=s.substring(i+1);
	if(s.endsWith(".gz")) s=s.substring(0,s.length()-3);
	if(s.endsWith(".fq")) s=s.substring(0,s.length()-3);
	if(s.endsWith(".fastq")) s=s.substring(0,s.length()-6);
	return s;
	}

process APPLY_FASTQC_01 {
    tag "${row.fastq}"
    afterScript "rm -rf TMP"
    input:
	val(meta)
        val(row)
    output:
	tuple val(row), path("OUT/${sampleOf(row.fastq)}_fastqc.zip") ,path("OUT/${sampleOf(row.fastq)}_fastqc_data.txt"),emit:output
        path("version.xml"),emit:version
    script:
	def sample = sampleOf(row.fastq.toString())
      """
      hostname 1>&2
      ${moduleLoad("fastqc")}
      set -o pipefail
      mkdir -p TMP
      mkdir -p OUT

      fastqc --dir TMP -o OUT \
		--noextract \
		${row.adapters.name.equals("NO_FILE")?"":"--adapters "+ row.adapters} \
		${row.contaminants.name.equals("NO_FILE")?"":"--contaminants "+ row.contaminants} \
		--quiet -f "fastq" ${row.fastq}


	unzip -p OUT/${sample}_fastqc.zip '${sample}_fastqc/fastqc_data.txt' > 'OUT/${sample}_fastqc_data.txt'

##################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">fastqc</entry>	
	<entry key="basename">${sample}</entry>
	<entry key="fastq">${row.fastq}</entry>
	<entry key="contaminants">${row.contaminants}</entry>
	<entry key="adapters">${row.adapters}</entry>
	<entry key="fastqc.version">\$( fastqc --version)</entry>
</properties>
EOF
"""
}

