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
include { isBlank } from '../../../modules/utils/functions.nf'

workflow READ_SAMPLESHEET {
take:
	meta
	filename
main:
	versions = Channel.empty()

	if(meta.arg_name==null) {
		log.warn("READ_SAMPLESHEET meta.arg_name undefined");
		}
	if(filename==null) {
		throw new IllegalArgumentException("READ_SAMPLESHEET  : --${meta.arg_name?:"sampleheet"} is not defined");
		}
	else if(filename.endsWith(".json")) {
		ch1 = Channel.fromPath(filename).parseJson()
		}
	else if(filename.endsWith(".csv")) {
		ch1 = Channel.fromPath(filename).splitCsv(header:true,sep:',')
		}
	else if(filename.endsWith(".tsv")) {
		ch1 = Channel.fromPath(filename).splitCsv(header:true,sep:'\t')
		}
	else if(filename.endsWith(".ped") || filename.endsWith(".pedigree")) {
		ch1 = Channel.fromPath(filename)
			.splitText()
			.filter{!(it.startsWith("#") || it.trim().isEmpty())}
			.map{it.split("\\s+")}
			.filter{it[0]!="FID"}
			/*
			.map{
				if(it.size()< ) throw new IllegalArgumentException("Read pedigree not enough columnsin "${it}"  of --${meta.arg_name?:"sampleheet"} = ${filename}");
				def h=[
					family : it[0],
					sample : it[1],
					id : it[1]
					]
				if(!(isBlank(it[2]) ||it[2]=="0")) h=h.plus(father:it[2]);
				if(!(isBlank(it[3]) ||it[3]=="0")) h=h.plus(mother:it[3]);
				if(it.size()>4 && !isBlank(it[4])) h=h.plus(sex:it[4]);
				if(it.size()>5 && !isBlank(it[5])) h=h.plus(status:it[5]);
				return h;
			}*/

		}
	else if(filename.endsWith(".vcf.gz") || filename.endsWith(".bcf") || filename.endsWith(".vcf")) {
		ch1 = Channel.of([vcf:filename])
		}
	else if(filename.endsWith(".list")) {
		ch1 = Channel.fromPath(filename)
			.splitText()
			.map{it.trim()}
			.filter{!(it.startsWith("#") || it.trim().isEmpty())}
			.toList()
			.map{array->
				def label="";
				if(array.every{s->s.endsWith(".bam") || s.endsWith(".cram")}) {
					label="bam";
					}
				else if ( array.every{s->s.endsWith(".vcf") || s.endsWith(".vcf.gz") || s.endsWith(".bcf")} ) {
					label="vcf";
					}
				else if ( array.every{s->s.endsWith(".ora")} ) {
					label="ora";
					}
				else	
					{
					throw new IllegalArgumentException("Cannot infer the content of --${meta.arg_name?:"sampleheet"} = ${filename}");
					}
				def L=[];
				for(int i=0;i< array.size();i++) {
					L.add([(label as String):array.get(i)]);
					}
				return L;
				}
			.flatMap()
		}
	else
		{
		throw new IllegalArgumentException("Cannot interprete file suffix of --${meta.arg_name?:"sampleheet"} = ${filename} . Valid are json,csv,tsv.");
		}

	ch1.count().filter{it==0}.map{log.warn("No data (N=$it) was found for --${meta.arg_name?:"sampleheet"} = ${filename}"); return it;}

emit:
	versions
	samplesheet  = ch1
}
