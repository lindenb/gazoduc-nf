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
nextflow.enable.dsl=2

//include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { HC_COMBINE1 } from '../../../modules/gatk/gatk4.combine.gvcfs.01.nf'
include { HC_COMBINE2 } from '../../../modules/gatk/gatk4.combine.gvcfs.02.nf'
include { HC_GENOTYPE } from '../../../modules/gatk/gatk4.genotype.gvcfs.01.nf'
include { HC_GENOMICDB_IMPORT} from '../../../modules/gatk/gatk4.genomicdb.import.01.nf'
include { HC_GENOMICDB_GENOTYPE} from '../../../modules/gatk/gatk4.genomicdb.genotype.01.nf'
include { SPLIT_BED; SPLIT_BED as SPLIT_BED2 } from './part.split.bed.nf'
include {VCF_STATS                           } from '../../../subworkflows/vcfstats'
include {COMPILE_VERSIONS                    } from '../../../modules/versions/main.nf'
include {MULTIQC                             } from '../../../modules/multiqc'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
//validateParameters()

// Print summary of supplied parameters
//log.info paramsSummaryLog(workflow)


List makeSQRT(def L1) {
	def key = L1.get(0);
	def L = L1.get(1).sort();
	int n = (int)Math.ceil(Math.sqrt(L.size()));
	if(n<25) n=25;
	def returnList = [];
	def currList = [];
	int i=0;
	for(;;) {
		if(i<L.size()) currList.add(L.get(i));
		if(i==L.size() || currList.size()==n) {
			if(!currList.isEmpty()) returnList.add([key,currList]);
			if(i==L.size()) break;
			currList=[];
			}
		i++;
		}
	return returnList;
	}


workflow {

	beds_ch = Channel.fromPath(params.beds).
		splitText().
		map{file(it.trim())}

	bams_ch = Channel.empty()
	versions = Channel.empty()
	multiqc = Channel.empty()

	ref_hash = [
		id:   file(params.fasta).simpleName,
		name: file(params.fasta).simpleName,
		ucsc_name : params.ucsc_name?:"undefined"
		]

	fasta = [ref_hash,file(params.fasta)]
	fai   = [ref_hash,file(params.fai)]
	dict  = [ref_hash,file(params.dict)]
	def gtf    = [ref_hash, file(params.gtf), file(params.gtf+".tbi") ]
	def gff3    = [ref_hash, file(params.gff3), file(params.gff3+".tbi") ]


        if(params.dbsnp!=null)
                {
                dbsnp =     [ref_hash, [] ]
                dbsnp_tbi = [ref_hash, [] ]
                }
        else
                {
                dbsnp =     [ ref_hash, file(params.dbsnp)]
                dbsnp_tbi = [ ref_hash, file(params.dbsnp+".tbi") ]
                }

	hc_ch = Channel.empty()
        if(!params.gvcfs.equals("NO_FILE")) {
		gvcfs_ch = Channel.fromPath(params.gvcfs).
			splitText().
			map{it.trim()}.
			map{[it,it+".tbi"]}.
			map{[file(it[0],checkIfExists:true),file(it[1],checkIfExists:true)]}
		hc_ch = hc_ch.mix(  beds_ch.combine(gvcfs_ch) )
		}


	if(params.samplesheet.contains(".")) {
		bams_ch = Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')
			.map{if(!it.containsKey("bam")) throw new IllegalArgumentException("${it} : bam missing"); return it;}
			.map{if(!(it.bam.endsWith(".cram") || it.bam.endsWith(".bam"))) throw new IllegalArgumentException("${it}.bam should end with bam or cram"); return it;}
			.map{it.bai?it: (it.bam.endsWith(".bam") ? it.plus(["bai":it.bam+".bai"]):  it.plus(["bai":it.bam+".crai"]))}
                        .map{!it.fasta || it.fasta.isEmpty()?it.plus(["fasta":params.fasta]):it}
                        .map{!it.fai || it.fai.isEmpty()?it.plus(["fai":it.fasta+".fai"]):it}
                        .map{!it.dict || it.dict.isEmpty()?it.plus(["dict":it.fasta.replaceAll("\\.(fasta|fa)\$","")+".dict"]):it}
			.map{[
					[
					"id":it.id?:file(it.bam).simpleName,
					"sex":(it.sample?:"undefined"),
					"father":(it.father?:""),
					"mother":(it.mother?:"")
					],
					file(it.bam,checkIfExists:true),
					file(it.bai,checkIfExists:true),
					file(it.fasta,checkIfExists:true),
					file(it.fai,checkIfExists:true),
					file(it.dict,checkIfExists:true)
			]}
		}




	HC_BAM_BED(
		fasta,
		fai,
		dict,
		dbsnp,
		dbsnp_tbi,
		bams_ch.combine(beds_ch)
		)
	versions = versions.mix(HC_BAM_BED.out.versions )
		
	hc_ch = hc_ch.mix( HC_BAM_BED.out.gvcf)
		

    if(params.combine_method.equalsIgnoreCase("combinegvcf")) {
		combine1_input_ch = hc_ch.map{[it[0].toRealPath(),[it[1],it[2]]]}.
			groupTuple().
			flatMap{makeSQRT(it)}.
			map{[it[0],it[1].flatten()]}

		hc1 = HC_COMBINE1(fasta,fai,dict, dbsnp, dbsnp_tbi, combine1_input_ch )

		combine2_input_ch = hc1.output.
			map{ [it[0].toRealPath(), it[1],it[2] ]}. /* duplicate bed file so we can extract distinct contig in first , keep bed content in second */
			groupTuple(). /* group by contig/bed */
			map{[it[0], it[1].flatten().plus(it[2].flatten())]}

		splitbed_ch = SPLIT_BED(combine2_input_ch)

		combine2_input_ch = splitbed_ch.output.
			map{T->(T[0] instanceof List?T:[[T[0]],T[1]])}.//only one bed ? convert to array for the flatMap below
			flatMap{T->T[0].collect{X->[X,T[1]]} }

		hc2 = HC_COMBINE2(fasta,fai,dict, dbsnp, dbsnp_tbi, combine2_input_ch )
		hg = HC_GENOTYPE(fasta,fai,dict, dbsnp, dbsnp_tbi, hc2.output )
        }
	else if(params.combine_method.equalsIgnoreCase("genomicsDB") ||
			params.combine_method.equalsIgnoreCase("genomicsDBAndGenotype") ||
			params.combine_method.equalsIgnoreCase("glnexus")) {
		bed_gvcfs = hc_ch.
			map{[it[0].toRealPath(), it[1], it[2] ]}.
			groupTuple().
			map{[it[0], it[1].flatten().plus(it[2].flatten())]}


		beds2_ch = SPLIT_BED2(bed_gvcfs).output.
			flatMap{L1->{
			def key = (L1.get(0) instanceof List? L1.get(0):[L1.get(0)]);
			def L = []
			for(i=0;i< key.size();i++) {
				L.add([key.get(i),L1.get(1)])
				}
			return L;
			}}
		beds2_ch.view{"beds2_ch $it"}

		if( params.combine_method.equalsIgnoreCase("genomicsDB")) {
			genomiddb_ch = HC_GENOMICDB_IMPORT(fasta,fai,dict,beds2_ch)
			hg = HC_GENOMICDB_GENOTYPE(fasta,fai,dict, dbsnp, dbsnp_tbi, genomiddb_ch.output)
			}
		else if( params.combine_method.equalsIgnoreCase("genomicsDBAndGenotype")){
			hg = HC_GENOMICDB_IMPORT_AND_GENOTYPE(
				fasta,fai,dict, 
				dbsnp,
				dbsnp_tbi,
				beds2_ch
				)
			}
		else  if( params.combine_method.equalsIgnoreCase("glnexus"))  {
			hg = GLNEXUS(
				fasta,fai,dict, 
				beds2_ch
				)
			}
		else {
			throw new IllegalArgumentException("${params.combine_method}");
			}
	} else {
                throw new IllegalArgumentException("${params.combine_method}");
                }
	
		
	concat0_ch = VCF_CONCAT1(hg.output.collate(10).map{it.flatten()})
	concat_ch  = VCF_CONCAT2(concat0_ch.output.flatten().collect())

	if(params.with_multiqc==true) {
		VCF_STATS(
			ref_hash,
			fasta,
			fai,
			dict,
			Channel.of(gtf),
			Channel.of(gff3),
			[[id:"noped"],[]],
			[[id:"nogroup2sample"],[]],
			[[id:"nobed"],[]],
			concat_ch.output
			)
		versions = versions.mix(VCF_STATS.out.versions)
		multiqc = versions.mix(VCF_STATS.out.multiqc)


        COMPILE_VERSIONS(versions.collect())
        multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

		MULTIQC(multiqc.collect().map{[[id:"gatk4"],it]})
		}

	}


process HC_GENOMICDB_SAMPLE_MAP {
label "process_single"
afterScript "rm -rf TMP"
input:
        path("VCFS/*")
output:
        path("sample.*.map"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP
set -x

find VCFS -name "*.g.vcf.gz" | sort > TMP/jeter.list
test -s TMP/jeter.list

MD5=`cat TMP/jeter.list | sha1sum | cut -d ' ' -f1`

touch TMP/sample.map

cat TMP/jeter.list | while read F
do
        gunzip -c "\${F}" | grep "^#CHROM" -m1 | cut -f 10 | tr "\\n" "\t" >> TMP/sample.map
        echo "\${F}" >> TMP/sample.map
done

LC_ALL=C sort -t '\t' -k1,1 -T TMP TMP/sample.map > "sample.\${MD5}.map"
test -s "sample.\${MD5}.map"

cat << END_VERSIONS > versions.yml
"${task.process}":
	sort: todo
END_VERSIONS
"""
}


process HC_BAM_BED {
tag "${bed.name} ${bam.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_single"
afterScript "rm -rf TMP"
input:
	tuple val(meta1), path("_reference.fa")
	tuple val(meta2), path("_reference.fa.fai")
	tuple val(meta3), path("_reference.dict")
	tuple val(meta4), path(dbsnp)
	tuple val(meta5), path(dbsnp_tbi)
	tuple val(meta),path(bam),path(bai),path(fasta2),path(fai2),path(dict2),path(bed)
output:
	tuple path(bed),path("*g.vcf.gz"),path("*.g.vcf.gz.tbi"),emit:gvcf
	path("versions.yml"),emit:versions
script:
	def reference="_reference"
	def prefix = task.ext.prefix?:bam.getBaseName()
	def mapq = ((task.ext.mapq?:-1) as int)
"""
hostname 1>&2
mkdir -p TMP
set -x

${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" |\\
	cut -f1-3 |\\
	sort -t '\t' -T TMP -k1,1 -k2,2n|\\
	bedtools merge  > TMP/jeter.bed

if ! cmp "${reference}.fa.fai" "${fai2}"
then
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -f "${fasta2}" --column 1 --convert SKIP TMP/jeter.bed > TMP/jeter2.bed
	mv -v  TMP/jeter2.bed TMP/jeter.bed
fi

test -s TMP/jeter.bed

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \\
     -L TMP/jeter.bed \\
     -R "${fasta2}" \\
     -I "${bam}" \\
     -ERC GVCF \\
     ${mapq <1?"":" --minimum-mapping-quality ${mapq}"}\\
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \\
     -O "TMP/jeter.g.vcf.gz"

if ! cmp "${reference}.fa.fai" "${fai2}"
then
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfsetdict \\
		-n SKIP \\
		-R "${reference}.fa" \\
		TMP/jeter.g.vcf.gz > TMP/jeter2.vcf
	
	bcftools sort  -T TMP/sort  --max-mem "${task.memory.giga}G" -O z -o TMP/jeter.g.vcf.gz TMP/jeter2.vcf
	bcftools index --threads ${task.cpus} --force --tbi TMP/jeter.g.vcf.gz
fi

mv TMP/jeter.g.vcf.gz "${prefix}.g.vcf.gz"
mv TMP/jeter.g.vcf.gz.tbi "${prefix}.g.vcf.gz.tbi"


cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: todo
	bcftools: todo
END_VERSIONS
"""
}


process VCF_CONCAT1 {
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path("VCFS/*")
output:
	tuple path("*.bcf"),path("*.bcf.csi"),emit:output
script:
	
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP
set -x
	
find ./VCFS/ -name "*.vcf.gz" -o -name "*.bcf" > TMP/jeter.list

MD5=`cat TMP/jeter.list | md5sum | cut -d ' ' -f1`

bcftools concat --threads ${task.cpus} --allow-overlaps --rm-dups exact --file-list TMP/jeter.list -O b -o  TMP/output.\${MD5}.bcf
bcftools index --threads ${task.cpus} -f TMP/output.\${MD5}.bcf
mv -v TMP/output.\${MD5}.bcf ./
mv -v TMP/output.\${MD5}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""
}


process VCF_CONCAT2 {
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_single"
afterScript "rm -rf TMP"
input:
	path("VCFS/*")
output:
	tuple path("output.bcf"),path("output.bcf.csi"),emit:output
script:
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP
set -x
	
find ./VCFS -name "*.bcf" > TMP/jeter.list

if [[ \$(wc -l < TMP/jeter.list) -eq 1 ]]
then
	cp -v VCFS/*.bcf TMP/output.bcf
	cp -v VCFS/*.bcf.csi TMP/output.bcf.csi
else
	bcftools concat --threads ${task.cpus} --allow-overlaps --rm-dups exact --file-list TMP/jeter.list -O b -o  TMP/output.bcf
	bcftools index --threads ${task.cpus} -f TMP/output.bcf
fi

mv -v TMP/output.bcf ./
mv -v TMP/output.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""
}



process HC_GENOMICDB_IMPORT_AND_GENOTYPE {
tag "${bed.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
		tuple val(meta4),path(dbsnp)
		tuple val(meta5),path(dbsnp_tbi)
        tuple path(bed),path("VCFS/*")
output:
	tuple path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:output
script:
	def batchSize=-1
	def maxAlternateAlleles=6
	def prefix=task.ext.prefix?:bed.getBaseName()
"""
hostname 1>&2
mkdir -p TMP
set -x

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

cat TMP/jeter.list | while read F
do
        bcftools query -l "\${F}" | cut -f 10 | tr "\\n" "\t" >> TMP/sample.map
        echo "\${F}" >> TMP/sample.map
done


# sort on sample name
LC_ALL=C sort -t '\t' -k1,1 -T TMP TMP/sample.map > "TMP/jeter.map"
test -s TMP/jeter.map
mv TMP/jeter.map TMP/sample.map



SQRT=`awk 'END{X=NR;if(X<10){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' TMP/sample.map`

gatk --java-options "-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP" GenomicsDBImport \\
    -R ${fasta} \\
    --batch-size ${(batchSize as Integer) <= 0 ? "\${SQRT}" : ""+batchSize} \\
    ${(task.cpus as Integer) > 1 ? "  --reader-threads " +task.cpus : "" } \\
    --sample-name-map TMP/sample.map \
    -L "${bed}" \
    --genomicsdb-workspace-path "TMP/database"



gatk --java-options "-Xmx${task.memory.giga}g -XX:-UsePerfData -Djava.io.tmpdir=TMP" GenotypeGVCFs \\
        -R "${fasta}"  \\
        -V "gendb://TMP/database" \\
        -L "${bed}" \\
         ${dbsnp?"--dbsnp ${dbsnp}":""} \\
       --max-alternate-alleles ${maxAlternateAlleles} \\
       --seconds-between-progress-updates 60 \\
       -G StandardAnnotation -G StandardHCAnnotation \\
       --verbosity INFO \\
       -O "TMP/jeter.vcf.gz"

mv TMP/jeter.vcf.gz "${prefix}.vcf.gz"
mv TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"
rm -rf TMP
"""
}

process GLNEXUS {
    tag "${bed.name}"
     array 100
    label 'process_single'
	conda "${moduleDir}/../../../conda/glnexus.yml"
    afterScript "rm -rf GLnexus.DB TMP"
    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
        tuple path(bed),path("GVCFS/*")
    output:
        tuple  path("*.bcf"),path("*.bcf.csi")   , emit: output
    when:
        task.ext.when == null || task.ext.when
    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${bed.baseName}"
        def config = task.ext.config?:"gatk"
    """
    mkdir -p TMP

    glnexus_cli \\
        --threads ${task.cpus} \\
        --mem-gbytes ${task.memory.giga} \\
        --bed ${bed} \\
        --config ${config} \\
        ${args} \\
        GVCFS/*vcf.gz > TMP/jeter.bcf

    bcftools +fill-tags --threads ${task.cpus} -O b9 -o TMP/jeter2.bcf TMP/jeter.bcf -- -t  AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
    bcftools index --threads ${task.cpus}  --force TMP/jeter2.bcf

    mv TMP/jeter2.bcf ./${prefix}.bcf
    mv TMP/jeter2.bcf.csi ./${prefix}.bcf.csi

    """
}
