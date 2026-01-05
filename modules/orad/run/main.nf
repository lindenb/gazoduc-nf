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
include {removeCommonSuffixes  } from '../../utils/functions.nf'
include {isBlank               } from '../../utils/functions.nf'
include {verify                } from '../../utils/functions.nf'

process RUN_ORAD {
label 'process_single'
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
        tuple val(meta1),path(oradir)
        tuple val(meta ),path(fastq_files)
output:
        tuple val(meta),path("*q.gz",arity:'1..*'),emit:fastq
        path("versions.yml"),emit:versions
script:
        def args1 =  task.ext.args1?:""
        def prefix = task.ext.prefix?:"${meta.prefix?:""}" //yes, prefix and NOT id
        def filenames = (fastq_files instanceof List?fastq_files:[fastq_files])
        if(isBlank(prefix)) {
                if((fastq_files instanceof List ) && fastq_files.size()==1) {
                        def fastq = fastq_files[0];
                        filenames = [fastq];
                        //prefix not used
                        prefix= removeCommonSuffixes(fastq.name);
                        }
                else if(fastq_files instanceof Path) {
                        def fastq = fastq_files;
                        filenames = [fastq];
                        //prefix not used
                        prefix= removeCommonSuffixes(fastq.name);
                        }
                else 
                        {
                        verify( (fastq_files instanceof List ) , "${task.process} expected a List of ora");
                        verify( !fastq_files.isEmpty() , "${task.process} expected a non empty List");
                        verify( !isBlank(prefix) , "${task.process}  prefix shouldn't be blank");
                        
                        filenames = fastq_files
                        }
                }
"""
mkdir -p TMP


${filenames.size()==1?"":" cat " + filenames.collect{it.name}.sort().join(" ")+" |\\"}
${oradir}/orad \\
        ${args1} \\
        --threads ${task.cpus} \\
        --path TMP \\
        --ora-reference "${oradir}/oradata" \\
        ${filenames.size()==1? filenames[0].name : " - "}

if test -f TMP/-R1.fastq.gz
then
        mv TMP/-R1.fastq.gz  "./${prefix}.R1.fastq.gz"
fi

if test -f TMP/-R2.fastq.gz
then
        mv TMP/-R2.fastq.gz  "./${prefix}.R2.fastq.gz"
fi


if test -f TMP/decomp_from_stdin.fastq.gz
then
	mv TMP/decomp_from_stdin.fastq.gz  "./${prefix}.fastq.gz"
fi

# check no file starts with R, case not handled for now
find TMP -name "-R*" > TMP/other.txt
cat TMP/other.txt 1>&2
test ! -s TMP/other.txt

mv -v TMP/*.gz ./ || true
    
cat <<-END_VERSIONS > versions.yml    
"${task.process}":
    orad: \$(${oradir}/orad --version | awk '(NR==1) {print \$3;}')
END_VERSIONS
"""

stub:
"""
touch ${meta.id}.R1.fq.gz
touch ${meta.id}.R2.fq.gz
touch versions.yml
"""
}
