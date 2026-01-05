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

/**

example ouptut:

Total number of sequences            :   1998      
Total number of nucleotides          :   139860    
Total bytes uncompressed             :   350836    
Total bytes compressed               :   42168     
Compression ratio                    :   8.320 x 
File compressed with Ora version     :   2.7.0
Checksum                             :   84859a5156942a55
Contains interleaved data            :   YES
Original file name R1                :   S1.R1.fq.gz
Original file name R2                :   S1.R2.fq.gz
Species name                         :   homo_sapiens
Reference checksum                   :   d1caec1cd374cdb2

**/


process ORAD_INFO {
label 'process_single'
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
        tuple val(meta1),path(oradir)
        tuple val(meta ),path(ora_file)
output:
        tuple val(meta),path("*.info.txt"),path(ora_file),emit:info
        path("versions.yml"),emit:versions
script:
        def args1 =  task.ext.args1?:""
        def prefix = task.ext.prefix?:removeCommonSuffixes(ora_file.name)
"""
mkdir -p TMP


${oradir}/orad \\
        ${args1} \\
        --ora-reference "${oradir}/oradata" \\
        --info  \\
        ${ora_file} > TMP/jeter.info

 cut -f1 -d ':' TMP/jeter.info | sed 's/[ ]*\$//' | tr " " "_"| paste -s -d ',' > TMP/jeter.a
 cut -f2 -d ':' TMP/jeter.info | sed 's/^[ ]*//'  | sed 's/[ ]*\$//' | paste -s -d ',' >> TMP/jeter.a

mv  TMP/jeter.a "${prefix}.info.txt"


cat <<-END_VERSIONS > versions.yml    
"${task.process}":
    orad: \$(${oradir}/orad --version | awk '(NR==1) {print \$3;}')
END_VERSIONS
"""

stub:
    def prefix = task.ext.prefix?:removeCommonSuffixes(ora_file.name)
"""
cat << __EOF__ > "${prefix}.info.txt"
Total_number_of_sequences,Total_number_of_nucleotides,Total_bytes_uncompressed,Total_bytes_compressed,Compression_ratio,File_compressed_with_Ora_version,Checksum,Contains_interleaved_data,Original_file_name_R1,Original_file_name_R2,Species_name,Reference_checksum
1998,139860,350836,42168,8.320 x,2.7.0,84859a5156942a55,YES,${meta.id}.R1.fq.gz,${meta.id}.R2.fq.gz,homo_sapiens,d1caec1cd374cdb2
__EOF__

touch versions.yml
"""
}
