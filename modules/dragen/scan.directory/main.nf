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

process DRAGEN_SCAN_DIRECTORY {
tag "${directory.name}"
label "process_single"
input:
    tuple val(meta),path(directory)
output:
    tuple val(meta),path("*.json"),emit:json
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def dir = directory.toRealPath().toString()
    def grep_cmd =  "grep -F '/germline_seq/' " 
"""
find "${dir}" -type f -name "dragenlicenceUsageAndCoverage.txt"  > jeter.txt

if test -s jeter.txt
then

find "${dir}" -type f -name "*.cram" |\\
    awk  -F '/' '{printf("{\\"type\\":\\"bam\\",\\"id\\":\\"%s\\",\\"bam\\":\\"%s\\",\\"bai\\":\\"%s.crai\\"}\\n",\$(NF - 1),\$0,\$0);}' > lines.txt

find "${dir}" -type f -name "*.bam" ||\\
    awk -F '/' '{printf("{\\"type\\":\\"bam\\",\\"id\\":\\"%s\\",\\"bam\\":\\"%s\\",\\"bai\\":\\"%s.bai\\"}\\n",\$(NF - 1),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.cnv_sv.vcf.gz" |\\
    awk -F '/' '{printf("{\\"type\\":\\"cnv_sv\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF - 1),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.sv.vcf.gz" |\\
    awk -F '/' '{printf("{\\"type\\":\\"sv\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF- 1 ),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.hard-filtered.vcf.gz" |\\
    awk -F '/' '{printf("{\\"type\\":\\"vcf\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF - 1 ),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.ploidy.vcf.gz" |\\
    awk -F '/' '{printf("{\\"type\\":\\"ploidy\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF- 1 ),\$0,\$0);}' >> lines.txt

else


find "${dir}" -type f -name "*.cram" | ${grep_cmd} | \\
    awk  -F '/' '{printf("{\\"type\\":\\"bam\\",\\"id\\":\\"%s\\",\\"bam\\":\\"%s\\",\\"bai\\":\\"%s.crai\\"}\\n",\$(NF - 2),\$0,\$0);}' > lines.txt

find "${dir}" -type f -name "*.bam" | ${grep_cmd} |\\
    awk -F '/' '{printf("{\\"type\\":\\"bam\\",\\"id\\":\\"%s\\",\\"bam\\":\\"%s\\",\\"bai\\":\\"%s.bai\\"}\\n",\$(NF - 2),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.cnv_sv.vcf.gz" | ${grep_cmd} |\\
    awk -F '/' '{printf("{\\"type\\":\\"cnv_sv\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF - 2),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.sv.vcf.gz" | ${grep_cmd} |\\
    awk -F '/' '{printf("{\\"type\\":\\"sv\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF - 2 ),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.hard-filtered.vcf.gz" | ${grep_cmd} |\\
    awk -F '/' '{printf("{\\"type\\":\\"vcf\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF - 2 ),\$0,\$0);}' >> lines.txt

find "${dir}" -type f -name "*.ploidy.vcf.gz" | ${grep_cmd} |\\
    awk -F '/' '{printf("{\\"type\\":\\"ploidy\\",\\"id\\":\\"%s\\",\\"vcf\\":\\"%s\\",\\"tbi\\":\\"%s.tbi\\"}\\n",\$(NF - 2 ),\$0,\$0);}' >> lines.txt


fi

echo '[' > "${prefix}.json"
cat lines.txt | paste -s -d ',' >> "${prefix}.json"
echo ']' >> "${prefix}.json"


rm jeter.txt
touch versions.yml
"""
stub:
    def prefix = task.ext.prefix?:""
"""
touch versions.yml
echo '[]' >  "${prefix}.json"
"""
}
