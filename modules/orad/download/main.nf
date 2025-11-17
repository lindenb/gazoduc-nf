/*

Copyright (c) 2025 Pierre Lindenbaum

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
process DOWNLOAD_ORAD {
label 'process_single'
afterScript "rm -rf jeter.tar.gz"
input:
        val(meta)
output:
        tuple val(meta),path("orad.*.linux"),emit:oradir
        path("versions.yml"),emit:versions
script:
        def version = task.ext.version?:"2.7.0"
	def local_dir = task.ext.local_dir?:"NO_DIR"
        def url = task.ext.url?:"https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.${version}.linux.tar.gz";
"""

if test -d "${local_dir}"
then

   ln -s "${local_dir}" "orad.local.linux"

else

   curl -L -o jeter.tar.gz "${url}"
   tar xvfz jeter.tar.gz
   rm jeter.tar.gz

fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    orad: "${version}"
    local: "${local_dir}"
    url: "${url}"
END_VERSIONS

"""
stub:
""" 
touch orad.stub.linux
touch versions.yml
"""
}
