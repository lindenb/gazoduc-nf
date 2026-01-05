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

include {getKeyValue} from '../utils/functions.nf'


process GRIDSS_SETUP_REFERENCE {
conda "${meta.conda}/gridss"
input:
	val(meta)
	path(reference)
output:
	path("SETUP/${reference.name}"),emit:preproc_reference
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
mkdir SETUP
ln -s ${reference.toRealPath()} \${PWD}/SETUP/

for X in amb ann bwt pac sa fai
do
	ln -s "${reference.toRealPath()}.\${X}" \${PWD}/SETUP/
done

ln -s "${file(reference.toRealPath()).getParent()}/${reference.getSimpleName()}.dict" \${PWD}/SETUP/


gridss \
  -r "SETUP/${reference.name}" \
  -s setupreference


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">preprocess reference for gridss</entry>
        <entry key="reference">${reference}</entry>
</properties>
EOF
"""
}
