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
include {PREPARE_REFERENCE } from '../prepare.ref'

/** just like PREPARE_REFERENCE but for one fasta sequence only */
workflow PREPARE_ONE_REFERENCE {
take:
	workflow_meta
	fasta
main:
	PREPARE_REFERENCE(workflow_meta,fasta)
emit:
	versions = PREPARE_REFERENCE.out.versions
	fasta = PREPARE_REFERENCE.out.fasta.first()
	fai = PREPARE_REFERENCE.out.fai.first()
	dict = PREPARE_REFERENCE.out.dict.first()
	bed = PREPARE_REFERENCE.out.bed.first()
	scatter_bed = PREPARE_REFERENCE.out.scatter_bed.first()
	complement_bed = PREPARE_REFERENCE.out.complement_bed.first()
}
