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

def gazoduc = gazoduc.Gazoduc.getInstance()

gazoduc.make("maxAlternateAlleles",6).
        description("Max number of GenotypeGVCF ALT alleles").
	setInt().
        put()


include {isBlank} from './../../modules/utils/functions.nf'

/**
  compute common options for gatk CombineGVCFs or GenotypeGVCF
*/
def gatkGetArguments(def meta) {
	def optDbsnp =  (isBlank(meta.dbsnp)?"":"--dbsnp \"${meta.dbsnp}\" ");
	return " -R \"${meta.reference}\" " + optDbsnp;
	}
	
/**
  compute common options for gatk CombineGVCFs
*/
def gatkGetArgumentsForCombineGVCFs(def meta) {
	return gatkGetArguments(meta) +" -G StandardAnnotation  -G AS_StandardAnnotation ";
	}

/**
  compute common options for gatk GenotypeGVCF
*/
def gatkGetArgumentsForGenotypeGVCF(def meta) {
	def maxAlternateAlleles = meta.maxAlternateAlleles?:6
	def pedigree = (meta.pedigree.name.equals("NO_FILE")?"":meta.pedigree.toRealPath())
	
	// BUG HAPLOID https://github.com/broadinstitute/gatk/issues/7304#issuecomment-1497966273
        def optPed = (isBlank(pedigree)?"":" -A PossibleDeNovo --pedigree \"${pedigree}\" ");
        
	return gatkGetArguments(meta) +" --max-alternate-alleles ${maxAlternateAlleles}  -G StandardAnnotation  -G AS_StandardAnnotation ${optPed} ";
	}
