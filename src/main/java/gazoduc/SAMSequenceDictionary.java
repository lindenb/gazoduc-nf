/*

Copyright (c) 2022 Pierre Lindenbaum

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
package gazoduc;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Small re-implementation of htsjdk.samtools.SAMSequenceDictionary
 */
public class SAMSequenceDictionary implements Iterable<SAMSequenceRecord> {
	private final Map<String,SAMSequenceRecord> id2ssr;
	private final List<SAMSequenceRecord> array;
	private final long len;
	SAMSequenceDictionary(final List<SAMSequenceRecord> L) {
		this.array = Collections.unmodifiableList(L);
		id2ssr = new HashMap<>(this.array.size());
		for(int i=0;i< this.array.size();i++) {
			final SAMSequenceRecord ssr = this.array.get(i);
			if(id2ssr.containsKey(ssr.getName())) throw new IllegalArgumentException("found sequence "+ssr+" twice in dictionary.");
			id2ssr.put(ssr.getName(),ssr);
			}
		this.len = this.array.stream().mapToLong(S->S.getLength()).sum();
		}
	
	/**
	@return a list of chromosome names for this dictionary
	*/
	public List<String> getChromosomes() {
		return this.array.stream().map(T->T.getName()).collect(Collectors.toList());
		}
	
	/**
	Alias of getChromosomes
	@return a list of chromosome names for this dictionary
	*/
	public final List<String> getContigs() {
		return getChromosomes();
		}

	/** @return the number of contigs in this dictionary */
	public int size() {
		return array.size();
		}
	/** @return the length of this genome */
	public long getReferenceLength() {
		return this.len;
		}
	/**
	@param index the index of the chromosome
	@return the index-th SAMSequenceRecord
	*/
	public SAMSequenceRecord getSequence(int index) {
		return array.get(index);
		}

	/**
	@param name the name of the chromosome
	@return SAMSequenceRecord for this chromosome name or null if not found
	*/
	public SAMSequenceRecord getSequence(final String name) {
		if(name==null) throw new IllegalArgumentException("name is null");
		return id2ssr.get(name);
		}
	@Override
	public Iterator<SAMSequenceRecord> iterator() {
		return this.array.iterator();
		}
	
	private boolean hasChromX(final String id,int expect) {
		SAMSequenceRecord rec = getSequence(id);
		if(rec==null)  rec = getSequence("chr"+id);
		if(rec!=null) {
			if(rec.getSequenceLength()==expect) return true;
			}
		return false;
		}
	
	/**
	@return an optional name for this genome for the UCSC database (hg19, hg38....)
	*/
	public Optional<String> getUcscName() {
		if(isHg19()) return Optional.of("hg19");
		if(isHg38()) return Optional.of("hg38");
		return Optional.empty();
		}
	/** @return an optional name for this genome */
	public Optional<String> getName() {
		return getUcscName();
		}
	
	/** @return true if dict looks like GRCh37  */
	public boolean isHg19() {
		return hasChromX("1",249_250_621);
		}
	/** @return true if dict looks like GRCh38  */
	public  boolean isHg38() {
		return hasChromX("1",248_956_422);
		}
	/** @return true if genome isHg19() or isHg38()  */
	public boolean isHuman() {
		return isHg19() || isHg38();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof SAMSequenceDictionary)) return false;
		final SAMSequenceDictionary o = SAMSequenceDictionary.class.cast(obj);
		return this.array.equals(o.array);
		}
	
	@Override
	public int hashCode() {
		return this.array.hashCode();
		}
	
	@Override
	public String toString() {
		return this.array.toString();
		}
	}
