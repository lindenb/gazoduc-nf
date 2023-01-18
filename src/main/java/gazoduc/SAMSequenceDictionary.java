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
	public int size() {
		return array.size();
		}
	public long getReferenceLength() {
		return this.len;
		}
	public SAMSequenceRecord getSequence(int index) {
		return array.get(index);
		}
	public SAMSequenceRecord getSequence(final String name) {
		return id2ssr.get(name);
		}
	@Override
	public Iterator<SAMSequenceRecord> iterator() {
		return this.array.iterator();
		}
	@Override
	public String toString() {
		return this.array.toString();
		}
	}
