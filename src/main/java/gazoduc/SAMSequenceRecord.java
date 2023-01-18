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

/**
 * Small re-implementation of htsjdk.samtools.SAMSequenceRecord
 */
public class SAMSequenceRecord {
private final String name;
private final int length;
SAMSequenceRecord(final String name,int length) {
	this.name = name;
	this.length = length;
	}

@Override
public boolean equals(final Object obj) {
	if(obj==this) return true;
	if(obj==null || !(obj instanceof SAMSequenceRecord)) return false;
	final SAMSequenceRecord other = SAMSequenceRecord.class.cast(obj);
	return getName().equals(other.getName()) && this.getLength()==other.getLength();
	}

@Override
public int hashCode() {
	return name.hashCode()*31 + Integer.hashCode(this.length);
	}

public String getName() {
	return this.name;
	}
public int getLength() {
	return this.length;
	}

@Override
public String toString() {
	return getName()+"("+this.getLength()+"bp)";
	}

}
