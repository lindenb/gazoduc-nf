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
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.*;

import org.w3c.dom.*;
import org.w3c.dom.Node;



/**
 * One genome definition stored in a genome file
 *
 */
public class Genome {
	private final Genomes genomes;
	private final String id;
	private final Element root;
	private SAMSequenceDictionary dict = null;

	Genome(final Genomes genomes,final Element root) {
		this.root = root;
		this.genomes = genomes;
		if(!root.hasAttribute("id")) throw new IllegalArgumentException("<genome> element is missing id attribute in "+ genomes.getFile());
		this.id = root.getAttribute("id").trim();
		if(this.id.isEmpty()) throw new IllegalArgumentException("<genome/@id> is empty in "+ genomes.getFile());
		}

	private List<Element> elements() {
		final List<Element> array = new ArrayList<>();
		 for(Node n1 = this.root.getFirstChild();n1!=null;n1=n1.getNextSibling()) {
            if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
			array.add(Element.class.cast(n1));
			}
		return array;
		}
	
	public final SAMSequenceDictionary getDict() {
		return getDictionary();
		}

	
	public SAMSequenceDictionary getDictionary() {
		if(this.dict==null) {
			final Path fai = Paths.get(getFasta()+".fai");
			try(BufferedReader br = Files.newBufferedReader(fai, StandardCharsets.UTF_8)) {
				this.dict = new SAMSequenceDictionary(
						br.lines().
						map(line->line.split("[\t]")).
						map(tokens->new SAMSequenceRecord(
							tokens[0],
							Integer.parseInt(tokens[1])
							)).
						collect(Collectors.toList())
						);
				}
			catch(IOException err) {
				throw new IllegalArgumentException("Cannot get dictionary for "+about(),err);
				}
			}
		return this.dict;
		}
	
	private Optional<String> get(final String nodeName) {
		final List<String> L = elements().
				stream().
				filter(N->N.getNodeName().equals(nodeName) && N.hasChildNodes()).
				map(E->E.getTextContent()).
				collect(Collectors.toList());
		if(L.isEmpty()) return Optional.empty();
		if(L.size()!=1) throw new IllegalArgumentException("Found "+L.size()+" nodes <"+nodeName+"> , but expected not more than one for "+ about());
		return Optional.of(L.get(0));
		}
	
	public boolean hasProperty(final String nodeName) {
		return elements().
				stream().
				anyMatch(N->N.getNodeName().equals(nodeName) && N.hasChildNodes())
				;
		}
	
	public String getProperty(final String nodeName,final String defaultValue) {
		return get(nodeName).orElse(defaultValue);
		}
	public String getRequiredProperty(final String nodeName) {
		return get(nodeName).orElseThrow(()->new IllegalArgumentException("Cannot find <"+nodeName+"> for "+ about()));
		}
	
	public Set<String> getAliases() {
		final Set<String> set = new HashSet<>();
		set.add(getId());
		elements().
				stream().
				filter(N->N.getNodeName().equals("alias") && N.hasChildNodes()).
				map(E->E.getTextContent().trim()).
				filter(S->!S.isEmpty()).
				forEach(S->set.add(S))
				;
		return set;
		}
	
	public Genomes getGenomes() {
		return this.genomes;
		}

	public final String getReference() {
		return get("fasta").orElseThrow(()->new IllegalArgumentException("Cannot find <fasta> for "+ about()));
		}
	
	public final String getFasta() {
		return getReference();
		}
	public final String getGff() {
		return get("gff").orElse(get("gff3").orElse(""));
		}
	public final String getGtf() {
		return get("gtf").orElse("");
		}

	public String getId() {
		return this.id;
		}

	public String getDescription() {
		return get("description").orElse(this.getId());
		}
	
	private String about() {
		return "<genome id=\""+getId()+"\"> in " + getGenomes().getFile() ;
		}
	
	/** test if dict looks like GRCh37  */
	public boolean isHg19() {
		return getDictionary().isHg19();
		}
	/** test if dict looks like GRCh38  */
	public  boolean isHg38() {
		return getDictionary().isHg38();
		}
	/** true if isHg19() || isHg38()  */
	public boolean isHuman() {
		return getDictionary().isHuman();
		}
	
	@Override
	public int hashCode() {
		return getId().hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj ==null || !(obj instanceof Genome )) return false;
		Genome other = Genome.class.cast(obj);
		return this.root == other.root;
		}
	
	@Override
	public String toString() {
		return this.getId();
		}
	}
