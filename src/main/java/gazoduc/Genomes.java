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
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;


/**
 * A list of genomes stored as an XML file
 * @author lindenb
 *
 */
public class Genomes {
    private static final Logger LOG = Logger.getLogger(Genomes.class.getSimpleName());
	private final File file;
	private final Map<String,Genome> id2genome = new HashMap<>();

	/** @return the path to the  XML source */
	public File getFile() {
		return file;
		}
	
	public Set<String> getReferences() {
		return id2genome.values().stream().map(G->G.getReference()).collect(Collectors.toSet());
		}
	
	public Genome getById(final String id) {
		Genome g = this.id2genome.get(id);
		if(g==null) {
			g = this.id2genome.values().stream().filter(G->G.getAliases().contains(id)).findFirst().orElse(null);
			if(g!=null) {
				LOG.warning("Cannot find genome id="+id+" but found "+g+" search its aliases.");
				}
			}

		if(g==null) throw new IllegalArgumentException("Cannot find genome id "+id+" in "+ getFile()+". Available are :" + String.join(" ; ", this.id2genome.keySet()));
		return g;
		}
	public Genome getByReference(final String R) {
		final List<Genome> L=this.id2genome.values().stream().filter(G->G.getFasta().equals(R)).collect(Collectors.toList());
		if(L.isEmpty()) throw new IllegalArgumentException("Cannot find genome with fasta \""+R+"\" in "+ getFile());
		if(L.size()==1) return L.get(0);
		throw new IllegalArgumentException("Found multiple genome with fasta \""+R+"\" in "+ getFile()+" : "+ L);
		}

	private Genomes(final File file) {
		this.file = file;
		}

	private Genomes loadXml() {
		final DocumentBuilderFactory f=DocumentBuilderFactory.newInstance();
		f.setCoalescing(true);
		f.setValidating(false);
		try {
			final DocumentBuilder docBuilder= f.newDocumentBuilder();
			final Document dom = docBuilder.parse(getFile());
			final Element root = dom.getDocumentElement();
			if(root==null) throw new IOException("not root in "+getFile());
			if(!root.getNodeName().equals("genomes")) {
				throw new IllegalArgumentException("root is not <genomes/> for "+ getFile() );
				}
			for(Node n1 = root.getFirstChild();n1!=null;n1=n1.getNextSibling()) {
				if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
				final Element e1 = (Element)n1;
				if(e1.getNodeName().equals("genome")) {
					final Genome g = new Genome(this,e1);
					for(Genome g2: this.id2genome.values()) {
						if(g.getAliases().stream().anyMatch(A->g2.getAliases().contains(A))) {
							throw new IllegalArgumentException("genome id="+g.getId()+" is present twice in "+this.file);
							}
						}
					this.id2genome.put(g.getId(),g);
					}
				}
			}
		catch (ParserConfigurationException | SAXException | IOException e) {
			LOG.severe("Cannot read XML from "+getFile());
			throw new RuntimeException("Cannot load genomes XML :" + getFile(), e);
			}
		return this;
		}
	
	@Override
	public String toString() {
		return this.file.toString();
		}

	public static Genomes load(final String path) throws Exception {
		return new Genomes(new File(path)).loadXml();
		}

	}
