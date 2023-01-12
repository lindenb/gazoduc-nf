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
import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.net.MalformedURLException;
import java.util.function.*;
import java.util.logging.Logger;
import java.util.stream.*;

import javax.xml.parsers.*;
import org.w3c.dom.*;



public class Genomes {
        private static final Logger LOG = Logger.getLogger(Genomes.class.getSimpleName());
	private final File file;
	private final Map<String,Genome> id2genome = new HashMap<>();


	private static List<Element> elements(final Node root) {
		final List<Element> array = new ArrayList<>();
		if(root==null) return array;
		 for(Node n1 = root.getFirstChild();n1!=null;n1=n1.getNextSibling()) {
                        if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
	                final Element e1 = (Element)n1;
			array.add(e1);
			}
		return array;
		}


	public class Genome {
		final String id;
		final Element root;

		Genome(final Element root) {
			this.root = root;
			if(!root.hasAttribute("id")) throw new IllegalArgumentException("<genome> element is missing id attribute in "+ Genomes.this.file);
			this.id = root.getAttribute("id").trim();
			if(this.id.isEmpty()) throw new IllegalArgumentException("<genome/@id> is empty in "+ Genomes.this.file);
			}

		
		public String getFasta() {
			return getRequired("fasta");
			}

		public final String getReference() {
			return getFasta();
			}

		public String getId() {
			return this.id;
			}

		private String about() {
			return "<genome id=\""+getId()+"\"> in " + Genomes.this.file ;
			}

		@Override
		public String toString() {
			return this.getId();
			}
		}

	private Genomes(final File file) {
		this.file = file;
		}

	private Genomes loadXml() {
		final DocumentBuilderFactory f=DocumentBuilderFactory.newInstance();
		f.setCoalescing(true);
		f.setValidating(false);
		final DocumentBuilder docBuilder= f.newDocumentBuilder();
		final Document dom = docBuilder.parse(this.file);
		final Element root = dom.getDocumentElement();
		for(Node n1 = root.getFirstChild();n1!=null;n1=n1.getNextSibling()) {
			if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
			final Element e1 = (Element)n1;
			if(e1.getNodeName().equals("genome")) {
				final Genome g = new Genome(e1);
				if( id2genome.containsKey(g.getId()) throw new IllegalArgumentException("genome id="+g.getId()+" is present twice in "+this.file);
				id2genome.put(g.getId(),g);
				}
			}
		return this;
		}
	
	@Override
	public String toString() {
		return this.file.toString();
		}

	public static Genomes load(final String path) throws Exception {
		return new Genomes(path).loadXml();
		}

	}
