/* NFConfigParser.java */
/* Generated By:JavaCC: Do not edit this line. NFConfigParser.java */
package gazoduc.config;
import java.util.*;
import java.io.*;
import java.nio.file.*;
import javax.xml.*;
import javax.xml.stream.*;


public class NFConfigParser implements NFConfigParserConstants {
private static abstract class Param {
        Token doc=null;
        void print(XMLStreamWriter out,int depth,String title) throws XMLStreamException{
                out.writeStartElement("dt");
                out.writeCharacters(title);

                if(doc!=null) {
                        String c="";
                        Token t = doc;
                        while(t!=null) {
                                c+=t.image;
                                t=t.next;
                                }
                        if(c.startsWith("/*")) {
                                c=c.substring(1);
                                while(c.startsWith("*")) c=c.substring(1);
                                }
                        while(c.endsWith("*/")) {
                                c=c.substring(0,c.length()-1);
                                while(c.endsWith("*"))c=c.substring(0,c.length()-1);
                                }
                        c=c.trim();
                        if(!c.isEmpty()) {
                                out.writeCharacters(" ");
                                out.writeStartElement("q");
                                out.writeCharacters(c);
                                out.writeEndElement();//qiote
                                }
                        }
                else
                        {
                        }

                out.writeEndElement();//dt
                }
        }

private static class ObjectParam extends Param {
        Map<String,Param> map = new LinkedHashMap<>();
        @Override
        public String toString() {
                return comment(doc) +" "+map.toString();
                }
        @Override
        void print(XMLStreamWriter out,int depth,String title) throws XMLStreamException {

                super.print(out,depth,title);
                out.writeStartElement("dd");
                out.writeStartElement("dl");
                for(final String k: keySet()) {
                        if(map.get(k)==null) {
                                System.err.println("NULLL?? "+k);
                                continue;
                                }
                        map.get(k).print(out,depth+1,k);
                        }
                out.writeEndElement();//dl
                out.writeEndElement();//dd
                }

        Set<String> keySet() { return map.keySet();}

        void put(String t,Param p) {
                if(p instanceof ValueParam) {
                        this.map.put(t,p);
                        }
                else
                        {
                        ObjectParam x = (ObjectParam)this.map.get(t);
                        if(x==null) {
                                x=new ObjectParam();
                                this.map.put(t,x);
                                }

                        final ObjectParam o  = ObjectParam.class.cast(p);
                        for(final String t2:o.keySet()) {
                                x.put(t2,o.map.get(t2));
                                }
                        }
                }

        }

private static class ValueParam extends Param {
        final String value;
        ValueParam(final String s) {
                this.value = s;
                }
        @Override
        public String toString() {
                return value;
                }
        @Override  void print(XMLStreamWriter out,int depth,String title) throws XMLStreamException {
                super.print(out,depth,title);
                out.writeStartElement("dd");
                out.writeStartElement("pre");
                out.writeCharacters(title+" = "+ value);
                out.writeEndElement();
                out.writeEndElement();
                }
        }

private Path path;
private ObjectParam params;
NFConfigParser(final ObjectParam params, final Path path) throws IOException {
        this(Files.newInputStream(path));
        System.err.println("parsing "+path);
        this.params= params;
        this.path = path;
        }
static String unescape(final String s) {
        if(s.startsWith("'") && s.endsWith("'")) return s.substring(1,s.length()-1);
        if(s.startsWith("\"") && s.endsWith("\"")) return s.substring(1,s.length()-1);
        return s;
        }


private static String comment(Token t) {
        List<String> L= new ArrayList<>();
        while(t!=null) {
                L.add(t.image);
                t=t.next;
                }
        return String.join(" ",L);
        }

public static void main(final String[] args) {
        try {
                final ObjectParam root=new ObjectParam();
                for(final String f: args) {
                        System.err.println(f);
                                final NFConfigParser p = new NFConfigParser(root,Paths.get(f));
                        p.input();
                        System.err.println(f+"="+root);
                        }


                if(root!=null) {
                        final XMLOutputFactory xof = XMLOutputFactory.newFactory();
                        final XMLStreamWriter w = xof.createXMLStreamWriter(System.out,"UTF-8");
                        w.writeStartElement("html");
                        w.writeStartElement("head");
                        w.writeStartElement("style");
                        w.writeCharacters("pre {background:lightgray;}\n");
                        w.writeCharacters("dt {color:blue;font-style: bold;}\n");
                        w.writeCharacters("q {color:blue;font-style: italic; color:black;}\n");
                        w.writeEndElement();
                        w.writeEndElement();
                        w.writeStartElement("body");
                        w.writeStartElement("dl");
                        root.print(w,0,"params");
                        w.writeEndElement();
                        w.writeEndElement();
                        w.writeEndElement();
                        w.close();
                        }
                } catch(Throwable err) {
                err.printStackTrace();
                System.exit(-1);
                }
        }

  final public void input() throws ParseException {
    label_1:
    while (true) {
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case INCLUDE_CONFIG:
      case IDENTIFIER:{
        ;
        break;
        }
      default:
        jj_la1[0] = jj_gen;
        break label_1;
      }
      any();
    }
    jj_consume_token(0);
System.err.println("DONE "+path+" "+this.params.keySet());
}

  final private void any() throws ParseException {Map.Entry<String,ObjectParam> h;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case IDENTIFIER:{
      h = map();
System.err.println("XXXXXXX adding "+h+" to params");
                for(final String k: h.getValue().keySet()) {
                        this.params.put(k,h.getValue().map.get(k));
                }
      break;
      }
    case INCLUDE_CONFIG:{
      includeConfig();
      break;
      }
    default:
      jj_la1[1] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
}

  final private void includeConfig() throws ParseException {String f;
    jj_consume_token(INCLUDE_CONFIG);
    f = string();
Path p = this.path.getParent().resolve(f);
                try {
                        final NFConfigParser c = new NFConfigParser(this.params,p);
                        c.input();
                        System.err.println("done "+p+" "+this.params+" ############### "+c.params);
                        }
                catch(IOException err) {
                        err.printStackTrace();
                        System.exit(-1);
                        }
}

  final private Map.Entry<String,ObjectParam> map() throws ParseException {String title;Token t;Param v;final ObjectParam h=new ObjectParam();Map.Entry<String,ObjectParam> h2=null;List<String> chain = new ArrayList<>();
    t = jj_consume_token(IDENTIFIER);
title=t.image;h.doc = t.specialToken;
    jj_consume_token(OBRACKET);
    label_2:
    while (true) {
      switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
      case IDENTIFIER:{
        ;
        break;
        }
      default:
        jj_la1[2] = jj_gen;
        break label_2;
      }
      if (jj_2_1(2)) {
        h2 = map();
h.put(h2.getKey(),h2.getValue());
      } else {
        switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
        case IDENTIFIER:{
          t = jj_consume_token(IDENTIFIER);
          jj_consume_token(EQ);
          v = value();
v.doc = t.specialToken;
                                        h.put(t.image,v);
          break;
          }
        default:
          jj_la1[3] = jj_gen;
          jj_consume_token(-1);
          throw new ParseException();
        }
      }
    }
    jj_consume_token(CBRACKET);
{if ("" != null) return new AbstractMap.SimpleEntry<String,ObjectParam>(title,h);}
    throw new Error("Missing return statement in function");
}

  final private Param value() throws ParseException {String s; Token t; ObjectParam m;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case SINGLE_QUOTE_LITERAL:
    case DOUBLE_QUOTE_LITERAL:{
      s = string();
{if ("" != null) return new ValueParam(s);}
      break;
      }
    case BOOLEAN:{
      t = jj_consume_token(BOOLEAN);
{if ("" != null) return new ValueParam(t.image);}
      break;
      }
    case INT:{
      t = jj_consume_token(INT);
{if ("" != null) return new ValueParam(t.image);}
      break;
      }
    case DOUBLE:{
      t = jj_consume_token(DOUBLE);
{if ("" != null) return new ValueParam(t.image);}
      break;
      }
    case IDENTIFIER:{
      t = jj_consume_token(IDENTIFIER);
      label_3:
      while (true) {
        switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
        case DOT:{
          ;
          break;
          }
        default:
          jj_la1[4] = jj_gen;
          break label_3;
        }
        jj_consume_token(DOT);
        jj_consume_token(IDENTIFIER);
      }
{if ("" != null) return new ValueParam("${"+t.image+"}");}
      break;
      }
    default:
      jj_la1[5] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    throw new Error("Missing return statement in function");
}

  final private String string() throws ParseException {Token t;
    switch ((jj_ntk==-1)?jj_ntk_f():jj_ntk) {
    case SINGLE_QUOTE_LITERAL:{
      t = jj_consume_token(SINGLE_QUOTE_LITERAL);
{if ("" != null) return unescape(t.image);}
      break;
      }
    case DOUBLE_QUOTE_LITERAL:{
      t = jj_consume_token(DOUBLE_QUOTE_LITERAL);
{if ("" != null) return unescape(t.image);}
      break;
      }
    default:
      jj_la1[6] = jj_gen;
      jj_consume_token(-1);
      throw new ParseException();
    }
    throw new Error("Missing return statement in function");
}

  private boolean jj_2_1(int xla)
 {
    jj_la = xla; jj_lastpos = jj_scanpos = token;
    try { return (!jj_3_1()); }
    catch(LookaheadSuccess ls) { return true; }
    finally { jj_save(0, xla); }
  }

  private boolean jj_3R_map_251_9_4()
 {
    if (jj_scan_token(IDENTIFIER)) return true;
    if (jj_scan_token(OBRACKET)) return true;
    return false;
  }

  private boolean jj_3_1()
 {
    if (jj_3R_map_251_9_4()) return true;
    return false;
  }

  /** Generated Token Manager. */
  public NFConfigParserTokenManager token_source;
  SimpleCharStream jj_input_stream;
  /** Current token. */
  public Token token;
  /** Next token. */
  public Token jj_nt;
  private int jj_ntk;
  private Token jj_scanpos, jj_lastpos;
  private int jj_la;
  private int jj_gen;
  final private int[] jj_la1 = new int[7];
  static private int[] jj_la1_0;
  static {
	   jj_la1_init_0();
	}
	private static void jj_la1_init_0() {
	   jj_la1_0 = new int[] {0x180000,0x180000,0x100000,0x100000,0x400,0x1320c0,0xc0,};
	}
  final private JJCalls[] jj_2_rtns = new JJCalls[1];
  private boolean jj_rescan = false;
  private int jj_gc = 0;

  /** Constructor with InputStream. */
  public NFConfigParser(java.io.InputStream stream) {
	  this(stream, null);
  }
  /** Constructor with InputStream and supplied encoding */
  public NFConfigParser(java.io.InputStream stream, String encoding) {
	 try { jj_input_stream = new SimpleCharStream(stream, encoding, 1, 1); } catch(java.io.UnsupportedEncodingException e) { throw new RuntimeException(e); }
	 token_source = new NFConfigParserTokenManager(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 7; i++) jj_la1[i] = -1;
	 for (int i = 0; i < jj_2_rtns.length; i++) jj_2_rtns[i] = new JJCalls();
  }

  /** Reinitialise. */
  public void ReInit(java.io.InputStream stream) {
	  ReInit(stream, null);
  }
  /** Reinitialise. */
  public void ReInit(java.io.InputStream stream, String encoding) {
	 try { jj_input_stream.ReInit(stream, encoding, 1, 1); } catch(java.io.UnsupportedEncodingException e) { throw new RuntimeException(e); }
	 token_source.ReInit(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 7; i++) jj_la1[i] = -1;
	 for (int i = 0; i < jj_2_rtns.length; i++) jj_2_rtns[i] = new JJCalls();
  }

  /** Constructor. */
  public NFConfigParser(java.io.Reader stream) {
	 jj_input_stream = new SimpleCharStream(stream, 1, 1);
	 token_source = new NFConfigParserTokenManager(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 7; i++) jj_la1[i] = -1;
	 for (int i = 0; i < jj_2_rtns.length; i++) jj_2_rtns[i] = new JJCalls();
  }

  /** Reinitialise. */
  public void ReInit(java.io.Reader stream) {
	if (jj_input_stream == null) {
	   jj_input_stream = new SimpleCharStream(stream, 1, 1);
	} else {
	   jj_input_stream.ReInit(stream, 1, 1);
	}
	if (token_source == null) {
 token_source = new NFConfigParserTokenManager(jj_input_stream);
	}

	 token_source.ReInit(jj_input_stream);
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 7; i++) jj_la1[i] = -1;
	 for (int i = 0; i < jj_2_rtns.length; i++) jj_2_rtns[i] = new JJCalls();
  }

  /** Constructor with generated Token Manager. */
  public NFConfigParser(NFConfigParserTokenManager tm) {
	 token_source = tm;
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 7; i++) jj_la1[i] = -1;
	 for (int i = 0; i < jj_2_rtns.length; i++) jj_2_rtns[i] = new JJCalls();
  }

  /** Reinitialise. */
  public void ReInit(NFConfigParserTokenManager tm) {
	 token_source = tm;
	 token = new Token();
	 jj_ntk = -1;
	 jj_gen = 0;
	 for (int i = 0; i < 7; i++) jj_la1[i] = -1;
	 for (int i = 0; i < jj_2_rtns.length; i++) jj_2_rtns[i] = new JJCalls();
  }

  private Token jj_consume_token(int kind) throws ParseException {
	 Token oldToken;
	 if ((oldToken = token).next != null) token = token.next;
	 else token = token.next = token_source.getNextToken();
	 jj_ntk = -1;
	 if (token.kind == kind) {
	   jj_gen++;
	   if (++jj_gc > 100) {
		 jj_gc = 0;
		 for (int i = 0; i < jj_2_rtns.length; i++) {
		   JJCalls c = jj_2_rtns[i];
		   while (c != null) {
			 if (c.gen < jj_gen) c.first = null;
			 c = c.next;
		   }
		 }
	   }
	   return token;
	 }
	 token = oldToken;
	 jj_kind = kind;
	 throw generateParseException();
  }

  @SuppressWarnings("serial")
  static private final class LookaheadSuccess extends java.lang.Error {
    @Override
    public Throwable fillInStackTrace() {
      return this;
    }
  }
  static private final LookaheadSuccess jj_ls = new LookaheadSuccess();
  private boolean jj_scan_token(int kind) {
	 if (jj_scanpos == jj_lastpos) {
	   jj_la--;
	   if (jj_scanpos.next == null) {
		 jj_lastpos = jj_scanpos = jj_scanpos.next = token_source.getNextToken();
	   } else {
		 jj_lastpos = jj_scanpos = jj_scanpos.next;
	   }
	 } else {
	   jj_scanpos = jj_scanpos.next;
	 }
	 if (jj_rescan) {
	   int i = 0; Token tok = token;
	   while (tok != null && tok != jj_scanpos) { i++; tok = tok.next; }
	   if (tok != null) jj_add_error_token(kind, i);
	 }
	 if (jj_scanpos.kind != kind) return true;
	 if (jj_la == 0 && jj_scanpos == jj_lastpos) throw jj_ls;
	 return false;
  }


/** Get the next Token. */
  final public Token getNextToken() {
	 if (token.next != null) token = token.next;
	 else token = token.next = token_source.getNextToken();
	 jj_ntk = -1;
	 jj_gen++;
	 return token;
  }

/** Get the specific Token. */
  final public Token getToken(int index) {
	 Token t = token;
	 for (int i = 0; i < index; i++) {
	   if (t.next != null) t = t.next;
	   else t = t.next = token_source.getNextToken();
	 }
	 return t;
  }

  private int jj_ntk_f() {
	 if ((jj_nt=token.next) == null)
	   return (jj_ntk = (token.next=token_source.getNextToken()).kind);
	 else
	   return (jj_ntk = jj_nt.kind);
  }

  private java.util.List<int[]> jj_expentries = new java.util.ArrayList<int[]>();
  private int[] jj_expentry;
  private int jj_kind = -1;
  private int[] jj_lasttokens = new int[100];
  private int jj_endpos;

  private void jj_add_error_token(int kind, int pos) {
	 if (pos >= 100) {
		return;
	 }

	 if (pos == jj_endpos + 1) {
	   jj_lasttokens[jj_endpos++] = kind;
	 } else if (jj_endpos != 0) {
	   jj_expentry = new int[jj_endpos];

	   for (int i = 0; i < jj_endpos; i++) {
		 jj_expentry[i] = jj_lasttokens[i];
	   }

	   for (int[] oldentry : jj_expentries) {
		 if (oldentry.length == jj_expentry.length) {
		   boolean isMatched = true;

		   for (int i = 0; i < jj_expentry.length; i++) {
			 if (oldentry[i] != jj_expentry[i]) {
			   isMatched = false;
			   break;
			 }

		   }
		   if (isMatched) {
			 jj_expentries.add(jj_expentry);
			 break;
		   }
		 }
	   }

	   if (pos != 0) {
		 jj_lasttokens[(jj_endpos = pos) - 1] = kind;
	   }
	 }
  }

  /** Generate ParseException. */
  public ParseException generateParseException() {
	 jj_expentries.clear();
	 boolean[] la1tokens = new boolean[21];
	 if (jj_kind >= 0) {
	   la1tokens[jj_kind] = true;
	   jj_kind = -1;
	 }
	 for (int i = 0; i < 7; i++) {
	   if (jj_la1[i] == jj_gen) {
		 for (int j = 0; j < 32; j++) {
		   if ((jj_la1_0[i] & (1<<j)) != 0) {
			 la1tokens[j] = true;
		   }
		 }
	   }
	 }
	 for (int i = 0; i < 21; i++) {
	   if (la1tokens[i]) {
		 jj_expentry = new int[1];
		 jj_expentry[0] = i;
		 jj_expentries.add(jj_expentry);
	   }
	 }
	 jj_endpos = 0;
	 jj_rescan_token();
	 jj_add_error_token(0, 0);
	 int[][] exptokseq = new int[jj_expentries.size()][];
	 for (int i = 0; i < jj_expentries.size(); i++) {
	   exptokseq[i] = jj_expentries.get(i);
	 }
	 return new ParseException(token, exptokseq, tokenImage);
  }

  private boolean trace_enabled;

/** Trace enabled. */
  final public boolean trace_enabled() {
	 return trace_enabled;
  }

  /** Enable tracing. */
  final public void enable_tracing() {
  }

  /** Disable tracing. */
  final public void disable_tracing() {
  }

  private void jj_rescan_token() {
	 jj_rescan = true;
	 for (int i = 0; i < 1; i++) {
	   try {
		 JJCalls p = jj_2_rtns[i];

		 do {
		   if (p.gen > jj_gen) {
			 jj_la = p.arg; jj_lastpos = jj_scanpos = p.first;
			 switch (i) {
			   case 0: jj_3_1(); break;
			 }
		   }
		   p = p.next;
		 } while (p != null);

		 } catch(LookaheadSuccess ls) { }
	 }
	 jj_rescan = false;
  }

  private void jj_save(int index, int xla) {
	 JJCalls p = jj_2_rtns[index];
	 while (p.gen > jj_gen) {
	   if (p.next == null) { p = p.next = new JJCalls(); break; }
	   p = p.next;
	 }

	 p.gen = jj_gen + xla - jj_la; 
	 p.first = token;
	 p.arg = xla;
  }

  static final class JJCalls {
	 int gen;
	 Token first;
	 int arg;
	 JJCalls next;
  }

}