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
include {moduleLoad;isBlank;getVersionCmd} from '../utils/functions.nf'

/* split input fastq into a set of N files */
process COMPILE_FASTQ_SPLIT2FILE {
input:
	val(meta)
output:
	path("split2file"),emit:executable
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools/1.15.1 gcc")} 

echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}" 1>&2

cat << EOF > jeter.cc
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <unistd.h>
#include <cassert>
#include <getopt.h>
#include <stdint.h>
#include <signal.h>
#include <zlib.h>
#include <string>
#include <vector>
#include "htslib/kseq.h"

using namespace std;

KSEQ_INIT(gzFile, gzread)

class SplitFile {
	public:
		string filename;
		gzFile out;
		long count;
	SplitFile(const char* prefix,int side,int id,level):filename(prefix),out(NULL),count(0L) {
	        char mode[5];
		char tmp[30];
		FILE* exists=NULL;
		if(side==1 || side==2) {
			filename.append(side==1?".R1":".R2");
			}		

		sprintf(tmp,".%05d.fastq.gz",(id+1));
		filename.append(tmp);
		
		exists = fopen(filename.c_str(),"rb");
		if(exists!=NULL) {
			fclose(exists);
			fprintf(stderr,"File already exists %s.\\n", filename.c_str() );
			exit(EXIT_FAILURE);
			}

		sprintf(mode,"wb%d",level);
		this->.out = gzopen(this->filename.c_str(),mode);
		if(out==NULL) {
			fprintf(stderr,"Cannot open %s for writing.(%s)\\n", this->filename.c_str(), strerror(errno));
			return EXIT_FAILURE;
			}

		}
	void close() {
		if(this->OUT==NULL) return;
		fprintf(stderr,"[LOG] closing \\"%s\\" N=%ld.\\n",this->filename.c_str(),this->count);
		gzclose(this->out);
		this->out = NULL;
		}
	~SplitFile() {
		close();
		}

#define GZWRITE(X) gzwrite(this->out,(void*)(X.s),X.l)

	void write(kseq_t* ks1) {
		gzputc(this->out,'@');
		GZWRITE(ks1->name);
		gzputc(this->out,'\\n');
		GZWRITE(ks1->seq);
		gzputc(this->out,'\\n');
		gzputc(this->out,'+');
		gzputc(this->out,'\\n');
		GZWRITE(ks1->qual);
		if(gzputc(this->out,'\\n')==-1) {
			fprintf(stderr,"[split2file]I/O error %s\\n",filename->c_str());
			exit(EXIT_FAILURE);
			}
		this->count++;
		}
#undef GZWRITE

	};

class Reader {
public:
        gzFile fp1;
	kseq_t* ks;
	string filename;
	Reader(const char* fname):ks(NULL),filename(fname==NULL?"<STDIN>":fname) {
		fp1 =  (fname == NULL || strcmp(fname,"-")==0)
			? gzdopen(fileno(stdin), "r") 
			: gzopen(fname], "r");
		if(fp1==NULL) {
			fprintf(stderr,"Cannot open %s .(%s)\\n",filename.c_str(),strerror(errno));
			exit(EXIT_FAILURE);
			}
		ks = kseq_init(fp1);
		if(ks==NULL) {
			fprintf(stderr,"Cannot initialize reader for %s\\n",filename.c_str()); 
			exit(EXIT_FAILURE);
			}
		}
	~Reader() {
		if(ks!=NULL) kseq_destroy(ks);
                gzclose(fp1);
		}
	bool next() {
		return kseq_read(ks) >= 0;
		}
	};


class AbstractEnd {
	public:
	AbstractEnd() {
		}
	virtual ~AbstractEnd() {
		}
	};

class SingleEnd: public AbstractEnd {
	public:
		SplitFile* sf;
		SingleEnd(const char* prefix,int id,level) {
			sf = new SplitFile(prefix,-1,id,level);
			}
		virtual ~SingleEnd() {
			delete sf;
			}
	};


class PairedEnd: public AbstractEnd {
	public:
		SplitFile* sf1;
		SplitFile* sf2;
		PairedEnd(const char* prefix,int id,level) {
			sf1 = new SplitFile(prefix,1,id,level);
			sf2 = new SplitFile(prefix,2,id,level);
			}
		virtual ~PairedEnd() {
			delete sf1;
			delete sf2;
			}
	};

static void usage(FILE* out) {
fprintf(out,"split2file. Pierre Lindenbaum 2020.\\n");
fprintf(out,"Usage:\\n");
fprintf(out,"  split2file [-s for single end] -C (compression-level 0-9) -n (num)  -o <PREFIX> <fastq|stdin> \\n");
fprintf(out,"\\n");
}

int main(int argc,char** argv) {
    int opt;
    int i;
    long nsplits=0;
    char* prefix = NULL;
    int level=5;
    int single  = 0 ;
    vector<AbstractEnd*> splitFiles;
    long nReads = 0L;

    while ((opt = getopt(argc, argv, "sho:n:C:")) != -1) {
	    switch (opt) {
	    case 'h':
		    usage(stdout);
		    return 0;
	    case 's':
		    single = 1;
		    break;
            case 'C':
		    level = atoi(optarg);
		    break;
	     case 'o':
		    prefix = optarg;
		    break;
	    case 'n':
		    nsplits =atol(optarg);
		    break;
	    case '?':
		    fprintf(stderr,"unknown option '%c'.\\n",(char)optopt);
		    return EXIT_FAILURE;
	    default: /* '?' */
		    fprintf(stderr,"unknown option\\n");
		    return EXIT_FAILURE;
	    }
	}
    if(prefix==NULL || strlen(prefix)==0UL) {
	fprintf(stderr,"empty prefix.\\n" );
	usage(stderr);
	return EXIT_FAILURE;
	}
    if(nsplits<=0)  {
	fprintf(stderr,"Bad value for nsplit : %ld.\\n",nsplits );
	usage(stderr);
	return EXIT_FAILURE;
        }
    if(level<0 || level>9)  {
	fprintf(stderr,"Bad compression level: %d.\\n",level );
	usage(stderr);
	return EXIT_FAILURE;
        }

    if(!(argc==optind || optind+1==argc))
	{
	fprintf(stderr,"Illegal number of arguments.\\n");
	usage(stderr);
	return EXIT_FAILURE;
	}

    if(single==0 && !(argc==optind || optind+2==argc)) {
	fprintf(stderr,"Illegal number of arguments (paired).\\n");
	usage(stderr);
	return EXIT_FAILURE;
	}


    for(i=0;i< nsplits;i++) {
	splitFiles.push_back(
		single==0
		? new PairedEnd(prefix,(i+1),level)
		: new SingleEnd(prefix,(i+1),level)
		);
	}
	

	if(single==0) /* paired */ {
		Reader* r1;
		Reader* r2;
		if(optind==argc) {
			r1 = new Reader(NULL);
			r2 = r1;
			}
		else
			{
			r1 = new Reader(argv[optind  ]);
			r2 = new Reader(argv[optind+1]);
			}
		while (r1->next())  {
			PairedEnd* out  = (PairedEnd*)splitFiles[nReads%nsplits];
			out->sf1->write(r1->ks);
			if(!r2->next()) {
				fprintf(stderr,"Paired end missing.\\n");
                                return EXIT_FAILURE;
				}
			out->sf2->write(r2->ks);
			nReads++;
			}
		if(r2->next()) {
			fprintf(stderr,"More reads in R2.\\n");
                        return EXIT_FAILURE;
			}

		delete r1;
		if(optind!=argc) {
			delete r2;
			}
		}
	else
		{
		Reader* r1 = new Reader(optind==argc?NULL:argv[optind]);
		while (r1->next())  {
			SingleEnd* out  = (SingleEnd*)splitFiles[nReads%nsplits];
			out->sf->write(r1->ks);
			nReads++;			
			}
		delete r1;
		}

	    

     
    for(i=0;i< nsplits;i++) {
	delete splitFiles[i];
	}

     return EXIT_SUCCESS;
     }
EOF

## FIX THIS FULL PATH!!
g++ -o split2file -Wall -L /CONDAS/users/ifb-tools/miniconda3/envs/samtools-1.15.1/lib -O3 jeter.cc -lz -lhts

rm jeter.cc

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Compile split2file</entry>
	<entry key="version">${getVersionCmd("gcc samtools")}</entry>
</properties>
EOF
"""
}
