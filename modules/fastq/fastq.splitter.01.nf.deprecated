/*

Copyright (c) 2026 Pierre Lindenbaum

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


process COMPILE_FASTQ_SPLITTER {
input:
	val(meta)
output:
	path("splitfastq"),emit:executable
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools/1.15.1 gcc")} 

echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}" 1>&2

cat << EOF > jeter.c
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include "signal.h"
#include <zlib.h>
#include "htslib/kseq.h"

KSEQ_INIT(gzFile, gzread)

static void usage(FILE* out) {
fprintf(out,"Usage:\\n");
fprintf(out,"  splitfastq -n (num) -m (modulo) <fastq1> \\n");
fprintf(out,"  splitfastq -n (num) -m (modulo) <fastq1> <fastq2> \\n");
fprintf(out,"other options:\\n");
fprintf(out,"  -B (int) buffer-size [%d]\\n",BUFSIZ);
fprintf(out,"\\n");
}

int main(int argc,char** argv) {
    int opt;
    int bufsize = BUFSIZ;
    long nsplits=0;
    long modulo = -1;
    char* buffer=NULL;
    while ((opt = getopt(argc, argv, "hn:m:B:")) != -1) {
	    switch (opt) {
	    case 'h':
		    usage(stdout);
		    return 0;
	    case 'n':
		    nsplits =atol(optarg);
		    break;
	    case 'm':
		    modulo =atol(optarg);
		    break;
	    case 'B':
		    bufsize = atoi(optarg);
		    break;
	    case '?':
		    fprintf(stderr,"unknown option '%c'.\\n",(char)optopt);
		    return EXIT_FAILURE;
	    default: /* '?' */
		    fprintf(stderr,"unknown option\\n");
		    return EXIT_FAILURE;
	    }
	}
    if(bufsize<1) {
	fprintf(stderr,"Bad value for buffer-size : %d.\\n",bufsize );
	usage(stderr);
	return EXIT_FAILURE;
	}
    if(nsplits<=0)  {
	fprintf(stderr,"Bad value for nsplit : %ld.\\n",nsplits );
	usage(stderr);
	return EXIT_FAILURE;
        }
    if(modulo<0 || modulo>=nsplits)  {
	fprintf(stderr,"Bad value for modulo : 0<=%ld<%ld \\n",modulo,nsplits );
	usage(stderr);
	return EXIT_FAILURE;
        }
     buffer=(char*)malloc(sizeof(char)*bufsize);
     if(buffer==NULL) {
	fprintf(stderr,"Cannot alloc buffer size: %d.\\n",bufsize );
	return EXIT_FAILURE;
	}
     setbuf(stdout, buffer);


#define OPENFQ(FP,KS,FNAME) \\
	FP = gzopen(FNAME, "r");\\
	if(FP==NULL) { fprintf(stderr,"Cannot open %s.(%s)\\n",FNAME,strerror(errno)); exit(EXIT_FAILURE);}\\
	KS = kseq_init(FP);\\
	if(KS==NULL) { fprintf(stderr,"Cannot initialize reader for \\"%s\\".\\n",FNAME); exit(EXIT_FAILURE);}

#define NEXT nReads++; if(nReads%nsplits!=modulo) continue;


#define KSWRITE(X) fwrite((void*)(X.s),sizeof(char),X.l,stdout)

#define WRITEFQ(ks) \\
	fputc('@',stdout);\\
	KSWRITE(ks->name);\\
	fputc('\\n',stdout);\\
	KSWRITE(ks->seq);\\
	fputc('\\n',stdout);\\
	fputc('+',stdout);\\
	fputc('\\n',stdout);\\
	KSWRITE(ks->qual);\\
	if(fputc('\\n',stdout)==EOF || ferror(stdout)) {fprintf(stderr,"[splitfastq]I/O error\\n");exit(EXIT_FAILURE);}


#define CLOSEFQ(FP,KS)  kseq_destroy(KS);gzclose(FP);

    long nReads = 0L;
    long nWrite = 0L;
    if(optind+1==argc) {
		gzFile fp1;
		kseq_t* ks1;
		OPENFQ(fp1,ks1,argv[optind])
		while (kseq_read(ks1) >= 0)  {
			NEXT
			WRITEFQ(ks1)
			nWrite++;
			}
		CLOSEFQ(fp1,ks1)
		fflush(stdout);
		fprintf(stderr,"[INFO] Read %ld, wrote %ld records.\\n",nReads,nWrite);
		}
    else if(optind+2==argc) {
		gzFile fp1;
		kseq_t* ks1;
		gzFile fp2;
		kseq_t* ks2;
		OPENFQ(fp1,ks1,argv[optind  ])
		OPENFQ(fp2,ks2,argv[optind+1])
		while (kseq_read(ks1) >= 0)  {
			if(kseq_read(ks2) <0 ) {
				fprintf(stderr,"Cannot more reads in  '%s' than in '%s'.\\n",argv[optind],argv[optind+1]);
				exit(EXIT_FAILURE);
				}
			NEXT
			WRITEFQ(ks1)
			WRITEFQ(ks2)
			nWrite++;
			}
		if(kseq_read(ks2) >=0 ) {
			fprintf(stderr,"Cannot more reads in  '%s' than in '%s'.\\n",argv[optind+1],argv[optind]);
			exit(EXIT_FAILURE);
			}
		
		CLOSEFQ(fp1,ks1)
		CLOSEFQ(fp2,ks2)
		fflush(stdout);
	   	fprintf(stderr,"[INFO] Read %ld, wrote %ld pairs.\\n",nReads,nWrite);
		}
    else
		{
		fprintf(stderr,"Illegal number of arguments.\\n");
		usage(stderr);
		return EXIT_FAILURE;
		}
     free(buffer);
     return EXIT_SUCCESS;
     }
EOF

## FIX THIS FULL PATH!!
g++ -o splitfastq -Wall -L /CONDAS/users/ifb-tools/miniconda3/envs/samtools-1.15.1/lib -O3 jeter.c -lz -lhts

rm jeter.c

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Compile splitfastq</entry>
	<entry key="version">${getVersionCmd("gcc samtools")}</entry>
</properties>
EOF
"""
}
