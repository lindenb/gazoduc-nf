/*

Copyright (c) 2025 Pierre Lindenbaum

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
include { isBlank     } from '../../../modules/utils/functions.nf'

List mapper(String id,List row){
    def L2=[];
    L2.add(row[0].plus(batch:id));
    L2.addAll(row.subList(1,row.size()));
    return L2;
    }

List groupBy(List L0) {
    def id2data = [:];
    def used = [:];
    L0.each{item->{
        if(!(item instanceof List)) throw new IllegalArgumentException("Not a LIST: ${item}}");
        if(isBlank(item[0].id))  throw new IllegalArgumentException("id missing in ${item}}");
        def sn = item[0].id;
        if(id2data[sn]!=null) throw new IllegalArgumentException("duplicate sample ${sn}");

        id2data[sn]=item;
        used[sn]=false;
        }}
    def L=[];

    L0.each{item->{
        def meta = item[0];
        def sn = meta.id;
        if(used[sn]==true) return;
        def father= (!isBlank(meta.father) && id2data[meta.father]!=null? id2data[meta.father]:null);
        def mother= (!isBlank(meta.mother) && id2data[meta.mother]!=null? id2data[meta.mother]:null);
        if(father!=null && mother!=null && used[father[0].id]==false && used[mother[0].id]==false) {
            def batch="trio_"+ sn;
            L.add(mapper(batch,item));
            L.add(mapper(batch,father));
            L.add(mapper(batch,mother));
            used[meta.id]=true;
            used[father[0].id]=true;
            used[mother[0].id]=true;
            }
        }}



    L0.each{item->{
        def meta = item[0];
        def sn = meta.id;
        if(used[sn]==true) return;
        def father= (!isBlank(meta.father) && id2data[meta.father]!=null? id2data[meta.father]:null);
        if(father!=null  && used[father[0].id]==false) {
            def batch="duo_"+ sn;
            L.add(mapper(batch,item));
            L.add(mapper(batch,father));
            used.put(meta.id,true);
            used.put(father[0].id,true);
            }
        }}

      L0.each{item->{
        def meta = item[0];
        def sn = meta.id;
        if(used[sn]==true) return;
        def mother= (!isBlank(meta.mother) && id2data[meta.mother]!=null? id2data[meta.mother]:null);
        if(mother!=null  && used[mother[0].id]==false) {
            def batch="duo_"+ sn;
            L.add(mapper(batch,item));
            L.add(mapper(batch,mother));
            used.put(meta.id,true);
            used.put(mother[0].id,true);
            }
        }}

     L0.each{item->{
        def sn = item[0].id;
        if(used[sn]==true) return;
        L.add(item);
        }}

    if(L.size()!=L0.size()) throw new IllegalStateException("GROUP_BY_TRIOS L.size()=${L.size()}!=L0.size()=${L0.size()}");
    return L;
    }

/** append meta.batch if it doesn't exists to trios */
workflow GROUP_BY_TRIOS {
    take:
        workflow_meta
        rows //[meta,any...] , first element must be a meta containing sample id/father/mother
    main:
        versions = Channel.empty()
        ch1 =  rows
                .map{
                        if(!(it instanceof List)) throw new IllegalStateException("GROUP_BY_TRIOS: input should be a list");
                        return it;
                }
		.branch{v->
	            	sporadic : isBlank(v[0].id) || !isBlank(v[0].batch)
        	    	family : true
            		}

        ch2 = ch1.family
            .collect(flat:false)
            .flatMap{groupBy(it)}

    emit:
        versions
        output = ch1.sporadic.mix(ch2)
}
