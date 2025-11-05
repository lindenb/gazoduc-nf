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
    };

List groupBy(List L0) {
    def id2data = [:];
    def used = [:];
    def i=0;
    for( i=0; i< L0.size();i++) {
        def item = L0[i];
        def sn = item[0].id;
        if(id2data.containsKey(sn)) throw new IllegalArgumentException("duplicate sample ${sn}");
        id2data.put(sn,item);
        used.put(sn,false);
        }
    def L=[];
    

    for( i=0; i< L0.size();i++) {
        def item = L0[i];
        def sn = item[0].id;
        if(used[sn]) continue;
        def father= (!isBlank(item[0].father) && id2data.containsKey(item[0].father)? id2data[item[0].father]:null);
        def mother= (!isBlank(item[0].mother) && id2data.containsKey(item[0].mother)? id2data[item[0].mother]:null);
        if(father!=null && mother!=null && !used[item.father] && !used[item.mother]) {
            def batch="trio_"+ sn;
            L.add(mapper(batch,item));
            L.add(mapper(batch,father));
            L.add(mapper(batch,mother));
            used.put(item.it,true);
            used.put(item.father,true);
            used.put(item.mother,true);
            }
        }

     for( i=0; i< L0.size();i++) {
        def item = L0[i];
        def sn = item[0].id;
        if(used[sn]) continue;
        def father= (!isBlank(item[0].father) && id2data.containsKey(item[0].father)? id2data[item[0].father]:null);
        if(father!=null  && !used[item.father]) {
            def batch="duo_"+ sn;
            L.add(mapper(batch,item));
            L.add(mapper(batch,father));
            used.put(item.it,true);
            used.put(item.father,true);
            }
        }

     for( i=0; i< L0.size();i++) {
        def item = L0[i];
        def sn = item[0].id;
        if(used[sn]) continue;
        def mother= (!isBlank(item[0].father) && id2data.containsKey(item[0].mother)? id2data[item[0].mother]:null);
        if(mother!=null  && !used[item.mother]) {
            def batch="duo_"+ sn;
            L.add(mapper(batch,item));
            L.add(mapper(batch,mother));
            used.put(item.it,true);
            used.put(item.mother,true);
            }
        }
        
     for( i=0; i< L0.size();i++) {
        def item = L0[i];
        def sn = item[0].id;
        if(used[sn]) continue;
        L.add(item);
        }
    if(L.size()!=L0.size()) throw new IllegalStateException("");
    return L;
    }

/** append meta.batch if it doesn't exists to trios */
workflow GROUP_BY_TRIOS {
    take:
        workflow_meta
        rows //[meta,any...] , first element must be a meta containing sample id/father/mother
    main:
        versions = Channel.empty()
        ch1 =  rows.branch{v->
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
