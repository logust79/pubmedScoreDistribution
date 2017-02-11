#!/usr/bin/env python
'''
draw distribution of pubmedscore on retinal related genes and else
'''
from __future__ import print_function, division
import sys
import json
sys.path.append('../phenopolis/offline_analysis/pubmedScore')
sys.path.append('../phenopolis/offline_analysis/commons')
sys.path.append('../BioTools')
# using BioTools/Genes for symbol translation
from Genes import *
from phenopolis_utils import *
from pubmedScore import *
import time
import sqlite3
import random

sqlitedb = sqlite3.connect('irdc.db')
retnetJSON = json.load(open('../phenopolis/offline_analysis/gene_list/retnet.json','r'))
class pubmedTest():
    def __init__(self,keywords=None,output='report'):
        keywords = keywords or "retina,retinal,retinitis,blindness,macula,macular,stargardt,pigmentosa,amaurosis"
        self.keywords = keywords.split(',')
        self.now = time.mktime(time.localtime())
        self.dbs = get_mongo_collections()
        self.output = output
        #self.life = 2592000 # simply search everything, without using cached result
        self.G = Genes(sqlitedb)
    
    @property
    def retnet_genes(self):
        if getattr(self, '_retnet', None) is None:
            myd = self.G.symbols_to_ensemblIds(retnetJSON.keys())
            for k in retnetJSON.keys():
                # change symbols to ensembl ids
                if k not in myd: continue
                retnetJSON[myd[k]] = retnetJSON.pop(k)
            mg = list(self.dbs['phenopolis_db'].genes.find({'gene_id':{'$in':retnetJSON.keys()}}))
            result = {}
            for k in retnetJSON:
                gene_names = [j['gene_name'] for j in mg if j['gene_id'] == k]
                if gene_names:
                    result[k] = gene_names[0]
            self._retnet = result
        return self._retnet

    def get_non_retina_genes(self,size=1000):
        # get all non retina gene names
        if getattr(self, '_non_retnet', None) is None:
            mg = self.dbs['phenopolis_db'].genes.aggregate([
                {'$match':{'gene_id':{'$nin':self.retnet_genes.keys()}}},
                {'$sample':{'size':size}}
            ])['result']
            self._non_retnet = {k['gene_id']:k['gene_name'] for k in mg}
        return self._non_retnet

    def pubmedbatch(self, retinal=True):
        if retinal:
            genes = self.retnet_genes
        else:
            genes = self.get_non_retina_genes()
        result = {}
        for k,v in genes.items():
            print(v)
            result[k] = pubmed(v, self.keywords, self.now)['score']
        return result #{k:pubmed(v, self.keywords, self.now)['total_score'] for k,v in genes.items()}

    def run(self):
            # output two files, one for retina, one for non-retina
            output = {
                'retinal':self.output+'_retinal',
                'non_retinal':self.output+'_nonretinal',
            }
            my_pubmed = {
                'retinal': self.pubmedbatch(retinal=True),
                'non_retinal': self.pubmedbatch(retinal=False),
            }
            header = ['gene_id','pubmed_score']
            for i in ['retinal', 'non_retinal']:
                print('processing %s' % i)
                with open(output[i],'w') as outf:
                    # write headder
                    outf.write('\t'.join(header)+'\n')
                    for k,v in my_pubmed[i].items():
                        outf.write('\t'.join([k,str(v)])+'\n')
                        
            print('job done')
            return(0)
    
    def test(self):
        print(len(self.retnet_genes))

if __name__ == '__main__':
    PUB = pubmedTest()
    #print(PUB.retnet)
    PUB.run()
