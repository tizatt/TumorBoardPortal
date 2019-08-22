from pymongo import MongoClient
import pymongo
client = MongoClient("labdb01.tgen.org", 25755)
annotateDb = client['AnnotateBy']
geneCollection = annotateDb['gene']

db = client['markers']
tumorCollection = db['tumor']


def getCancerCensus():
    cancerCensusGenes = geneCollection.find({'cancerCensus': {"$exists": "true"}},{"gene":1})
    ccGenes = []
    for ccGene in cancerCensusGenes:
        ccGenes.append(ccGene['gene'])
    return ccGenes


ccGenes = getCancerCensus()

studyResult = tumorCollection.find({"variants.study": "C024", 'variants.filter': 'PASS',
                                    'aberration.aberration_type2': {"$ne": "Synonymous"},
                                    "$or": [{'aberration.aberration_type': 'SNV'},
                                            {'aberration.aberration_type': 'FUSED_GENE'},
                                            {'aberration.aberration_type': 'CNV_LOSS'},
                                            {'aberration.aberration_type': 'CNV_GAIN'},
                                            {'aberration.aberration_type': 'INDEL'}
                                            ]
                                    }
                                   )

for result in studyResult:
    variants = result['variants']
    gene = result['gene']
    alteration_type = result['aberration']['aberration_type2']
    origin = "DNA"
    patient = result['variants'][0]['studyPatient']
    effect = result['aberration']['aberration_value']
    tissue = result['variants'][0]['tissue']

    if gene in ccGenes:
        print str(patient) + " " + str(gene) + " " + str(alteration_type) + " " + origin + " " + str(effect) + " " + str(tissue)

db.close()
annotateDb.close()