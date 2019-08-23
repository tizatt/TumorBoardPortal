from pymongo import MongoClient
import pymongo
import pprint


client = MongoClient("labdb01.tgen.org", 25755)
annotateDb = client['AnnotateBy']
geneCollection = annotateDb['gene']

db = client['markers']
tumorCollection = db['tumor']
study = "C024"
origin = "DNA"


def getCancerCensus():
    cancerCensusGenes = geneCollection.find({'cancerCensus': {"$exists": "true"}},{"gene":1})
    ccGenes = []
    for ccGene in cancerCensusGenes:
        ccGenes.append(ccGene['gene'])
    return ccGenes


def printBiomarker( patient, variants ):
    biomarkers = {"source" : "TGen",
                  "drug_rules_version": "1.3",
                  "order": {"tumor_collection_date":"01/01/2000", "primary_diagnosis": "", "specimen_type" : "", "specimen_sit": ""},
                  "biomarkers" : variants}
    file = open(patient + ".json", "w")
    pp = pprint.PrettyPrinter(indent=4)
    with open( patient + ".json", "w") as jout:
        jout.write(pprint.pformat(vars(pprint)))
    pp.pprint( biomarkers)

def getVariants():
    ccGenes = getCancerCensus()
    studyResult = tumorCollection.find({"variants.study" : study, 'variants.filter' : 'PASS', "gene" : { "$in": ccGenes},
                                    'aberration.aberration_type2': {"$ne": "Synonymous"},
                                    "$or": [{'aberration.aberration_type': 'SNV'},
                                            {'aberration.aberration_type': 'FUSED_GENE'},
                                            {'aberration.aberration_type': 'CNV_LOSS'},
                                            {'aberration.aberration_type': 'CNV_GAIN'},
                                            {'aberration.aberration_type': 'INDEL'}
                                            ]
                                    }
                                   ).sort("variants.studyPatient")
    return studyResult

def parseVariants(studyResult):
    variants = []
    patients = []
    currentPatient = ""
    for result in studyResult:
        gene = result['gene']
        alteration_type = result['aberration']['aberration_type2']
        origin = "DNA"
        patient = result['variants'][0]['studyPatient']
        effect = result['aberration']['aberration_value']
        tissue = result['variants'][0]['tissue']

        variant = { "gene": gene, "alteration_type": alteration_type, "origin": origin, "effect" : effect}
        variants.append( variant )

        if patient != currentPatient:
            currentPatient = patient
            printBiomarker( patient, variants )
            #print( patient + " " + gene + " " + alteration_type + " " + origin + " " + effect + " " + tissue)



queryResult = getVariants()
parseVariants( queryResult )


client.close()