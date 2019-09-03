from pymongo import MongoClient
import pymongo
import json
import pprint

client = MongoClient("labdb01.tgen.org", 25755)
annotateDb = client['AnnotateBy']
geneCollection = annotateDb['gene']

statsDb = client['statsTemp']
componentsCollection = statsDb['components']

db = client['markers']
tumorCollection = db['tumor']
study = "C024"
origin = "DNA"


def get_patient_id(study_patient_tissue):
    convert_plan_c024 = []
    with open('PLAN_C024.txt') as f:
        convert_plan_c024 = dict(x.rstrip().split(None, 1) for x in f)

    patient = study_patient_tissue[:-3]
    if patient in convert_plan_c024:
        return convert_plan_c024[patient]
    else:
        return study_patient_tissue


def get_cancer_census():
    cancer_census_genes = geneCollection.find({'cancerCensus': {"$exists": "true"}, 'cancerCensus.Tier': '1'},
                                              {"gene": 1})
    cc_genes = []
    for ccGene in cancer_census_genes:
        cc_genes.append(ccGene['gene'])
    return cc_genes


def print_biomarker(patient, variants):
    biomarkers = {
        "api_version": "2.0.1",
        "source": "TGen",
        "drug_rules_version": "Cosmic v82",
        "biomarkers": variants
    }
    biomarkers_json = json.dumps(biomarkers)
    patient_id = get_patient_id(patient)
    file = open("patients/" + patient_id + ".json", "w")
    length = len(variants)
    print("")
    print(str(length) + " " + patient_id)
    print("")
    for row in variants:
        print(row['gene'] + " " + row['effect'] + " " + row['alteration_type'])
    file.write(biomarkers_json)
    file.close()


def getVariants():
    ccGenes = get_cancer_census()
    # {'variants.filename': {'$regex':'snpeff.final.vcf'}}
    # {'variants.filename': {'$regex':'rna.final.seurat.vcf'}}
    # {'variants.filename': {'$regex': 'Freebayes.filt.snpEff.vcf'}}
    studyResult = tumorCollection.find({"variants.study": study, 'variants.filter': 'PASS', "gene": {"$in": ccGenes},
                                        '$and': [{'aberration.aberration_type2': {"$ne": "Synonymous"}},
                                                 {'aberration.aberration_type2': {"$ne": "Intron Variant"}},
                                                 {'aberration.aberration_type2': {"$ne": "UTR"}}],
                                        '$or': [{'variants.filename': {'$regex': 'snpeff.final.vcf'}},
                                                {'variants.filename': {'$regex': 'Freebayes.filt.snpEff.vcf'}},
                                                {'variants.filename': {'$regex': 'rna.final.seurat.vcf'}},
                                                {'variants.filename': {'$regex': 'cna.seg.vcf'}}
                                                ]
                                        }
                                       ).sort([("variants.studyPatientTissue", pymongo.ASCENDING)])
    return studyResult


def parseVariants(studyResult):
    variants = []

    current_patient = ""

    for result in studyResult:

        gene = result['gene']
        alteration_type = result['aberration']['aberration_type2']
        origin = "DNA"
        studyPatientTissue = result['variants'][0]['studyPatientTissue']

        effect = ""
        if not current_patient:
            print("changing current patient to : " + studyPatientTissue)
            current_patient = studyPatientTissue

        if 'AminoAcidChange' in result['snpeff']:
            effect = result['snpeff']['AminoAcidChange']
        else:
            effect = result['aberration']['aberration_value']

        if alteration_type.startswith("Focal Copy Number Gain"):
            effect = "Amplification"
        elif alteration_type.startswith("Focal Copy Number Loss"):
            effect = "Deletion"
        elif alteration_type.startswith("High Or Moderate Variant"):
            alteration_type = "Missense"
        elif alteration_type.startswith("Splic"):
            effect = alteration_type

        if studyPatientTissue != current_patient:
            # id = getPatientId(studyPatientTissue)
            print_biomarker(current_patient, variants)
            variants = []
            current_patient = studyPatientTissue


        variant = {"gene": gene, "alteration_type": alteration_type, "origin": origin, "effect": effect}
        variants.append(variant)


queryResult = getVariants()
parseVariants(queryResult)

client.close()
