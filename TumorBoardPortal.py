from pymongo import MongoClient
import pymongo
import json
import pprint

client = MongoClient("labdb01.tgen.org", 25755)
annotateDb = client['AnnotateBy']
geneCollection = annotateDb['gene']
out_dir = "/ngd-data/prodCentralDB/C024/"
statsDb = client['statsTemp']
componentsCollection = statsDb['components']

db = client['markers']
tumorCollection = db['tumor']
study = "C024"
origin = "DNA"


# disease = visit.diagnosis
#

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


def get_clinical_info( project_run ):
    

def print_biomarker(patient, project_run, variants):
    
    clinical_info = get_clinical_info( project_run )
    biomarkers = {
        "api_version": "2.0.1",
        "source": "TGen",
        "drug_rules_version": "Cosmic v82",
        "biomarkers": variants
    }
    biomarkers_json = json.dumps(biomarkers)
    patient_id = get_patient_id(patient)
    file = open(out_dir + "/patients/" + patient_id + ".json", "w")
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


def parseVariants(study_result):
    variants = []
    project_run = ""
    current_patient = ""
    origin = "DNA"

    for result in study_result:

        gene = result['gene']
        alteration_type = result['aberration']['aberration_type2']
        dna_allele_freq = ""
        project_run = result['variants'[0]['projectRun']
	study_patient_tissue = result['variants'][0]['studyPatientTissue']

        if 'SEURAT_AR_TUMOR' in result['variants'][0]:
            dna_allele_freq = result['variants'][0]['SEURAT_AR_TUMOR']
        elif 'AR2' in result['variants'][0]:
            dna_allele_freq = result['variants'][0]['AR2']

        effect = ""
        if not current_patient:
            print("changing current patient to : " + study_patient_tissue)
            current_patient = study_patient_tissue

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

        # Switching patients when we get to new patient
        if study_patient_tissue != current_patient:
            # id = getPatientId(studyPatientTissue)
            print_biomarker(current_patient, variants)
            variants = []
            current_patient = study_patient_tissue
    
        variant = {"gene": gene, "alteration_type": alteration_type, "origin": origin, "effect": effect, "dna_allele_frequency": dna_allele_freq}
        variants.append(variant)

    print_biomarker("C024_0008_T1", project_run, variants)

queryResult = getVariants()
parseVariants(queryResult)

client.close()
