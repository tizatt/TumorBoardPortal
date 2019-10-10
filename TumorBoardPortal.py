from pymongo import MongoClient
import pymongo
import json
import pprint

# Database setup
client = MongoClient("labdb01.tgen.org", 25755)
annotateDb = client['AnnotateBy']
geneCollection = annotateDb['gene']
db = client['markers']
tumorCollection = db['tumor']
statsDb = client['stats']
componentsCollection = statsDb['components']

# files
plan_conversion_file = "PLAN_C024.txt"
out_dir = "/ngd-data/prodCentralDB/C024/patients/"
convert_plan_c024 = []
# constants
study = "C024"
origin = "DNA"
api_version = "2.0.1"


# convert C024 ids to PLAN ids
def load_plan_conversion_file():
    global convert_plan_c024
    with open(plan_conversion_file) as f:
        rows = (line.split(' ') for line in f)
        convert_plan_c024 = {row[0]: row[1:] for row in rows}

    # convert_plan_c024 = dict(x.rstrip().split(None, 1) for x in f)

    # patient = study_patient_tissue[:-3]
    # if patient in convert_plan_c024:
    #    return convert_plan_c024[patient][0]
    # else:
    #    return study_patient_tissue


# Get a list of genes in the Cosmic cancer census
def get_cancer_census():
    cancer_census_genes = geneCollection.find({'cancerCensus': {"$exists": "true"}, 'cancerCensus.Tier': '1'},
                                              {"gene": 1})
    cc_genes = []
    for ccGene in cancer_census_genes:
        cc_genes.append(ccGene['gene'])
    return cc_genes


# Query the clinical db for information about the tissue
def get_clinical_info(project_run):
    kb_result = componentsCollection.find({"projectRun": project_run})
    for result in kb_result:
        if 'kb' in result:
            diagnosis = result['kb']['visit']['diagnosis']
            specimen_site = ""
            if 'samples' in result['kb']:
                specimen_site = result['kb']['samples']['sampleSource']
            tumor_collection_date = result['kb']['visit']['CollectionDate']
            order = {
                "primary_diagnosis": diagnosis,
            }
            return order


# Print out the json records to separate patient json files
def print_biomarker(study_patient_tissue, assay, variants):
    # order = get_clinical_info( project_run )
    global convert_plan_c024
    patient = study_patient_tissue[:-3]
    diagnosis = ""
    patient_id = study_patient_tissue
    print ("PATIENT ID= " + patient_id)
    if patient in convert_plan_c024:
        patient_list = convert_plan_c024[patient]
        if len(patient_list) > 0:
            patient_id = patient_list[0]
        if len(patient_list) > 1:
            diagnosis = patient_list[1].rstrip()

    print ("PATIENT ID= " + patient_id)
    order = {
        "primary_diagnosis": diagnosis
    }
    biomarkers = {
        "api_version": api_version,
        "source": "TGen",
        "order": order,
        "drug_rules_version": assay,
        "biomarkers": variants
    }

    biomarkers_json = json.dumps(biomarkers)
    # patient_id = get_patient_id(patient)

    file = open(out_dir + patient_id + ".json", "w")
    print(file)
    file.write(biomarkers_json)
    file.close()


# Query the mongo db for "Pass" variants within certain filetypes
def get_variants():
    ccGenes = get_cancer_census()
    print(ccGenes)
    studyResult = tumorCollection.find({"variants.study": study, 'variants.filter': 'PASS', "gene": {"$in": ccGenes},
                                        '$and': [{'aberration.aberration_type2': {"$ne": "Synonymous"}},
                                                 {'aberration.aberration_type2': {"$ne": "Intron Variant"}},
                                                 {'aberration.aberration_type2': {"$ne": "UTR"}}],
                                        '$or': [{'variants.filename': {'$regex': 'snpeff.final.vcf'}},
                                                {'variants.filename': {'$regex': 'Freebayes.filt.snpEff.vcf'}},
                                                {'variants.filename': {'$regex': 'rna.final.seurat.vcf'}},
                                                {'variants.filename': {'$regex': 'cna.seg.vcf'}},
                                                {'variants.filename': {'$regex': 'trn.vcf'}},
                                                {'variants.filename': {'$regex': 'trn2.vcf'}},
                                                {'variants.filename': {'$regex': 'copy_number.vcf'}},
                                                ]
                                        }
                                       ).sort([("variants.studyPatientTissue", pymongo.ASCENDING)])
    return studyResult


# Parse variant records from mongodb.  Extract relevant fields for the TumorBoardPortal json format
def parse_variants(study_result):
    variants = []
    current_patient = ""
    origin = "DNA"
    assay = ""

    for result in study_result:

        gene = result['gene']
        alteration_type = result['aberration']['aberration_type2']
        dna_allele_freq = ""
        project_run = result['variants'][0]['projectRun']
        assay = result['variants'][0]['assay']
        study_patient_tissue = result['variants'][0]['studyPatientTissue']

        if assay.endswith("STX"):
            assay = "Strexome"
        elif "ID" in assay:
            assay = "GEMExTra"

        # Convert variant fields to tumor board portal syntax
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
            print_biomarker(current_patient, assay, variants)
            variants = []
            current_patient = study_patient_tissue

        variant = {"gene": gene, "alteration_type": alteration_type, "origin": origin, "effect": effect,
                   "dna_allele_frequency": dna_allele_freq}
        variants.append(variant)

    print_biomarker(current_patient, assay, variants)


load_plan_conversion_file()
queryResult = get_variants()
parse_variants(queryResult)

client.close()
