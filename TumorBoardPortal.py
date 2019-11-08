from pymongo import MongoClient
import pymongo
import json
import pprint
import os.path

# Database setup
host="labdb01.tgen.org"
port=25755

# files
plan_conversion_file = "PLAN_C024.txt"
other_vcf = "other.vcf"
out_dir = "/Volumes/ngd-data/prodCentralDB/C024/patients/"

# constants
study = "C024"
origin = "DNA"
api_version = "2.0.1"

# global variable
convert_plan_c024 = []

client = MongoClient(host, port)
annotateDb = client['AnnotateBy']
geneCollection = annotateDb['gene']
db = client['markers']
tumorCollection = db['tumor']
statsDb = client['stats']
componentsCollection = statsDb['components']


# convert ids to PLAN ids
def load_plan_conversion_file():
    global convert_plan_c024
    with open(plan_conversion_file) as f:
        rows = (line.split(' ') for line in f)
        convert_plan_c024 = {row[0]: row[1:] for row in rows}


def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

    # Parse the "other.vcf" file that Ashion provides.  This has the TMB / MSI values


def parse_other_vcf():
    tmb_msi_record = []
    tmb_category = ""
    tmb_value = ""
    msi_category = ""
    with open(other_vcf) as fp:
        for line in fp:
            if line.startswith("##"):
                a = 5
            else:
                vcf_list = line.split("\t")
                info_group = vcf_list[7]
                if info_group.startswith("TMB"):
                    tmb_line = info_group.split(";")
                    tmb_value_float = float(tmb_line[0].split("=")[1])
                    tmb_value = truncate(tmb_value_float, 2)
                    tmb_category = tmb_line[1].split("=")[1]
                if info_group.startswith("MSI"):
                    msi_category = info_group.split("=")[1]

    tmb_record = {
        "gene": "TMB",
        "alteration_type": "TMB " + tmb_category.strip(),
        "effect": str(tmb_value) + " mut/Mb"
    }

    msi_record = {
        "gene": "MSI",
        "alteration_type": "MSI: " + msi_category.strip()
    }
    tmb_msi_record.append(tmb_record)
    tmb_msi_record.append(msi_record)

    fp.close()

    return tmb_msi_record

def get_specimen_information(study_id):
    # Query the clinical db for information about the studyID
    specimen_type = ""
    kb_result = componentsCollection.find({"projectRun": {'$regex': study_id }})
    for result in kb_result:
        if 'kb' in result:
            diagnosis = result['kb']['visit']['diagnosis']
            specimen_site = ""
            if 'samples' in result['kb']:
                specimen_site = result['kb']['samples']['sampleSource']
                specimen_type = result['kb']['samples']['sampleType']
            tumor_collection_date = result['kb']['visit']['CollectionDate']
            order = {
                "primary_diagnosis": diagnosis,
                "specimen_site" : specimen_site,
                "specimen_type" : specimen_type,
                "tumor_collection_date" : tumor_collection_date
            }
            return order





# Get a list of genes in the Cosmic cancer census
def get_cancer_census():
    cancer_census_genes = geneCollection.find({'cancerCensus': {"$exists": "true"}, 'cancerCensus.Tier': '1'},
                                              {"gene": 1})
    cc_genes = []
    for ccGene in cancer_census_genes:
        cc_genes.append(ccGene['gene'])
    return cc_genes



# Print out the json records to separate patient json files
def print_biomarker(study_patient_tissue, assay, variants):

    global convert_plan_c024
    patient = study_patient_tissue[:-3]
    patient_id = study_patient_tissue
    print ("PATIENT ID= " + patient_id)
    if patient in convert_plan_c024:
        patient_list = convert_plan_c024[patient]
        if len(patient_list) > 0:
            patient_id = patient_list[0]
        if len(patient_list) > 1:
            diagnosis = patient_list[1].rstrip()
    print(patient_id)
    order = get_specimen_information( patient )

    biomarkers = {
        "api_version": api_version,
        "source": "TGen",
        "order": order,
        "drug_rules_version": assay,
        "biomarkers": variants
    }

    biomarkers_json = json.dumps(biomarkers)
    # patient_id = get_patient_id(patient)
    filename = out_dir + patient_id + ".json"
    if not os.path.exists( filename ):
        file = open(filename,'w')
        print(file)
        print(json.dumps(biomarkers_json,indent=4))
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
