from pymongo import MongoClient
import pymongo
import json
import portal_convert
import os.path

# Database setup
host = "labdb01.tgen.org"
port = 25755

study_dir = "/ngd-data/reports/C024/"

# files
# plan_conversion_file = "PLAN_C024.txt"
other_vcf = "other.vcf"
out_dir = "/NMTRC/TumorBoardPortal/C024/SGR/"

# constants

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



def find_patients(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


def get_specimen_information(study_id):
    # Query the clinical db for information about the studyID
    specimen_type = ""
    kb_result = componentsCollection.find({"projectRun": {'$regex': study_id}})
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
                "specimen_site": specimen_site,
                "specimen_type": specimen_type,
                "tumor_collection_date": tumor_collection_date
            }
            return order


def get_patient_id(study_id):

    kb_result = componentsCollection.find({"projectRun": {'$regex': study_id}})
    for result in kb_result:
        if 'kb' in result:
            patient_id = result['kb']['patients']['studyPatientID']
            return patient_id


# Print out the json records to separate patient json files
def print_biomarker(study_patient_tissue, assay, variants):
    patient = study_patient_tissue[:-3]
    patient_id = study_patient_tissue

    if patient in convert_plan_c024:
        patient_list = convert_plan_c024[patient]
        if len(patient_list) > 0:
            patient_id = patient_list[0]
        if len(patient_list) > 1:
            diagnosis = patient_list[1].rstrip()
    order = get_specimen_information(patient)

    plan_id = get_patient_id(patient)
    if plan_id is not None:
        patient_id = plan_id

    biomarkers = {
        "api_version": api_version,
        "source": "TGen",
        "order": order,
        "drug_rules_version": assay,
        "biomarkers": variants
    }

    biomarkers_json = json.dumps(biomarkers)

    # portal_convert.print_to_csv(biomarkers_json)
    # portal_convert.parse_other_vcf()

    filename = out_dir + patient_id + ".json"
    print(filename)
    if not os.path.exists(filename):
        csv_output = portal_convert.print_to_csv(biomarkers)
        file = open(filename, 'w')
        csv_file = open(filename.replace(".json","_SGR.csv"), 'w')
        csv_file.write(csv_output)
        file.write(biomarkers_json)
        file.close()
        csv_file.close()


# Need to change this to run on single patients so that if patient variant count is zero, it will still generate a report
#  Start by finding the "other.vcf" files and from there if that studyPatient is new then generate a report. This would give
#  the studypatient info needed, which is currently missing.
# Query the mongo db for "Pass" variants within certain filetypes

def get_patients():
    study_patient_tissue = []
    other_vcfs = []
    results = find_patients("other.vcf", study_dir)
    for result in results:
        result_arr = result.split("/")
        patient = str(result_arr[4])
        tissue = str(result_arr[6])
        study = str(result_arr[3])
        study_patient_tissue.append(study + "_" + patient + "_" + tissue)
        other_vcfs.append(result)

    return other_vcfs, study_patient_tissue


def get_variants(study_patient_tissue, cc_genes):
    print(study_patient_tissue)

    study_result = tumorCollection.find(
        {"variants.studyPatientTissue": study_patient_tissue, 'variants.filter': 'PASS', "gene": {"$in": cc_genes},
         '$and': [{'aberration.aberration_type2': {"$ne": "Synonymous"}},
                  {'aberration.aberration_type2': {"$ne": "Intron Variant"}},
                  {'aberration.aberration_type2': {"$ne": "UTR"}}],
         '$or': [{'variants.filename': {'$regex': 'snpeff.final.vcf'}},
                 {'variants.filename': {'$regex': 'Freebayes.filt.snpEff.vcf'}},
                 {'variants.filename': {'$regex': 'rna.final.seurat.vcf'}},
                 {'variants.filename': {'$regex': 'cna.seg.vcf'}},
                 {'variants.filename': {'$regex': 'trn.vcf'}},
                 {'variants.filename': {'$regex': 'trn2.vcf'}},
                 {'variants.filename': {'$regex': 'copy_number.vcf'}}
                 ]
         }
    ).sort([("variants.studyPatientTissue", pymongo.ASCENDING)])
    return study_result


# Parse variant records from mongodb.  Extract relevant fields for the TumorBoardPortal json format
def parse_variants():
    other_vcfs, study_patient_tissues = get_patients()
    cc_genes = get_cancer_census()
    variants = []
    assay = ""

    for (study_patient_tissue, other_vcf) in zip(study_patient_tissues, other_vcfs):
        study_result = get_variants(study_patient_tissue, cc_genes)

        for result in study_result:

            gene = result['gene']
            alteration_type = result['aberration']['aberration_type2']
            dna_allele_freq = ""
            assay = result['variants'][0]['assay']

            if assay.endswith("STX"):
                assay = "Strexome"
            elif "ID" in assay:
                assay = "GEMExTra"

            # Convert variant fields to tumor board portal syntax
            if 'SEURAT_AR_TUMOR' in result['variants'][0]:
                dna_allele_freq = result['variants'][0]['SEURAT_AR_TUMOR']
            elif 'AR2' in result['variants'][0]:
                dna_allele_freq = result['variants'][0]['AR2']

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

            variant = {"gene": gene, "alteration_type": alteration_type, "origin": origin, "effect": effect,
                       "dna_allele_frequency": dna_allele_freq}
            variants.append(variant)

        tmb_msi_records = portal_convert.parse_other_vcf( other_vcf)
        for tm_record in tmb_msi_records:
            variants.append(tm_record)
        print_biomarker(study_patient_tissue, assay, variants)
        variants = []


# Get a list of genes in the Cosmic cancer census
def get_cancer_census():
    cancer_census_genes = geneCollection.find({'cancerCensus': {"$exists": "true"},
                                               "$or": [{"cancerCensus.Tier": "1"},
                                                       {"cancerCensus.Tier": "2"}]},
                                              {'gene': 1})
    cc_genes = []
    for ccGene in cancer_census_genes:
        cc_genes.append(ccGene['gene'])
    return cc_genes


parse_variants()

client.close()
