import xlrd
from pymongo import MongoClient
import json

host = "labdb01.tgen.org"
port = 25755

StudyID = 2
Gene = 3
Drug = 4
DrugInfo = 6
ZScore = 9
TumorPercent = 11

client = MongoClient(host, port)
statsDb = client['stats']
componentsCollection = statsDb['components']

bcc_file = '/Users/tizatt/Desktop/PLAN-094-09.xlsx'
other_vcf = '/Volumes/ngd-data/prodCentralDB/TumorBoardPortal/other.vcf'


def read_xls_file():
    wb = xlrd.open_workbook(bcc_file)
    bcc_data = wb.sheet_by_name('nmtrc_2')

    return bcc_data


def get_specimen_information(studyId):
    # Query the clinical db for information about the studyID

    kb_result = componentsCollection.find({"kb.patients.studyPatientID": studyId})
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
        "alteration_type": tmb_category.strip(),
        "effect": str(tmb_value) + " mut/Mb"
    }

    msi_record = {
        "gene": "MSI",
        "alteration_type": msi_category.strip()
    }
    tmb_msi_record.append(tmb_record)
    tmb_msi_record.append(msi_record)

    fp.close()

    return tmb_msi_record


def print_to_csv(json):
    header = "source,api version,drug rules version," \
             "specimen collection date,specimen type," \
             "specimen site,disease type,gene,origin," \
             "alteration type,effect,tmb level," \
             "tmb quantity,msi,recommended medications"
    csv_output = ""
    csv_first_row = json["source"] + "," + json["api_version"] + "," + json["drug_rules_version"] \
                    + "," + json["order"]["tumor_collection_date"] + "," + json["order"]["specimen_type"] \
                    + "," + json["order"]["specimen_site"] + "," + json["order"]["primary_diagnosis"]
    first_row_drugs = ""
    tmb_quantity = ""
    tmb_level = ""
    msi = ""

    counter = 0
    for record in json["biomarkers"]:

        if record["gene"] == "MSI" or record["gene"] == "TMB":
            a =1
        else:
            csv_record = record["gene"] + "," + record["origin"] + "," + record["alteration_type"] \
                     + "," + "\""+record["effect"]+"\""

        if record["gene"] == "MSI" or record["gene"] == "TMB":
            if record["gene"] == "MSI":
                msi = record["alteration_type"]
            else:
                tmb_level = record["alteration_type"]
                tmb_quantity = record["effect"]
        elif counter == 0:
            csv_first_row = csv_first_row + "," + csv_record
            for drug in record["drugs"]:
                first_row_drugs = first_row_drugs + "," + drug["drug_name"]
        else:
            csv_record = csv_record + ",,,,"
            for drug in record["drugs"]:
                csv_record = "BeatCC,,,,,,,"+ csv_record + drug["drug_name"]
            csv_output = csv_output + csv_record + "\n"

        counter += 1

    csv_first_row = csv_first_row + "," + tmb_level + "," + tmb_quantity + "," + msi + first_row_drugs
    csv_output = header + "\n" + csv_first_row + "\n" + csv_output
    print(csv_output.strip())


def parse_beatcc_rna(bcc_data):
    header = bcc_data.col(0)
    row_counter = 0
    aberration_records = []
    drug_name = "%"
    study_id = ""
    report_version = ""

    for row in header:
        row_val = str(row)
        if "Study ID" in row_val:
            study_id = bcc_data.cell_value(row_counter, StudyID)
        if "Report Version" in row_val:
            report_version = bcc_data.cell_value(row_counter, StudyID)[1:]
        if "Drug -" in row_val:
            drug_name = bcc_data.cell_value(row_counter, Drug)
        if drug_name in row_val:
            gene = bcc_data.cell_value(row_counter, Gene)
            drug_info = bcc_data.cell_value(row_counter, DrugInfo)

            expression_type = bcc_data.cell_value(row_counter, ZScore)

            zscore_val = truncate(bcc_data.cell_value(row_counter + 1, ZScore), 2)

            tumor_percent_val = truncate(bcc_data.cell_value(row_counter + 1, TumorPercent), 2)
            drugs = [{"drug_name": drug_name.lower(), "drug_status": "Predicted Beneficial"}]

            aberration_record = dict(
                gene=gene,
                alteration_type=expression_type + "expression",
                origin="RNA",
                effect=expression_type + ", NRZ=" + str(zscore_val) + ", CRC=" + str(tumor_percent_val),
                drugs=drugs
            )

            aberration_records.append(aberration_record)

        row_counter += 1

    tmb_msi_records = parse_other_vcf()
    for tm_record in tmb_msi_records:
        aberration_records.append(tm_record)

    specimen_info = get_specimen_information(study_id)
    bcc_json = {"api_version": "1", "drug_rules_version": report_version, "source": "BeatCC",
                "order": specimen_info,
                "biomarkers": aberration_records
                }
    print_to_csv(bcc_json)


# print(json.dumps(bcc_json, indent=4))

#   for rownum in range(bcc_data.nrows):
#      print(bcc_data.row_values(rownum))
#      print(bcc_data.cell(1))


bcc_data = read_xls_file()
parse_beatcc_rna(bcc_data)
