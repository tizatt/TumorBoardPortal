
def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def print_to_csv(json):
    header = "source,api version,drug rules version," \
             "specimen collection date,specimen type," \
             "specimen site,disease type,gene,origin," \
             "alteration type,effect,tmb level," \
             "tmb quantity,msi,recommended medications"
    csv_output = ""
    print(json)
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
                tmb_level = record["effect"]
                tmb_quantity = record["alteration_type"]
        elif counter == 0:
            csv_first_row = csv_first_row + "," + csv_record
            if "drugs" in record:
                for drug in record["drugs"]:
                    first_row_drugs = first_row_drugs + "," + drug["drug_name"]

        else:
            csv_record = csv_record + ",,,,"
            csv_record = "TGen,,,,,,,"+ csv_record

            if "drugs" in record:
                for i in range(len(record["drugs"])):
                    if i > 0:
                        csv_record=csv_record + ";"
                    csv_record = csv_record + record["drugs"][i]["drug_name"]

            csv_output = csv_output + csv_record + "\n"

        counter += 1

    csv_first_row = csv_first_row + "," + tmb_level + "," + tmb_quantity + "," + msi + first_row_drugs
    csv_output = header + "\n" + csv_first_row + "\n" + csv_output
    print("csv output")
    return csv_output.strip()




def parse_other_vcf( other_vcf ):
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
                    tmb_value = round(tmb_value_float, 2)
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

