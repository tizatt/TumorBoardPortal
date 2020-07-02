import os.path
####
# Constants
####

# Database setup
mongo_settings = []
mongo_settings['host'] = "labdb01.tgen.org"
mongo_settings['port'] = 25755
mongo_settings['annotate_db'] = "AnnotateBy"
mongo_settings['gene_coll'] = "gene"
mongo_settings['markers_db'] = "markers"
mongo_settings['tumor_coll'] = "tumor"
mongo_settings["stats_db"] = "stats"
mongo_settings["components_coll"] = "components"

# other_vcf
other_vcf = []
other_vcf['filename'] = "other_vcf"
other_vcf['search_dir'] = "/labs/ngd-data/reports/C024/"

# beatcc_dna settings
beatcc_dna = []
beatcc_dna['out_dir'] = "/labs/NMTRC/TumorBoardPortal/C024/SGR/"
beatcc_dna['api_version'] = "2.0.1"
beatcc_dna['origin'] = "DNA"
beatcc_dna['source'] = "TGen"
beatcc_dna["outfile_extension"] = "_SGR.csv"

# beatcc_rna settings
beatcc_rna = []
beatcc_rna['worksheet'] = "nmtrc_2"
beatcc_rna['api_version'] = "1"
beatcc_rna['source'] = "BeatCC"