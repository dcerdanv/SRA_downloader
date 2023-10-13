import pandas as pd
from snakemake.utils import validate

#### GLOBAL PARAMETERS ####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config['outdir']
LOGDIR = config['logdir']
TMPDIR = config['tmpdir']


compress = config['parameters']['compress']
remove_sra = config['parameters']['remove_sra']

#### LOAD SAMPLES TABLE ###
samples = pd.read_table(config["samples"]).set_index("sample_ID", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")


#### GLOBAL scope functions ####
def get_resource(rule,resource) -> int:
	'''
	Attempt to parse config.yaml to retrieve resources available for a given
	rule. It will revert to default if a key error is found. Returns an int.
	with the allocated resources available for said rule. Ex: "threads": 1
	'''

	try:
		return config['resources'][rule][resource]
	except KeyError: # TODO: LOG THIS
		print(f'Failed to resolve resource for {rule}/{resource}: using default parameters')
		return config["resources"]['default'][resource]


#### Load rules ####
include: 'rules/sra_downloader.smk'


### RULES ###
rule all:
	input:
		get_all_input

