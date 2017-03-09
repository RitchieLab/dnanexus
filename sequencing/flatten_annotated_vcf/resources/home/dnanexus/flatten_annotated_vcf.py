#!/usr/bin/python

import argparse, gzip, re, json,copy
from collections import OrderedDict

def is_number(S):
	try:
		float(S)
		return True
	except ValueError:
		return False

def parseArguments():

	parser = argparse.ArgumentParser()

	parser.add_argument('-F', action='store',
						dest='frequency',
						type=float,
						help='Filter Variants to those with a frequency less than input.  VEP must be set to filter by public database frequencies.' )

	parser.add_argument('-i', action='store', dest='VEP_Impact',
						help='Highest VEP Impact to filter too: 1=Modifier, 2=Low, 3=Moderate,4=High',
						type=int)

	parser.add_argument('--b_snp', action='store_true', default=False,
						dest='b_snp',
						help='Restrict to BiAllelic SNPs')

	parser.add_argument('--Cannonical', action='store_true', default=False,
						dest='Cannonical',
						help='Do restrict VEP IMPACT Filtering to Cannonical transcripts. Defualt is select highest overal IMPACT.')

	parser.add_argument('--Ensembl', action='store_true', default=False,
						dest='Ensembl',
						help='Use Transcripts from Ensembl for VEP filtration.  Default is to select Transcripts from RefSeq.')

	parser.add_argument('--RefSeq', action='store_true', default=False,
						dest='RefSeq',
						help='Use Transcripts from RefSeq for VEP filtration.  Default is to select Transcripts from RefSeq.')

	parser.add_argument('-c', action='store', dest='clinVarStar',
						help='ClinVar Star Level to Filter to.  Must be set to select ClinVar variants.  Set to 0 to retrive all',
						type=int)

	parser.add_argument('-p', action='store', dest='ClinicalSignificance',
						help='Highest ClinVar Significance to filter too: 8="Pathogenic", 7="Likely pathogenic", 6="Drug Response", 5="Protective",' \
						'4="Risk Factor" or "association" or "Affects", 3="Uncertain Significance" or "not provided" 2="Likely benign" 1="Benign"',
						type=int)

	parser.add_argument('--HGMD', action='store', type=int,
						dest='HGMD',
						help='HIGHEST Include variants from HGMD. 6=DM, 5=DM?, 4=DP, 3=DFP, 2=FTV, 1=FP')

	parser.add_argument('--vcf_file', action='store', type=str,
						dest='vcf_file',
						help='VCF.gz File to Filter')

	parser.add_argument('--gene_list', action='store', type=str,
						dest='gene_list',
						help='File List of Genes to Filter To. One gene per line. Note: VEP_Impact must be set.')

	parser.add_argument('--sample_list', action='store', type=str,
						dest='sample_list',
						help='File List of Samples to Filter To. One gene per line. Must be in format GHS_PT###_####')

	parser.add_argument('--out_file', action='store', type=str,
						dest='out_file',
						help='Prefix of output files')

	return parser.parse_args()

def getStar(ReviewStatus):
	ReviewStatus=ReviewStatus.strip().lower()
	if ReviewStatus == 'classified_by_single_submitter':
		return 1
	if ReviewStatus == 'criteria_provided,_single_submitter':
		return 1
	if ReviewStatus == 'classified_by_multiple_submitters':
		return 2
	if ReviewStatus == 'criteria_provided,_multiple_submitters,_no_conflicts':
		return 2
	if ReviewStatus == 'reviewed_by_expert_panel':
		return 3
	if ReviewStatus == 'reviewed_by_professional_society':
		return 4
	if ReviewStatus == 'practice_guideline':
		return 4
	return 0

def getimpact_level(IMPACT):
	IMPACT=IMPACT.strip().upper()
	if IMPACT=="HIGH":
		return 4
	if IMPACT=="MODERATE":
		return 3
	if IMPACT=="LOW":
		return 2
	if IMPACT=="MODIFIER":
		return 1
	return 0

def getHGMD_level(CLASS):
	CLASS=CLASS.strip().upper()
	if CLASS=="DM":
		return 6
	if CLASS=="DM?":
		return 5
	if CLASS=="DP":
		return 4
	if CLASS=="DFP":
		return 3
	if CLASS=="DFP":
		return 2
	if CLASS=="DFP":
		return 1
	return 0

def getClinVar_level(ClinicalSignificance):
	ClinicalSignificance=ClinicalSignificance.strip()
	if "Pathogenic" in ClinicalSignificance:
		return 8
	if "pathogenic" in ClinicalSignificance:
		return 7
	if "drug" in ClinicalSignificance:
		return 6
	if "protect" in ClinicalSignificance:
		return 5
	if "risk" in ClinicalSignificance or "Affect" in ClinicalSignificance or "association" in ClinicalSignificance:
		return 4
	if "certain" in ClinicalSignificance or "provided" in ClinicalSignificance:
		return 3
	if "benign" in ClinicalSignificance:
		return 2
	if "Benign" in ClinicalSignificance:
		return 1
	return 0

def fill_dict(data_list, ordered_data_dict):
	copyDict=copy.deepcopy(ordered_data_dict)
	for i, d in enumerate(ordered_data_dict.keys()):
		if i>=len(data_list):
			continue
		if data_list[i] is None:
			continue
		elif data_list[i].strip()=="":
			continue
		else:
			copyDict[d]=data_list[i]
	return copyDict

def parse_INFO(line):
	fields = line.split(':')[-1].replace('"',"").replace("'","").replace("'","").replace(">","").strip().split("|")
	if line.startswith("##INFO=<ID=CSQ"):
		VEP_Fields=OrderedDict()
		VEP_Fields["MAF"]=set()
		for i,f in enumerate(fields):
			VEP_Fields[f]=i
		for i,f in enumerate(fields):
			if f.endswith("MAF"):
				VEP_Fields["MAF"].add(i)
		return VEP_Fields
	elif line.startswith("##INFO=<ID=ClinVar"):
		ClinVar_fields=OrderedDict()
		for i,f in enumerate(fields):
			ClinVar_fields[f]=i
		return ClinVar_fields


def filterLine(all_fields,cli_arguments,VEP_Fields,ClinVar_fields,geneSet,IDs):
	info_field=[]
	alleles=[all_fields[3]]
	returnDict=dict()

	returnDict["CHROM"]=all_fields[0].strip()
	returnDict["POS"]=all_fields[1].strip()
	returnDict["ID"]=all_fields[2].strip()
	returnDict["REF"]=all_fields[3].strip()
	returnDict["ALT"]=all_fields[4].strip()
	returnDict["UNIQ_IDs"]=set()

	if cli_arguments.b_snp:
		pattern = re.compile("^[AGCT]$")
		if "," in fields[4]:
			return None
		elif  pattern.match(all_fields[4]) and pattern.match(all_fields[3]):
			info_field = all_fields[7].split(";")
			alleles.extend(all_fields[4].strip())
			returnDict["BI_ALLELIC"]=True
		else:
			return None
	else:

		info_field = all_fields[7].strip().split(";")
		if "," in all_fields[4]:
			alleles.extend(all_fields[4].strip().split(","))

			returnDict["BI_ALLELIC"]=False


		else:
			alleles.extend([all_fields[4].strip()])
			returnDict["BI_ALLELIC"]=True
	returnDict["ALL_ALLELES"]=alleles
	returnDict["PASS_ALLELES"]=set()

	geneList=[]
	frequency=[]
	annotations=[]

	for field in info_field:
		if field.startswith("AF"):
			AFs = field.replace("AF=","").split(",")
			AFs.insert(0,0)
			for frq_index, f in enumerate(AFs):
				if is_number(f):
					if f=="0":
						frequency.append(False)
					elif float(f)>cli_arguments.frequency:
						frequency.append(False)
					else:
						frequency.append(True)
						alleleDict=dict()
						alleleDict["ALLELE_IDs"]=set()
						alleleDict["FRQ"]=[float(f)]
						returnDict[alleles[frq_index]]=alleleDict
				else:
					return None
			if len(frequency)>0:
				if not any(frequency):
					return None
		elif field.startswith("CSQ"):
			CSQs =[x.strip() for x in field.replace("CSQ=","").strip().split(",")]
			if cli_arguments.VEP_Impact is not None:
				for i,CSQ in enumerate(CSQs):
					Cs=CSQ.split("|")
					if Cs[0] not in returnDict:
						continue
					alleleDict=returnDict[Cs[0]]
					if "VEP" in alleleDict:
						continue
					for f_field in VEP_Fields["MAF"]:
						#print VEP_Fields["MAF"]
						#print f_field
						Cs[f_field]=Cs[f_field].strip("&").strip(":").replace("&&","&").replace("&&&","&").replace("&&","&")
						if Cs[f_field]=="":
							continue
						else:
							if "&" in Cs[f_field] or ":" in Cs[f_field]:
								pubMAFs=dict()
								if "&" in Cs[f_field]:
									pubMAFs=dict(item.split(":") for item in Cs[f_field].split("&"))
								else:
									pubMAFs={Cs[f_field].split(":")[0]:Cs[f_field].split(":")[1]}
								for a in pubMAFs.keys():
									if a==Cs[0]:
										if float(pubMAFs[a])==0:
											frequency.append(False)
										elif float(pubMAFs[a])>cli_arguments.frequency:
											frequency.append(False)
										else:
											frequency.append(True)
											alleleDict["FRQ"].append(float(pubMAFs[a]))
							elif is_number(Cs[f_field]):
								if float(Cs[f_field])==0:
									frequency.append(False)
								elif float(Cs[f_field])>cli_arguments.frequency:
									frequency.append(False)
								else:
									frequency.append(True)
									alleleDict["FRQ"].append(float(Cs[f_field]))
					if len(frequency)>0:
						if not any(frequency):
							return None
					if cli_arguments.Ensembl:
						if not Cs[VEP_Fields["Feature"]].startswith("E"):
							continue
						if cli_arguments.Cannonical and Cs[VEP_Fields["CANONICAL"]]=="YES":
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
								returnDict["PASS_ALLELES"].add(alleles.index(Cs[0]))
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
						else:
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
								returnDict["PASS_ALLELES"].add(alleles.index(Cs[0]))
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
					elif cli_arguments.RefSeq:
						if Cs[VEP_Fields["Feature"]].startswith("E"):
							continue
						if cli_arguments.Cannonical and Cs[VEP_Fields["CANONICAL"]]=="YES":
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
								returnDict["PASS_ALLELES"].add(alleles.index(Cs[0]))
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
						else:
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
								returnDict["PASS_ALLELES"].add(alleles.index(Cs[0]))
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
					if cli_arguments.gene_list is not None:
						if not any(geneList):
							return None
					alleleDict["VEP"]=fill_dict(Cs, VEP_Fields)
					returnDict[Cs[0]]=alleleDict
			else:
				for i,CSQ in enumerate(CSQs):
					Cs=CSQ.split("|")
					if Cs[0] not in returnDict:
						continue
					alleleDict=returnDict[Cs[0]]
					if "VEP" in alleleDict:
						continue
					if Cs[VEP_Fields["Feature"]].startswith("E"):
						continue
					if cli_arguments.Cannonical and Cs[VEP_Fields["CANONICAL"]]=="YES":
						if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
							annotations.append(True)
							returnDict["PASS_ALLELES"].add(alleles.index(Cs[0]))
						else:
							annotations.append(False)
						if cli_arguments.gene_list is not None:
							if Cs[VEP_Fields["SYMBOL"]] in geneSet:
								geneList.append(True)
				alleleDict["VEP"]=fill_dict(Cs, VEP_Fields)
				returnDict[Cs[0]]=alleleDict
		elif field.startswith("ClinVar.TSV.Jan2017="):
			clinvar = field.replace("ClinVar.TSV.Jan2017=","").split("|")
			if cli_arguments.clinVarStar is not None:
				if len(clinvar)<ClinVar_fields["ReviewStatus"]:
					annotations.append(False)
				if clinvar[0] not in returnDict:
					continue
				alleleDict=returnDict[clinvar[0]]
				if getStar(clinvar[ClinVar_fields["ReviewStatus"]])>=cli_arguments.clinVarStar:
					if getClinVar_level(clinvar[ClinVar_fields["ClinicalSignificance"]])>=cli_arguments.ClinicalSignificance:
						annotations.append(True)
						returnDict["PASS_ALLELES"].add(alleles.index(Cs[0]))
					else:
						annotations.append(False)
				else:
					annotations.append(False)
			alleleDict["ClinVar"]=fill_dict(clinvar, ClinVar_fields)
			returnDict[clinvar[0]]=alleleDict
		elif field.startswith("CLASS="):
			if cli_arguments.HGMD is not None:
				if getHGMD_level(field.replace("CLASS=",""))>=cli_arguments.HGMD:
					annotations.append(True)
					returnDict["PASS_ALLELES"].add(alleles.index(Cs[0]))
				else:
					annotations.append(False)
			returnDict["HGMD_CLASS"]=field.replace("CLASS=","")
		elif field.startswith("GENE="):
			returnDict["HGMD_GENE"]=field.replace("GENE=","")
		elif field.startswith("STRAND="):
			returnDict["HGMD_STRAND"]=field.replace("STRAND=","")
		elif field.startswith("PROT="):
			returnDict["HGMD_PROT"]=field.replace("PROT=","")
		elif field.startswith("DB="):
			returnDict["HGMD_DB"]=field.replace("DB=","")
		elif field.startswith("PHEN="):
			returnDict["HGMD_PHEN"]=field.replace("PHEN=","")

	if not any(annotations):
		if (cli_arguments.VEP_Impact is not None or cli_arguments.HGMD is not None or cli_arguments.clinVarStar is not None):
			return None
		elif cli_arguments.gene_list is not None:
			if not any(geneList):
				return None
	if len(frequency)>0:
		if not any(frequency):
			return None
	formatFields=all_fields[8].split(":")
	GT=-1
	for i, f in enumerate(formatFields):
		if f=="GT":
			GT=i
			break
	if len(list(returnDict["PASS_ALLELES"]))==0:
		return None
	for i, s in enumerate(all_fields):
		if i>8:
			if i in IDs:
				GT_fields=s.split(":")[GT].split("/")
				if GT_fields[0]=="." or GT_fields[1]==".":
					continue
				if GT_fields[0]!="0" and GT_fields[0] in returnDict["PASS_ALLELES"]:
					returnDict[alleles[GT_fields[0]]]["ALLELE_IDs"].add(IDs[i])
					returnDict["UNIQ_IDs"].add(IDs[i])
				if GT_fields[1]!="0" and GT_fields[1] in returnDict["PASS_ALLELES"]:
					returnDict[alleles[GT_fields[1]]]["ALLELE_IDs"].add(IDs[i])
					returnDict["UNIQ_IDs"].add(IDs[i])
	returnDict["PASS_ALLELES"]=list(returnDict["PASS_ALLELES"])
	if len(list(returnDict["UNIQ_IDs"]))==0:
		return None
	returnDict["UNIQ_IDs"]=list(returnDict["UNIQ_IDs"])
	for a in alleles:
		if a in returnDict:
			if len(list(returnDict[a]["ALLELE_IDs"]))==0:
				returnDict.pop(a,None)
			else:
				returnDict[a]["ALLELE_IDs"]=list(returnDict[a]["ALLELE_IDs"])
	return returnDict


def main():

	cli_arguments =  parseArguments()

	outDict = dict()
	outDict["UNIQ_IDs"]=set()

	outDict["CLI_ARGUMENTS"]=str(cli_arguments)

	if not cli_arguments.vcf_file.endswith("vcf.gz"):
		print "Not a bgziped vcf file"
		return
	geneSet=set()
	if cli_arguments.gene_list is not None:
		#print cli_arguments.gene_list
		geneList = open(cli_arguments.gene_list,'r')
		for line in geneList:
			line=line.strip().upper()
			if line=="":
				continue
			else:
				geneSet.add(line)
	outDict["INPUT_GENE_LIST"]=list(geneSet)
	sampleSet=set()
	if cli_arguments.sample_list is not None:
		#lsprint cli_arguments.sample_list
		sampleList = open(cli_arguments.sample_list,'r')
		for line in sampleList:
			line=line.strip().upper()
			if line=="":
				continue
			else:
				sampleSet.add(line)

	vcf_file = gzip.open(str(cli_arguments.vcf_file),'r')
	VEP_Fields=dict()
	ClinVar_fields=dict()
	IDs=dict()

	for line in vcf_file:
		if line.startswith("#"):
			#print line.strip()
			if line.startswith("##INFO=<ID=CSQ"):
				VEP_Fields=parse_INFO(line.strip())
				outDict["VEP_FIELDS"]=list(VEP_Fields.keys())
			elif line.startswith("##INFO=<ID=ClinVar"):
				ClinVar_fields=parse_INFO(line.strip())
				outDict["CLINVAR_FIELDS"]=list(ClinVar_fields.keys())
			elif line.startswith("#CHROM"):
				fields = [x.strip() for x in line.strip().split("\t")]
				for i,f in enumerate(fields):
					if f in sampleSet:
						IDs[i]=f
				if not IDs:
					#print "No Samples found"
					return
				else:
					outDict["INPUT_SAMPLES"]=list(sampleSet)
					outDict["VCF_SAMPLES"]=list(IDs.values())
					outDict["HGMD_CLASS"]="Mutation Category, https://portal.biobase-international.com/hgmd/pro/global.php#cats"
					outDict["HGMD_GENE"]="Gene symbol"
					outDict["HGMD_STRAND"]="Gene strand"
					outDict["HGMD_PROT"]="Protein annotation"
					outDict["HGMD_DB"]="dbSNP identifier, build 137"
					outDict["HGMD_PHEN"]="HGMD primary phenotype"
		else:
			all_fields=line.strip().split("\t")
			lineDict=filterLine(all_fields,cli_arguments,VEP_Fields,ClinVar_fields,geneSet,IDs)
			if lineDict is not None:
				outDict["UNIQ_IDs"].update(lineDict["UNIQ_IDs"])
				outDict[all_fields[0]+":"+all_fields[1]+":"+all_fields[3]+":"+all_fields[4]]=lineDict
	vcf_file.close()
	outDict["UNIQ_IDs"]=list(outDict["UNIQ_IDs"])
	out_file_json = gzip.open(cli_arguments.vcf_file.replace(".vcf.",".filtered.json."),'w')
	out_file_json.write(json.dumps(outDict, sort_keys=True,indent=4, separators=(',', ': ')))
	return

if __name__ == '__main__':
	main()
	exit()
