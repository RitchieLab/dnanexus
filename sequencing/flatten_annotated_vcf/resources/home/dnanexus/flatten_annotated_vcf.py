#!/usr/bin/python

import argparse, gzip, re, json,copy
from collections import OrderedDict
import pprint

def is_number(S):
	try:
		float(S)
		return True
	except ValueError:
		return False

def stripAllelesLIST(id_list):
	pos = int(id_list[1])
	ref = id_list[2]
	alt = id_list[3]

	if ref=="-":
		ref="."
	if alt=="-":
		alt="."
	if ref=="." or alt == ".":
		pass
	if len(ref)==1 and len(ref)==len(alt):
		pass
	elif len(ref)>=len(alt):
		numMatchChar = 0
		for idx, char in enumerate(alt):
			if alt[idx]==ref[idx]:
				numMatchChar+=1
			else:
				break

		if numMatchChar>0:
			ref=ref[numMatchChar:]
			pos+=numMatchChar
			if numMatchChar<len(alt):
				alt=alt[numMatchChar:]
			else:
				alt="."

	elif len(ref)<=len(alt):
		numMatchChar = 0
		for idx, char in enumerate(ref):
			if alt[idx]==ref[idx]:
				numMatchChar+=1
			else:
				break

		if numMatchChar>0:
			alt=alt[numMatchChar:]
			pos+=numMatchChar
			if numMatchChar<len(ref):
				ref=ref[numMatchChar:]
			else:
				ref="."

	if len(ref) == len(alt) and ref != "." and alt != "."  and len(ref)>1 and ref[-1]==alt[-1]:
		en_ref = ref
		for idx, char in reversed(list(enumerate(en_ref))):
			if char == alt[idx] and len(alt)>1:
				ref=ref[:-1]
				alt=alt[:-1]
			else:
				break

	id_list[1]=str(pos)
	id_list[2]=ref
	id_list[3]=alt

	return id_list

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
	copyDict.pop('MAF', None)
	newDict=dict()

	for i, d in enumerate(copyDict.keys()):
		if i>=len(data_list):
			continue
		if data_list[i] is None:
			continue
		elif data_list[i].strip()=="":
			continue
		else:
			newDict[d]=data_list[i]
	#print ordered_data_dict
	#print copyDict
	#print data_list
	return newDict

def parse_INFO(line):
	fields = line.split(':')[-1].replace('"',"").replace("'","").replace("'","").replace(">","").strip().split("|")
	if line.startswith("##INFO=<ID=CSQ"):
		VEP_Fields=OrderedDict()
		VEP_Fields["MAF"]=set()
		for i,f in enumerate(fields):
			VEP_Fields[f.strip()]=i
		for i,f in enumerate(fields):
			if f.endswith("MAF"):
				VEP_Fields["MAF"].add(i)
		return VEP_Fields
	elif line.startswith("##INFO=<ID=ClinVar"):
		ClinVar_fields=OrderedDict()
		for i,f in enumerate(fields):
			ClinVar_fields[f.strip()]=i
		return ClinVar_fields

def filterLine(all_fields,cli_arguments,VEP_Fields,ClinVar_fields,geneSet,IDs):
	pp = pprint.PrettyPrinter(indent=4)

	info_field=[]
	alleles=[all_fields[3].strip()]
	returnDict=dict()

	returnDict["CHROM"]=all_fields[0].strip()
	returnDict["POS"]=all_fields[1].strip()
	if all_fields[2].strip()!=".":
		returnDict["ID"]=all_fields[2].strip()
	returnDict["REF"]=all_fields[3].strip()
	returnDict["ALT"]=all_fields[4].strip()
	returnDict["UNIQ_IDs"]=set()
	returnDict["STRIPPED_ALLELS"]=dict()

	if cli_arguments.b_snp:
		pattern = re.compile("^[AGCT]$")
		if "," in all_fields[4]:
			return None
		elif  pattern.match(all_fields[4]) and pattern.match(all_fields[3]):
			info_field = all_fields[7].split(";")
			alleles.extend(all_fields[4].strip())
			returnDict["BI_ALLELIC"]=True
			returnDict["STRIPPED_ALLELS"][stripAllelesLIST([all_fields[0].strip(),all_fields[1].strip(),all_fields[3].strip(),all_fields[4].strip()])[-1]]=1
		else:
			return None
	else:
		info_field = all_fields[7].strip().split(";")
		if "," in all_fields[4]:
			alleles.extend(all_fields[4].strip().split(","))
			returnDict["BI_ALLELIC"]=False
			for i,a in enumerate(all_fields[4].strip().split(",")):
				if a == "*":
					continue
				else:
					returnDict["STRIPPED_ALLELS"][stripAllelesLIST([all_fields[0].strip(),all_fields[1].strip(),all_fields[3].strip(),a])[-1]]=(i+1)
		else:
			alleles.extend([all_fields[4].strip()])
			returnDict["BI_ALLELIC"]=True
			returnDict["STRIPPED_ALLELS"][stripAllelesLIST([all_fields[0].strip(),all_fields[1].strip(),all_fields[3].strip(),all_fields[4].strip()])[-1]]=1
	returnDict["ALL_ALLELES"]=alleles
	returnDict["PASS_ALLELES"]=set()

	geneList=[]
	frequency=[]
	annotations=[]


	for field in info_field:
		if field.startswith("AF="):
			AFs = field.replace("AF=","").split(",")
			for o_f in info_field:
				if o_f.startswith("AF_Orig="):
					AFs = o_f.replace("AF_Orig=","").split(",")
			AFs.insert(0,'0')
			for frq_index, f in enumerate(AFs):
				if is_number(f):
					if f.strip()=="0" or frq_index==0 or float(f)==0:
						frequency.append(False)
					elif float(f)>cli_arguments.frequency:
						frequency.append(False)
					else:
						frequency.append(True)
						alleleDict=dict()
						alleleDict["ALLELE_IDs"]=set()
						alleleDict["HOM_SAMPLES"]=set()
						alleleDict["HET_SAMPLES"]=set()
						alleleDict["FRQ"]=[float(f)]
						alleleDict["AC"]=-1
						for a_c in info_field:
							if a_c.startswith("AC="):
								ACs =  a_c.replace("AC=","").split(",")
								for o_c in info_field:
									if o_c.startswith("AC_Orig="):
										ACs =  o_c.replace("AC_Orig=","").split(",")
								ACs.insert(0,'0')
								for ac_index, ac in enumerate(ACs):
									if ac_index==0:
										continue
									alleleDict["AC"]=int(ac)
						returnDict[alleles[frq_index]]=alleleDict
						returnDict["AN"]=-1


				else:
					return None
			if len(frequency)>0:
				if not any(frequency):
					return None
		elif field.startswith("AN="):
			AN = int(field.replace("AN=",""))
			for o_n in info_field:
				if o_n.startswith("AN_Orig="):
					AN = int(o_n.replace("AN_Orig=",""))
			returnDict["AN"]=AN
		elif field.startswith("CSQ="):
			CSQs =[x.strip() for x in field.replace("CSQ=","").strip().split(",")]
			if cli_arguments.VEP_Impact is not None:
				for i,CSQ in enumerate(CSQs):
					Cs=CSQ.split("|")
					if not Cs[0] in returnDict:
						if Cs[0] == "-":
							Cs[0] ="."
						if Cs[0] in returnDict["STRIPPED_ALLELS"]:
							Cs[0]=alleles[returnDict["STRIPPED_ALLELS"][Cs[0]]]
						if not Cs[0] in returnDict:
							continue
					alleleDict=returnDict[Cs[0]]
					if "VEP" in alleleDict:
						continue
					for f_field in VEP_Fields["MAF"]:
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
						if cli_arguments.Cannonical:
							if Cs[VEP_Fields["CANONICAL"]]=="YES":
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
						else:
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
					elif cli_arguments.RefSeq:
						if Cs[VEP_Fields["Feature"]].startswith("E"):
							continue
						if cli_arguments.Cannonical:
							if Cs[VEP_Fields["CANONICAL"]]=="YES":
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

						else:
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
			if cli_arguments.gene_list is not None:
				if not any(geneList):
					return None
		elif field.startswith("ClinVar.TSV.Jan2017="):
			clinvar = field.replace("ClinVar.TSV.Jan2017=","").split("|")
			if not cli_arguments.clinVarStar is None:
				if len(clinvar)<ClinVar_fields["ReviewStatus"]:
					annotations.append(False)
<<<<<<< HEAD
					continue 
=======
>>>>>>> parent of f80d591... asdf
				if clinvar[0] not in returnDict:
					if clinvar[0] == "-":
						clinvar[0] ="."
					if clinvar[0] in returnDict["STRIPPED_ALLELS"]:
						clinvar[0]=alleles[returnDict["STRIPPED_ALLELS"][clinvar[0]]]
					if not clinvar[0] in returnDict:
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
			elif clinvar[0] in returnDict:
				alleleDict=returnDict[clinvar[0]]
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
				if GT_fields[0]=="0" and GT_fields[1]=="0":
					continue
				if GT_fields[0]==GT_fields[1] and int(GT_fields[0]) in returnDict["PASS_ALLELES"]:
						returnDict["UNIQ_IDs"].add(IDs[i])
						returnDict[alleles[int(GT_fields[0])]]["ALLELE_IDs"].add(IDs[i])
						returnDict[alleles[int(GT_fields[0])]]["HOM_SAMPLES"].add(IDs[i])
				else:
					if GT_fields[0]!="0" and int(GT_fields[0]) in returnDict["PASS_ALLELES"]:
						returnDict[alleles[int(GT_fields[0])]]["ALLELE_IDs"].add(IDs[i])
						returnDict["UNIQ_IDs"].add(IDs[i])
						returnDict[alleles[int(GT_fields[0])]]["HET_SAMPLES"].add(IDs[i])
					if GT_fields[1]!="0" and int(GT_fields[1]) in returnDict["PASS_ALLELES"]:
						returnDict[alleles[int(GT_fields[1])]]["ALLELE_IDs"].add(IDs[i])
						returnDict["UNIQ_IDs"].add(IDs[i])
						returnDict[alleles[int(GT_fields[1])]]["HET_SAMPLES"].add(IDs[i])
	if len(list(returnDict["UNIQ_IDs"]))==0:
		return None
	returnDict["UNIQ_IDs"]=list(returnDict["UNIQ_IDs"])
	returnDict["PASS_ALLELES"]=list(returnDict["PASS_ALLELES"])
	for a in alleles:
		if a in returnDict:
			if len(list(returnDict[a]["ALLELE_IDs"]))==0:
				returnDict.pop(a,None)
			else:
				returnDict[a]["ALLELE_IDs"]=list(returnDict[a]["ALLELE_IDs"])
				if "HET_SAMPLES" in returnDict[a]:
					if len(returnDict[a]["HET_SAMPLES"])>0:
						returnDict[a]["HET_SAMPLES"]=list(returnDict[a]["HET_SAMPLES"])
					else:
						returnDict[a].pop("HET_SAMPLES",None)
				if "HOM_SAMPLES" in returnDict[a]:
					if len(returnDict[a]["HOM_SAMPLES"])>0:
						returnDict[a]["HOM_SAMPLES"]=list(returnDict[a]["HOM_SAMPLES"])
					else:
						returnDict[a].pop("HOM_SAMPLES",None)

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
					elif not sampleSet:
						IDs[i]=f
				if not IDs:
					print "NO IDS!"
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
	#out_file_json.write(str(outDict))
	out_file_json.close()
	print "json printed"
	out_tsv(cli_arguments,outDict)
	print "tsv printed"
	return

def out_tsv(cli_arguments,outDict):
	out_file_tsv = gzip.open(cli_arguments.vcf_file.replace(".vcf.",".filtered.tsv."),'w')
	out_file_tsv.write("SAMPLE_COUNT\tHOM_COUNT\tHET_Count\tBI_ALLELIC\tCHROM\tPOS\tREF\tALT\tMAF\tAN\tAC\t")
	outDict["VEP_FIELDS"].remove('MAF')
	out_file_tsv.write("VEP:"+"\tVEP:".join(outDict["VEP_FIELDS"])+"\t")
	out_file_tsv.write("CLINVAR:"+"\tCLINVAR:".join(outDict["CLINVAR_FIELDS"])+"\t")
	out_file_tsv.write("ID\tHGMD_CLASS\tHGMD_DB\tHGMD_GENE\tHGMD_PHEN\tHGMD_PROT\tHGMD_STRAND\tALL_SAMPLES\tHOM_SAMPLES\tHET_SAMPLES\n")

	for key in outDict.keys():
		if ":" in key:
			lineDict=outDict[key]
			for a in lineDict["PASS_ALLELES"]:
				if not lineDict["ALL_ALLELES"][a] in lineDict:
					print lineDict
					alleleDict=lineDict[lineDict["ALL_ALLELES"][a]]
					exit()
				alleleDict=lineDict[lineDict["ALL_ALLELES"][a]]
				out_file_tsv.write(str(len(alleleDict["ALLELE_IDs"]))+"\t")
				if "HOM_SAMPLES" in alleleDict:
					out_file_tsv.write(str(len(alleleDict["HOM_SAMPLES"]))+"\t")
				else:
					out_file_tsv.write("0\t")
				if "HET_SAMPLES" in alleleDict:
					out_file_tsv.write(str(len(alleleDict["HET_SAMPLES"]))+"\t")
				else:
					out_file_tsv.write("0\t")
				out_file_tsv.write(str(lineDict["BI_ALLELIC"])+"\t"+lineDict["CHROM"]+"\t"+lineDict["POS"]+"\t"+lineDict["REF"]+"\t"+lineDict["ALL_ALLELES"][a]+"\t"+str(alleleDict["FRQ"][0])+"\t")
				out_file_tsv.write(str(lineDict["AN"])+"\t"+str(alleleDict["AC"])+"\t")
				for v in outDict["VEP_FIELDS"]:
					if "VEP" in alleleDict:
						if v in alleleDict["VEP"]:
							out_file_tsv.write(alleleDict["VEP"][v]+"\t")
						else:
							out_file_tsv.write(".\t")
					else:
						out_file_tsv.write(".\t")

				for c in outDict["CLINVAR_FIELDS"]:
					if "ClinVar" in alleleDict:
						if c in alleleDict["ClinVar"]:
							out_file_tsv.write(alleleDict["ClinVar"][c]+"\t")
						else:
							out_file_tsv.write(".\t")
					else:
						out_file_tsv.write(".\t")
				if "ID" in lineDict:
					out_file_tsv.write(lineDict["ID"]+"\t")
				else:
					out_file_tsv.write(".\t")
				if "HGMD_CLASS" in lineDict:
					out_file_tsv.write(lineDict["HGMD_CLASS"]+"\t")
				else:
					out_file_tsv.write(".\t")
				if "HGMD_DB" in lineDict:
					out_file_tsv.write(lineDict["HGMD_DB"]+"\t")
				else:
					out_file_tsv.write(".\t")
				if "HGMD_GENE" in lineDict:
					out_file_tsv.write(lineDict["HGMD_GENE"]+"\t")
				else:
					out_file_tsv.write(".\t")
				if "HGMD_PHEN" in lineDict:
					out_file_tsv.write(lineDict["HGMD_PHEN"]+"\t")
				else:
					out_file_tsv.write(".\t")
				if "HGMD_PROT" in lineDict:
					out_file_tsv.write(lineDict["HGMD_PROT"]+"\t")
				else:
					out_file_tsv.write(".\t")
				if "HGMD_STRAND" in lineDict:
					out_file_tsv.write(lineDict["HGMD_STRAND"]+"\t")
				else:
					out_file_tsv.write(".\t")
				out_file_tsv.write(",".join(alleleDict["ALLELE_IDs"])+"\t")
				if "HOM_SAMPLES" in alleleDict:
					out_file_tsv.write(",".join(alleleDict["HOM_SAMPLES"])+"\t")
				else:
					out_file_tsv.write(".\t")
				if "HET_SAMPLES" in alleleDict:
					out_file_tsv.write(",".join(alleleDict["HET_SAMPLES"])+"\n")
				else:
					out_file_tsv.write(".\n")

	out_file_tsv.write("UniqueSamples:\t"+str(len(outDict["UNIQ_IDs"]))+"\n")
	out_file_tsv.close()
	return


if __name__ == '__main__':
	main()
	exit()
