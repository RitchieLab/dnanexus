#!/usr/bin/python

import argparse, gzip, re

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

def parse_INFO(line):
	fields = line.split(':')[-1].replace('"',"").replace("'","").replace("'","").replace(">","").strip().split("|")
	if line.startswith("##INFO=<ID=CSQ"):
		VEP_Fields=dict()
		VEP_Fields["MAF"]=set()
		for i,f in enumerate(fields):
			if f=="IMPACT":
				VEP_Fields["IMPACT"]=i
			elif f=="SYMBOL":
				VEP_Fields["SYMBOL"]=i
			elif f=="Feature":
				VEP_Fields["Feature"]=i
			elif f=="CANONICAL":
				VEP_Fields["CANONICAL"]=i
			elif f.endswith("MAF"):
				VEP_Fields["MAF"].add(i)
		return VEP_Fields
	elif line.startswith("##INFO=<ID=ClinVar"):
		ClinVar_fields=dict()
		for i,f in enumerate(fields):
			if f=="ClinicalSignificance":
				ClinVar_fields["ClinicalSignificance"]=i
			elif f=="ReviewStatus":
				ClinVar_fields["ReviewStatus"]=i
		return ClinVar_fields

def filterLine(line,cli_arguments,VEP_Fields,ClinVar_fields,geneSet):
	info_field=[]
	if cli_arguments.b_snp:
		fields = line.split('\t')
		pattern = re.compile("^[AGCT]$")
		if "," in fields[4]:
			return False
		elif  pattern.match(fields[4]) and pattern.match(fields[3]):
			info_field = fields[7].split(";")
		else:
			return False
	else:
		info_field = line.split('\t')[7].split(";")
	geneList=[]
	frequency=[]
	annotations=[]
	for field in info_field:
		if field.startswith("AF"):
			AFs = field.replace("AF=","").split(",")
			for f in AFs:
				if is_number(f):
					if f=="0":
						frequency.append(False)
					elif float(f)>cli_arguments.frequency:
						frequency.append(False)
					else:
						frequency.append(True)
				else:
					return False
			if len(frequency)>0:
				if not any(frequency):
					return False
		elif field.startswith("CSQ"):
			if cli_arguments.VEP_Impact is not None:
				CSQs =[x.strip() for x in field.replace("CSQ=","").strip().split(",")]
				for i,CSQ in enumerate(CSQs):
					Cs=CSQ.split("|")
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
										if pubMAFs[a]=="0":
											frequency.append(False)
										elif float(pubMAFs[a])>cli_arguments.frequency:
											frequency.append(False)
										else:
											frequency.append(True)
							elif is_number(Cs[f_field]):
								if float(Cs[f_field])>cli_arguments.frequency:
									frequency.append(False)
								else:
									frequency.append(True)
					if len(frequency)>0:
						if not any(frequency):
							return False
					if cli_arguments.Ensembl:
						if not Cs[VEP_Fields["Feature"]].startswith("E"):
							continue
						if cli_arguments.Cannonical and Cs[VEP_Fields["CANONICAL"]]=="YES":
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
						else:
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
					if cli_arguments.RefSeq:
						if Cs[VEP_Fields["Feature"]].startswith("E"):
							continue
						if cli_arguments.Cannonical and Cs[VEP_Fields["CANONICAL"]]=="YES":
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
						else:
							if getimpact_level(Cs[VEP_Fields["IMPACT"]])>=cli_arguments.VEP_Impact:
								annotations.append(True)
							else:
								annotations.append(False)
							if cli_arguments.gene_list is not None:
								if Cs[VEP_Fields["SYMBOL"]] in geneSet:
									geneList.append(True)
					if cli_arguments.gene_list is not None:
						if not any(geneList):
							return False
		elif field.startswith("ClinVar.TSV.Jan2017="):
			if cli_arguments.clinVarStar is not None:
				clinvar = field.replace("ClinVar.TSV.Jan2017=","").split("|")
				if len(clinvar)<ClinVar_fields["ReviewStatus"]:
					annotations.append(False)
				elif getStar(clinvar[ClinVar_fields["ReviewStatus"]])>=cli_arguments.clinVarStar:
					if getClinVar_level(clinvar[ClinVar_fields["ClinicalSignificance"]])>=cli_arguments.ClinicalSignificance:
						annotations.append(True)
					else:
						annotations.append(False)
				else:
					annotations.append(False)
		elif field.startswith("CLASS="):
			if cli_arguments.HGMD is not None:
				if getHGMD_level(field.replace("CLASS=",""))>=cli_arguments.HGMD:
					annotations.append(True)
				else:
					annotations.append(False)
	if not any(annotations):
		return False
	elif cli_arguments.gene_list is not None:
		if not any(geneList):
			return False
	return True


def main():

	cli_arguments =  parseArguments()

	if not cli_arguments.vcf_file.endswith("vcf.gz"):
		print "Not a bgziped vcf file"
		return
	geneSet=set()
	if cli_arguments.gene_list is not None:
		print cli_arguments.gene_list
		geneList = open(cli_arguments.gene_list,'r')
		for line in geneList:
			line=line.strip().upper()
			if line=="":
				continue
			else:
				geneSet.add(line)

	vcf_file = gzip.open(str(cli_arguments.vcf_file),'r')
	VEP_Fields=dict()
	ClinVar_fields=dict()
	IDs=dict()

	for line in vcf_file:
		if line.startswith("#"):
			print line.strip()
			if line.startswith("##INFO=<ID=CSQ"):
				VEP_Fields=parse_INFO(line.strip())
			elif line.startswith("##INFO=<ID=ClinVar"):
				ClinVar_fields=parse_INFO(line.strip())
			elif line.startswith("#CHROM"):
				fields = [x.strip() for x in line.strip().split("\t")]
				for i,f in enumerate(fields):
					IDs[i]=f
		elif filterLine(line,cli_arguments,VEP_Fields,ClinVar_fields,geneSet):
			print line.strip()
			pass
	vcf_file.close()
	return

if __name__ == '__main__':
	main()
	exit()
