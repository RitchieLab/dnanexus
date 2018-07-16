#!/usr/bin/env python
"""Script to reformat the ClinVar TSV file into two VCF files, one for b37 and
one for b38.
"""
__author__ = "Thomas Nate Person"
__license__ = "GPLv3"
__email__ = "thomas.n.person@gmail.com"

import re, sys, time, gzip, csv


def write_vcf(assembly, new_vcf, vcf_csv_writer):

    for site_id, site in new_vcf[assembly].items():

        outsite = site_id.split("\t") + [",".join(site[site_id]), ".", "."]
        site.pop(site_id, None)
        line = outsite + [new_vcf["source"] + ",".join(site.values())]
        vcf_csv_writer.writerows([line])


dbFile = gzip.open(sys.argv[1], "r")
with open(sys.argv[1].replace("txt.gz", "") + "b37.vcf", "w") as db_b37_vcf:
    with open(sys.argv[1].replace("txt.gz", "") + "b38.vcf", "w") as db_b38_vcf:

        # grab header fields
        header = dbFile.next()
        headerFields = "|".join(header.lstrip("#").strip().split("\t"))
        headerFields = "Allele|" + headerFields

        # get dates
        MonthYear = time.strftime("%b%Y")
        DayMonthYear = time.strftime("%d%b%Y")

        # initialize csv writers
        db_b37_vcf_writer = csv.writer(
            db_b37_vcf,
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            quotechar="",
            escapechar="\\",
            lineterminator="\n",
        )
        db_b38_vcf_writer = csv.writer(
            db_b38_vcf,
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            quotechar="",
            escapechar="\\",
            lineterminator="\n",
        )

        # write b37 header
        db_b37_vcf_writer.writerow(["##fileformat=VCFv4.2"])
        db_b37_vcf_writer.writerow(
            [
                '##INFO=<ID=ClinVar.TSV.{},Number=.,Type=String,Description="ClinVar.TSV.{}:'.format(
                    DayMonthYear, DayMonthYear
                )
                + " '{}'\">".format(headerFields)
            ]
        )
        db_b37_vcf_writer.writerow(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        )

        # write b38 header
        db_b38_vcf_writer.writerow(["##fileformat=VCFv4.2"])
        db_b38_vcf_writer.writerow(
            [
                '##INFO=<ID=ClinVar.TSV.{},Number=.,Type=String,Description="ClinVar.TSV.{}:'.format(
                    DayMonthYear, DayMonthYear
                )
                + " '{}'\">".format(headerFields)
            ]
        )
        db_b38_vcf_writer.writerow(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        )

        CHROM, POS, REF, ALT, ASSEMBLY = -1, -1, -1, -1, -1
        new_vcf, new_vcf["GRCh37"], new_vcf["GRCh38"] = {}, {}, {}
        new_vcf["source"] = "ClinVar.TSV.{}=".format(DayMonthYear)
        dbFile.seek(0)
        header = dbFile.next()
        fields = [x.strip() for x in header.split("\t")]
        for i, f in enumerate(fields):

            if f == "Chromosome":
                CHROM = i
            elif f == "Start":
                POS = i
            elif f == "ReferenceAllele":
                REF = i
            elif f == "AlternateAllele":
                ALT = i
            elif f == "Assembly":
                ASSEMBLY = i

        for line in dbFile:

            fields = [x.strip() for x in line.split("\t")]
            if fields[REF] == fields[ALT]:
                continue
            elif fields[ASSEMBLY] == "NCBI36":
                continue

            if bool(re.search("^[AGCT]+$", fields[REF])) and bool(
                re.search("^[AGCT]+$", fields[ALT])
            ):

                site_id = fields[CHROM] + "\t" + fields[POS] + "\t.\t" + fields[REF]
                variant_id = (
                    fields[CHROM]
                    + "\t"
                    + fields[POS]
                    + "\t"
                    + fields[REF]
                    + "\t"
                    + fields[ALT]
                )

                INFO = fields[ALT] + "|" + "|".join(fields)
                INFO = INFO.replace(";", "/")
                INFO = INFO.replace(",", "/")
                INFO = INFO.replace(" ", "_")

                if site_id in new_vcf[fields[ASSEMBLY]]:

                    site = new_vcf[fields[ASSEMBLY]][site_id]
                    site[site_id].add(fields[ALT])
                    if variant_id in site:
                        site[variant_id] = site[variant_id] + "," + INFO
                    else:
                        site[variant_id] = INFO

                    new_vcf[fields[ASSEMBLY]][site_id] = site

                else:

                    new_vcf[fields[ASSEMBLY]][site_id] = {
                        variant_id: INFO,
                        site_id: {fields[ALT]},
                    }

        write_vcf("GRCh37", new_vcf, db_b37_vcf_writer)
        write_vcf("GRCh38", new_vcf, db_b38_vcf_writer)


exit()
