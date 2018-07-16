#!/usr/bin/env python
"""Script to reformat the ClinVar TSV file into two VCF files, one for b37 and
one for b38.
"""
__author__ = "Thomas Nate Person"
__license__ = "GPLv3"
__email__ = "thomas.n.person@gmail.com"

import re, sys, time, gzip
from string import maketrans, translate


def write_vcf_fields(assembly, new_vcf, file_obj):
    """ writes non-header lines to a VCF file """
    for site_id, site in new_vcf[assembly].items():

        outsite = site_id + ",".join(site[site_id]) + "\t.\t.\t" + new_vcf["source"]
        site.pop(site_id, None)
        outsite = outsite + ",".join(site.values()) + "\n"
        file_obj.write(outsite)


def write_vcf_header(DayMonthYear, headerFields, file_obj):
    """ writes header lines to a VCF file """
    file_obj.write("##fileformat=VCFv4.2\n")
    file_obj.write(
        "##INFO=<ID=ClinVar.TSV."
        + DayMonthYear
        + ',Number=.,Type=String,Description="ClinVar.TSV.'
        + DayMonthYear
        + ": '"
        + headerFields
        + "'\">\n"
    )
    file_obj.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def main():

    # intialize database file containing ClinVar variants
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

            # write header fields
            write_vcf_header(DayMonthYear, headerFields, db_b37_vcf)
            write_vcf_header(DayMonthYear, headerFields, db_b38_vcf)

            CHROM, POS, REF, ALT, ASSEMBLY = -1, -1, -1, -1, -1
            new_vcf, new_vcf["GRCh37"], new_vcf["GRCh38"] = {}, {}, {}
            new_vcf["source"] = "ClinVar.TSV.{}=".format(DayMonthYear)
            dbFile.seek(0)
            header = dbFile.next()
            fields = [x.strip() for x in header.split("\t")]
            info_trans_table = maketrans(";, ", "//_")
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

                    site_id = (
                        fields[CHROM]
                        + "\t"
                        + fields[POS]
                        + "\t.\t"
                        + fields[REF]
                        + "\t"
                    )
                    variant_id = (
                        fields[CHROM]
                        + "\t"
                        + fields[POS]
                        + "\t"
                        + fields[REF]
                        + "\t"
                        + fields[ALT]
                    )

                    # translate info fiels with trans table
                    INFO = fields[ALT] + "|" + "|".join(fields)
                    INFO = translate(INFO, info_trans_table)

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

            # write non-header VCF fields
            write_vcf_fields("GRCh37", new_vcf, db_b37_vcf)
            write_vcf_fields("GRCh38", new_vcf, db_b38_vcf)


if __name__ == "__main__":
    main()
