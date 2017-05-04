import os
import sys
import subprocess
#print os.environ["filter_options"]

subprocess.check_call("perl /usr/share/ensembl-tools-release-87/scripts/variant_effect_predictor/filter_vep.pl --gz --format vcf -i " + sys.argv[1] +" " + os.environ["filter_options"], shell=True)
