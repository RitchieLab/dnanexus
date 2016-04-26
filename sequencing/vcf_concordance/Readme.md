<!-- dx-header -->
# VCF Concordance (DNAnexus Platform App)

Generate concordance metrics from VCF file(s) using [GATK GenotypeConcordance](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeConcordance.php).  Not recommended for large files.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->

This app generates a variety of concordance metrics from the given VCF file(s).  It can be run in one of 3 ways, outlined as follows:

----

* If 2 VCF files are provided, will calculate the concordance metrics of all IDs in common between the two files.
* If one VCF file is provided along with a mapping of sample ID -> subject ID, this will calculate the internal consistency of the VCF file for all subjects with more than one sample.  In the case of three or more samples per subject, two samples will be chosen at random to be compared.  Note that the sample in the "truth" set and the sample in the "eval" set is random.
* If only one VCF file is provided without a subject-sample mapping file, the app will assume the mapping scheme of "GHS\_PT[0-9]\*\_[0-9]\*" -> "PT[0-9]\*" and proceed as above.

----
The app defines concordance according to the GATK GenotypeConcordance tool, which generates a concordance matrix with the following definitions:

----

* Homozygous Referent (HR)

  0/0
* Heterozygous (HET)

  0/1
* Homozygous Variant (HV)

  1/1

* No Call (NC)

  ./. OR filtered at genotype level (i.e., 0/1:BadCall)
* Unavailable (U)

  Site not present in VCF OR filtered at variant level (i.e., QDFilter)
  
----

This results in the following table of values, numbered here for reference when defining the specific concordance metrics below.

<!--

<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;}
.tg .tg-s6z2{text-align:center}
.tg .tg-8o5d{background-color:#34ff34}
.tg .tg-8e65{background-color:#f8a102}
.tg .tg-5fb6{background-color:#fe0000}
.tg .tg-v88f{background-color:#f8ff00}
</style>

<table class="tg">
  <tr>
    <th class="tg-031e" colspan="2" rowspan="2"></th>
    <th class="tg-s6z2" colspan="5">Truth ("comp") set<br></th>
  </tr>
  <tr>
    <td class="tg-031e">HR</td>
    <td class="tg-031e">HET</td>
    <td class="tg-031e">HV</td>
    <td class="tg-031e">NC</td>
    <td class="tg-031e">U</td>
  </tr>
  <tr>
    <td class="tg-s6z2" rowspan="5">Test<br>("eval")<br>set<br></td>
    <td class="tg-031e">HR<br></td>
    <td class="tg-8o5d">1</td>
    <td class="tg-8e65">2</td>
    <td class="tg-5fb6">3</td>
    <td class="tg-v88f">4</td>
    <td class="tg-v88f">5</td>
  </tr>
  <tr>
    <td class="tg-031e">HET</td>
    <td class="tg-8e65">6</td>
    <td class="tg-8o5d">7</td>
    <td class="tg-8e65">8</td>
    <td class="tg-v88f">9</td>
    <td class="tg-v88f">10</td>
  </tr>
  <tr>
    <td class="tg-031e">HV</td>
    <td class="tg-5fb6">11</td>
    <td class="tg-8e65">12</td>
    <td class="tg-8o5d">13</td>
    <td class="tg-v88f">14</td>
    <td class="tg-v88f">15</td>
  </tr>
  <tr>
    <td class="tg-031e">NC</td>
    <td class="tg-v88f">16</td>
    <td class="tg-v88f">17</td>
    <td class="tg-v88f">18</td>
    <td class="tg-8o5d">19</td>
    <td class="tg-8o5d">20</td>
  </tr>
  <tr>
    <td class="tg-031e">U</td>
    <td class="tg-v88f">21</td>
    <td class="tg-v88f">22</td>
    <td class="tg-v88f">23</td>
    <td class="tg-8o5d">24</td>
    <td class="tg-8o5d">25</td>
  </tr>
</table>
-->

|                   |     | Truth ("comp") set |     |    |    |    |
|-------------------|-----|:------------------:|-----|----|----|----|
|                   |     | HR                 | HET | HV | NC | U  |
| Test ("eval") set | HR  | 1                  | 2   | 3  | 4  | 5  |
|                   | HET | 6                  | 7   | 8  | 9  | 10 |
|                   | HV  | 11                 | 12  | 13 | 14 | 15 |
|                   | NC  | 16                 | 17  | 18 | 19 | 20 |
|                   | U   | 21                 | 22  | 23 | 24 | 25 |


The output of this app is the GATK output for each of the following comparisons:

----
* raw, all sites
* filtered, all sites
* raw, SNPs only 
* filtered, SNPs only
* raw, INDELs only
* filtered, INDELs only
----

In the case of the raw concordance, any annotated filters in either the FILTER field or the FT genotype annotation are completely ignored when calculating the concordance metrics.  In the case fo filtered concordance, these genotype calls are set to "Unavailable" or "No Call", as described above.

In addition to the GATK output, we also provide a summary of the concordance for each of the above comparisons, both overall and per subject.  Our summary includes the following metrics, defined by the calculations using the cell numebrs in the above table:

----

* Non-Reference Sensitivity (NRS)

  (7+8+12+13) / (2+3+7+8+12+13+17+18+22+23)
* Non-Reference Sensitivity, Reversed (NRS-R)

  (7+8+12+13) / (6+7+8+9+10+11+12+13+14+15)
* Overall Genotype Concordance (OGC)

  AKA: Concordance, where both called
  
  (1+7+13) / (1+2+3+6+7+8+11+12+13)
* Overall Genotype Concordance w/ Missing (OGCM)

  (1+7+13+19) / (1+2+3+4+6+7+8+9+11+12+13+14+16+17+18+19)
* Overall genotype Concordance w/ Missing and Unavailable (OGCMU)

  AKA: Concordance, missing as discordant

  (1+7+13+19+20+24+25) / SUM(1:25)
* Heteroxygous Concordance (HHC)

  7 / (2+6+7+8+12)
* Heterozygous Concordance w/ Missing (HHC-M)

  7 / (2+6+7+8+9+10+12+17+22)
* Non-Reference Concordance (NRC)

  AKA: Precision
  
  (7+13)/(6+7+8+11+12+13)
* Non-Reference Concordance, Reversed (NRC-R)

  (7+13)/(2+3+7+8+12+13)
* Recall (REC)

  (7+13)/(2+3+7+8+12+13+17+18+22+23)
* Recall, Reversed (REC-R)

  (7+13)/(6+7+8+9+10+11+12+13+14+15)

----

<!--
TODO: This app directory was automatically generated by dx-app-wizard;
please edit this Readme.md file to include essential documentation about
your app that would be helpful to users. (Also see the
Readme.developer.md.) Once you're done, you can remove these TODO
comments.

For more info, see https://wiki.dnanexus.com/Developer-Portal.
-->
