######### Quality Control #################;
# Investigate missingness per individual and per SNP and make histograms.;
..\plink.exe --bfile HEBON --missing --geno 0.02 --mind 0.02 --check-sex --make-bed --out temp1;

############################################################;
# Delete individuals with sex discrepancy.;
grep "PROBLEM" temp1.sexcheck| gawk '{print$1,$2}'> sex_discrepancy.txt;
# This command generates a list of individuals with the status “PROBLEM”.;
..\plink.exe --bfile temp1 --remove sex_discrepancy.txt --make-bed --out temp2 ;

############################################################
# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).;
gawk "{ if ($1 >= 1 && $1 <= 22) print $2 }" temp2.bim > snp_1_22.txt;
..\plink.exe --bfile temp2 --extract snp_1_22.txt --make-bed --out temp3;

..\plink.exe --bfile temp3 --maf 0.05 --make-bed --out temp4;

############################################################
# By default the --hwe option in plink only filters for controls.
# Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data.
..\plink.exe --bfile temp4 --hwe 1e-6 --make-bed --out hwe_filter_step1;

# The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE. 
# This second HWE step only focusses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed
..\plink.exe --bfile hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out temp5;

############################################################
# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions (inversion.txt [High LD regions]) and prune the SNPs using the command --indep-pairwise’.
# The parameters ‘50 5 0.2’ stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.

..\plink.exe --bfile temp5 --indep-pairwise 50 5 0.2 --out indepSNP;
..\plink.exe --bfile temp5 --extract indepSNP.prune.in --het --out R_check;
# This file contains your pruned data set.

# Plot of the heterozygosity rate distribution
G:\R\R-4.0.4\bin\Rscript.exe --no-save check_heterozygosity_rate.R;

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
# For data manipulation we recommend using UNIX. However, when performing statistical calculations R might be more convenient, hence the use of the Rscript for this step:
G:\R\R-4.0.4\bin\Rscript.exe --no-save heterozygosity_outliers_list.R;

# Output of the command above: fail-het-qc.txt .
# When using our example data/the HapMap data this list contains 2 individuals (i.e., two individuals have a heterozygosity rate deviating more than 3 SD's from the mean).
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | gawk '{print$1, $2}'> het_fail_ind.txt;

# Remove heterozygosity rate outliers.
..\plink.exe --bfile temp5 --remove het_fail_ind.txt --make-bed --out temp6 ;

############################################################
# It is essential to check datasets you analyse for cryptic relatedness.
# Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2 in this tutorial.

# Check for relationships between individuals with a pihat > 0.2.
..\plink.exe --bfile temp6 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2;

# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome;

# The generated plots show a considerable amount of related individuals (explentation plot; PO = parent-offspring, UN = unrelated individuals) in the Hapmap data, this is expected since the dataset was constructed as such.
# Normally, family based data should be analyzed using specific family based methods. In this tutorial, for demonstrative purposes, we treat the relatedness as cryptic relatedness in a random population sample.
# In this tutorial, we aim to remove all 'relatedness' from our dataset.
# To demonstrate that the majority of the relatedness was due to parent-offspring we only include founders (individuals without parents in the dataset).

..\plink.exe --bfile temp6 --filter-founders --make-bed --out temp7;

############################################################
# Now we will look again for individuals with a pihat >0.2.
..\plink.exe --bfile temp7 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders;
# The file 'pihat_min0.2_in_founders.genome' shows that, after exclusion of all non-founders, only 1 individual pair with a pihat greater than 0.2 remains in the HapMap data.
# This is likely to be a full sib or DZ twin pair based on the Z values. Noteworthy, they were not given the same family identity (FID) in the HapMap data.

# For each pair of 'related' individuals with a pihat > 0.2, we recommend to remove the individual with the lowest call rate. 
..\plink.exe --bfile temp7 --missing;
# Use an UNIX text editor (e.g., vi(m) ) to check which individual has the highest call rate in the 'related pair'. 

# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
..\plink.exe --bfile temp7 --remove 0.2_low_call_rate_pihat.txt --make-bed --out QC_HEBON;

rm temp1.bed temp1.bim temp1.fam temp2.bed temp2.bim temp2.fam ;
rm temp3.bed temp3.bim temp3.fam temp4.bed temp4.bim temp4.fam;
rm temp5.bed temp5.bim temp5.fam temp6.bed temp6.bim temp6.fam;
..\plink.exe --bfile QC_HEBON --update-map ..\Illumina_GSA_rsids.txt --update-name --make-bed --out HEB.QC.mapped;

######### Alignment #################
# Home DIR
filename="HEBON.QC.mapped";
path="HEBON/";

java -Xmx40g -jar GenotypeHarmonizer/GenotypeHarmonizer.jar --inputType PLINK_BED --input "$path$filename" --update-id --outputType PLINK_BED --output "$path"alignment/all_chrs --refType PLINK_BED --ref 1000G;

## Linux
cd "$path";
mkdir unp;
for i in {1..22} do ../plink --bfile alignment/all_chrs --chr ${i} --recode vcf --out unp/unphased_chr${i};
bcftools sort unp/unphased_chr${i}.vcf -Oz -o unp/sorted_chr${i}.vcf.gz;
rm unp/unphased_chr${i}.vcf;
done