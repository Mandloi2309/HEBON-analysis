touch all_files.txt;
for i in {1..22}
do 
7z x chr_${i}.zip -p "";
echo "imputed_chr${i}.bed imputed_chr${i}.bim imputed_chr${i}.fam" >> all_files.txt;
rm chr_${i}.zip;
bcftools norm -d both chr${i}.dose.vcf.gz -O z --threads 4 -o chr${i}.vcf.gz;
gunzip chr${i}.info.gz;
rm chr${i}.dose.vcf.gz ;
python3 DelDups.py chr${i}.info ${i};
../../plink --vcf chr${i}.vcf.gz --exclude snps-excluded.txt --make-bed --out tempchr;
../../plink --bfile tempchr --list-duplicate-vars;
../../plink --bfile tempchr --exclude plink.dupvar --make-bed --out tempchr2;
../../plink --bfile tempchr2 --qual-scores chr${i}_p.info 7 1 '#'  --qual-threshold 0.8 --make-bed --out tempchr3;
../../plink --bfile tempchr3 --hwe 1e-6 --maf 0.01 --mind 0.98 --make-bed --out imputed_chr${i};
rm tempchr.bed tempchr.bim tempchr.fam tempchr2.bed tempchr2.bim tempchr2.fam tempchr3.bed tempchr3.bim tempchr3.fam ;
done;

filename='HEBON';
../../plink --bfile imputed_chr1 --merge-list all_files.txt --make-bed --out $filename ;
../../plink --bfile $filename --list-duplicate-vars;
../../plink --bfile $filename --exclude plink.dupvar --make-bed --out $filename.imps;
../../plink --bfile $filename.imps --update-sex $filename.clinical.txt --make-bed --out temp;
../../plink --bfile temp --pheno $filename.clinical.txt --mpheno 2 --make-bed --out $filename.HRC.merged;
rm temp.bim temp.bed temp.fam $filename.imps.bim $filename.imps.bed $filename.imps.fam;

## Check gender and phenotype once
../plink --bfile $filename.HRC.merged --exclude plink.dupvar  --hwe 1e-6 --maf 0.01 --mind 0.98 --make-bed --out $filename.merged;

sumFile="$filename.Sum"
## Summary 
../../plink --bfile $filename.merged --assoc counts --out $sumFile;

R
Sum = read.table("HEBON.Sum.assoc", header = TRUE, row.names = NULL);
Sum$row.names = NULL;
write.table(Sum, file = "HEBON.merged.sumstats", sep='\t', quote = FALSE, row.names = FALSE);
quit();
n;

cut -f2 HEBON.merged.sumstats|tail -n+2|awk '{n=split($0,a,/[:_]/); print a[1]"\t"a[2]"\t.\t"a[3]"\t"a[4]"\t.\t.\t."}'|java -jar ../snpEff/SnpSift.jar annotate -id ../00-All.vcf.gz|cut -f 1,2,3,4,5 > rsId_list.txt;

R
RSIDS = read.table("rsId_list.txt", sep='\t');
RSIDS = `colnames<-`(RSIDS, c('CHR','BP','DBSNP','A1','A2'));
RSIDS$SNPID = paste(RSIDS$CHR,RSIDS$BP,RSIDS$A1, RSIDS$A2, sep=':');
write.table(RSIDS, file = "HEBON.RSID.txt", sep='\t', quote = FALSE, row.names = FALSE);
quit();
n;

R
Sum = read.table("HEBON.sumstats", header = TRUE, row.names = NULL);
Sum$row.names = NULL;
RSIDS = read.table("rsId_list.txt", sep='\t');
RSIDS = `colnames<-`(RSIDS, c('CHR','BP','DBSNP','A1','A2'));
RSIDS$SNPID = paste(RSIDS$CHR,RSIDS$BP,RSIDS$A1, RSIDS$A2, sep=':');
newdf = merge(x = Sum, y = RSIDS, by.x = "SNP", by.y="SNPID");
newdf$DBSNP = gsub(";.*", "", newdf$DBSNP);
newdf[grep("\\.", newdf$DBSNP), 13] = newdf[grep("\\.", newdf$DBSNP), 1];
newdf = newdf[, c(2,13,3,4,7,5,6,8,9,10,1)];
newdf = `colnames<-`(newdf, c('CHR','SNP', 'BP','A1','A2','C_A', 'C_U', 'CHISQ', 'P', 'OR', 'name'));
newdf = newdf[with(newdf, order(CHR, BP)), ];
write.table(newdf, file = "Meta.1KG.RSID.sumstats", sep='\t', quote = FALSE, row.names = FALSE);
library(qqman);
tiff(filename = "Meta_Manhattan.tiff", width = 3200, height = 2400, units = "px", pointsize = 75);
manhattan(Sum, chr = "CHR",bp = "BP",p = "P", logp = T, ylim= c(0,30), , cex = 0.45, cex.axis = 0.45);
dev.off();
quit();
n;

conda activate imlabtools;
./Predict.py --model_db_path ../GTEx-V6p-1KG/TW_Breast_Mammary_Tissue_0.5_1KG.db --vcf_genotypes ../../Merged/1KG/HO.RS.Xcan.vcf --vcf_mode genotyped --prediction_output ../../Merged/1KG/HO.predict.Breast.txt --prediction_summary_output ../../Merged/1KG/HO.predict.txt;
./PrediXcanAssociation.py --hdf5_expression_file ../../Merged/1KG/HO.predict.Breast.h5 --mode linear --output ../../Merged/1KG/GTEx.HO.txt --input_phenos_file ../../Merged/1KG/HO.clinical.txt --input_phenos_column 'Pheno';

### PRS
Rscript PRSice/PRSice.R --dir ../. --prsice PRSice/PRSice_linux --base mmc3.txt.gz --snp 0 --chr 1 --bp 2 --A1 3 --A2 4 --stat 11 --pvalue 12 --index --target HEBON/HEBON.HRC.merged --pheno HEBON/CLIN.txt --perm 10000 --thread 4 --beta --binary-target F --out HEBON_PRS_313;