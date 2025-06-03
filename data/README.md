Original VCFs are provided by Severus people ([variant_calls_and_benchmarks.tar.gz](https://zenodo.org/records/14541057)). To run truvari bench, we had to remove BND from some callsets (sniffles2, cuteSV, and severus):
```
bcftools view -Oz --exclude "INFO/SVTYPE='BND'" input.vcf.gz > output.vcf.gz
tabix -p vcf output.vcf.gz
```
