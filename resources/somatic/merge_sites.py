import sys, collections
WD, FAI = sys.argv[1], sys.argv[2]
sites = collections.OrderedDict()
for line in open(f"{WD}/sites_raw.tsv"):
    c, p, r, a, d = line.rstrip("\n").split("\t")
    sites.setdefault((c, int(p), r, a), set()).add(d)
order = {ln.split("\t")[0]: i for i, ln in enumerate(open(FAI))}
keys = sorted(sites, key=lambda k: (order.get(k[0], 10**9), k[1]))
contigs = [ln.split("\t")[:2] for ln in open(FAI)]
edit = {("A", "G"), ("T", "C")}
with open(f"{WD}/merged_somatic.vcf", "w") as vcf, open(f"{WD}/site_donor_map.tsv", "w") as mp:
    vcf.write("##fileformat=VCFv4.2\n")
    for cn, cl in contigs:
        vcf.write(f"##contig=<ID={cn},length={cl}>\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    mp.write("chrom\tpos\tref\talt\tdonors\tediting_prone\n")
    for (c, p, r, a) in keys:
        vcf.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n")
        mp.write(f"{c}\t{p}\t{r}\t{a}\t{','.join(sorted(sites[(c,p,r,a)]))}\t{int((r,a) in edit)}\n")
edit_n = sum(1 for k in keys if (k[2], k[3]) in edit)
shared = sum(1 for k in keys if len(sites[k]) > 1)
print(f"union somatic sites: {len(keys)}  (A>G/T>C editing-prone: {edit_n}, {edit_n/len(keys):.1%}; shared by >1 donor: {shared})")
