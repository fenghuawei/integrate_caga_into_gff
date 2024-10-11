import os
f = open("GCF_000004515.6_Glycine_max_v4.0_genomic_chr-name-new.gff3",'r')

# sort the TSS file
cmd = "python sort_TSS_file.py Tag_clusters_all.txt  seqnames  dominant_ctss.start  Tag_clusters_all_sort.txt"
os.system(cmd)




f2 = open("Tag_clusters_all_sort.txt",'r')
f2.readline()
line_caga = f2.readline()
m = f.readlines()
n = f2.readlines()
gff_chr = []
TSS_chr = []


#extract the gff chr list
for line in m:
    if line.startswith("#"):
        continue
    else:
        line_ = line.strip("\n").split("\t")
        chr_value = line_[0]
        if not(chr_value in gff_chr):
            gff_chr.append(chr_value)
#extract the TSS chr list
for line in n:
    line_ = line.strip("\n").split("\t")
    chr_value = line_[0]
    if not(chr_value in TSS_chr):
        TSS_chr.append(chr_value)
f.close()
f2.close()


#split the gff  file and TSS file into many files by chr
def generate_file(gff,TSS,key):
    gff_id = key + 'gff.txt'
    TSS_id = key + 'TSS.txt'
    bajie = open(gff_id,'w')
    shaseng = open(TSS_id,'w')
    for line in m:
        line_ = line.strip('\n').split("\t")
        if line_[0] == key:
            bajie.write(line)
        else:
            continue
    bajie.close()
    for line in n:
        line_ = line.strip('\n').split("\t")
        if line_[0] == key:
            shaseng.write(line)
        else:
            continue
    shaseng.close()


for key in TSS_chr:
    file_gff_name = key + 'gff.txt'
    file_TSS_name = key + 'TSS.txt'
    file_gff_new_name = key + 'gff_new.txt'
    generate_file(file_gff_name,file_TSS_name,key)
    # change the start pos by caga data by chr
    cmd = "python sort_gff.py  " + file_gff_name + "  "+  file_gff_new_name + "  " + file_TSS_name
    #python  sort_gff.py  chr3gff.txt  chr3gff_new.txt  chr3TSS.txt
    print(cmd)
    os.system(cmd)
f1 = open("GCF_000004515.6_Glycine_max_v4.0_genomic_chr-name-new_with_CGAR.gff3",'w')

# combine the generated chr file to the final new gff file
def add_gff(new_file):
    wukong = open(new_file,"r")
    for line in wukong:
        f1.write(line)
    wukong.close()
niu = []
for line in m:
    if line.startswith("#"):
        f1.write(line)
    else:

        line_ = line.strip("\n").split("\t")
        if line_[0] in TSS_chr:
            if not(line_[0] in niu):
                file_gff_new_name = line_[0] + 'gff_new.txt'
                add_gff(file_gff_new_name)
                niu.append(line_[0])
            else:
                continue
        else:
            f1.write(line)



