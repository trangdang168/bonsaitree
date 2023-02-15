# always use the same pedigree and genotyped files!
AMISH_GENOTYPED_PEDIGREE = "/homes/thdang/trang_amish/ped-aware_pipeline/output/amish1_chr21_amr_geno.ped"
AMISH_FULL_PEDIGREE = "/homes/stan1/sam_amish/amish_sim_pipeline/input/extended_pedigree_final.fam"

def find_genotyped_individuals():
    genotyped_individuals = set()
    with open(AMISH_GENOTYPED_PEDIGREE, 'r') as geno_ped_file:
        for line in geno_ped_file:
            line = line.split()
            genotyped_individuals.add(line[1])
    return genotyped_individuals


GENOTYPED_INDIVIDUALS = find_genotyped_individuals()

def amish_id_to_bonsai_id(id: str):
    """
    translate original amish id (numbers) to bonsai id 
    (negative numbers for ungenotyped, positive for genotyped)
    """
    assert(type(id) == str)
    if id not in GENOTYPED_INDIVIDUALS:
        return "-" + id
    else:
        return id

def amish_to_bonsai_ped():
    """
    Translate ped file (amish) to a list of individuals (bonsai)

    fam file row structure: Individual ID, Paternal ID, Maternal ID, 
        Sex (1=male; 2=female; other=unknown)
    bonsai list item structure: [sex, age, father, mother] 
        (1 for male, 0 for female)
    """
    
    bonsai_ped = {}
    with open(AMISH_FULL_PEDIGREE, 'r') as full_ped_file:
        for line in full_ped_file:
            if "FATHER" in line:
                continue

            line = line.split()
            individual_id = int(amish_id_to_bonsai_id(line[0]))
            father = int(amish_id_to_bonsai_id(line[1]))
            mother = int(amish_id_to_bonsai_id(line[2]))
            sex = line[3]

            if sex == "1": # male from 1(ped) to 1 (bonsai)
                sex = 1
            else: # female from 2 (ped) to 0 (bonsai)
                sex = 2
            
            bonsai_indv = [sex, None, father, mother]
            bonsai_ped[individual_id] = bonsai_indv

    return bonsai_ped

def match_file_to_segments(match_file):
    """
    Match file: [1 103.0 1 468.1 21      10867286 16173940       rs1 rs1 4096    3.227   cM      0       1       1]
    Bonsai: ['1', '2', '21', 10867286, 11158632, False, .911326],
    """
    bonsai_ibds = []
    with open(match_file, 'r') as file:
        
        for line in file:
            line = line.split()
            hap1 = line[1]
            hap2 = line[3]
            chrom = line[4]
            start = int(line[5])
            end = int(line[6])
            length = float(line[10])

            indv1 = int(amish_id_to_bonsai_id(hap1.split(".")[0]))
            indv2 = int(amish_id_to_bonsai_id(hap2.split(".")[0]))

            segment = [indv1, indv2, chrom, start, end, False, length]
            bonsai_ibds.append(segment)
    
    return bonsai_ibds


def tests():
    bonsai_ibds = match_file_to_segments("/homes/thdang/trang_amish/ped-aware_pipeline/output/amish1_chr21_amr_geno_germline.match")
    print(bonsai_ibds[:5])