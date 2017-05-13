import sys
import decimal
import re
import matplotlib.pyplot as plt

mass_h20 = decimal.Decimal('18.01528')
mass_h = decimal.Decimal('1.007276466')
mass_co = decimal.Decimal('28.0101')
mass_oh = decimal.Decimal('17.008')
mass_nh2 = decimal.Decimal('16.02258')

masses = open("aaMass.txt")
#generate dictionary of amino acid residue mass
#average (not monoisotopic) masses
massDict = {}
for m in masses:
    m = m.split('\t')
    massDict[m[0]] = decimal.Decimal(m[4])
#check if masses are correctly assigned
#for m in massDict.keys():
#    print(m + ": " + str(massDict[m]))  

def calc_tmass(protein):
    prot_mass = decimal.Decimal(0)
    #add mass of each residue
    for r in protein.upper():
        prot_mass = prot_mass + massDict[r]
    #add mass of water
    prot_mass = prot_mass + mass_h20
    return prot_mass

def sum_res(peptide):
    pep_mass = decimal.Decimal(0)
    for r in peptide.upper():
        pep_mass = pep_mass + massDict[r]
    return pep_mass

def calc_mzmass():
    mToz = decimal.Decimal("777.25439453")
    z = decimal.Decimal("7")
    parent_mass = (mToz * z) - (z * mass_h)
    return parent_mass

def gen_peptides(protein):
    if protein[-1] == '\n':
        protein = protein[:-1]
    #a is length of peptide
    pep_list = []
    for a in range(1, len(protein)):
        #b is start
        for b in range(0, len(protein)-a+1):
            pep_list.append(protein[b:b+a])
    return pep_list

def checkAgainstAccurate(acc_mass, peptides):
    err = decimal.Decimal(".00005")
    #err is for m/z 
    #multiply by 7 to account for charge
    mass_lower = acc_mass - (acc_mass * err)*7
    mass_upper = acc_mass + (acc_mass * err)*7

    out = []
    for pep in peptides:
        tmass = calc_tmass(pep)
        if mass_lower <= tmass and tmass <= mass_upper:
            out.append(pep)
    return out

def gen_spectrum(peptide):
    #dictionary mapping m/z value to peptide seq
    spectrum = {}
    for a in range(1, len(peptide)):
        n_term = peptide[:a]
        c_term = peptide[a:]
        #generate prefix ions (a,b,c)
        for q in range(1, min(len(n_term)+1,7)):
            Q = d(q)
            res = sum_res(n_term)
            a_ion = (res - mass_co + (mass_h*2) + (mass_h*Q))/Q
            b_ion = (res + (mass_h*Q))/Q
            c_ion = (res + mass_nh2 + (mass_h*Q))/Q
            for ion in [a_ion, b_ion, c_ion]:
                spectrum[str(ion)] = n_term 
        for q in range(1, min(len(c_term)+1,7)):
            res = sum_res(c_term)
            x_ion = (res + mass_co + mass_oh + (mass_h*(Q-1)))/Q
            y_ion = (res + mass_oh + mass_h + (mass_h*Q))/Q
            z_ion = (res - mass_nh2 + mass_oh + mass_h + (mass_h*(Q-1)))/Q
            for ion in [x_ion, y_ion, z_ion]:
                spectrum[str(ion)] = c_term
    return spectrum

def score_match(spec_list, mz_list, inten_list):
    for cand_list in spec_list:
        score = decimal.Decimal(0)
        sc1 = d('0')
        sc2 = d('0')
        hits = 0
        for a in range(0, len(mz_list)):
            mz = mz_list[a]
            err = decimal.Decimal(50) / decimal.Decimal(10**6)
            mz_lower = mz - d(2)*(err*mz)
            mz_upper = mz + d(2)*(err*mz)
            #print(mz_lower)
            #print(mz_upper)
            #print()
            b = 0
            while b<len(cand_list):
                if mz_lower <= d(cand_list[b]) <= mz_upper:
                    score += mz * inten_list[a]
                    sc1 += mz
                    sc2 += inten_list[a]
                    b = len(cand_list)
                    hits += 1
                b += 1
        print(score)
        print(sc1)
        print(sc2)
        print(hits/len(mz_list))
        print()

def d(string):
    return decimal.Decimal(string)

def gen_acc_spec(infile):
    data = re.findall('\n[0-9\.]+\s[0-9\.]+', infile)
    mz_list = []
    inten_list = []
    for d in data:
        mz_list.append(decimal.Decimal(d.split(' ')[0][1:]))
        inten_list.append(decimal.Decimal(d.split(' ')[1]))
    return [mz_list, inten_list]

def norm(inten_list):
    tot = sum(inten_list)
    out = []
    for i in range(0, len(inten_list)):
        out.append(inten_list[i] / tot)
    return out

prots = open(sys.argv[1])

accurate_mass = calc_mzmass()
print("Parent mass: " + str(accurate_mass) + '\n')

candidate_pep = []

for prot in prots:
    if '>' not in prot and prot != '\n':
        peps = gen_peptides(prot)
        candidate_pep += checkAgainstAccurate(accurate_mass, peps)


#get list of m/z values and intensities
data = open(sys.argv[2]).read()
mz_list, inten_list = gen_acc_spec(data)
inten_list_normed = norm(inten_list)

#list to hold ideal spectrums of candidate peptides 
spec_list = []

cand_spec = open("idealSpectrums.txt", "w")

for c in candidate_pep:
    spec = gen_spectrum(c)
    specKeys = list(spec.keys())
    specKeys.sort(key=decimal.Decimal)
    spec_list.append(specKeys)
    cand_spec.write(c + '\n')
    for s in specKeys:
        cand_spec.write(s + '\n')
    cand_spec.write('\n')
    print(c)
    score_match([specKeys], mz_list, inten_list_normed)

mz_float = [float(m) for m in mz_list]
inten_float = [float(i) for i in inten_list]
plt.bar(mz_float, inten_float)
plt.xlabel("m/z")
plt.ylabel("Intensity")
plt.show()

