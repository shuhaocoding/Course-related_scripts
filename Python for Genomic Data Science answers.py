# final exam for Python for Genomic Data Science
f = open('/Users/shuhao/Desktop/dna2.fasta')

count = 0
d1= {}
d2= {}
l= {}
orf= {}
name = ''
seq = ''
import re
def findorf(s):
    return max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',s), key = len)

for line in f:
    if line.startswith('>'):
        count += 1
        if name != '':
            d1[name] = seq
            d2[seq] = name
            l[name] = len(seq)
            orf[name] = findorf(seq)
        seq = ''
        name = line[1:]
    else:
        seq += line.strip()
else:
    d1[name] = seq
    d2[seq] = name
    l[name] = len(seq)
print(count)

for i in d1.values():
    if len(i) == max(l.values()):
        print('longest sequence is', d2[i], '\n', i)
        print(len(i))

for i in d1.values():
    if len(i) == min(l.values()):
        print('shortest sequence is', d2[i], '\n', i)
        print(len(i))

maxorf = max(orf.values(), key=len)
maxorflen = len(maxorf)
for n in d1.keys():
    if orf[n] == maxorf:
        maxorfname = n

print('longest orf is', maxorflen, 'nt long from', maxorfname)


n2 =
print('longest orf for this seq is', orf[n2])
print('starting at', n2.index(orf[n2]))

n = 12
nmers = []
repeats = []
for seq in d1.values():
    for i in range(len(seq)-n+1):
        if seq[i:i+n] in nmers and seq[i:i+n] not in repeats:
            repeats.append(seq[i:i+n])
        if seq[i:i+n] not in nmers:
            nmers.append(seq[i:i+n])

allstring = ' '.join(d1.values())
rc = {}
for i in repeats:
    count2 = 0
    for j in range(len(allstring)-n+1):
        if allstring[j:j+n] == i:
            count2+=1
    rc[i] = count2
Max = max(rc.values())
for i in repeats:
    if rc[i] == Max:
        print('the most frequent repeat is', i)
print(Max)
