import os, re, csv
def qmax(cif):
    q=[]
    for l in open(cif,encoding='utf-8',errors='ignore'):
        s=l.split()
        if len(s)==8 and re.match(r'^[A-Z][a-z]?$',s[0]) and re.match(r'^[A-Z][a-z]?\d',s[1]):
            try: q.append(abs(float(s[7])))
            except: pass
    return max(q) if q else None
out=open(os.path.expanduser('~/qmax_7820.csv'),'w',newline='')
w=csv.writer(out); w.writerow(['nama','qmax'])
n=0
for f in sorted(os.listdir(os.path.expanduser('~/cif_semua'))):
    if not f.endswith('.cif'): continue
    v=qmax(os.path.expanduser('~/cif_semua/'+f))
    if v is not None:
        w.writerow([f[:-4], round(v,4)]); n+=1
out.close()
print('ditulis', n)
