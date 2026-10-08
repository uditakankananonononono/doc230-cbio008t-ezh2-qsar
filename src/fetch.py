import numpy as np
import pandas as pd,requests,time,re
d=pd.read_csv('/tmp/ezh2/conc.csv',low_memory=False)
e=d[(d['Target GeneID']==2146)&d.CID.notna()&d['Activity Value [uM]'].notna()&(d['Activity Name']=='IC50')].copy()
e['CID']=e.CID.astype(int)
bad=e['Assay Name'].str.contains(r'\bEED\b|cell|prolifer|viability|proliferation|Karpas|WSU|Pfeiffer',case=False,regex=True)
e=e[~bad]
e['patent']=e['Assay Name'].str.extract(r'(US\d{7,8})')[0]
e=e[e.patent.notna()]
e=e[(e['Activity Value [uM]']>0)]
e['pIC50']=-np.log10(e['Activity Value [uM]']*1e-6) if False else 6-np.log10(e['Activity Value [uM]'])
cids=sorted(e.CID.unique()); smi={}
for i in range(0,len(cids),200):
    b=cids[i:i+200]
    for _ in range(3):
        r=requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/IsomericSMILES/JSON',data={'cid':','.join(map(str,b))},timeout=60)
        if r.status_code==200: break
        time.sleep(3)
    for p in r.json()['PropertyTable']['Properties']: smi[p['CID']]=p.get('IsomericSMILES') or p.get('SMILES')
e['smiles']=e.CID.map(smi); e=e[e.smiles.notna()]
e[['AID','CID','patent','Assay Name','pIC50','smiles']].to_csv('/tmp/ezh2/ezh2_ic50.csv',index=False)
print(len(e),e.CID.nunique(),e.patent.nunique(),e.AID.nunique())
