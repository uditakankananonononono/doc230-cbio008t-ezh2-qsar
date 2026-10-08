import pandas as pd,numpy as np,json,requests
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import AllChem,MACCSkeys,rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import KFold,GroupKFold
from sklearn.metrics import roc_auc_score,mean_squared_error
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
raw=pd.read_csv('/tmp/ezh2/ezh2_ic50.csv'); R={}
raw['mol']=raw.smiles.map(Chem.MolFromSmiles); raw=raw[raw.mol.notna()]
raw=raw[(raw.pIC50>3)&(raw.pIC50<11.5)]
g=raw.groupby('CID'); df=pd.DataFrame({'pIC50':g.pIC50.median(),'patent':g.patent.first(),'smiles':g.smiles.first(),'n_meas':g.size(),'n_assays':g.AID.nunique()}).reset_index()
df['mol']=df.smiles.map(Chem.MolFromSmiles)
R['d1']=dict(n_compounds=len(df),n_measurements=len(raw),n_patents=int(df.patent.nunique()),pIC50_mean=round(float(df.pIC50.mean()),2),pIC50_sd=round(float(df.pIC50.std()),2),top_patent_share=round(float(df.patent.value_counts(normalize=True).iloc[0]),3))
def fp(kind):
    out=[]
    for m in df.mol:
        if kind=='m2': v=AllChem.GetMorganFingerprintAsBitVect(m,2,2048); a=np.zeros(2048); DataStructs.ConvertToNumpyArray(v,a)
        elif kind=='m3c':
            v=AllChem.GetHashedMorganFingerprint(m,3,2048); a=np.zeros(2048)
            for k,c in v.GetNonzeroElements().items(): a[k]=c
        elif kind=='maccs': v=MACCSkeys.GenMACCSKeys(m); a=np.zeros(167); DataStructs.ConvertToNumpyArray(v,a)
        elif kind=='rdk': v=Chem.RDKFingerprint(m,fpSize=2048); a=np.zeros(2048); DataStructs.ConvertToNumpyArray(v,a)
        out.append(a)
    return np.array(out)
X=fp('m2'); y=df.pIC50.values; grp=df.patent.fillna('NA').values
def mk(n):
    return {'ridge':lambda:Ridge(alpha=10),'rf':lambda:RandomForestRegressor(300,max_features=0.33,n_jobs=2,random_state=0),'mlp':lambda:MLPRegressor(hidden_layer_sizes=(512,128),early_stopping=True,random_state=0,max_iter=300)}[n]()
def cvrun(X,y,splits,names=('ridge','rf','mlp')):
    P={n:np.zeros(len(y)) for n in names}
    for tr,te in splits:
        for n in names: P[n][te]=mk(n).fit(X[tr],y[tr]).predict(X[te])
    return P
def met(y,p): return dict(spearman=round(float(spearmanr(y,p)[0]),3),rmse=round(float(mean_squared_error(y,p)**0.5),3))
rs=list(KFold(5,shuffle=True,random_state=0).split(X)); P=cvrun(X,y,rs); R['d2']={n:met(y,P[n]) for n in P}
ls=list(GroupKFold(5).split(X,y,grp)); PL=cvrun(X,y,ls); R['d3']={n:met(y,PL[n]) for n in PL}
best=max(PL,key=lambda n:R['d3'][n]['spearman']); R['d3_best']=best; R['d3_gate']=bool(R['d3'][best]['spearman']>=0.5)
R['d4_gate_mlp_beats_rf']=bool(R['d3']['mlp']['spearman']-R['d3']['rf']['spearman']>=0.03)
# d5
fps=[AllChem.GetMorganFingerprintAsBitVect(m,2,2048) for m in df.mol]; mx=np.zeros(len(y))
for tr,te in ls:
    trf=[fps[i] for i in tr]
    for i in te: mx[i]=max(DataStructs.BulkTanimotoSimilarity(fps[i],trf))
bins={'>=0.7':mx>=0.7,'0.5-0.7':(mx>=0.5)&(mx<0.7),'<0.5':mx<0.5}
R['d5']={k:dict(n=int(v.sum()),**met(y[v],PL['rf'][v])) if v.sum()>10 else dict(n=int(v.sum())) for k,v in bins.items()}
# d6
act=(y>=7).astype(int); R['d6']={n:dict(auroc=round(float(roc_auc_score(act,PL[n])),3),top10_recall=round(float(act[np.argsort(-PL[n])[:len(y)//10]].sum()/act.sum()),3)) for n in PL}; R['d6']['base_rate']=round(float(act.mean()),3)
# d7
R['d7']={}
for k in ['m2','m3c','maccs','rdk']:
    Xk=X if k=='m2' else fp(k); R['d7'][k]=met(y,cvrun(Xk,y,ls,names=('rf',))['rf'])
# d8
mult=raw.groupby(['CID','AID']).pIC50.median().reset_index(); piv=mult.groupby('CID').filter(lambda d:d.AID.nunique()>=2)
pairs=[]
for cid,d in piv.groupby('CID'):
    v=d.pIC50.values; pairs.append((v[0],v[1]))
if len(pairs)>10:
    a,b=np.array(pairs).T; R['d8']=dict(n_compounds=len(pairs),spearman=round(float(spearmanr(a,b)[0]),3),mae=round(float(np.mean(abs(a-b))),3),rmse_between_assays=round(float(np.sqrt(np.mean((a-b)**2))),3))
# d9
names=['tazemetostat','GSK126','EPZ005687','CPI-1205']; R['d9']={}
for nm in names:
    try:
        r=requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{nm}/cids/JSON',timeout=30).json()['IdentifierList']['CID'][0]
    except Exception as e: R['d9'][nm]='lookup failed'; continue
    i=np.where(df.CID.values==r)[0]
    R['d9'][nm]=dict(CID=int(r),in_data=bool(len(i)),measured=round(float(y[i[0]]),2) if len(i) else None,LPO_pred_rf=round(float(PL['rf'][i[0]]),2) if len(i) else None)
# d10
sc=df.mol.map(lambda m:MurckoScaffold.MurckoScaffoldSmiles(mol=m)); sg=pd.factorize(sc)[0]
R['d10_n_scaffolds']=int(len(set(sg))); ss=list(GroupKFold(5).split(X,y,sg)); PS=cvrun(X,y,ss,names=('rf',)); R['d10_scaffold_rf']=met(y,PS['rf'])
json.dump(R,open('results/results.json','w'),indent=1); print(json.dumps(R,indent=1))
