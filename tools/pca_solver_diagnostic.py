"""Standalone software regression for CI's crashing numerical dependency versions."""
import faulthandler,json
faulthandler.enable()
from vdjtools.cores import numerical_kernels,kernel_threads
import numpy as np,scipy,sklearn
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_info
print(json.dumps({'numpy':np.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'pools':threadpool_info()}),flush=True)
for repeat in range(4):
 for shape,k in [((200,700),6),((45,80),2),((1050,32),5),((80,550),70)]:
  x=np.random.RandomState(46).normal(size=shape)
  x[np.random.RandomState(47).uniform(size=shape)<.6]=0
  budget=1 if repeat%2 else kernel_threads(0)
  print(json.dumps({'repeat':repeat,'shape':shape,'threads':budget}),flush=True)
  with numerical_kernels(budget):
   p=make_pipeline(StandardScaler(),PCA(k,random_state=23)).fit(x)
   assert np.isfinite(p.transform(x)).all()
print('PCA_RUNTIME_PASS',flush=True)
