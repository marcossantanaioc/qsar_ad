# qsar_ad
> Calculate the applicability domain.


## How to use

You can use any function to generate descriptors for the molecules in the dataset. For instance, we could use morgan fingerprints from RDkit to generate a vector of 2048 bits for each molecule. 

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from qsar_ad.applicability_domain import kNNDomain
from nbdev.showdoc import show_doc
```

```python
mol_list = ['CCCCCCC','c1ccccc1', 'c1ccccn1', 'c1cc(O)ccn1', 'c1cc(Cl)c(C(=O))cn1']
```

```python
Xtrain = np.array([AllChem.GetMorganFingerprintAsBitVect(x, radius=2048) for x in list(map(Chem.MolFromSmiles, mol_list))])[1:]
```

```python
Xtrain.shape
```




    (4, 2048)



```python
Xsample = Xtrain[0][None,:]
Xsample.shape
```




    (1, 2048)



## Calculating AD

```python
ad_domain = kNNDomain(Xtrain)
```

```python
avg_distance = ad_domain.calculate_applicability_domain(Xsample)
```

```python
ad_domain.ad_threshold # Threshold of the applicability domain. A compound with Dc higher than the threshold is considered out of the domain.
```




    4.364247092576795



```python
avg_distance
```




    (array([3.06075884]), array([ True]))


