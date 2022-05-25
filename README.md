# qsar_ad
> Calculate the applicability domain.


## How to use

You can use any function to generate descriptors for the molecules in the dataset. For instance, we could use morgan fingerprints from RDkit to generate a vector of 2048 bits for each molecule. 

```
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from qsar_ad.applicability_domain import kNNDomain
```


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    <ipython-input-1-50c2e0bb6a11> in <module>
          3 import numpy as np
          4 import pandas as pd
    ----> 5 from qsar_ad.applicability_domain import kNNDomain
    

    ModuleNotFoundError: No module named 'qsar_ad'


```
mol_list = ['CCCCCCC','c1ccccc1', 'c1ccccn1', 'c1cc(O)ccn1', 'c1cc(Cl)c(C(=O))cn1']
```

```
Xtrain = np.array([AllChem.GetMorganFingerprintAsBitVect(x, radius=2048) for x in list(map(Chem.MolFromSmiles, mol_list))])[1:]
```

```
Xtrain.shape
```

```
Xsample = Xtrain[0][None,:]
Xsample.shape
```

```
kNNDomain??
```

## Calculating AD

```
ad_domain = kNNDomain(Xtrain)
```

```
avg_distance = ad_domain.calculate_applicability_domain(Xsample)
```

```
ad_domain.ad_threshold # Threshold of the applicability domain. A compound with Dc higher than the threshold is considered out of the domain.
```




    4.364247092576795



```
avg_distance
```




    (3.060758843123258, True)


