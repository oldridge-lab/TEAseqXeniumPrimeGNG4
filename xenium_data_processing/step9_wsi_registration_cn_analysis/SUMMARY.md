Categories	|   tasks	        |   env	     |   notebooks/scripts

CN analysis	| CN identification	|  neighbor	 |   neighborhood_identification.ipynb
            | GC merge	        |  neighbor	 |   neighborhood_merge_GC.ipynb
			
				
Xenium-hae  | registration	    |   wsireg	 |   wsireg.py
	        | validation	    |  stardist	 |   registration_validation.ipynb
			
			
* yml files of all env used are saved in ./env		

Main packages used in each env/scripts:

1. env: neighbor
neighborhood_identification.ipynb
neighborhood_merge_GC.ipynb
- pandas=2.1.4
- numpy=1.26.4
- scikit-learn==0.22.1
- seaborn==0.9.0
- matplotlib==3.5.3
- scikit-image=0.19.3
- scipy=1.7.3

2. env: wsireg
wsireg.py
- pandas==2.2.3
- wsireg==0.3.10

3. env: stardist
registration_validation.ipynb
- numpy=1.26.4
- pandas=2.1.4
- scipy=1.10.1
- tifffile=2025.5.10
- rasterio=1.3.11
- matplotlib==3.10.7