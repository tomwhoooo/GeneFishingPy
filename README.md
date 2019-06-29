# GeneFishingPy

Please download .RData into /data folder from https://www.synapse.org/#!Synapse:syn18434549/files/
and then use data_processing.R to process the multiple data into data/csv_data (the directory will be created automatically).

You can run it in the following way:

```python GeneFishingPy.py /data/csv_data /data/c1.csv /CFR [Options]```

The first argument is the *directory* where all csv files are stored. The second argument is the path to your bait csv. The third one is where the output will be automatically saved to. They will be named as ```[tissue_name]_CFR.csv```. One sample output file has been provided.

If you would like to tune parameters, please use ```python GeneFishingPy.py /data/csv_data /data/c1.csv /CFR -h``` to find out what parameters you may adjust, what do they mean and their default values.

An example would be:

```python GeneFishingPy.py /data/csv_data /data/c1.csv /CFR --number_of_iteration 100 --way_of_cluster K-means``` 
to decrease the number of iterations from 1000 to 100 and change the clustering method as pure k-means in order to save computational time.
