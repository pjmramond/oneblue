# oneblue
Processing metagenomics and metabarcoding data from the project ONEBLUE (https://one-blue.eu/).

## 1/ Import on MARBITS
For now I created a project folder in my <emm1> account on MARBITs:
```
/mnt/smart/scratch/emm1/users/pramond/oneblue
```

Metagenomics data were processed by MACROGEN. A first batch of samples (seawater samples, 26) was received June 18th, 2026.
The transfer was performed directly from the hard drive that MACROGEN sent us to MARBITS ("only" 500GB).
We checked that the transfer was made correctly with their md5sum:
```
cd /mnt/smart/scratch/emm1/users/pramond/oneblue/MACROGEN/HN00276367/0.RawData
md5sum *.fastq.gz > marbits_md5
```

