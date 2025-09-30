# The code for the genome assembly and analysis of cointree

## Use the WGDI prescribed workflow to invoke the configuration file.
### 1.Generate an inter-species dot plot using the WGDI pipeline with the -d configuration option.
```
wgdi -d 1.d.conf
```
### 2.Perform collinearity calculation using the -icl parameter.
```
wgdi -icl 2.icl.conf
```
### 3.-bi is used for the identification of syntenic blocks.
```
wgdi -bi 3.bi.conf
```
### 4.-c parameter is used to filter the results of syntenic blocks.
```
wgdi -c 4.c.conf
```
### 5.-ak parameter is used to specify the ancestral karyotype file for chromosome evolution inference.
```
wgdi -ak 5.ak.conf
```
### 6.

