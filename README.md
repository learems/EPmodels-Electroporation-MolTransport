# EPmodels-Electroporation-MolTransport
Models coupling electroporation and molecular transport at the single cell level. 

These are the models presented in [Miklavcic and Towhidi 2010](https://doi.org/10.2478/v10019-010-0002-3), [Li and Lin 2011](https://doi.org/10.1016/j.bioelechem.2011.04.006), and [Smith 2011](https://dspace.mit.edu/handle/1721.1/63085), which have been used in the study of [Scuderi et al.](http://dx.doi.org/10.1016/j.bioelechem.2022.108216). 

MT2010 model is implemented directly in Comsol. See the file MT2010_nosolutions.mph. 

S2011 and LL2011  models need to be created and run from Comsol with Matlab. 
For S2011 and LL2011 models adapt the Main.m script to run the simulations. 
After completing the simulations, see examples in the ExtractData.m script to extract the data from the models via Comsol with Matlab. 

If you use these models, please cite the following paper:
Scuderi M, Dermol-Černe J, Amaral da Silva C, Muralidharan A, Boukany PE, Rems L. **Models of electroporation and the associated transmembrane molecular transport should be revisited**. Bioelectrochemistry 147: 108216, 2022. doi: [10.1016/j.bioelechem.2022.108216](http://dx.doi.org/10.1016/j.bioelechem.2022.108216)

For additional info write to lea.rems@fe.uni-lj.si. 
