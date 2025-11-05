# FriendlyClearMap-scripts

This is the collection of scripts to build your own Dockerfiles for ClearMap1 and ClearMap2's CellMap component, as outlined in our 2022 GigaScience publication: 
Negwer, M., Bosch, B., Bormann, M., Hesen, R., Lütje, L., Aarts, L., Rossing, C., Nadif Kasri, N., & Schubert, D. (2022). 
FriendlyClearMap: an optimized toolkit for mouse brain mapping and analysis. GigaScience, 12. 
https://doi.org/10.1093/gigascience/giad035 

ClearMap 1 is as originally described in: 
Renier et al. Cell 2016. Mapping of Brain Activity by Automated Volume Analysis of Immediate Early Genes. Cell, 165(7), 1789–1802. https://doi.org/10.1016/j.cell.2016.05.007

 Code from: https://github.com/ChristophKirst/ClearMap and adapted as follows:
  - ported to Python 3, adapted to include ImageJ's BigWarp landmark output (see here for an earlier stage of this work: https://github.com/BrammBosch/clearmap) 
  - Adapted for Dockerization 
  
  
ClearMap 2 as originally described in:
Kirst et al., Cell 2020. Mapping the Fine-Scale Organization and Plasticity of the Brain Vasculature. Cell, 180(4), 780–795.e25. https://doi.org/10.1016/j.cell.2020.01.028
  
Code from: https://github.com/ChristophKirst/ClearMap2 and adapted the  CellMap portion (the TubeMap portion is not included) as follows: 
  - Adapted the Elastix registration to include ImageJ's BigWarp landmark output
  - Adapted for Dockerization
  - Included an adapted version of ClearMap1's Ilastik Cell Detection mode 
  
Please see the Appendices 1-4 of our Gigascience publication for detailed instructions on how to use the pipelines. 
dx.doi.org/10.17504/protocols.io.eq2lynnkrvx9/v2
dx.doi.org/10.17504/protocols.io.yxmvmn9pbg3p/v2
dx.doi.org/10.17504/protocols.io.dm6gpbdwdlzp/v1 
dx.doi.org/10.17504/protocols.io.36wgq77m5vk5/v1

All raw and processed data is available on GigaDB, linked to the GigaScience article. Please see here: https://doi.org/10.5524/102385

Specifically, find the docker containers here: 
https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102385/docker_containers.tar.gz

And the data here: 
Parvalbumin P14: https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102385/dataset_01-1_PV_P14.tar.gz
Parvalbumin P56: https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102385/dataset_01-2_PV_P56.tar.gz
Somatostatin P28: https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102385/dataset_02_SST_P28.tar.gz
VIP P28: https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102385/dataset_03_VIP_P28.tar.gz

---
Other work used in this article: 

- Ilastik: Berg et al. ilastik: interactive machine learning for (bio)image analysis. Nat Methods 16, 1226–1232 (2019). https://doi.org/10.1038/s41592-019-0582-9

- BigWarp: Bogovic et al. (2016). Robust registration of calcium images by learned contrast synthesis. In 2016 IEEE 13th International Symposium on Biomedical Imaging (ISBI). IEEE. https://doi.org/10.1109/isbi.2016.7493463

- BrainRender: Claudi et al. (2021) Visualizing anatomically registered data with brainrender eLife 10:e65751 https://doi.org/10.7554/eLife.65751

- Young Mouse Brain Atlases: Newmaster et al. 2020. Quantitative cellular-resolution map of the oxytocin receptor in postnatally developing mouse brains. Nat Commun 11, 1885 (2020). https://doi.org/10.1038/s41467-020-15659-1

- Elastix: 
  Klein et al., 2010. elastix: a toolbox for intensity based medical image registration, IEEE Transactions on Medical Imaging, vol. 29, no. 1, pp. 196 - 205, January 2010. http://dx.doi.org/10.1109/TMI.2009.2035616
  Shamonin et al., 2014. Fast Parallel Image Registration on CPU and GPU for Diagnostic Classification of Alzheimer's Disease, Frontiers in Neuroinformatics, vol. 7, no. 50, pp. 1-15, January 2014. http://dx.doi.org/10.3389/fninf.2013.00050
  
  - Cortical Flatmap projections: Wang et al. 2020: The Allen Mouse Brain Common Coordinate Framework: A 3D Reference Atlas. Cell, Volume 181, Issue 4, 936 - 953.e20. https://doi.org/10.1016/j.cell.2020.04.007 
https://ccf-streamlines.readthedocs.io/en/latest/guide.html
https://github.com/AllenInstitute/ccf_streamlines
