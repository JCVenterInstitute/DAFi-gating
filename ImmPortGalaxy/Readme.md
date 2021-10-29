# Reproducible Analytics for Secondary Analyses of ImmPort Vaccination-Related Cytometry Data

The project is focused on implementation of cutting edge computational analysis methods and applications of these methods for re-analysis of publicaly avaialble flow cytometry data in ImmPort influenza studies (https://www.immport.org).

This GitHub folder contains the Galaxy workflows we created for automated FLOCK and DAFi analysis of ImmPort flow cytometry data. Note that these workflows include ImmPortGalaxy accessory tools such as FCSTrans, RearrangeHeader, and GenerateOverview. They can only be run on an ImmPortGalaxy instance.

### Main contributors:
Yu "Max" Qian, Aishwarya Mandava, Mehmet Kuscuoglu, Yun "Renee" Zhang, and Richard H. Scheuermann.

Informatics Department, J. Craig Venter Institute

### Copyright:
Authors and J. Craig Venter Institute

### Funding source:
US NIAID UH2AI132342 (PI: Yu "Max" Qian, PhD; Co-Investigator: Richard H. Scheuermann, PhD)

Project period: July 2018 to June 2021

### Example data analysis results and Galaxy tools implemented
https://am794.github.io/ImmPort_DataAnalysis/index.html

dafi_intel.tgz: XML file of DAFi tool for Galaxy installation

Galaxy-Workflow-DAFI_FLOCK_TXT.ga: Galaxy installation file for DAFi+FLOCK workflow with TXT files as input

Galaxy-Workflow-DAFI_FLOCK_COMP.ga: Galaxy installation file for DAFi+FLOCK workflow with support for external compensation

Galaxy-Workflow-DAFI_FLOCK.ga: Galaxy installation file for DAFi+FLOCKv1 workflow

Galaxy-Workflow-DAFI_FLOCK_2.ga: Galaxy installation file for DAFi+FLOCKv2 workflow

UH2_ImmPort_Influenza_Studies_Metadata.xlsx: metadata of 55 ImmPort influenza vaccine studies, organized based on studies, reagent panels, and cell types

### Manuscript under review
https://www.medrxiv.org/content/10.1101/2021.09.14.21263182v1 

For questions, feel free to contact: mqian@jcvi.org or qianyu.cs@gmail.com
