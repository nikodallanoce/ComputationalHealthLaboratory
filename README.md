# ComputationalHealthLaboratory
Computational Health Laboratory project for the a.y. 2021/2022.
## Group Members
- [Niko Dalla Noce](https://github.com/nikodallanoce)
- [Alessandro Ristori](https://github.com/RistoAle97)
- [Andrea Zuppolini](https://github.com/AndreZupp)
## Gene Set Pathway Merging and Analysis
Starting from one or more genes, extract from interaction databases the genes they interact with. Using the expanded gene set, perform pathway analysis and obtain all disease pathways in which the genes appear. Merge the pathways to obtain a larger graph. Perform further network analysis to extract central biomarkers and communities beyond pathways. Compute a distance between the initial gene set and the various pathways (diseases).
## Repository structure
```bash
📂ComputationalHealthLaboratory
├── 💼0_Pathway_Enrichment.ipynb  # Pathway gene dataset expansion and pathway enrichment
├── 💼1_Network_Analysis.ipynb  # Network building and analysis
├── 💼2_Community_Analysis.ipynb  # Community detection and analysis
├── 💼3_Plots.ipynb  # Methods to plot the protein, disease and community graphs
├── 💼4_Project_CHL.ipynb  # Entire project, the previous four notebooks combined
├── 📄config_example.yml  # Replace this with your customized configuration file
├── 📄config.py  # Method to retrieve data from BioGRID
├── 📂datasets  # Datasets used by the project
│   ├── 🗃️biomarkers.csv  # Central nodes
│   ├── 🗃️communities.csv  # Communities of the protein-to-protein graph
│   ├── 🗃️communities_metrics.csv
│   ├── 🗃️community_gene_metrics.csv
│   ├── 🗃️diseases_pathways.csv  # Disease pathways retrieved from DisGeNET
│   ├── 🗃️diseases_scores.csv  # Disease pathways with their metrics
│   ├── 🗃️genes.csv  # Expanded gene dataset
│   ├── 🗃️geneset.csv  # Starting gene interactions, retrieved by BioGRID
│   ├── 🗃️interactions.csv  # Expanded gene interactions dataset
│   ├── 🗃️mean_distances.csv
│   ├── 🗃️modularities.csv
│   └── 🗃️protein_graph.gpickle  # Protein-to-protein graph
├── 📄README.md
├── 📄requirements.txt
└── 📂src  # Project methods
    ├── 📄communities.py
    ├── 📄disease.py
    ├── 📄plot_graphs.py
    ├── 📄protein_to_protein_graph.py
    └── 📄utilities.py
```
## How to run the project
First, clone the repo
```
git clone https://github.com/nikodallanoce/ComputationalHealthLaboratory
```
Install all the required packages
```
pip install -r requirements.txt
```
Then you can work with the notebooks and our package, for a deeper understanding of our work, use **4_Project_CHL.ipynb** to run the entire project.
## Resources
- [BioGRID](https://thebiogrid.com/)
- [DisGeNET](https://www.disgenet.org/)
- [GSEApy](https://github.com/zqfang/GSEApy)
- [networkx](https://github.com/networkx/networkx)
