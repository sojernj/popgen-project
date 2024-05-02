# popgen-project
Repo for the exam project in PopulationGenetics

# Project Plan: Exploring Non-African Archaic Genomic Segments

## 1. Data Collection and Setup
- **Data Sources**: 
  - `ArchaicSegments.txt`: Contains information on segments of archaic origin, including details about the population, region, start and end positions, length, and probability of archaic content.
  - `SNP.txt`: Lists non-African SNPs, indicating whether they are shared with high-coverage archaic genomes (Altai or Vindija Neanderthals, Denisova), or if they are classified as "human".
- **Environment Setup**: Ensure you have the necessary tools and libraries for data analysis. Python with Pandas and visualization libraries like Matplotlib or Seaborn could be helpful.

## 2. Data Preprocessing and Quality Control
- **Data Import**: Load both data sets into your analysis environment. Convert them into a suitable format if needed (e.g., Pandas data frames).
- **Data Cleaning**: Check for and address missing or inconsistent data in both files. Ensure consistent column naming and structure.

## 3. Analyzing Archaic Segments
- **Distribution of Archaic Segments**: 
  - Examine the distribution of archaic segments across chromosomes, individuals, populations, and regions. Calculate summary statistics like mean segment length and number of segments per individual.
  - Create visualizations (e.g., density plots, scatter plots) to depict the distribution of these segments.
- **Shared Archaic Content**: 
  - Analyze the shared content with Altai Neanderthals, Denisovans, and Vindija Neanderthals. Determine if certain populations or regions have higher sharing with specific archaic genomes.
  - Investigate patterns in "Shared_with_Altai", "Shared_with_Denisova", and "Shared_with_Vindija" columns to identify notable trends.

## 4. Investigating SNP Data
- **Archaic SNP Analysis**: 
  - Examine the `SNP.txt` data to identify which SNPs are shared with archaic genomes and which are uniquely human. Determine their distribution across the chromosome and population groups.
  - Investigate correlations between archaic SNPs and their corresponding archaic segments.
- **Adaptive Introgression**: Explore whether specific regions or SNPs suggest evidence of adaptive introgression. Focus on those with high probability of archaic origin or associated with known adaptive traits (e.g., EPAS1).

## 5. Exploring Specific Regions and Genes
- **EPAS1 Gene Region**: 
  - Extract the data for the 1 Mb window surrounding the EPAS1 gene. Compare this to other regions to determine if it has unique archaic segment characteristics.
  - Analyze the total length of archaic segments in this region and correlate it with other genomic regions.
- **Other Genes/Regions of Interest**: Identify additional genes or regions with high archaic content. Investigate their significance and potential connection to adaptive traits.

## 6. Interpretation and Reporting
- **Results Interpretation**: Summarize the key findings from the analyses. Compare them with existing studies on archaic introgression and human evolutionary history.
- **Additional Analyses**: Based on your findings, consider further analyses to explore new questions or unexpected trends.
- **Visualizations and Plots**: Create visual representations of your results to aid in interpretation and communication. Consider bar plots, scatter plots, heatmaps, etc.
- **Report Writing**: Compile the results, interpretations, and additional analyses into a comprehensive report. Include references to relevant studies and previous research.
- **Prepare Presentation**: Develop a presentation to communicate your findings to peers, colleagues, or at conferences.
